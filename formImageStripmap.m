clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stripmap SAR image formation from range-compressed pulses and        %
% INS/GPS navigation (lat, lon, alt, roll, pitch, yaw, NED velocity).  %
%                                                                      %
%  1. navToPulses:      nav records -> antenna phase center position,  %
%                       velocity and beam pointing at every pulse      %
%                       (Hermite interpolation using the velocity, so  %
%                       the nav rate can be slower than the PRF)       %
%  2. makeGeoImageGrid: ground grid along the flight track, on the     %
%                       WGS-84 ellipsoid or a DEM                      %
%  3. bpNav:            backprojection, each pulse gated to the beam   %
%                       footprint given by the measured attitude       %
%                                                                      %
% Set use_simulated = 0 and fill in section 1 to process real data.    %
% Requires MATLAB R2016b+ or GNU Octave; no toolboxes.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(fullfile(fileparts(mfilename('fullpath')), 'navutils'));
c = 299792458;

%% 1. DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_simulated = 1;

if use_simulated
    p.mode = 'stripmap';
    p.ref = struct('lat', 39.9990, 'lon', -83.0150, 'alt', 220);  % scene center
    p.fc = 9.6e9;  p.B = 200e6;  p.fs = 250e6;  p.Nr = 512;
    p.PRF = 500;   p.speed = 120;  p.heading = 35;   % 0.24 m pulse spacing
    p.H = 3000;    p.Yc = 4000;                      % AGL height, ground range
    p.nav_rate = 50;                                 % nav interval 20 ms vs PRT 2 ms
    p.motion = struct('cross_amp', 2, 'cross_T', 6, 'up_amp', 1, 'up_T', 4);
    p.att = struct('roll_amp', 2, 'roll_T', 3, 'pitch_amp', 1, 'pitch_T', 5, ...
                   'yaw_amp', 0.5, 'yaw_T', 7, 'crab', 0.5);
    p.antenna = struct('lever_arm', [1.5; 0.3; 0.8], 'squint', 0, 'D', 1.0);
    p.targets = [   0    0 0 1.0;      % [along cross height amp], m
                 0.75    0 0 0.7;      % pair 0.75 m apart along track
                  -50  -20 0 1.0;
                   50   20 0 1.0;
                  -30   25 0 1.0;
                   60  -25 0 1.0 ];
    p.SNR_dB = 0;
    % Residual INS position error the processor does not know about
    % (what autofocus must remove); set amp to [0 0 0] for perfect nav
    p.nav_err = struct('amp', [0.03 0.03 0.02], 'T', [0.7 1.6 3.3]);
    p.nav_single = true;          % INS reports lat/lon in single precision
                                  % (position jitter of ~0.2 m)
    p.nav_time_offset = 0.02;     % Nav time tags 20 ms late vs. radar clock
    [nav, radar, antenna, truth] = simulateNavSARData(p);
else
    % ---- LOAD YOUR DATA HERE --------------------------------------------
    % nav  (vectors, nav clock):  t (s), lat, lon (deg, WGS-84),
    %       alt (m above ellipsoid), roll, pitch, yaw (deg, Z-Y-X, yaw from
    %       true north), vN, vE, vD (m/s)
    % radar: rc (Nr x Np range-compressed complex pulses),
    %        t_pulse (1 x Np transmit times on the nav clock, s),
    %        fs (range sample rate, Hz), tau0 (two-way delay of sample 1,
    %        s, scalar or 1 x Np), fc (Hz), B (Hz),
    %        phase_sign (+1 if a target at delay tau has phase
    %        exp(-j*2*pi*fc*tau))
    % antenna: lever_arm (INS->antenna phase center, body x fwd / y right /
    %        z down, m), look ('right'|'left'), depression (deg),
    %        squint (deg, + toward nose), D (antenna length, m)
    load('my_stripmap_data.mat');     % <- replace with your loader
end

%% 2. NAVIGATION -> PULSE GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ground_h = 220;               % Terrain height above ellipsoid (m), or a
                              % DEM handle @(lat,lon) ... (see section 3)
nav_opts.position_source = 'velocity';  % Trajectory from integrated INS
                              % velocity (smooth); 'nav' = interpolate the
                              % reported positions (use only if they are
                              % smooth to a few mm)
nav_opts.anchor = 'linear';   % Tie integrated track to nav positions:
                              % 'first' | 'mean' | 'linear' (fits drift)
nav_opts.interp = 'hermite';  % 'hermite' (uses velocity) or 'linear'
nav_opts.time_offset = 0;     % Radar-to-nav clock offset (s): pulse n is
                              % at nav time t_pulse(n) + time_offset

% Local ENU origin: middle of the nav track, at ground height
mid = round(numel(nav.t)/2);
ref = struct('lat', nav.lat(mid), 'lon', nav.lon(mid), 'alt', ground_h);
geom = navToPulses(nav, radar.t_pulse, ref, antenna, nav_opts);

%% 3. IMAGE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scene.center  = [];           % struct('lat',..,'lon',..); [] = beam center
                              % on the ground at mid-aperture
scene.heading = [];           % deg from north; [] = mean flight track
scene.along_vec = -80:0.25:80;   % Along-track pixel positions (m)
scene.cross_vec = -35:0.25:35;   % Cross-track pixel positions (m), + = right
scene.height = ground_h;      % scalar, Ny x Nx matrix, or @(lat,lon) DEM

if isempty(scene.heading)
    v = mean(geom.vel, 2);
    scene.heading = atan2(v(1), v(2))*180/pi;
end
if isempty(scene.center)
    n = round(size(geom.apc,2)/2);
    if isa(scene.height,'function_handle'), h0 = ground_h; else h0 = scene.height(1); end
    s = (h0 - ref.alt - geom.apc(3,n)) / geom.bore(3,n);
    [la, lo] = nav_enu2geodetic(geom.apc(:,n) + s*geom.bore(:,n), ref);
    scene.center = struct('lat', la, 'lon', lo);
end
if use_simulated, scene.center = p.ref; end   % so target positions are known

grid = makeGeoImageGrid(ref, scene.center, scene.heading, ...
    scene.along_vec, scene.cross_vec, scene.height);
fprintf('Scene center %.6f, %.6f; grid heading %.2f deg; %d x %d pixels\n', ...
    scene.center.lat, scene.center.lon, scene.heading, numel(scene.cross_vec), numel(scene.along_vec));

%% 4. BACKPROJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = c/radar.fc;
bp_opts.mode = 'stripmap';
bp_opts.bw_az = lambda/antenna.D;   % Processed beamwidth (rad) -> ~D/2 resolution
bp_opts.bw_el = [];                 % Elevation gate (rad), [] = off
bp_opts.D = antenna.D;
bp_opts.ant_weight = true;          % Azimuth matched filter (antenna pattern)
bp_opts.upsample = 8;
bp_opts.stop_and_hop = false;       % Account for motion during echo flight
fprintf('Expected resolution: %.2f m slant range, %.2f m azimuth\n', ...
    c/(2*radar.B), lambda/(2*bp_opts.bw_az));

% Radar/nav clock offset estimated from the data by maximizing focus
% (estimateNavTimeOffset).  Only reliable when the timing error dominates
% the other motion errors (see its header); otherwise set
% nav_opts.time_offset from a calibration and let autofocus remove the
% remaining phase error.
timing_flag = 0;
t_opts.range = 0.05;          % Search +/- this around nav_opts.time_offset (s)
if timing_flag
    [nav_opts.time_offset, t_est] = estimateNavTimeOffset(nav, radar, ref, ...
        antenna, nav_opts, grid, bp_opts, t_opts);
    geom = navToPulses(nav, radar.t_pulse, ref, antenna, nav_opts);
end
img = bpNav(radar, geom, grid, bp_opts);

%% 4b. AUTOFOCUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrects residual motion errors (INS drift, lever-arm/timing errors)
% directly on the time-domain backprojection output.  See autofocusBP.
af_flag = 1;
af_opts.block = [64 16 4];    % Pulses per phase estimate, coarse to fine;
                              % last stage shorter than the fastest residual
                              % motion period
af_opts.n_pixels = 5000;      % Bright pixels used for the sharpness metric
af_opts.n_iter = 20;          % Max coordinate-descent passes per stage

img_noaf = img;
if af_flag
    [img, af] = autofocusBP(radar, geom, grid, bp_opts, img_noaf, af_opts);
end

%% 5. DISPLAY / SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dyn_range = 50;
out_file = '';                % e.g. 'stripmap_image.mat' to save

im_dB = 20*log10(abs(img.im) + eps);  im_dB = im_dB - max(im_dB(:));
figure
imagesc(grid.along_vec, grid.cross_vec, im_dB);
colormap(green14); axis image xy; caxis([-dyn_range 0]); colorbar
xlabel('Along track (m)'); ylabel('Cross track, + away from radar (m)');
title('Stripmap backprojection, autofocused (dB)');

figure
imagesc(grid.along_vec, grid.cross_vec, img.hits); axis image xy; colorbar
xlabel('Along track (m)'); ylabel('Cross track (m)');
title('Pulses integrated per pixel');

if af_flag
    nf_dB = 20*log10(abs(img_noaf.im) + eps);  nf_dB = nf_dB - max(nf_dB(:));
    figure
    imagesc(grid.along_vec, grid.cross_vec, nf_dB);
    colormap(green14); axis image xy; caxis([-dyn_range 0]); colorbar
    xlabel('Along track (m)'); ylabel('Cross track (m)');
    title('Before autofocus (dB)');

    figure
    tp = geom.t - geom.t(1);
    plot(tp, af.phi_pulse, 'b'); hold on
    leg = {'autofocus estimate'};
    if use_simulated && isfield(truth, 'nav_err_pulse')
        % True correction: 4*pi/lambda * (line of sight . nav error).
        % Compare modulo 2*pi after removing the best-fit constant + linear
        % phase, which only shift the image and cannot be estimated.
        los = repmat(grid.center_enu, 1, numel(tp)) - geom.apc;
        los = los ./ repmat(sqrt(sum(los.^2,1)), 3, 1);
        ph_true = 4*pi/lambda * sum(los .* truth.nav_err_pulse, 1);
        d = exp(1i*(af.phi_pulse - ph_true));
        M = 16*numel(d); [~, m] = max(abs(fft(d, M)));
        d = d .* exp(-2i*pi*(m-1)*(0:numel(d)-1)/M);
        resid = angle(d * exp(-1i*angle(sum(d))));
        plot(tp, af.phi_pulse - resid, 'r--');
        leg{end+1} = 'truth (aligned)';
        fprintf('Autofocus phase vs. truth: %.2f rad RMS residual; nav-induced phase error was %.2f rad RMS\n', ...
            sqrt(mean(resid.^2)), std(ph_true));
    end
    xlabel('Time (s)'); ylabel('Phase correction (rad)'); legend(leg);
    title('Autofocus phase correction');
end

if use_simulated
    if af_flag
        fprintf('\nBefore autofocus:'); reportPointTargets(img_noaf, grid, truth, 12);
        fprintf('\nAfter autofocus:'); reportPointTargets(img, grid, truth, 12);
    else
        reportPointTargets(img, grid, truth);
    end
end
if ~isempty(out_file)
    image = img.im; lat = grid.lat; lon = grid.lon; h = grid.h; %#ok<NASGU>
    save(out_file, 'image', 'lat', 'lon', 'h', 'scene', 'bp_opts');
    fprintf('Saved %s\n', out_file);
end
