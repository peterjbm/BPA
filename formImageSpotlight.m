clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spotlight SAR image formation from range-compressed pulses and       %
% INS/GPS navigation (lat, lon, alt, roll, pitch, yaw, NED velocity).  %
%                                                                      %
%  1. navToPulses:      nav records -> antenna phase center position   %
%                       and velocity at every pulse (Hermite           %
%                       interpolation using the velocity; roll, pitch  %
%                       and yaw rotate the INS-to-antenna lever arm)   %
%  2. makeGeoImageGrid: ground grid around the spotlight aim point,    %
%                       on the WGS-84 ellipsoid or a DEM               %
%  3. bpNav:            backprojection of every pulse to every pixel   %
%                       (the antenna is steered at the scene)          %
%                                                                      %
% Set use_simulated = 0 and fill in section 1 to process real data.    %
% For the AFRL GOTCHA phase-history datasets use formImage3DSAR.       %
% Requires MATLAB R2016b+ or GNU Octave; no toolboxes.                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(fullfile(fileparts(mfilename('fullpath')), 'navutils'));
c = 299792458;

%% 1. DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_simulated = 1;

if use_simulated
    p.mode = 'spotlight';
    p.ref = struct('lat', 39.9990, 'lon', -83.0150, 'alt', 220);  % aim point
    p.fc = 9.6e9;  p.B = 600e6;  p.fs = 750e6;  p.Nr = 1024;
    p.PRF = 500;   p.speed = 120;  p.heading = 35;
    p.H = 3000;    p.Yc = 4000;
    p.azSpan = 3.5;                  % Aperture angle (deg) -> ~0.25 m cross-range
    p.nav_rate = 50;                 % nav interval 20 ms vs PRT 2 ms
    p.motion = struct('cross_amp', 2, 'cross_T', 6, 'up_amp', 1, 'up_T', 4);
    p.att = struct('roll_amp', 2, 'roll_T', 3, 'pitch_amp', 1, 'pitch_T', 5, ...
                   'yaw_amp', 0.5, 'yaw_T', 7, 'crab', 0.5);
    p.antenna = struct('lever_arm', [1.5; 0.3; 0.8], 'squint', 0, 'D', 1.0);
    p.targets = [   0     0 0 1.0;     % [cross-range down-range height amp], m
                 0.75     0 0 0.7;     % pair 0.75 m apart in cross-range
                    0  0.75 0 0.7;     % pair 0.75 m apart in range
                  -15    10 0 1.0;
                   15   -10 0 1.0;
                   20    20 0 1.0;
                  -20   -20 0 1.0 ];
    p.SNR_dB = 0;
    % Residual INS position error the processor does not know about
    % (what autofocus must remove); set amp to [0 0 0] for perfect nav
    p.nav_err = struct('amp', [0.03 0.03 0.02], 'T', [0.7 1.6 3.3]);
    [nav, radar, antenna, truth] = simulateNavSARData(p);
    aim = p.ref;
else
    % ---- LOAD YOUR DATA HERE --------------------------------------------
    % nav, radar and antenna as described in formImageStripmap.m, plus
    % the spotlight aim point (scene center):
    %   aim = struct('lat', .., 'lon', .., 'alt', ..)   % alt above ellipsoid
    load('my_spotlight_data.mat');    % <- replace with your loader
end

%% 2. NAVIGATION -> PULSE GEOMETRY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nav_opts.interp = 'hermite';  % 'hermite' (uses velocity) or 'linear'
nav_opts.time_offset = 0;     % Radar-to-nav clock offset (s)

ref = aim;                    % Local ENU origin at the aim point
geom = navToPulses(nav, radar.t_pulse, ref, antenna, nav_opts);

%% 3. IMAGE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scene.heading = [];           % deg from north; [] = mean flight track, so
                              % "along" = cross-range, "cross" = range
scene.along_vec = -30:0.125:30;  % Cross-range pixel positions (m)
scene.cross_vec = -30:0.125:30;  % Down-range pixel positions (m), + = away
scene.height = aim.alt;       % scalar, Ny x Nx matrix, or @(lat,lon) DEM

if isempty(scene.heading)
    v = mean(geom.vel, 2);
    scene.heading = atan2(v(1), v(2))*180/pi;
end
grid = makeGeoImageGrid(ref, aim, scene.heading, ...
    scene.along_vec, scene.cross_vec, scene.height);

% Aperture diagnostics: look-angle span and alias-free scene size
los = geom.apc ./ repmat(sqrt(sum(geom.apc.^2,1)), 3, 1);
span = acos(max(min(los(:,1)'*los(:,end), 1), -1));
dth = span/(size(geom.apc,2)-1);
lambda = c/radar.fc;
fprintf(['Aperture %.2f deg: expected resolution %.2f m range, %.2f m cross-range; ' ...
    'alias-free cross-range extent %.0f m\n'], span*180/pi, c/(2*radar.B), ...
    lambda/(2*span), lambda/(2*dth));

%% 4. BACKPROJECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bp_opts.mode = 'spotlight';
bp_opts.upsample = 8;
bp_opts.stop_and_hop = false;
taper_flag = 1;               % Hanning taper across pulses (sidelobe control)

if taper_flag
    Np = size(radar.rc, 2);
    radar.rc = radar.rc .* repmat(0.5 - 0.5*cos(2*pi*(1:Np)/(Np+1)), size(radar.rc,1), 1);
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
out_file = '';                % e.g. 'spotlight_image.mat' to save

im_dB = 20*log10(abs(img.im) + eps);  im_dB = im_dB - max(im_dB(:));
figure
imagesc(grid.along_vec, grid.cross_vec, im_dB);
colormap(green14); axis image xy; caxis([-dyn_range 0]); colorbar
xlabel('Cross-range (m)'); ylabel('Down-range, + away from radar (m)');
title('Spotlight backprojection, autofocused (dB)');

if af_flag
    nf_dB = 20*log10(abs(img_noaf.im) + eps);  nf_dB = nf_dB - max(nf_dB(:));
    figure
    imagesc(grid.along_vec, grid.cross_vec, nf_dB);
    colormap(green14); axis image xy; caxis([-dyn_range 0]); colorbar
    xlabel('Cross-range (m)'); ylabel('Down-range (m)');
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
        fprintf('\nBefore autofocus:'); reportPointTargets(img_noaf, grid, truth, 5);
        fprintf('\nAfter autofocus:'); reportPointTargets(img, grid, truth, 5);
    else
        reportPointTargets(img, grid, truth);
    end
end
if ~isempty(out_file)
    image = img.im; lat = grid.lat; lon = grid.lon; h = grid.h; %#ok<NASGU>
    save(out_file, 'image', 'lat', 'lon', 'h', 'scene', 'bp_opts');
    fprintf('Saved %s\n', out_file);
end
