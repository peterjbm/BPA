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
nav_opts.interp = 'hermite';  % 'hermite' (uses velocity) or 'linear'
nav_opts.time_offset = 0;     % Radar-to-nav clock offset (s)

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

img = bpNav(radar, geom, grid, bp_opts);

%% 5. DISPLAY / SAVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dyn_range = 50;
out_file = '';                % e.g. 'stripmap_image.mat' to save

im_dB = 20*log10(abs(img.im) + eps);  im_dB = im_dB - max(im_dB(:));
figure
imagesc(grid.along_vec, grid.cross_vec, im_dB);
colormap(green14); axis image xy; caxis([-dyn_range 0]); colorbar
xlabel('Along track (m)'); ylabel('Cross track, + away from radar (m)');
title('Stripmap backprojection (dB)');

figure
imagesc(grid.along_vec, grid.cross_vec, img.hits); axis image xy; colorbar
xlabel('Along track (m)'); ylabel('Cross track (m)');
title('Pulses integrated per pixel');

if use_simulated
    reportPointTargets(img, grid, truth);
end
if ~isempty(out_file)
    image = img.im; lat = grid.lat; lon = grid.lon; h = grid.h; %#ok<NASGU>
    save(out_file, 'image', 'lat', 'lon', 'h', 'scene', 'bp_opts');
    fprintf('Saved %s\n', out_file);
end
