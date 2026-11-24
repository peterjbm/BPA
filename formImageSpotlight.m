clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spotlight SAR demo: simulates point targets with                     %
% generateDataSpotlight, then forms the image with bpBasic (AFRL       %
% spotlight backprojection).  For real GOTCHA data see formImage3DSAR. %
% Set compare_flag = 1 to also form a matched-filter (mfBasic) image   %
% on a small patch around the scene center as an exact reference.     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RADAR PARAMETERS %
data.fc = 10e9;          % Center frequency (Hz)
data.B  = 600e6;         % Bandwidth (Hz) -> 0.25 m range resolution
data.K  = 512;           % Frequency samples per pulse (sets range extent)
data.D  = 0.3;           % Antenna length (m); Inf = isotropic

% COLLECTION GEOMETRY %
data.path   = 'circular'; % 'circular' (GOTCHA-like) or 'linear'
data.Rc     = 5000;      % Slant range to scene center (m)
data.elev   = 30;        % Elevation angle (deg)
data.az0    = 90;        % Azimuth of aperture center (deg)
data.azSpan = 3.5;       % Total look-angle change (deg) -> ~0.25 m cross-range
data.Np     = 512;       % Number of pulses (sets cross-range extent)
data.SNR_dB = 0;         % Per-sample SNR of a unit target (Inf = no noise)

% POINT TARGETS: [x y z amplitude], scene center at the origin %
data.targets = [   0    0 0 1.0;
                0.75    0 0 0.7;   % pair 0.75 m apart in x (cross-range)
                   0 0.75 0 0.7;   % pair 0.75 m apart in y (range)
                 -15   10 0 1.0;
                  15  -10 0 1.0;
                  20   20 0 1.0;
                 -20  -20 0 1.0 ];

% IMAGE GRID (ground plane, z = 0) %
data.x_vec = linspace(-25, 25, 401);
data.y_vec = linspace(-25, 25, 401);

% PROCESSING %
data.Nfft = 8*data.K;    % Zero-padded FFT size for range interpolation
taper_flag = 1;          % 2-D Hanning taper for sidelobe control
dyn_range  = 50;         % Displayed dynamic range (dB)
compare_flag = 0;        % 1 = also run mfBasic on a small patch

% INPUT PARAMETERS END HERE %

data = generateDataSpotlight(data);
fprintf('%s path, %d pulses, %d x %d pixels\n', data.path, data.Np, ...
    numel(data.y_vec), numel(data.x_vec));

% Taper across frequency (rows) and pulses (columns)
if taper_flag
    hw = @(n) 0.5 - 0.5*cos(2*pi*(1:n)'/(n+1));
    data.phdata = data.phdata .* (hw(data.K) * hw(data.Np)');
end

t0 = tic; data = bpBasic(data); tBP = toc(t0);
fprintf('bpBasic: %.1f s\n', tBP);

im_dB = 20*log10(abs(data.im_final) + eps);
im_dB = im_dB - max(im_dB(:));

figure
imagesc(data.x_vec, data.y_vec, im_dB);
colormap(green14); axis image xy; caxis([-dyn_range 0]); colorbar
xlabel('x (m)'); ylabel('y (m)');
title('Spotlight backprojection (dB)');

if compare_flag
    % Matched filter is exact but costs K times more; use a small patch
    mf = data;
    mf.x_vec = linspace(-2, 2, 33);
    mf.y_vec = linspace(-2, 2, 33);
    [mf.x_mat, mf.y_mat] = meshgrid(mf.x_vec, mf.y_vec);
    mf.z_mat = zeros(size(mf.x_mat));
    t0 = tic; mf = mfBasic(mf); fprintf('mfBasic: %.1f s\n', toc(t0));
    mf_dB = 20*log10(abs(mf.im_final) + eps);
    figure
    imagesc(mf.x_vec, mf.y_vec, mf_dB - max(mf_dB(:)));
    colormap(green14); axis image xy; caxis([-dyn_range 0]); colorbar
    title('Matched filter, scene center patch (dB)');
end
