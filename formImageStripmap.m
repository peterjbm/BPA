clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stripmap SAR demo: simulates point targets along a strip with        %
% generateDataStripmap, then forms the image with bpStripmap.          %
% Set compare_flag = 1 to also run full-aperture bpBasic (spotlight    %
% assumption) for comparison.                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RADAR PARAMETERS %
data.fc = 10e9;          % Center frequency (Hz)
data.B  = 300e6;         % Bandwidth (Hz) -> 0.5 m range resolution
data.K  = 256;           % Frequency samples per pulse (sets swath depth)
data.D  = 1.0;           % Antenna length (m) -> ~0.5 m azimuth resolution
data.squint = 0;         % Beam squint from broadside (rad), + forward

% GEOMETRY %
data.H  = 1000;          % Platform altitude (m)
data.Yc = 2000;          % Ground range to swath center (m)
data.x_start = -140;     % Flight start, along track (m)
data.x_stop  =  140;     % Flight end, along track (m)
data.dxPulse = 0.25;     % Pulse spacing (m); keep below D/2
data.SNR_dB  = 10;       % Per-sample SNR of a unit target (Inf = no noise)

% POINT TARGETS: [x y z amplitude] %
data.targets = [ -80  1990 0 1.0;
                 -40  2010 0 1.0;
                   0  2000 0 1.0;
                   1  2000 0 0.5;   % close pair to show azimuth resolution
                  40  1985 0 1.0;
                  80  2015 0 1.0 ];

% IMAGE GRID (ground plane, z = 0) %
data.x_vec = -100:0.25:100;
data.y_vec = 1975:0.25:2025;

% PROCESSING %
data.Nfft = 4*data.K;        % Zero-padded FFT size for range interpolation
data.ant_weight_flag = 1;    % Weight pulses by the antenna pattern
data.taper_flag = 1;         % Hanning window across frequency
data.dyn_range = 50;         % Displayed dynamic range (dB)
compare_flag = 0;            % 1 = also form a full-aperture bpBasic image

% INPUT PARAMETERS END HERE %

data = generateDataStripmap(data);
fprintf('%d pulses, %d x %d pixels\n', data.Np, numel(data.y_vec), numel(data.x_vec));

t0 = tic; data = bpStripmap(data); tStrip = toc(t0);
fprintf('bpStripmap: %.1f s\n', tStrip);

im_dB = 20*log10(abs(data.im_final) + eps);
im_dB = im_dB - max(im_dB(:));

figure
imagesc(data.x_vec, data.y_vec, im_dB);
colormap(green14); axis image xy; caxis([-data.dyn_range 0]); colorbar
xlabel('Along track x (m)'); ylabel('Ground range y (m)');
title('Stripmap backprojection (dB)');

figure
imagesc(data.x_vec, data.y_vec, data.hit_count);
axis image xy; colorbar
xlabel('Along track x (m)'); ylabel('Ground range y (m)');
title('Pulses integrated per pixel');

if compare_flag
    dataBP = data;
    t0 = tic; dataBP = bpBasic(dataBP); tBasic = toc(t0);
    fprintf('bpBasic (full aperture): %.1f s\n', tBasic);
    imB = 20*log10(abs(dataBP.im_final) + eps);
    figure
    imagesc(data.x_vec, data.y_vec, imB - max(imB(:)));
    colormap(green14); axis image xy; caxis([-data.dyn_range 0]); colorbar
    title('bpBasic, full aperture (dB)');
end
