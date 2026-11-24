function data = generateDataStripmap(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates stepped-frequency stripmap SAR phase history for a set of  %
% point targets, in the same format the AFRL loaders produce (so the   %
% output can be fed to bpStripmap, bpBasic or mfBasic).                %
%                                                                      %
% The platform flies a straight line along +x at height H, offset      %
% from the scene in y.  A real-aperture antenna of length D, squinted  %
% by data.squint, illuminates a strip of ground; each target is only   %
% seen while it lies inside the beam, weighted by the two-way pattern. %
%                                                                      %
% Required fields:                                                     %
% data.fc:        Center frequency (Hz)                                %
% data.B:         Bandwidth (Hz)                                       %
% data.K:         Number of frequency samples per pulse                %
% data.D:         Antenna length along track (m)                       %
% data.H:         Platform altitude (m)                                %
% data.Yc:        Ground-range (y) to swath center (m)                 %
% data.x_start, data.x_stop:  Along-track extent of the flight (m)     %
% data.dxPulse:   Along-track spacing between pulses (m)               %
% data.targets:   Nt x 4 matrix [x y z amplitude] of point targets     %
% data.x_vec, data.y_vec:  Image grid vectors (m), z = 0               %
%                                                                      %
% Optional fields:                                                     %
% data.squint:    Beam squint from broadside, + is forward (rad) [0]   %
% data.SNR_dB:    Per-sample SNR of a unit target; Inf = no noise [Inf]%
%                                                                      %
% Output fields: phdata, freq, minF, deltaF, AntX/Y/Z, R0, Np,         %
% x_mat, y_mat, z_mat, lambda, azBeamwidth                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 299792458;

if ~isfield(data,'squint'), data.squint = 0;   end
if ~isfield(data,'SNR_dB'), data.SNR_dB = Inf; end

% Frequency samples
data.freq   = data.fc - data.B/2 + (0:data.K-1)' * data.B/data.K;
data.deltaF = data.freq(2) - data.freq(1);
data.lambda = c / data.fc;

% Nominal azimuth beamwidth (rad), first-null half width = lambda/D
data.azBeamwidth = data.lambda / data.D;

% Platform track
data.AntX = data.x_start:data.dxPulse:data.x_stop;
data.Np   = numel(data.AntX);
data.AntY = zeros(1, data.Np);
data.AntZ = data.H * ones(1, data.Np);
data.minF = min(data.freq) * ones(1, data.Np);

% Reference (motion-compensation) range for each pulse: slant range to
% the swath-center line along the beam axis.  Phase history is
% "dechirped" relative to this range, as in the GOTCHA data.
data.R0 = sqrt((data.Yc - data.AntY).^2 + data.AntZ.^2) ./ cos(data.squint);

% Check that the swath fits in the unambiguous range window
maxWr = c / (2*data.deltaF);
dR_t = sqrt(data.targets(:,2).^2 + data.H^2) - data.R0(1);
if any(abs(dR_t) >= maxWr/2)
    warning('Some targets lie outside the unambiguous range window (+/- %.1f m)', maxWr/2);
end

% Simulate returns
data.phdata = zeros(data.K, data.Np);
for t = 1:size(data.targets,1)
    tx = data.targets(t,1); ty = data.targets(t,2);
    tz = data.targets(t,3); ta = data.targets(t,4);

    dx = tx - data.AntX; dy = ty - data.AntY; dz = tz - data.AntZ;
    R  = sqrt(dx.^2 + dy.^2 + dz.^2);

    % Angle off broadside (flight direction is +x)
    theta = asin(dx ./ R);

    % Two-way amplitude pattern of a uniform aperture: sinc^2
    u = data.D/data.lambda * sin(theta - data.squint);
    pat = ones(size(u));
    nz  = u ~= 0;
    pat(nz) = (sin(pi*u(nz)) ./ (pi*u(nz))).^2;

    data.phdata = data.phdata + ta * repmat(pat, data.K, 1) .* ...
        exp(-1i*4*pi/c * data.freq * (R - data.R0));
end

% Additive complex white Gaussian noise
if isfinite(data.SNR_dB)
    sigma = 10^(-data.SNR_dB/20);
    data.phdata = data.phdata + sigma/sqrt(2) * ...
        (randn(size(data.phdata)) + 1i*randn(size(data.phdata)));
end

% Image grid (z = 0 ground plane)
[data.x_mat, data.y_mat] = meshgrid(data.x_vec, data.y_vec);
data.z_mat = zeros(size(data.x_mat));

return
