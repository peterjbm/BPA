function data = generateDataSpotlight(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates stepped-frequency spotlight SAR phase history for a set of %
% point targets, in the AFRL/GOTCHA data format, so the output can be  %
% fed straight to bpBasic or mfBasic.                                  %
%                                                                      %
% In spotlight mode the antenna is steered to keep the scene center    %
% (the origin) in the beam for every pulse, and each pulse is          %
% motion-compensated to the scene center: R0 is the range from the     %
% antenna to the origin.  Cross-range resolution comes from the total  %
% change in look angle, c/(2*totalAz*fmin).                            %
%                                                                      %
% Required fields:                                                     %
% data.fc:       Center frequency (Hz)                                 %
% data.B:        Bandwidth (Hz)                                        %
% data.K:        Number of frequency samples per pulse                 %
% data.Np:       Number of pulses                                      %
% data.Rc:       Slant range from aperture center to scene center (m)  %
% data.elev:     Elevation (grazing) angle at aperture center (deg)    %
% data.az0:      Azimuth of aperture center (deg, from +x toward +y)   %
% data.azSpan:   Total change in look angle over the aperture (deg)    %
% data.targets:  Nt x 4 matrix [x y z amplitude] of point targets      %
% data.x_vec, data.y_vec:  Image grid vectors (m), z = 0               %
%                                                                      %
% Optional fields:                                                     %
% data.path:     'circular' (constant range and elevation, GOTCHA-     %
%                like) or 'linear' (straight flight, beam steered)     %
%                ['circular']                                          %
% data.D:        Antenna length (m) for the two-way beam pattern;      %
%                Inf = isotropic [Inf]                                 %
% data.SNR_dB:   Per-sample SNR of a unit target; Inf = no noise [Inf] %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 299792458;

if ~isfield(data,'path'),   data.path = 'circular'; end
if ~isfield(data,'D'),      data.D = Inf;           end
if ~isfield(data,'SNR_dB'), data.SNR_dB = Inf;      end

% Frequency samples
data.freq   = data.fc - data.B/2 + (0:data.K-1)' * data.B/data.K;
data.deltaF = data.freq(2) - data.freq(1);
data.lambda = c / data.fc;
data.minF   = min(data.freq) * ones(1, data.Np);

el  = data.elev * pi/180;
az0 = data.az0 * pi/180;
span = data.azSpan * pi/180;

switch lower(data.path)
    case 'circular'
        % Constant slant range and elevation, azimuth sweeps the span
        az = az0 + linspace(-span/2, span/2, data.Np);
        data.AntX = data.Rc*cos(el)*cos(az);
        data.AntY = data.Rc*cos(el)*sin(az);
        data.AntZ = data.Rc*sin(el) * ones(1, data.Np);
    case 'linear'
        % Straight, level line perpendicular to the aperture-center look
        % direction, long enough to subtend azSpan at the scene center
        Rg = data.Rc*cos(el);                  % ground range at center
        L  = 2*Rg*tan(span/2);                 % aperture length
        s  = linspace(-L/2, L/2, data.Np);
        data.AntX = Rg*cos(az0) - s*sin(az0);
        data.AntY = Rg*sin(az0) + s*cos(az0);
        data.AntZ = data.Rc*sin(el) * ones(1, data.Np);
    otherwise
        error('generateDataSpotlight: unknown path ''%s''', data.path);
end

% Motion compensation to scene center
data.R0 = sqrt(data.AntX.^2 + data.AntY.^2 + data.AntZ.^2);

% Check the scene against the alias-free extents bpBasic will report
maxWr = c/(2*data.deltaF);
dAz = abs(mean(diff(unwrap(atan2(data.AntY, data.AntX)))));
maxWx = c/(2*dAz*min(data.freq));
ext = max(max(data.x_vec)-min(data.x_vec), max(data.y_vec)-min(data.y_vec)) * sqrt(2);
if ext > min(maxWr, maxWx)
    warning(['Image diagonal (%.1f m) exceeds the alias-free scene size ' ...
        '(%.1f m range, %.1f m cross-range): increase K or Np'], ext, maxWr, maxWx);
end

% Simulate returns
data.phdata = zeros(data.K, data.Np);
for t = 1:size(data.targets,1)
    tx = data.targets(t,1); ty = data.targets(t,2);
    tz = data.targets(t,3); ta = data.targets(t,4);

    R = sqrt((tx-data.AntX).^2 + (ty-data.AntY).^2 + (tz-data.AntZ).^2);

    % Two-way pattern of an antenna steered at the scene center
    pat = ones(1, data.Np);
    if isfinite(data.D)
        cosA = ((tx-data.AntX).*(-data.AntX) + (ty-data.AntY).*(-data.AntY) + ...
            (tz-data.AntZ).*(-data.AntZ)) ./ (R .* data.R0);
        u = data.D/data.lambda * sin(acos(min(cosA, 1)));
        nz = u ~= 0;
        pat(nz) = (sin(pi*u(nz)) ./ (pi*u(nz))).^2;
    end

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
