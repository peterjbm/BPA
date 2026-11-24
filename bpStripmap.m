function data = bpStripmap(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stripmap backprojection.  Same inputs and conventions as bpBasic,    %
% but each pulse is only backprojected onto the pixels that lie inside %
% the antenna's azimuth beam at that pulse.  In spotlight mode (which  %
% bpBasic assumes) every pulse illuminates the whole scene; in         %
% stripmap mode the beam sweeps along the ground, so each pixel has    %
% its own synthetic aperture.  Restricting to the beam:                %
%   - avoids smearing returns from pulses that never saw the pixel,    %
%   - sets the azimuth resolution to about lambda/(2*beamwidth),       %
%   - and skips most pixel/pulse pairs over a long strip.              %
%                                                                      %
% Required fields (as for bpBasic):                                    %
% data.Nfft, data.deltaF, data.minF, data.x_mat, data.y_mat,           %
% data.z_mat, data.AntX, data.AntY, data.AntZ, data.R0, data.phdata    %
%                                                                      %
% Stripmap fields:                                                     %
% data.procBeamwidth: Full azimuth beamwidth used to gate pixels (rad).%
%                     Defaults to data.azBeamwidth, else lambda/D.     %
% data.squint:        Beam squint from broadside, + forward (rad) [0]  %
% data.D:             Antenna length (m); needed for ant_weight_flag   %
%                     and the default beamwidth                        %
% data.ant_weight_flag: 1 = weight each pulse by the two-way antenna   %
%                     pattern (azimuth matched filter) [0]             %
% data.taper_flag:    1 = Hanning window across frequency [0]          %
%                                                                      %
% The along-track direction at each pulse is taken from the platform   %
% track (AntX/Y/Z), so the flight path need not lie along an axis.     %
%                                                                      %
% Outputs:                                                             %
% data.im_final:   Complex image                                       %
% data.hit_count:  Number of pulses that contributed to each pixel     %
%                                                                      %
% Built on bpBasic by LeRoy Gorham, AFRL (Gorham & Moore, SPIE 2010).  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 299792458;

data.K  = size(data.phdata,1);
data.Np = size(data.phdata,2);

if ~isfield(data,'squint'),          data.squint = 0;          end
if ~isfield(data,'ant_weight_flag'), data.ant_weight_flag = 0; end
if ~isfield(data,'taper_flag'),      data.taper_flag = 0;      end

lambda = c / (mean(data.minF) + data.deltaF*data.K/2);
if ~isfield(data,'procBeamwidth')
    if isfield(data,'azBeamwidth')
        data.procBeamwidth = data.azBeamwidth;
    elseif isfield(data,'D')
        data.procBeamwidth = lambda / data.D;
    else
        error('bpStripmap: set data.procBeamwidth (or data.D)');
    end
end
if data.ant_weight_flag && ~isfield(data,'D')
    error('bpStripmap: ant_weight_flag needs data.D');
end
halfBW = data.procBeamwidth / 2;

% Scene size and resolution
data.maxWr = c/(2*data.deltaF);
data.dr = c/(2*data.deltaF*data.K);
data.dx = lambda/(2*data.procBeamwidth);
fprintf('Maximum range swath:  %.2f m\n', data.maxWr);
fprintf('Resolution:  %.2f m range, %.2f m azimuth (beamwidth %.2f deg)\n', ...
    data.dr, data.dx, data.procBeamwidth*180/pi);

data.r_vec = linspace(-data.Nfft/2,data.Nfft/2-1,data.Nfft)*data.maxWr/data.Nfft;

% Optional frequency taper
if data.taper_flag
    win = 0.5 - 0.5*cos(2*pi*(1:data.K)'/(data.K+1));
else
    win = ones(data.K,1);
end

% Unit along-track (velocity) direction for every pulse
pos = [data.AntX(:) data.AntY(:) data.AntZ(:)];
vel = zeros(size(pos));
if data.Np > 1
    vel(2:end-1,:) = pos(3:end,:) - pos(1:end-2,:);
    vel(1,:)   = pos(2,:)   - pos(1,:);
    vel(end,:) = pos(end,:) - pos(end-1,:);
else
    vel(1,:) = [1 0 0];
end
vhat = vel ./ repmat(sqrt(sum(vel.^2,2)), 1, 3);

data.im_final  = zeros(size(data.x_mat));
data.hit_count = zeros(size(data.x_mat));

t = zeros(1,data.Np);
for ii = 1:data.Np

    if ii > 1 && mod(ii,100)==0
        t_sofar = sum(t(1:(ii-1)));
        t_est = (t_sofar*data.Np/(ii-1)-t_sofar)/60;
        fprintf('Pulse %d of %d, %.02f minutes remaining\n',ii,data.Np,t_est);
    end
    tic

    % Range profile for this pulse
    rc = fftshift(ifft(win .* data.phdata(:,ii), data.Nfft));

    % Pixel geometry relative to the antenna
    dx = data.x_mat - data.AntX(ii);
    dy = data.y_mat - data.AntY(ii);
    dz = data.z_mat - data.AntZ(ii);
    R  = sqrt(dx.^2 + dy.^2 + dz.^2);
    dR = R - data.R0(ii);

    % Angle of each pixel off broadside, then off the beam axis
    sinTheta = (dx*vhat(ii,1) + dy*vhat(ii,2) + dz*vhat(ii,3)) ./ R;
    offBeam  = asin(sinTheta) - data.squint;

    % Pixels inside both the beam and the range window
    I = find(abs(offBeam) <= halfBW & dR > min(data.r_vec) & dR < max(data.r_vec));
    if isempty(I), t(ii) = toc; continue; end

    phCorr = exp(1i*4*pi*data.minF(ii)/c*dR(I));
    val = interp1(data.r_vec, rc, dR(I), 'linear') .* phCorr;

    if data.ant_weight_flag
        u = data.D/lambda * sin(offBeam(I));
        w = ones(size(u)); nz = u ~= 0;
        w(nz) = (sin(pi*u(nz)) ./ (pi*u(nz))).^2;
        val = val .* w;
    end

    data.im_final(I)  = data.im_final(I) + val;
    data.hit_count(I) = data.hit_count(I) + 1;

    t(ii) = toc;
end

return
