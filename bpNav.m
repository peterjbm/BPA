function img = bpNav(radar, geom, grid, opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-domain backprojection of range-compressed pulses using measured %
% navigation (from navToPulses) onto a geographic grid (from           %
% makeGeoImageGrid).  Handles stripmap and spotlight collections.      %
%                                                                      %
% radar:                                                               %
%   .rc        Range-compressed complex pulses, Nr x Np (fast time in  %
%              rows), compressed so a point at two-way delay tau peaks %
%              at sample (tau - tau0)*fs + 1                           %
%   .fs        Range sample rate of .rc (Hz)                           %
%   .tau0      Two-way delay of the first sample (s): scalar, or 1xNp  %
%              if the receive window moves                             %
%   .fc        Carrier frequency (Hz)                                  %
%   .phase_sign  +1 if a target at delay tau has phase                 %
%              exp(-j*2*pi*fc*tau) (demodulation by exp(-j*2*pi*fc*t)),%
%              -1 for the opposite convention                 [+1]     %
%                                                                      %
% opts:                                                                %
%   .mode      'stripmap' (gate each pulse to the antenna beam using   %
%              the measured attitude) or 'spotlight' (antenna steered  %
%              to the scene: every pulse to every pixel)               %
%   .bw_az     Full azimuth beamwidth used to gate pixels (rad)        %
%              [stripmap only]                                         %
%   .bw_el     Full elevation beamwidth (rad); [] = no elevation gate  %
%   .D         Antenna length (m); if given with .ant_weight, pulses   %
%              are weighted by the two-way pattern (matched filter)    %
%   .ant_weight  Weight by antenna pattern                  [false]    %
%   .upsample  FFT upsampling of each pulse before linear interp [8]   %
%   .stop_and_hop  true = ignore platform motion during the echo's     %
%              flight time; false = receive position = transmit        %
%              position + velocity * delay                  [false]    %
%   .verbose   Progress messages                            [true]     %
%   .store_pixels  Linear pixel indices at which to keep per-block     %
%              contributions (for autofocusBP)                  [[]]   %
%   .block     Pulses per block for .store_pixels               [1]    %
%                                                                      %
% img.im      Complex image (size of grid.X)                           %
% img.hits    Number of pulses integrated per pixel                    %
% img.contrib Nblocks x numel(store_pixels) contributions (if asked)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 299792458;

if ~isfield(radar,'phase_sign'),  radar.phase_sign = 1;   end
if ~isfield(opts,'mode'),         opts.mode = 'stripmap'; end
if ~isfield(opts,'bw_el'),        opts.bw_el = [];        end
if ~isfield(opts,'ant_weight'),   opts.ant_weight = false; end
if ~isfield(opts,'upsample'),     opts.upsample = 8;      end
if ~isfield(opts,'stop_and_hop'), opts.stop_and_hop = false; end
if ~isfield(opts,'verbose'),      opts.verbose = true;    end
if ~isfield(opts,'store_pixels'), opts.store_pixels = []; end
if ~isfield(opts,'block'),        opts.block = 1;         end

strip = strcmpi(opts.mode, 'stripmap');
if strip && ~isfield(opts,'bw_az')
    error('bpNav: stripmap mode needs opts.bw_az');
end
if opts.ant_weight && ~isfield(opts,'D')
    error('bpNav: ant_weight needs opts.D');
end

[Nr, Np] = size(radar.rc);
if size(geom.apc,2) ~= Np
    error('bpNav: %d pulses in radar.rc but %d in geom', Np, size(geom.apc,2));
end
tau0 = radar.tau0(:)'; if isscalar(tau0), tau0 = tau0*ones(1,Np); end
lambda = c/radar.fc;

% Upsample range profiles (band-limited FFT interpolation)
U = opts.upsample;
X = fft(radar.rc, [], 1);
NrU = Nr*U;
Y = zeros(NrU, Np);
nh = ceil(Nr/2); nl = Nr - nh;
Y(1:nh,:) = X(1:nh,:);
Y(end-nl+1:end,:) = X(end-nl+1:end,:);
rcu = ifft(Y, [], 1) * U;
fsU = radar.fs * U;

[Ny, Nx] = size(grid.X);
img.im   = zeros(Ny, Nx);
img.hits = zeros(Ny, Nx);

% Optional: keep each pulse block's time-domain contribution at selected
% pixels (used by autofocusBP).  contrib(k,j) = sum over pulses in block k
% of the backprojected value at pixel store_pixels(j).
store = ~isempty(opts.store_pixels);
if store
    slot = zeros(Ny, Nx);
    slot(opts.store_pixels) = 1:numel(opts.store_pixels);
    img.block_of_pulse = ceil((1:Np)/opts.block);
    img.contrib = zeros(img.block_of_pulse(end), numel(opts.store_pixels));
end

% Along-track pixel coordinate of each column (for the stripmap prefilter)
if strip
    col_along = grid.along_vec;
    if isfield(grid,'dpx'), dpx = grid.dpx; else dpx = abs(mean(diff(col_along))); end
    if isfield(grid,'corners')         % sub-grids carry the full grid's corners
        corners = grid.corners;
    else
        corners = [grid.X([1 end 1 end]); grid.Y([1 end 1 end]); grid.Z([1 end 1 end])];
    end
end

if opts.verbose
    fprintf('bpNav (%s): %d pulses, %d x %d pixels\n', opts.mode, Np, Ny, Nx);
end
tstart = tic;
for n = 1:Np
    if opts.verbose && mod(n, 500) == 0
        el = toc(tstart);
        fprintf('  pulse %d of %d, %.1f s elapsed, %.1f s remaining\n', n, Np, el, el*(Np-n)/n);
    end

    apc = geom.apc(:,n);

    % Stripmap: restrict to columns that can be inside the beam
    if strip
        b  = geom.bore(:,n);
        s_n = (apc - grid.center_enu)' * grid.u_along;
        ba = b' * grid.u_along;
        Rmax = sqrt(max(sum((corners - repmat(apc,1,4)).^2, 1)));
        lo = s_n + min(0, Rmax*ba) - Rmax*sin(0.75*opts.bw_az) - 2*dpx;
        hi = s_n + max(0, Rmax*ba) + Rmax*sin(0.75*opts.bw_az) + 2*dpx;
        cols = find(col_along >= lo & col_along <= hi);
        if isempty(cols), continue; end
    else
        cols = 1:Nx;
    end

    dx = grid.X(:,cols) - apc(1);
    dy = grid.Y(:,cols) - apc(2);
    dz = grid.Z(:,cols) - apc(3);
    Rtx = sqrt(dx.^2 + dy.^2 + dz.^2);
    if opts.stop_and_hop
        Rrx = Rtx;
    else
        v = geom.vel(:,n);
        dtv = 2*Rtx/c;
        Rrx = sqrt((dx - v(1)*dtv).^2 + (dy - v(2)*dtv).^2 + (dz - v(3)*dtv).^2);
    end
    tau = (Rtx + Rrx)/c;
    fi  = (tau - tau0(n))*fsU + 1;
    mask = fi >= 1 & fi < NrU;

    if strip || opts.ant_weight
        ux = dx./Rtx; uy = dy./Rtx; uz = dz./Rtx;
        b = geom.bore(:,n); a = geom.az_axis(:,n);
        cb = ux*b(1) + uy*b(2) + uz*b(3);
        th = atan2(ux*a(1) + uy*a(2) + uz*a(3), cb);
        if strip
            mask = mask & cb > 0 & abs(th) <= opts.bw_az/2;
            if ~isempty(opts.bw_el)
                e = geom.el_axis(:,n);
                the = atan2(ux*e(1) + uy*e(2) + uz*e(3), cb);
                mask = mask & abs(the) <= opts.bw_el/2;
            end
        end
    end

    I = find(mask(:));                      % column vectors throughout, so
    if isempty(I), continue; end            % single-row grids work too
    fI = fi(:); fI = fI(I);
    i0 = floor(fI); fr = fI - i0;
    prof = rcu(:,n);
    tI = tau(:);
    val = prof(i0).*(1-fr) + prof(i0+1).*fr;
    val = val .* exp(1i*radar.phase_sign*2*pi*radar.fc*tI(I));
    if opts.ant_weight
        thI = th(:);
        val = val .* nav_sinc(opts.D/lambda*sin(thI(I))).^2;
    end

    sub = img.im(:,cols);   sub(I) = sub(I) + reshape(val, size(sub(I)));   img.im(:,cols) = sub;
    hs  = img.hits(:,cols); hs(I) = hs(I) + 1;                             img.hits(:,cols) = hs;

    if store
        sl = slot(:,cols); sl = sl(:); sl = sl(I); hit = sl > 0;
        if any(hit)
            k = img.block_of_pulse(n);
            img.contrib(k, sl(hit)) = img.contrib(k, sl(hit)) + val(hit).';
        end
    end
end
if opts.verbose, fprintf('bpNav done in %.1f s\n', toc(tstart)); end

end
