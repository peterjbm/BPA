function [img_af, af, radar_af] = autofocusBP(radar, geom, grid, bp_opts, img0, af_opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autofocus for time-domain backprojection (BPA) images.               %
%                                                                      %
% Residual motion errors the navigation does not capture (INS drift,   %
% lever-arm or timing errors, a few cm) put an unknown phase on each   %
% pulse.  A backprojected image is a sum of per-pulse images,          %
%     I(x) = sum_n exp(j*phi_n) * b_n(x),                              %
% so the phases can be found directly in the image domain - no polar   %
% reformatting or FFT to a phase-history domain, which is why this     %
% works for stripmap (space-variant apertures) as well as spotlight.   %
%                                                                      %
% Method (after Ash, "An autofocus method for backprojection imagery   %
% in synthetic aperture radar", IEEE GRSL 9(1), 2012): maximize image  %
% sharpness S = sum |I(x)|^4 over a set of bright pixels by coordinate %
% descent.  For one block of pulses with contribution b and the rest   %
% of the image A = I - exp(j*phi)*b,                                   %
%     S(phi) = const + Re(P exp(j phi)) + Re(Q exp(j 2 phi)),          %
%     P = 4 sum (|A|^2+|b|^2) conj(A) b,   Q = 2 sum (conj(A) b)^2,    %
% which is maximized exactly per block (grid search + Newton).         %
%                                                                      %
% Pulses are grouped into blocks sharing one phase: fewer unknowns,    %
% less noise, less memory.  Constant and linear phase across the       %
% aperture do not change sharpness but shift the image, so they are    %
% removed to keep the navigation's geolocation.                        %
%                                                                      %
% Inputs: radar, geom, grid, bp_opts as for bpNav; img0 = bpNav output %
% before autofocus ([] to compute it).                                 %
% af_opts:                                                             %
%   .block     Pulses per phase block.  A vector runs coarse-to-fine   %
%              stages (e.g. [64 16 4]): long blocks first catch large, %
%              slow errors quickly; short blocks then resolve faster   %
%              motion.  Pixels are re-selected from the improved image %
%              between stages.  The last block should be shorter than %
%              the fastest residual motion period.        [[64 16 4]]  %
%   .n_pixels  Bright pixels used for the metric              [5000]   %
%   .n_tiles   Along-track tiles the pixels are spread over            %
%              [16 for stripmap, 1 for spotlight]                      %
%   .n_iter    Max coordinate-descent passes                    [10]   %
%   .tol       Stop when sharpness gain per pass is below this [1e-4]  %
%   .detrend   Remove constant + linear phase                 [true]   %
%   .verbose                                                  [true]   %
%                                                                      %
% Outputs:                                                             %
%   img_af     bpNav output for the phase-corrected pulses             %
%   af.phi_pulse   Phase correction applied to every pulse (rad)       %
%   af.sharpness (per pass), af.pixels (last selection)               %
%   radar_af   radar with corrected pulses                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6, af_opts = struct(); end
strip = strcmpi(bp_opts.mode, 'stripmap');
if ~isfield(af_opts,'block'),    af_opts.block = [64 16 4]; end
if ~isfield(af_opts,'n_pixels'), af_opts.n_pixels = 5000; end
if ~isfield(af_opts,'n_tiles'),  af_opts.n_tiles = 1 + 15*strip; end
if ~isfield(af_opts,'n_iter'),   af_opts.n_iter = 20;     end
if ~isfield(af_opts,'tol'),      af_opts.tol = 1e-4;      end
if ~isfield(af_opts,'detrend'),  af_opts.detrend = true;  end
if ~isfield(af_opts,'verbose'),  af_opts.verbose = true;  end

if isempty(img0), img0 = bpNav(radar, geom, grid, bp_opts); end
Np = size(radar.rc, 2);
phi_pulse = zeros(1, Np);
radar_af = radar;
img_cur = img0;
af.sharpness = [];

for outer = 1:numel(af_opts.block)
    bsz = af_opts.block(outer);
    blk = ceil((1:Np)/bsz);
    Nb = blk(end);
    % ---- Select bright pixels, spread along track ----------------------
    sel = selectPixels(abs(img_cur.im).^2, img_cur.hits, af_opts.n_pixels, af_opts.n_tiles);

    % ---- Per-block contributions at those pixels (sub-grid backprojection)
    [~, col] = ind2sub(size(grid.X), sel);
    sg = grid;
    sg.X = grid.X(sel)'; sg.Y = grid.Y(sel)'; sg.Z = grid.Z(sel)';
    sg.along_vec = grid.along_vec(col);
    sg.dpx = abs(mean(diff(grid.along_vec)));
    sg.corners = [grid.X([1 end 1 end]); grid.Y([1 end 1 end]); grid.Z([1 end 1 end])];
    o = bp_opts; o.verbose = false;
    o.store_pixels = 1:numel(sel); o.block = bsz;
    r = bpNav(radar_af, geom, sg, o);
    B = r.contrib;                                    % Nb x K
    B = B / max(abs(r.im(:)));                        % numeric scaling

    % ---- Coordinate descent on block phases ----------------------------
    phi = zeros(Nb, 1);
    I = sum(B, 1);
    energy = sum(abs(B).^2, 2);
    active = find(energy > 1e-8*max(energy));
    S = sum(abs(I).^4) / sum(abs(I).^2)^2;
    S0 = S;
    grid_phi = (0:255)/256*2*pi;
    for it = 1:af_opts.n_iter
        for k = active'
            b = B(k,:);
            A = I - exp(1i*phi(k))*b;
            z = conj(A).*b;
            P = 4*sum((abs(A).^2 + abs(b).^2).*z);
            Q = 2*sum(z.^2);
            f = real(P*exp(1i*grid_phi)) + real(Q*exp(2i*grid_phi));
            [~, m] = max(f); ph = grid_phi(m);
            for nt = 1:3                               % Newton refinement
                d1 = -imag(P*exp(1i*ph)) - 2*imag(Q*exp(2i*ph));
                d2 = -real(P*exp(1i*ph)) - 4*real(Q*exp(2i*ph));
                if d2 < 0, ph = ph - d1/d2; end
            end
            phi(k) = ph;
            I = A + exp(1i*ph)*b;
        end
        Snew = sum(abs(I).^4) / sum(abs(I).^2)^2;
        af.sharpness(end+1) = Snew;
        if af_opts.verbose
            fprintf('Autofocus stage %d (block %d), pass %d: normalized sharpness %.4g (stage start %.4g)\n', ...
                outer, bsz, it, Snew, S0);
        end
        if (Snew - S)/S < af_opts.tol, S = Snew; break; end
        S = Snew;
    end

    % ---- Fill inactive blocks, remove constant + linear phase ----------
    phi = unwrap(phi);
    if numel(active) < Nb && numel(active) > 1
        phi = interp1(active, phi(active), (1:Nb)', 'linear', 'extrap');
        phi(1:active(1)-1) = phi(active(1));
        phi(active(end)+1:end) = phi(active(end));
    end
    if af_opts.detrend && numel(active) > 1
        w = energy(active);
        Xd = [ones(numel(active),1) active(:)];
        cf = (Xd' * (Xd .* repmat(w,1,2))) \ (Xd' * (w .* phi(active)));
        phi = phi - (cf(1) + cf(2)*(1:Nb)');
    end

    % ---- Apply to the pulses -------------------------------------------
    dphi = phi(blk)';
    radar_af.rc = radar_af.rc .* repmat(exp(1i*dphi), size(radar.rc,1), 1);
    phi_pulse = phi_pulse + dphi;

    if af_opts.verbose, fprintf('Autofocus stage %d: re-forming image\n', outer); end
    img_cur = bpNav(radar_af, geom, grid, bp_opts);
end

img_af = img_cur;
af.phi_pulse = phi_pulse;
af.pixels = sel;
af.block = af_opts.block;

end


function sel = selectPixels(pw, hits, K, ntiles)
% Brightest pixels, the same number from each along-track tile
[~, Nx] = size(pw);
pw(hits == 0) = 0;
edges = round(linspace(0, Nx, ntiles+1));
per = ceil(K/ntiles);
sel = [];
for t = 1:ntiles
    cols = edges(t)+1:edges(t+1);
    tile = pw(:, cols);
    [~, ord] = sort(tile(:), 'descend');
    ord = ord(1:min(per, nnz(tile)));
    [r, c] = ind2sub(size(tile), ord);
    sel = [sel; sub2ind(size(pw), r, cols(c)')]; %#ok<AGROW>
end
end
