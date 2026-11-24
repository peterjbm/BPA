function diag = diagnoseResidualPhase(af, geom, grid, radar, opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explains the autofocus phase correction in terms of navigation       %
% errors, so they can be fixed at the source instead of only being     %
% absorbed by autofocus.                                               %
%                                                                      %
% The autofocus correction phi(t) (af.phi_pulse) is compared with the  %
% phase signature each candidate error would produce:                  %
%   time offset dt:  4*pi/lambda * (u . v(t)) * dt                     %
%                    (u = line of sight to the scene, v = velocity)    %
%   lever arm dL:    4*pi/lambda * u . C(t) dL   (C = body -> ENU)     %
% after removing a constant and linear term (which only shift the      %
% image).  Each candidate is tested on its own (value, % explained),   %
% and the phase is analysed by frequency: ripple at the nav            %
% rate (interpolation / quantization), ~1 Hz steps (GPS updates in the %
% INS solution), broadband (velocity noise, pulse time-tag jitter).    %
%                                                                      %
% Over a short aperture with smooth motion, different causes produce   %
% nearly the same phase and one pass cannot separate them; candidates  %
% that look alike are flagged.  Treat results as hypotheses: apply the %
% best one, re-form without autofocus, and confirm across passes.      %
%                                                                      %
% Inputs: af from autofocusBP, geom from navToPulses, grid, radar.     %
% opts.terms  any of {'time','lever'}         [both]                   %
% opts.verbose                                [true]                   %
%                                                                      %
% diag.time_offset_correction (s): ADD to nav_opts.time_offset if the  %
%   timing hypothesis holds; diag.lever_arm_correction (m, body x/y/z) %
%   likewise, each estimated as the only error                         %
% diag.explained: fraction of the phase each hypothesis explains       %
% diag.signature_corr, diag.freq, diag.psd, diag.peaks                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5, opts = struct(); end
if ~isfield(opts,'terms'),   opts.terms = {'time','lever'}; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

c = 299792458; k = 4*pi*radar.fc/c;
phi = af.phi_pulse(:);
Np  = numel(phi);
t   = geom.t(:) - mean(geom.t);
T   = max(abs(t));

u = repmat(grid.center_enu, 1, Np) - geom.apc;
u = u ./ repmat(sqrt(sum(u.^2,1)), 3, 1);

% ---- Physical signatures --------------------------------------------------
S = []; names = {};
if any(strcmp(opts.terms,'time'))
    S = [S, k*sum(u .* geom.vel, 1)'];
    names{end+1} = 'time offset';
end
if any(strcmp(opts.terms,'lever'))
    if ~isfield(geom,'C'), error('diagnoseResidualPhase: geom.C missing (update navToPulses)'); end
    L = zeros(Np,3);
    for n = 1:Np, L(n,:) = k * (u(:,n)' * geom.C(:,:,n)); end
    S = [S, L];
    names = [names, {'lever x','lever y','lever z'}];
end

% Remove constant + linear from everything (they only shift the image)
N = [ones(Np,1) t];
P = @(x) x - N*(N\x);
phd = P(phi);
Sd  = zeros(size(S));
for j = 1:size(S,2), Sd(:,j) = P(S(:,j)); end

% ---- Single-cause hypotheses -------------------------------------------
% For each candidate error alone: the value that best explains the phase
% and the fraction of the phase it explains.  Several candidates can fit
% about equally well when their signatures look alike over this
% aperture - then the data cannot tell them apart.
nP = numel(names);
val = zeros(nP,1); expl = zeros(nP,1);
for j = 1:nP
    nj = Sd(:,j)'*Sd(:,j);
    if nj > 0
        val(j) = (Sd(:,j)'*phd)/nj;
        expl(j) = 1 - sum((phd - Sd(:,j)*val(j)).^2)/sum(phd.^2);
    end
end
Cn = Sd ./ repmat(sqrt(sum(Sd.^2,1)) + eps, Np, 1);
R  = abs(Cn'*Cn);

res = phd;      % spectrum of the full phase (constant + linear removed)
diag.rms.total = std(phd);
diag.names = names; diag.value = -val; diag.explained = expl; diag.signature_corr = R;
it = find(strcmp(names,'time offset'));
il = find(strncmp(names,'lever',5));
diag.time_offset_correction = []; diag.lever_arm_correction = [];
if ~isempty(it), diag.time_offset_correction = -val(it); end
if ~isempty(il), diag.lever_arm_correction = -val(il); end

pr = @(varargin) opts.verbose && fprintf(varargin{:});
pr('\nResidual phase diagnosis: autofocus correction %.2f rad RMS (constant + linear removed)\n', diag.rms.total);
pr('  If the error were ONLY ...           correction to apply       explains\n');
[~, ord] = sort(expl, 'descend');
for j = ord(:)'
    if strcmp(names{j}, 'time offset')
        txt = sprintf('nav_opts.time_offset %+8.2f ms', -1e3*val(j));
    else
        txt = sprintf('antenna.lever_arm(%s) %+7.3f m', names{j}(end), -val(j));
    end
    alike = setdiff(find(R(j,:) > 0.9), j);
    note = '';
    if ~isempty(alike), note = ['  (looks like: ' strjoin(names(alike), ', ') ')']; end
    pr('    %-12s                  %-28s %4.0f%%%s\n', names{j}, txt, 100*max(expl(j),0), note);
end
pr(['  A cause is likely if it explains most of the phase with a physically sensible value.\n' ...
    '  Confirm by applying it and re-forming WITHOUT autofocus: the image should sharpen, and the\n' ...
    '  value should repeat across passes.  Fix timing/lever arm from a calibration where possible.\n']);

% ---- What is left, by frequency -------------------------------------------
dtp = median(diff(geom.t));
if isfield(af,'block') && ~isempty(af.block)
    fmax = 1/(2*dtp*af.block(end));       % autofocus sees up to its block rate
else
    fmax = 1/(2*dtp);
end
w = 0.5 - 0.5*cos(2*pi*(0:Np-1)'/(Np-1));
Nf = 2^nextpow2(4*Np);
X = fft(res.*w, Nf);
f = (0:Nf/2)'/(Nf*dtp);
psd = abs(X(1:Nf/2+1)).^2;
keep = f <= fmax;
f = f(keep); psd = psd(keep);
bands = [0 0.5; 0.5 5; 5 fmax];
tot = sum(psd);
pr('  phase by frequency:');
for b = 1:size(bands,1)
    in = f >= bands(b,1) & f < bands(b,2);
    pr('  %g-%g Hz %3.0f%%', bands(b,1), round(bands(b,2)*10)/10, 100*sum(psd(in))/max(tot,eps));
end
pr('\n');
% strongest peaks above 0.5 Hz
hi = find(f > 0.5);
[~, ord] = sort(psd(hi), 'descend');
pk = [];
for i = ord(:)'
    fi = f(hi(i));
    if all(abs(pk - fi) > 0.5), pk(end+1) = fi; end %#ok<AGROW>
    if numel(pk) == 3, break; end
end
if ~isempty(pk)
    pr('  strongest components above 0.5 Hz: %s Hz\n', sprintf('%.1f ', pk));
    pr(['  (peaks at the nav rate or its multiples -> nav interpolation/quantization; ~1 Hz -> GPS ' ...
        'updates in the INS solution;\n   broadband -> velocity noise or pulse time-tag jitter)\n']);
end

diag.residual = res'; diag.freq = f; diag.psd = psd; diag.peaks = pk;

end
