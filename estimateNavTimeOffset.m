function [dt, est] = estimateNavTimeOffset(nav, radar, ref, antenna, nav_opts, grid, bp_opts, t_opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimates the radar-to-navigation clock offset from the radar data:  %
% the offset dt (s) such that pulse n was transmitted at nav time      %
% t_pulse(n) + dt.  Use it as nav_opts.time_offset.                    %
%                                                                      %
% With the wrong offset, every antenna position is off by roughly      %
% v(t)*dt.  The constant-velocity part of that only shifts the image   %
% along track; the changing part (weave, climb, turns) defocuses it.   %
% For each trial offset this function removes the shift (mean          %
% velocity * trial offset) so the image stays in place, backprojects   %
% onto a set of bright pixels, and measures sharpness                  %
% sum|I|^4 / (sum|I|^2)^2.  The sharpest image gives the offset.       %
%                                                                      %
% The offset is only observable through velocity changes during the    %
% aperture: on a perfectly steady straight line it cannot be found     %
% from focus (and does not defocus the image either - it only shifts   %
% it; use a ground control point for that).                            %
%                                                                      %
% RELIABILITY: focus cannot tell a timing error from other smooth      %
% motion errors (INS drift) over a short aperture - a wrong offset can %
% partly cancel them and look sharper.  Trust the result only when     %
% timing dominates: a strongly maneuvering pass, a long strip, or after%
% checking the estimate is stable across several data takes.  The     %
% offset is a property of the system, so calibrate it once and reuse   %
% it.  The most robust calibration is a surveyed corner reflector:     %
% along-track geolocation error / ground speed = time offset.          %
%                                                                      %
% Inputs: as navToPulses / bpNav; grid from makeGeoImageGrid.          %
% t_opts:                                                              %
%   .range     Search +/- this around nav_opts.time_offset (s) [0.05]  %
%   .step      Coarse search step (s); [] = automatic from how fast    %
%              the line-of-sight velocity changes            [[]]      %
%   .tol       Final precision (s)                         [1e-4]      %
%   .n_pixels  Bright pixels used for the metric           [3000]      %
%   .img0      Image at the starting offset ([] = form it)  [[]]       %
%   .verbose                                               [true]      %
%                                                                      %
% Outputs: dt  estimated offset (s)                                    %
%   est.dt_curve, est.sharpness  coarse scan                           %
%   est.shift  along-track image shift relative to the start (m)       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 8, t_opts = struct(); end
if ~isfield(t_opts,'range'),    t_opts.range = 0.05;    end
if ~isfield(t_opts,'step'),     t_opts.step = [];       end
if ~isfield(t_opts,'tol'),      t_opts.tol = 1e-4;      end
if ~isfield(t_opts,'n_pixels'), t_opts.n_pixels = 3000; end
if ~isfield(t_opts,'img0'),     t_opts.img0 = [];       end
if ~isfield(t_opts,'verbose'),  t_opts.verbose = true;  end
if ~isfield(nav_opts,'time_offset'), nav_opts.time_offset = 0; end

c = 299792458; lambda = c/radar.fc;
dt0 = nav_opts.time_offset;
no = nav_opts; no.verbose = false;
g0 = navToPulses(nav, radar.t_pulse, ref, antenna, no);
vbar = mean(g0.vel, 2);

% Automatic coarse step: change in offset that moves the antenna by
% ~lambda/8 along the line of sight relative to the steady-flight shift
if isempty(t_opts.step)
    los = repmat(grid.center_enu, 1, size(g0.apc,2)) - g0.apc;
    los = los ./ repmat(sqrt(sum(los.^2,1)), 3, 1);
    dv = max(abs(sum((g0.vel - repmat(vbar,1,size(g0.vel,2))) .* los, 1)));
    t_opts.step = min(max(lambda/(8*max(dv, eps)), t_opts.tol), t_opts.range/5);
end

% Bright pixels from the starting image
img0 = t_opts.img0;
if isempty(img0)
    o = bp_opts; o.verbose = false;
    img0 = bpNav(radar, g0, grid, o);
end
strip = strcmpi(bp_opts.mode, 'stripmap');
sel = selectBrightPixels(abs(img0.im).^2, img0.hits, t_opts.n_pixels, 1 + 15*strip);
sg = makeSubGrid(grid, sel);
o = bp_opts; o.verbose = false;

    function S = sharp(d)
        no2 = no; no2.time_offset = d;
        g = navToPulses(nav, radar.t_pulse, ref, antenna, no2);
        g.apc = g.apc - vbar*(d - dt0);          % keep the image in place
        r = bpNav(radar, g, sg, o);
        I2 = abs(r.im).^2;
        S = sum(I2.^2) / sum(I2)^2;
    end

% ---- Coarse scan ---------------------------------------------------------
cand = dt0 + (-t_opts.range:t_opts.step:t_opts.range);
if t_opts.verbose
    fprintf('Timing offset search: %d trial offsets, %.2f ms apart, +/- %.0f ms\n', ...
        numel(cand), 1e3*t_opts.step, 1e3*t_opts.range);
end
Sc = zeros(size(cand));
for i = 1:numel(cand)
    Sc(i) = sharp(cand(i));
end
[~, ib] = max(Sc);

% ---- Golden-section refinement ------------------------------------------
a = cand(max(ib-1,1)); b = cand(min(ib+1,numel(cand)));
gr = (sqrt(5)-1)/2;
x1 = b - gr*(b-a); x2 = a + gr*(b-a); f1 = sharp(x1); f2 = sharp(x2);
while (b - a) > t_opts.tol
    if f1 > f2
        b = x2; x2 = x1; f2 = f1; x1 = b - gr*(b-a); f1 = sharp(x1);
    else
        a = x1; x1 = x2; f1 = f2; x2 = a + gr*(b-a); f2 = sharp(x2);
    end
end
dt = (a + b)/2;

est.dt_curve = cand; est.sharpness = Sc;
est.shift = norm(vbar(1:2)) * (dt - dt0);
est.S_start = sharp(dt0); est.S_best = sharp(dt);
if t_opts.verbose
    fprintf(['Estimated radar->nav time offset %.2f ms (start %.2f ms): sharpness x%.2f, ' ...
        'image moves %.2f m along track\n'], 1e3*dt, 1e3*dt0, est.S_best/est.S_start, est.shift);
    if ib == 1 || ib == numel(cand)
        warning('estimateNavTimeOffset: best offset is at the edge of the search; widen t_opts.range');
    end
end

end
