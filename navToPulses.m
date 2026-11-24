function geom = navToPulses(nav, t_pulse, ref, antenna, opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Converts INS/GPS navigation records to antenna phase center (APC)    %
% position, velocity and antenna pointing at every radar pulse, in a   %
% local East-North-Up (ENU) frame centered on ref.                     %
%                                                                      %
% The navigation and radar clocks can run at any rates; each pulse is  %
% interpolated from the nav records around it.  Position uses cubic    %
% Hermite interpolation that matches both the nav position AND the nav %
% velocity at each record, so a slow nav rate (nav interval longer     %
% than the PRT) still gives smooth, accurate pulse positions.          %
%                                                                      %
% nav (all vectors of length Nn, same clock as t_pulse):               %
%   .t           Time (s)                                              %
%   .lat, .lon   WGS-84 latitude/longitude of the INS (deg)            %
%   .alt         Height above the WGS-84 ellipsoid (m)                 %
%   .roll, .pitch, .yaw   Z-Y-X Euler angles (deg); yaw from true      %
%                north; body x forward, y right, z down                %
%   .vN, .vE, .vD  Velocity in local North-East-Down (m/s)             %
%                                                                      %
% t_pulse:       Transmit time of every pulse (s), 1 x Np              %
% ref:           ENU origin, .lat .lon (deg) .alt (m)                  %
% antenna:                                                             %
%   .lever_arm   INS -> APC offset in body axes [x;y;z] (m) [0;0;0]    %
%   .look        'right' or 'left'                         ['right']   %
%   .depression  Boresight angle below body x-y plane (deg)            %
%   .squint      Boresight rotation toward body +x (deg)       [0]     %
% opts:                                                                %
%   .interp      'hermite' (uses velocity) or 'linear'   ['hermite']   %
%   .time_offset Added to t_pulse to put it on the nav clock (s) [0]   %
%   .verbose     Print diagnostics                          [true]     %
%                                                                      %
% geom (per pulse, ENU):                                               %
%   .apc  3xNp APC position (m)      .vel  3xNp velocity (m/s)         %
%   .bore 3xNp antenna boresight     .az_axis 3xNp antenna azimuth     %
%   .el_axis 3xNp antenna elevation axis (all unit vectors)            %
%   .ins  3xNp INS position          .t    pulse times on nav clock    %
%   .roll/.pitch/.yaw  interpolated attitude (deg)                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5, opts = struct(); end
if ~isfield(opts,'interp'),      opts.interp = 'hermite'; end
if ~isfield(opts,'time_offset'), opts.time_offset = 0;    end
if ~isfield(opts,'verbose'),     opts.verbose = true;     end
if ~isfield(antenna,'lever_arm'), antenna.lever_arm = [0;0;0]; end
if ~isfield(antenna,'look'),      antenna.look = 'right';      end
if ~isfield(antenna,'squint'),    antenna.squint = 0;          end

tn = nav.t(:)';
Nn = numel(tn);
if any(diff(tn) <= 0), error('navToPulses: nav.t must be strictly increasing'); end

% ---- Nav positions and velocities in ENU --------------------------------
Pn = nav_geodetic2enu(nav.lat, nav.lon, nav.alt, ref);          % 3 x Nn
Renu = nav_ecef2enu_rot(ref.lat, ref.lon);
Vn = zeros(3, Nn);
for k = 1:Nn
    Vn(:,k) = Renu * nav_ned2ecef_rot(nav.lat(k), nav.lon(k)) * ...
        [nav.vN(k); nav.vE(k); nav.vD(k)];
end

% Consistency check: reported velocity vs. differenced positions.  A large
% mismatch usually means a wrong velocity frame/sign or a time-tag error.
if Nn > 2
    dt = diff(tn);
    vfd = diff(Pn,1,2) ./ repmat(dt, 3, 1);
    vmid = (Vn(:,1:end-1) + Vn(:,2:end)) / 2;
    vres = sqrt(mean(sum((vfd - vmid).^2, 1)));
    if opts.verbose
        fprintf('Nav: %d records, interval %.4f s (%.1f Hz), velocity vs. position RMS mismatch %.3f m/s\n', ...
            Nn, median(dt), 1/median(dt), vres);
    end
    if vres > 0.05*sqrt(mean(sum(vmid.^2,1))) + 0.2
        warning(['navToPulses: nav velocity disagrees with differenced positions ' ...
            '(%.2f m/s RMS). Check the NED convention and time tags.'], vres);
    end
end

% ---- Interpolate to pulse times -----------------------------------------
t = t_pulse(:)' + opts.time_offset;
Np = numel(t);
if t(1) < tn(1) || t(end) > tn(end)
    error('navToPulses: pulses span %.4f..%.4f s but nav covers only %.4f..%.4f s', ...
        t(1), t(end), tn(1), tn(end));
end

switch lower(opts.interp)
    case 'hermite'
        fi = interp1(tn, 1:Nn, t);
        k  = min(floor(fi), Nn-1);
        h  = tn(k+1) - tn(k);
        s  = (t - tn(k)) ./ h;
        h00 = 2*s.^3 - 3*s.^2 + 1;   h10 = s.^3 - 2*s.^2 + s;
        h01 = -2*s.^3 + 3*s.^2;      h11 = s.^3 - s.^2;
        d00 = (6*s.^2 - 6*s)./h;     d10 = 3*s.^2 - 4*s + 1;
        d01 = (-6*s.^2 + 6*s)./h;    d11 = 3*s.^2 - 2*s;
        Pins = zeros(3,Np); Vins = zeros(3,Np);
        for ax = 1:3
            p0 = Pn(ax,k); p1 = Pn(ax,k+1); v0 = Vn(ax,k); v1 = Vn(ax,k+1);
            Pins(ax,:) = h00.*p0 + h10.*h.*v0 + h01.*p1 + h11.*h.*v1;
            Vins(ax,:) = d00.*p0 + d10.*v0 + d01.*p1 + d11.*v1;
        end
    case 'linear'
        Pins = interp1(tn, Pn', t)';
        Vins = interp1(tn, Vn', t)';
    otherwise
        error('navToPulses: unknown interp ''%s''', opts.interp);
end

% Attitude: interpolate unwrapped angles
d2r = pi/180;
roll  = interp1(tn, unwrap(nav.roll(:)*d2r),  t);
pitch = interp1(tn, unwrap(nav.pitch(:)*d2r), t);
yaw   = interp1(tn, unwrap(nav.yaw(:)*d2r),   t);

% Local NED orientation at each pulse
[plat, plon] = nav_enu2geodetic(Pins, ref);

% ---- Antenna geometry in body axes --------------------------------------
switch lower(antenna.look)
    case 'right', sgn = 1;
    case 'left',  sgn = -1;
    otherwise, error('navToPulses: antenna.look must be ''right'' or ''left''');
end
dep = antenna.depression*d2r; sq = antenna.squint*d2r;
b_body = [cos(dep)*sin(sq); sgn*cos(dep)*cos(sq); sin(dep)];   % boresight
a_body = [1;0;0] - b_body(1)*b_body;                            % azimuth axis
a_body = a_body / norm(a_body);
e_body = cross(b_body, a_body);                                 % elevation axis

geom.apc = zeros(3,Np); geom.bore = zeros(3,Np);
geom.az_axis = zeros(3,Np); geom.el_axis = zeros(3,Np);
lever = antenna.lever_arm(:);
for n = 1:Np
    C = Renu * nav_ned2ecef_rot(plat(n), plon(n)) * nav_euler2dcm(roll(n), pitch(n), yaw(n));
    geom.apc(:,n)     = Pins(:,n) + C*lever;
    geom.bore(:,n)    = C*b_body;
    geom.az_axis(:,n) = C*a_body;
    geom.el_axis(:,n) = C*e_body;
end

geom.ins = Pins;
geom.vel = Vins;
geom.t   = t;
geom.roll = roll/d2r; geom.pitch = pitch/d2r; geom.yaw = yaw/d2r;
geom.ref = ref;

if opts.verbose
    spd = sqrt(sum(Vins.^2,1));
    trk = atan2(mean(Vins(1,:)), mean(Vins(2,:)))/d2r;
    fprintf(['Pulses: %d, PRF %.1f Hz, mean speed %.1f m/s, track %.1f deg, ' ...
        'mean pulse spacing %.3f m\n'], Np, 1/median(diff(t)), mean(spd), mod(trk,360), ...
        mean(sqrt(sum(diff(geom.apc,1,2).^2,1))));
end

end
