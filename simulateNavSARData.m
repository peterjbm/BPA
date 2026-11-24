function [nav, radar, antenna, truth] = simulateNavSARData(p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates an airborne stripmap or spotlight collection in the same   %
% form as real data: INS/GPS records (lat, lon, alt, roll, pitch, yaw, %
% NED velocity) at the nav rate, and range-compressed pulses at the    %
% PRF.  The aircraft weaves and its attitude oscillates; the antenna   %
% sits on a lever arm from the INS.  Used to test formImageStripmap    %
% and formImageSpotlight before loading real data.                     %
%                                                                      %
% p fields (see the two formImage scripts for typical values):         %
%  .mode 'stripmap'|'spotlight'   .ref (scene center lat/lon/alt)      %
%  .fc .B .fs .PRF .Nr            .speed .heading(deg) .H .Yc          %
%  .nav_rate (Hz)                 .targets [along cross dh amp]        %
%  .motion: .cross_amp .cross_T .up_amp .up_T  (m, s)                  %
%  .att:    .roll_amp .roll_T .pitch_amp .pitch_T .yaw_amp .yaw_T      %
%           .crab  (deg, s)                                            %
%  .antenna: .lever_arm .squint .D  (look is right)                    %
%  .azSpan (deg, spotlight aperture)   .along_extent [s0 s1] (optional)%
%  .SNR_dB  .nav_pos_noise (m)  .nav_att_noise (deg)                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = 299792458;
d2r = pi/180;
lambda = c/p.fc;
ref = p.ref;
hd = p.heading*d2r;
ua = [sin(hd); cos(hd); 0];       % along track (ENU)
uc = [cos(hd); -sin(hd); 0];      % right of track
uz = [0;0;1];
Rg = p.Yc; Rc = sqrt(p.Yc^2 + p.H^2);
if ~isfield(p,'SNR_dB'),        p.SNR_dB = Inf;      end
if ~isfield(p,'nav_pos_noise'), p.nav_pos_noise = 0; end
if ~isfield(p,'nav_att_noise'), p.nav_att_noise = 0; end

antenna = p.antenna;
antenna.look = 'right';
antenna.depression = atan2(p.H, p.Yc)/d2r;    % aim at swath center

% ---- Pulse along-track extent --------------------------------------------
if isfield(p,'along_extent') && ~isempty(p.along_extent)
    ext = p.along_extent;
elseif strcmpi(p.mode, 'spotlight')
    L = 2*Rc*tan(p.azSpan/2*d2r);       % subtends azSpan at the scene (slant)
    ext = [-L/2 L/2];
else
    foot = Rc*lambda/antenna.D;
    shift = Rc*abs(sin((p.att.crab + antenna.squint)*d2r));
    ext = [min(p.targets(:,1)) - foot/2 - shift - 10, ...
           max(p.targets(:,1)) + foot/2 + shift + 10];
end
t0 = 1.0;                                    % first pulse, s on nav clock
T  = diff(ext)/p.speed;
t_pulse = t0 + (0:1/p.PRF:T);
Np = numel(t_pulse);

% ---- Truth motion model (functions of time) ------------------------------
m = p.motion; a = p.att;
pos = @(t) ua*(ext(1) + p.speed*(t - t0)) + ...
    uc*(-p.Yc + m.cross_amp*sin(2*pi*t/m.cross_T)) + ...
    uz*(p.H + m.up_amp*sin(2*pi*t/m.up_T + 1));
vel = @(t) ua*p.speed*ones(size(t)) + ...
    uc*(m.cross_amp*2*pi/m.cross_T*cos(2*pi*t/m.cross_T)) + ...
    uz*(m.up_amp*2*pi/m.up_T*cos(2*pi*t/m.up_T + 1));
att = @(t) [a.roll_amp*sin(2*pi*t/a.roll_T); ...
            a.pitch_amp*sin(2*pi*t/a.pitch_T + 0.5); ...
            p.heading + a.crab + a.yaw_amp*sin(2*pi*t/a.yaw_T + 2)] * d2r;

Renu = nav_ecef2enu_rot(ref.lat, ref.lon);
dep = antenna.depression*d2r; sq = antenna.squint*d2r;
b_body = [cos(dep)*sin(sq); cos(dep)*cos(sq); sin(dep)];
a_body = [1;0;0] - b_body(1)*b_body; a_body = a_body/norm(a_body);

    function [apc, bore, azax] = antenna_at(t)
        P = pos(t); A = att(t);
        [la, lo] = nav_enu2geodetic(P, ref);
        apc = zeros(3,numel(t)); bore = apc; azax = apc;
        for i = 1:numel(t)
            C = Renu*nav_ned2ecef_rot(la(i), lo(i))*nav_euler2dcm(A(1,i), A(2,i), A(3,i));
            apc(:,i)  = P(:,i) + C*antenna.lever_arm(:);
            bore(:,i) = C*b_body; azax(:,i) = C*a_body;
        end
    end

% ---- Targets on the ellipsoid --------------------------------------------
Nt = size(p.targets,1);
Ttan = ua*p.targets(:,1)' + uc*p.targets(:,2)';
[tla, tlo] = nav_enu2geodetic(Ttan, ref);
th = ref.alt + p.targets(:,3)';
Tenu = nav_geodetic2enu(tla, tlo, th, ref);

% ---- Range-compressed pulses ---------------------------------------------
tau0 = 2*(Rc - p.Nr/2*c/(2*p.fs))/c;
tau_k = tau0 + (0:p.Nr-1)'/p.fs;
[apc, bore, azax] = antenna_at(t_pulse);
rc = zeros(p.Nr, Np);
for it = 1:Nt
    d = repmat(Tenu(:,it),1,Np) - apc;
    Rtx = sqrt(sum(d.^2,1));
    apc_rx = antenna_at(t_pulse + 2*Rtx/c);          % platform moves during echo
    Rrx = sqrt(sum((repmat(Tenu(:,it),1,Np) - apc_rx).^2,1));
    tau = (Rtx + Rrx)/c;
    u = d ./ repmat(Rtx,3,1);
    if strcmpi(p.mode, 'spotlight')
        % antenna steered at scene center (ENU origin)
        uc0 = -apc ./ repmat(sqrt(sum(apc.^2,1)),3,1);
        off = acos(min(sum(u.*uc0,1), 1));
    else
        off = atan2(sum(u.*azax,1), sum(u.*bore,1));
    end
    pat = nav_sinc(antenna.D/lambda*sin(off)).^2;
    rc = rc + p.targets(it,4) * nav_sinc(p.B*(repmat(tau_k,1,Np) - repmat(tau,p.Nr,1))) .* ...
        repmat(pat .* exp(-1i*2*pi*p.fc*tau), p.Nr, 1);
end
if isfinite(p.SNR_dB)
    sigma = 10^(-p.SNR_dB/20);
    rc = rc + sigma/sqrt(2)*(randn(size(rc)) + 1i*randn(size(rc)));
end

radar.rc = rc; radar.fs = p.fs; radar.tau0 = tau0; radar.fc = p.fc;
radar.B = p.B; radar.t_pulse = t_pulse; radar.phase_sign = 1;

% ---- Navigation records (INS point, not the APC) -------------------------
tn = (t_pulse(1) - 0.5):(1/p.nav_rate):(t_pulse(end) + 0.5);
Pn = pos(tn); Vn = vel(tn); An = att(tn);
Pn = Pn + p.nav_pos_noise*randn(size(Pn));
[nav.lat, nav.lon, nav.alt] = nav_enu2geodetic(Pn, ref);
nav.vN = zeros(size(tn)); nav.vE = nav.vN; nav.vD = nav.vN;
for k = 1:numel(tn)
    vned = nav_ned2ecef_rot(nav.lat(k), nav.lon(k))' * Renu' * Vn(:,k);
    nav.vN(k) = vned(1); nav.vE(k) = vned(2); nav.vD(k) = vned(3);
end
An = An/d2r + p.nav_att_noise*randn(size(An));
nav.t = tn; nav.roll = An(1,:); nav.pitch = An(2,:); nav.yaw = mod(An(3,:), 360);

truth.targets_enu = Tenu; truth.targets_lat = tla; truth.targets_lon = tlo;
truth.targets_h = th; truth.apc = apc; truth.scene_center = ref;

fprintf(['Simulated %s: %d pulses at %.0f Hz (PRT %.2f ms), nav at %.0f Hz ' ...
    '(interval %.2f ms), %d targets\n'], p.mode, Np, p.PRF, 1e3/p.PRF, ...
    p.nav_rate, 1e3/p.nav_rate, Nt);

end
