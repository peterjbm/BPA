function [lat, lon, h] = nav_ecef2geodetic(P)
% ECEF (3xN, m) -> WGS-84 geodetic (deg, deg, m), iterative (mm accuracy)
a = 6378137; f = 1/298.257223563; e2 = f*(2-f);
x = P(1,:); y = P(2,:); z = P(3,:);
lon = atan2(y, x);
p = sqrt(x.^2 + y.^2);
lat = atan2(z, p*(1-e2));
for it = 1:6
    N = a ./ sqrt(1 - e2*sin(lat).^2);
    h = p./cos(lat) - N;
    lat = atan2(z, p.*(1 - e2*N./(N+h)));
end
N = a ./ sqrt(1 - e2*sin(lat).^2);
h = p./cos(lat) - N;
lat = lat*180/pi; lon = lon*180/pi;
end
