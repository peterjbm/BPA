function P = nav_geodetic2ecef(lat, lon, h)
% WGS-84 geodetic (deg, deg, m above ellipsoid) -> ECEF (3xN, m)
a = 6378137; f = 1/298.257223563; e2 = f*(2-f);
lat = lat(:)'*pi/180; lon = lon(:)'*pi/180; h = h(:)';
N = a ./ sqrt(1 - e2*sin(lat).^2);
P = [(N+h).*cos(lat).*cos(lon); (N+h).*cos(lat).*sin(lon); (N*(1-e2)+h).*sin(lat)];
end
