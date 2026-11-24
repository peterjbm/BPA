function [lat, lon, h] = nav_enu2geodetic(E, ref)
% Local ENU (3xN, m) about ref.lat/lon/alt -> geodetic (deg, deg, m)
P0 = nav_geodetic2ecef(ref.lat, ref.lon, ref.alt);
P = nav_ecef2enu_rot(ref.lat, ref.lon)' * E + P0;
[lat, lon, h] = nav_ecef2geodetic(P);
end
