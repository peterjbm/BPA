function E = nav_geodetic2enu(lat, lon, h, ref)
% Geodetic (deg, deg, m) -> local ENU (3xN, m) about ref.lat/lon/alt
P0 = nav_geodetic2ecef(ref.lat, ref.lon, ref.alt);
E = nav_ecef2enu_rot(ref.lat, ref.lon) * (nav_geodetic2ecef(lat, lon, h) - P0);
end
