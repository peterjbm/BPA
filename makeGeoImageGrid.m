function grid = makeGeoImageGrid(ref, center, heading, along_vec, cross_vec, height)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Builds an image grid on the Earth's surface, aligned with a heading. %
%                                                                      %
% ref:        ENU origin, .lat .lon (deg) .alt (m)                     %
% center:     Grid center, .lat .lon (deg)                             %
% heading:    Direction of the grid's "along" axis, deg from north     %
%             (use the flight track); "cross" is 90 deg to its right   %
% along_vec:  Pixel coordinates along the heading (m), 1 x Nx          %
% cross_vec:  Pixel coordinates to the right of heading (m), 1 x Ny    %
% height:     Surface height above the WGS-84 ellipsoid (m): a scalar, %
%             an Ny x Nx matrix, or a function handle h = f(lat, lon)  %
%             (e.g. a DEM lookup)                                      %
%                                                                      %
% Pixels are laid out on the tangent plane at the center, then moved   %
% onto the ellipsoid/terrain height, so Earth curvature is included.   %
%                                                                      %
% grid.X, .Y, .Z   ENU pixel positions (m), Ny x Nx (rows = cross)     %
% grid.lat, .lon, .h   Geodetic pixel positions                        %
% grid.u_along, .u_cross, .center_enu, .along_vec, .cross_vec, .ref    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hd = heading*pi/180;
u_along = [sin(hd); cos(hd); 0];
u_cross = [cos(hd); -sin(hd); 0];

% Center on the tangent plane at the center's surface height
if isa(height, 'function_handle')
    hc = height(center.lat, center.lon);
elseif isscalar(height)
    hc = height;
else
    hc = height(ceil(end/2), ceil(end/2));
end
c_enu = nav_geodetic2enu(center.lat, center.lon, hc, ref);
Renu  = nav_ecef2enu_rot(ref.lat, ref.lon);
Rc    = nav_ecef2enu_rot(center.lat, center.lon);
% Center's local east/north/up, expressed in the reference ENU frame
Ec = Renu*Rc'*[1;0;0]; Nc = Renu*Rc'*[0;1;0];
ua = u_along(1)*Ec + u_along(2)*Nc;
uc = u_cross(1)*Ec + u_cross(2)*Nc;

[A, Cr] = meshgrid(along_vec(:)', cross_vec(:)');
Ptan = repmat(c_enu, 1, numel(A)) + ua*A(:)' + uc*Cr(:)';

% Move onto the surface
[lat, lon] = nav_enu2geodetic(Ptan, ref);
if isa(height, 'function_handle')
    h = height(lat, lon);
elseif isscalar(height)
    h = height*ones(size(lat));
else
    h = height(:)';
end
P = nav_geodetic2enu(lat, lon, h, ref);

sz = size(A);
grid.X = reshape(P(1,:), sz); grid.Y = reshape(P(2,:), sz); grid.Z = reshape(P(3,:), sz);
grid.lat = reshape(lat, sz); grid.lon = reshape(lon, sz); grid.h = reshape(h, sz);
grid.u_along = ua; grid.u_cross = uc; grid.center_enu = c_enu;
grid.along_vec = along_vec(:)'; grid.cross_vec = cross_vec(:)';
grid.ref = ref;

end
