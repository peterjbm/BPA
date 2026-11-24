function R = nav_ecef2enu_rot(lat, lon)
% Rotation matrix taking ECEF vectors to local East-North-Up at (lat, lon) deg
p = lat*pi/180; l = lon*pi/180;
R = [ -sin(l)          cos(l)         0;
      -sin(p)*cos(l)  -sin(p)*sin(l)  cos(p);
       cos(p)*cos(l)   cos(p)*sin(l)  sin(p) ];
end
