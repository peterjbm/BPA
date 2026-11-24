function C = nav_ned2ecef_rot(lat, lon)
% Rotation matrix taking local North-East-Down vectors at (lat, lon) deg to ECEF
p = lat*pi/180; l = lon*pi/180;
C = [ -sin(p)*cos(l)  -sin(l)  -cos(p)*cos(l);
      -sin(p)*sin(l)   cos(l)  -cos(p)*sin(l);
       cos(p)          0       -sin(p) ];
end
