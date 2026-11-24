function C = nav_euler2dcm(roll, pitch, yaw)
% Body -> NED direction cosine matrix from aerospace Z-Y-X Euler angles (rad).
% Body axes: x forward, y right, z down.  yaw from true north.
cr = cos(roll);  sr = sin(roll);
cp = cos(pitch); sp = sin(pitch);
cy = cos(yaw);   sy = sin(yaw);
C = [ cp*cy   sr*sp*cy - cr*sy   cr*sp*cy + sr*sy;
      cp*sy   sr*sp*sy + cr*cy   cr*sp*sy - sr*cy;
     -sp      sr*cp              cr*cp ];
end
