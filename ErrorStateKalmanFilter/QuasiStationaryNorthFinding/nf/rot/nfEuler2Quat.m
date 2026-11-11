function q = nfEuler2Quat(roll, pitch, yaw)
%NFEULER2QUAT Quaternion q_b^n of ZYX Euler angles [rad]: q_z(yaw) (x) q_y(pitch) (x) q_x(roll), scalar first.
%   nfQuat2Dcm(nfEuler2Quat(roll, pitch, yaw)) = nfEuler2Dcm(roll, pitch, yaw).

cr = cos(roll/2);  sr = sin(roll/2);
cp = cos(pitch/2); sp = sin(pitch/2);
cy = cos(yaw/2);   sy = sin(yaw/2);
q = [cr*cp*cy + sr*sp*sy;
     sr*cp*cy - cr*sp*sy;
     cr*sp*cy + sr*cp*sy;
     cr*cp*sy - sr*sp*cy];
end
