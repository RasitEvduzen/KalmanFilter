function C = nfEuler2Dcm(roll, pitch, yaw)
%NFEULER2DCM Body-to-navigation DCM from ZYX Euler angles [rad]: C_b^n = Rz(yaw)*Ry(pitch)*Rx(roll).
%   C = nfEuler2Dcm(roll, pitch, yaw) maps a vector in body axes to navigation axes,
%   v_n = C*v_b (Groves, 2nd ed., eq. 2.22 transposed). yaw is the heading of body x.
%   Inverse: nfDcm2Euler.

cr = cos(roll);  sr = sin(roll);
cp = cos(pitch); sp = sin(pitch);
cy = cos(yaw);   sy = sin(yaw);

C = [cp*cy, -cr*sy + sr*sp*cy,  sr*sy + cr*sp*cy;
     cp*sy,  cr*cy + sr*sp*sy, -sr*cy + cr*sp*sy;
       -sp,             sr*cp,             cr*cp];
end
