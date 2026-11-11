function [att, v, r, C] = nfMechanize(att, v, r, f_b, w_b, Ts, w_ie, g_D, df)
%NFMECHANIZE One quasi-static strapdown step in NED: attitude, velocity and displacement.
%   [att, v, r, C] = nfMechanize(att, v, r, f_b, w_b, Ts, w_ie, g_D, df)
%   att is the body-to-NED attitude, either the quaternion q_b^n (4x1, scalar first) or the DCM
%   C_b^n (3x3); C returns the DCM of the propagated attitude. v, r velocity [m/s] and
%   displacement from the start [m] (NED); f_b [m/s^2] and w_b [rad/s] bias-corrected body
%   samples; Ts step [s]; w_ie Earth rate in NED [rad/s]; g_D gravity, down [m/s^2]; df NED
%   specific-force correction (disturbance estimate) [m/s^2].
%       quaternion:  q+ = q(-w_ie Ts) (x) q (x) q(w_b Ts), normalised (nfRotVec2Quat, nfQuatMultiply)
%       DCM:         C+ = orth(C (I + [w_b x] Ts) - [w_ie x] C Ts)
%       v+ = v + (0.5 (C+ + C) f_b - df) Ts + [0; 0; g_D] Ts,    r+ = r + Ts/2 (v + v+)
%   Coriolis and transport rate are neglected: the IMU stays at one place, and at its residual
%   speed (~1e-4 m/s) they are below 2e-8 m/s^2 and 1e-5 deg/h. The Earth rate itself is kept.

if numel(att) == 4
    Cm = nfQuat2Dcm(att);
    att = nfQuatMultiply(nfQuatMultiply(nfRotVec2Quat(-w_ie * Ts), att), nfRotVec2Quat(w_b * Ts));
    att = att / norm(att);
    C = nfQuat2Dcm(att);
else
    Cm = att;
    C = att*(eye(3) + nfSkew(w_b * Ts)) - (nfSkew(w_ie) * att) * Ts;
    C = nfOrthonormalize(C);
    att = C;
end
f_n = 0.5*(C + Cm)*f_b - df;
vn = v + f_n*Ts + [0; 0; g_D]*Ts;
r = r + (Ts/2)*(v + vn);
v = vn;
end
