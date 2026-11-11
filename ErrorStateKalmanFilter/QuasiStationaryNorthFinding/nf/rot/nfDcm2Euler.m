function [roll, pitch, yaw] = nfDcm2Euler(C)
%NFDCM2EULER ZYX Euler angles [rad] of a body-to-navigation DCM C = Rz(yaw)*Ry(pitch)*Rx(roll).
%   roll in (-pi, pi], pitch in [-pi/2, pi/2], yaw in (-pi, pi] (heading of body x).

roll  = atan2(C(3,2), C(3,3));
pitch = -asin(max(-1, min(1, C(3,1))));
yaw   = atan2(C(2,1), C(1,1));
end
