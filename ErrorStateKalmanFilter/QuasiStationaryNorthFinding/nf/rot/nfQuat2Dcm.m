function C = nfQuat2Dcm(q)
%NFQUAT2DCM Body-to-navigation DCM of a unit quaternion q_b^n = [q0; q1; q2; q3] (scalar first, Hamilton).

q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
C = [q0^2 + q1^2 - q2^2 - q3^2, 2*(q1*q2 - q0*q3),         2*(q1*q3 + q0*q2);
     2*(q1*q2 + q0*q3),         q0^2 - q1^2 + q2^2 - q3^2, 2*(q2*q3 - q0*q1);
     2*(q1*q3 - q0*q2),         2*(q2*q3 + q0*q1),         q0^2 - q1^2 - q2^2 + q3^2];
end
