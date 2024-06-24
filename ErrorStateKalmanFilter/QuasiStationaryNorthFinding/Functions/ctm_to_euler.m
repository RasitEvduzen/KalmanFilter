function [roll,pitch,yaw] = ctm_to_euler(ctm)
% Coordinate Transform Matrix (Rotation Matrix) To Euler Angles
% ---------------- Input --------------
% Rotx * Roty * Rotz  -> 3x3 Rotation Matrix
% 
% ---------------- Output --------------

% roll  = atan2(ctm(2,3),ctm(3,3)); % rad
% pitch = -asin(ctm(1,3));          % rad
% yaw   = atan2(ctm(1,2),ctm(1,1)); % rad

tmp_angle = rotm2eul(ctm,'XYZ');
roll = tmp_angle(1);
pitch = tmp_angle(2);
yaw = tmp_angle(3);
end




