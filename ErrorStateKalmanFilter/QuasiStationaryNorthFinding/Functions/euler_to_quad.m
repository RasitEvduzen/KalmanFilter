function [q] = euler_to_quad(roll,pitch,yaw)
% Coordinate Transform Matrix(Rotation Matrix) to Quaternion
% RotX * RotY * RotZ -> CTM
% ----------------------------- Input -------------------------------------
% Roll  angle  rotation about  X axis Degree
% Pitch angle  rotation about  Y axis Degree
% Yaw   angle  rotation about  Z axis Degree
% ----------------------------- Output ------------------------------------
% Hamiltonian Quaternion

% ------------ Compute Quaternion -----------------------------------------


q(1,1) = cosd(roll/2)*cosd(pitch/2)*cosd(yaw/2)+sind(roll/2)*sind(pitch/2)*sind(yaw/2);
q(2,1) = sind(roll/2)*cosd(pitch/2)*cosd(yaw/2)-cosd(roll/2)*sind(pitch/2)*sind(yaw/2);
q(3,1) = cosd(roll/2)*sind(pitch/2)*cosd(yaw/2)+sind(roll/2)*cosd(pitch/2)*sind(yaw/2);
q(4,1) = cosd(roll/2)*cosd(pitch/2)*sind(yaw/2)-sind(roll/2)*sind(pitch/2)*cosd(yaw/2);
end