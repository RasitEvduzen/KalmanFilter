function [q] = ctm_to_quad(ctm)
% Coordinate Transform Matrix (Rotation Matrix) To Quaternion
% ---------------- Input --------------
% Rotx * Roty * Rotz  -> 3x3 Rotation Matrix

% ---------------- Output --------------
% Coordinate Transformation Matrix  Paul D Groves Principle of GNSS 2.
% edition page[41]

q(1,1) = 0.5 * sqrt(1 + ctm(1,1) + ctm(2,2) + ctm(3,3));
q(2,1) = (ctm(3,2) - ctm(2,3)) / (4 * q(1,1));
q(3,1) = (ctm(1,3) - ctm(3,1)) / (4 * q(1,1));
q(4,1) = (ctm(2,1) - ctm(1,2)) / (4 * q(1,1));


end