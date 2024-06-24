function [R_N,R_E]= Rad_Curv(L)
%   L ->  geodetic latitude (rad) [Pg 78, Fig 2.19]
%   R_N   meridian radius of curvature (m)
%   R_E   transverse radius of curvature (m)
global Eq_rad Eccen

% Calculate meridian radius of curvature [Eq 2.105]
Denominator = 1 - (Eccen * sin(L))^2; 
R_N = Eq_rad * (1 - Eccen^2) / Denominator^1.5;

% Calculate transverse radius of curvature  [Eq 2.105]
R_E = Eq_rad / sqrt(Denominator);

end