function g = Gravity_NED(L_b,h_b)
% Inputs:
%   L_b           latitude (rad)
%   h_b           height (m)
% Outputs:
%   g       Acceleration due to gravity (m/s^2)
global W_ie Eq_rad Pol_rad Flat Eccen mu

% Parameters
% Eq_rad = 6378137; %WGS84 Equatorial radius in meters
% Pol_rad = 6356752.31425; %WGS84 Polar radius in meters
% Eccen = 0.0818191908425; %WGS84 eccentricity
% Flat = 1 / 298.257223563; %WGS84 flattening
% mu = 3.986004418E14; %WGS84 Earth gravitational constant (m^3 s^-2)
% W_ie = 7.292115E-5;  % Earth rotation rate (rad/s)

% Calculate surface gravity using the Somigliana model, (2.134)
sinsqL = sin(L_b)^2;
g_0 = 9.7803253359 * (1 + 0.001931853 * sinsqL) / sqrt(1 - Eccen^2 * sinsqL);

% Calculate north gravity using (2.140)
g(1,1) = -8.08E-9 * h_b * sin(2 * L_b);

% East gravity is zero
g(2,1) = 0;

% Calculate down gravity using (2.139)
g(3,1) = g_0 * (1 - (2 / Eq_rad) * (1 + Flat * (1 - 2 * sinsqL) +...
    (W_ie^2 * Eq_rad^2 * Pol_rad / mu)) * h_b + (3 * h_b^2 / Eq_rad^2));

end