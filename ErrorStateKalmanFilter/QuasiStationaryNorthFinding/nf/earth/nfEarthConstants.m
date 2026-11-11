function k = nfEarthConstants()
%NFEARTHCONSTANTS Earth model (WGS-84) and unit conversion constants (SI units).
%   W_ie    Earth rotation rate, one turn per sidereal day of 86164 s [rad/s]
%   R_0     equatorial radius [m]          R_P   polar radius [m]
%   f       flattening [-]                 e     eccentricity [-]
%   mu      gravitational constant [m^3/s^2]
%   g_std   standard gravity [m/s^2] (unit conversions only)
%   DPH     deg/h -> rad/s                 UG    micro-g -> m/s^2

k.W_ie  = (360 / 86164) * (pi / 180);
k.R_0   = 6378137;
k.R_P   = 6356752.31425;
k.f     = 1 / 298.257223563;
k.e     = 0.0818191908425;
k.mu    = 3.986004418e14;
k.g_std = 9.80665;
k.DPH   = (pi / 180) / 3600;
k.UG    = 9.80665e-6;
end
