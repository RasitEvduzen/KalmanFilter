function g_D = nfGravity(L, h)
%NFGRAVITY Down component of gravity [m/s^2] at latitude L [rad] and height h [m].
%   Somigliana model with the height correction (Groves, 2nd ed., eqs. 2.134 and 2.139). The
%   north component (-8e-9 h sin 2L) is neglected.

k = nfEarthConstants();
sinsqL = sin(L)^2;
g_0 = 9.7803253359 * (1 + 0.001931853*sinsqL) / sqrt(1 - k.e^2*sinsqL);
g_D = g_0 * (1 - (2/k.R_0) * (1 + k.f*(1 - 2*sinsqL) + (k.W_ie^2 * k.R_0^2 * k.R_P / k.mu)) * h ...
      + (3*h^2/k.R_0^2));
end
