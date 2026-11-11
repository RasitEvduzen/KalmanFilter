function a = nfGyrocompass(t, W_b, F_b, window, cfg)
%NFGYROCOMPASS Leveling and gyrocompassing of a static IMU on a time window (Groves eqs. 5.101, 5.105).
%   a = nfGyrocompass(t, W_b, F_b, window, cfg)
%   t (1xN) record time [s]; W_b (3xN) angular rate [rad/s] and F_b (3xN) specific force
%   [m/s^2] in body axes (x forward, y right, z down); window [t0 t1) [s] (t1 = Inf: to the
%   end). With f and w the window means:
%       roll  = atan2(-f_y, -f_z),      pitch = atan(f_x / sqrt(f_y^2 + f_z^2))
%       W cosL sin(psi) = -w_y cos(roll) + w_z sin(roll)
%       W cosL cos(psi) =  w_x cos(pitch) + w_y sin(roll) sin(pitch) + w_z cos(roll) sin(pitch)
%   psi = heading of body x (exact for any roll and |pitch| < 90 deg). 1-sigma accuracy:
%       sigma^2 = (N_E / (W cosL sqrt(T)))^2 + (sigma_bg / (W cosL))^2 + (sigma_ba tanL / g)^2
%   N_E = gyro noise density projected on east (cfg.sensor.arw_dph), sigma_bg = residual gyro
%   bias (cfg.sensor.bg_dph), sigma_ba = accelerometer bias (cfg.sensor.ba_mg), T = span
%   length, L = cfg.site.lat_deg.
%
%   a: roll, pitch, heading, sigma, sigma_noise, sigma_bias, sigma_accel [rad], C (C_b^n),
%   f_mean, w_mean (3x1), n (samples), T [s], window.

D2R = pi / 180;
k = nfEarthConstants();
m = (t >= window(1) - 1e-9) & (t < window(2) - 1e-9);
if ~any(m)
    error('nf:align:window', 'No samples in the window [%g, %g) s; the record covers [%g, %g] s.', ...
        window(1), window(2), t(1), t(end));
end
f = mean(F_b(:, m), 2);
w = mean(W_b(:, m), 2);

roll  = atan2(-f(2), -f(3));
pitch = atan(f(1) / sqrt(f(2)^2 + f(3)^2));
heading = atan2(-w(2)*cos(roll) + w(3)*sin(roll), ...
    w(1)*cos(pitch) + w(2)*sin(roll)*sin(pitch) + w(3)*cos(roll)*sin(pitch));

C = nfEuler2Dcm(roll, pitch, heading);
N_b = (cfg.sensor.arw_dph(:)' * D2R / 60) .* ones(1, 3);      % rad/sqrt(s), body x y z
N_E = sqrt(sum((C(2, :).^2) .* (N_b.^2)));
OcL = k.W_ie * cos(cfg.site.lat_deg * D2R);
n = nnz(m);
T = n * median(diff(t));

a.roll = roll;
a.pitch = pitch;
a.heading = heading;
a.sigma_noise = N_E / (OcL * sqrt(T));
a.sigma_bias = cfg.sensor.bg_dph * k.DPH / OcL;
L = cfg.site.lat_deg * D2R;
a.sigma_accel = tan(L) * cfg.sensor.ba_mg * 1e-3 * k.g_std / nfGravity(L, cfg.site.h_m);
a.sigma = sqrt(a.sigma_noise^2 + a.sigma_bias^2 + a.sigma_accel^2);
a.C = C;
a.f_mean = f;
a.w_mean = w;
a.n = n;
a.T = T;
a.window = window;
end
