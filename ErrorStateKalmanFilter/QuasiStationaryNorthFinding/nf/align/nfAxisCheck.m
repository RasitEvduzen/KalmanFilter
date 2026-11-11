function c = nfAxisCheck(a, cfg)
%NFAXISCHECK Physical consistency of the sensor-to-body axis maps on a static span.
%   c = nfAxisCheck(a, cfg) uses the window means of nfGyrocompass:
%     gravity     |f| = g and f along -z_b (body z down: tilt < 30 deg)
%     Earth rate  vertical component w.f/|f| = +W sin(L) (an axis map with the wrong output
%                 convention flips its sign: implied latitude ~ -L) and horizontal component
%                 W cos(L)
%   c: g, g_model [m/s^2], tilt_deg, w_up_dph, w_up_model_dph, w_h_dph, w_h_model_dph [deg/h],
%   lat_implied_deg, ok (all checks passed), msg (what failed, '' if ok).

k = nfEarthConstants();
L = cfg.site.lat_deg * pi / 180;
f = a.f_mean;
w = a.w_mean;
up = f / norm(f);                                % the specific force points up at rest

c.g = norm(f);
c.g_model = nfGravity(L, cfg.site.h_m);
c.tilt_deg = acosd(max(-1, min(1, -up(3))));
c.w_up_dph = (w' * up) / k.DPH;
c.w_up_model_dph = k.W_ie * sin(L) / k.DPH;
c.w_h_dph = sqrt(max(w' * w - (w' * up)^2, 0)) / k.DPH;
c.w_h_model_dph = k.W_ie * cos(L) / k.DPH;
c.lat_implied_deg = asind(max(-1, min(1, (w' * up) / k.W_ie)));

msg = {};
if abs(c.g / c.g_model - 1) > 0.01
    msg{end + 1} = sprintf('|f| = %.4f m/s^2 instead of %.4f (units, dead channel or motion?)', c.g, c.g_model);
end
if c.tilt_deg > 30
    msg{end + 1} = sprintf('body z is %.1f deg from the down direction (accelAxes?)', c.tilt_deg);
end
if abs(c.w_up_dph - c.w_up_model_dph) > max(0.5 * abs(c.w_up_model_dph), 2)
    msg{end + 1} = sprintf(['vertical Earth rate %+.2f deg/h instead of %+.2f (implied latitude %.2f deg): ' ...
        'gyroAxes and accelAxes disagree in sign'], c.w_up_dph, c.w_up_model_dph, c.lat_implied_deg);
end
c.ok = isempty(msg);
c.msg = strjoin(msg, '; ');
end
