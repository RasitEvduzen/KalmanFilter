function nfPrintResults(res)
%NFPRINTRESULTS Print the north-finding result (res from nfSummary) to the command window.

s = res.summary;
c = res.check;
cfg = res.cfg;
fprintf('\n=== North finding: %s ===\n', res.imuReport.source);
fprintf('site: latitude %.6f deg, longitude %.6f deg, height %.2f m\n', cfg.site.lat_deg, cfg.site.lon_deg, cfg.site.h_m);
fprintf('axis check: |f| %.4f m/s^2 (model %.4f), vertical Earth rate %+.3f deg/h (model %+.3f), horizontal %.3f (model %.3f)', ...
    c.g, c.g_model, c.w_up_dph, c.w_up_model_dph, c.w_h_dph, c.w_h_model_dph);
if c.ok
    fprintf(' -> OK\n');
else
    fprintf('\n  WARNING: %s\n', c.msg);
end
w = cfg.align.window;
fprintf('leveling + gyrocompass %g-%s s: roll %.4f, pitch %.4f, heading %.3f +- %.2f deg\n', w(1), endLabel(w(2)), ...
    s.roll0_deg, s.pitch0_deg, s.coarse_heading_deg, s.coarse_sigma_deg);
for i = 1:numel(res.gc)
    w = res.gc(i).window;
    fprintf('static gyrocompass %g-%s s: %.3f +- %.2f deg\n', w(1), endLabel(w(2)), s.gc_heading_deg(i), s.gc_sigma_deg(i));
end
fprintf('ESKF (%s attitude) %.2f-%.2f s: %d samples, %d ZUPT updates, mean NIS %.3f (3 dof; conservative ZUPT sigma)\n', ...
    cfg.filter.attitude, s.t_start, s.t_end, s.n_samples, s.n_updates, s.nis_mean);
fprintf('  >>> HEADING of the case -Y axis (body x), true north, clockwise: %.3f +- %.3f deg\n', s.heading_deg, s.sigma_deg);
fprintf('  gyro bias  [%s] deg/h (1 sigma [%s])\n', sprintf(' %+.4f', s.b_g_dph), sprintf(' %.4f', s.b_g_sigma_dph));
fprintf('  accel bias [%s] mg    (1 sigma [%s]; x, y not observable at one position)\n', ...
    sprintf(' %+.4f', s.b_a_mg), sprintf(' %.4f', s.b_a_sigma_mg));
fprintf('  the sigma assumes a residual gyro bias of %.2f deg/h (%.2f deg of heading)\n', cfg.sensor.bg_dph, ...
    res.align.sigma_bias * 180 / pi);
end

% =====================================================================
function s = endLabel(t1)
if isinf(t1)
    s = 'end';
else
    s = sprintf('%g', t1);
end
end
