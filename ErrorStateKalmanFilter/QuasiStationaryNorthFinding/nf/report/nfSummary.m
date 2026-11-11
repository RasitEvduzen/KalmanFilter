function res = nfSummary(res, alignment, staticGc, axesCheck, imu, cfg)
%NFSUMMARY Attach the alignment results to the nfEskf output and collect the key numbers.
%   res = nfSummary(res, alignment, staticGc, axesCheck, imu, cfg) adds res.align (nfGyrocompass
%   on cfg.align.window), res.check (nfAxisCheck), res.gc (static gyrocompass, one per
%   cfg.align.gcWindows row), res.imuReport, res.cfg and res.summary:
%     heading_deg, sigma_deg            ESKF heading of body x and 1-sigma (last state)
%     coarse_heading_deg, coarse_sigma_deg, roll0_deg, pitch0_deg   initial alignment
%     gc_heading_deg, gc_sigma_deg      static gyrocompass per report window
%     t_start, t_end, n_samples, n_updates, nis_mean
%     b_a_mg, b_a_sigma_mg, b_g_dph, b_g_sigma_dph   final bias estimates (body)
%     axes_ok                           result of the axis-map check

k = nfEarthConstants();
R2D = 180 / pi;
s.heading_deg = mod(res.heading_deg(end), 360);
s.sigma_deg = res.heading_sigma_deg(end);
s.coarse_heading_deg = mod(alignment.heading * R2D, 360);
s.coarse_sigma_deg = alignment.sigma * R2D;
s.roll0_deg = alignment.roll * R2D;
s.pitch0_deg = alignment.pitch * R2D;
s.gc_heading_deg = mod([staticGc.heading] * R2D, 360);
s.gc_sigma_deg = [staticGc.sigma] * R2D;
s.t_start = res.time(1);
s.t_end = res.time(end);
s.n_samples = res.N;
s.n_updates = numel(res.nis);
s.nis_mean = mean(res.nis);
s.b_a_mg = res.b_a(:, end) / (1e-3 * k.g_std);
s.b_a_sigma_mg = res.b_a_sigma(:, end) / (1e-3 * k.g_std);
s.b_g_dph = res.b_g(:, end) / k.DPH;
s.b_g_sigma_dph = res.b_g_sigma(:, end) / k.DPH;
s.axes_ok = axesCheck.ok;

res.align = alignment;
res.check = axesCheck;
res.gc = staticGc;
res.imuReport = imu.report;
res.cfg = cfg;
res.summary = s;
end
