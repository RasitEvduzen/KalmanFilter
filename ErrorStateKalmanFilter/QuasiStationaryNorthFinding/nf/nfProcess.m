function res = nfProcess(cfg, imu)
%NFPROCESS The steps of main.m as one function, for batch runs (e.g. nfMonteCarlo).
%   res = nfProcess(cfg)        reads cfg.dataFile
%   res = nfProcess(cfg, imu)   uses an IMU record struct (nfReadImu or nfStaticImu)
%   No printing and no figures; failing axis maps give the warning 'nf:axes'.
%   res: the nfEskf output completed by nfSummary.

if nargin < 2
    imu = nfReadImu(cfg.dataFile, cfg.sensor.fs);
end
W_b = cfg.sensor.gyroAxes * imu.gyro;
F_b = cfg.sensor.accelAxes * imu.accel;

alignment = nfGyrocompass(imu.t, W_b, F_b, cfg.align.window, cfg);
axesCheck = nfAxisCheck(alignment, cfg);
if ~axesCheck.ok
    warning('nf:axes', 'Axis maps fail the static check: %s.', axesCheck.msg);
end
staticGc = repmat(alignment, 1, 0);
for iw = 1:size(cfg.align.gcWindows, 1)
    staticGc(iw) = nfGyrocompass(imu.t, W_b, F_b, cfg.align.gcWindows(iw, :), cfg);
end
res = nfEskf(imu.t, W_b, F_b, alignment, cfg);
res = nfSummary(res, alignment, staticGc, axesCheck, imu, cfg);
end
