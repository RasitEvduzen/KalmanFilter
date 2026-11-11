%% Quasi-stationary north finding with the IFOS-500 FOG IMU (RawData.csv) - main script
% Run with F5 (or type main). Every step leaves its result in the workspace:
%   imu         cleaned record in sensor axes (t [s], gyro [rad/s], accel [m/s^2])
%   W_b, F_b    angular rate and specific force in body axes
%   alignment   leveling + gyrocompassing on the stable span (initial attitude of the filter)
%   axesCheck   physical check of the sensor-to-body axis maps
%   staticGc    static gyrocompass headings on the report windows
%   res         ESKF result; res.summary holds the key numbers
% Settings: nf/nfConfig.m. Reference: P. D. Groves, Principles of GNSS, Inertial, and
% Multisensor Integrated Navigation Systems, 2nd ed.

projectDir = fileparts(mfilename('fullpath'));
if isempty(projectDir)                          % sections run one by one from the project folder
    projectDir = pwd;
end
addpath(genpath(fullfile(projectDir, 'nf')));

%% Settings
nfCfg = nfConfig();
% nfCfg = nfConfig('filter.attitude', 'dcm');                      % DCM instead of quaternion attitude
% nfCfg = nfConfig('sensor.bg_dph', 0.3, 'filter.t_start', 20);   % example: other settings

%% 1. IMU record: read, clean, sensor axes -> body axes (x = -Y_s, y = -Z_s, z = +X_s down)
imu = nfReadImu(nfCfg.dataFile, nfCfg.sensor.fs);
W_b = nfCfg.sensor.gyroAxes * imu.gyro;         % angular rate, body [rad/s]
F_b = nfCfg.sensor.accelAxes * imu.accel;       % specific force, body [m/s^2] (accelerometers output -f)

%% 2. Leveling + gyrocompassing on the stable span (initial attitude) and axis-map check
alignment = nfGyrocompass(imu.t, W_b, F_b, nfCfg.align.window, nfCfg);
axesCheck = nfAxisCheck(alignment, nfCfg);
if ~axesCheck.ok
    warning('nf:axes', 'Axis maps fail the static check: %s.', axesCheck.msg);
end

%% 3. Static gyrocompass on the report windows
staticGc = repmat(alignment, 1, 0);
for iw = 1:size(nfCfg.align.gcWindows, 1)
    staticGc(iw) = nfGyrocompass(imu.t, W_b, F_b, nfCfg.align.gcWindows(iw, :), nfCfg);
end

%% 4. Error-state Kalman filter with ZUPT
res = nfEskf(imu.t, W_b, F_b, alignment, nfCfg);

%% 5. Results
res = nfSummary(res, alignment, staticGc, axesCheck, imu, nfCfg);
nfPrintResults(res);
