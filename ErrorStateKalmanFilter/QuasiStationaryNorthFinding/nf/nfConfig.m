function cfg = nfConfig(varargin)
%NFCONFIG Configuration of the north-finding pipeline (SI units unless the name says otherwise).
%   cfg = nfConfig()                            defaults: RawData.csv, IFOS-500, case X axis down
%   cfg = nfConfig('group.field', value, ...)   overrides, e.g. nfConfig('sensor.bg_dph', 0.3)
%   cfg = nfConfig(S, ...)                      overrides from a (partial) struct of the same layout
%                                                (e.g. a saved res.cfg), then name-value pairs
%
%   dataFile          IMU record, rows [index, GX GY GZ (deg/s), AX AY AZ (m/s^2)], sensor axes
%   site              lat_deg, lon_deg [deg], h_m [m] of the 2023 recording
%   sensor.fs         sampling rate [Hz]
%   sensor.gyroAxes   3x3 signed permutation: w_b = gyroAxes * [GX; GY; GZ]
%   sensor.accelAxes  3x3 signed permutation: f_b = accelAxes * [AX; AY; AZ]
%                     During the recording the case +X axis pointed down and the accelerometers
%                     output -f (+g on an axis pointing down). Body frame (z down): x = -Y_s
%                     (the heading reference), y = -Z_s, z = +X_s. accelAxes = -gyroAxes.
%   sensor.arw_dph    gyro angle random walk, body x y z [deg/sqrt(h)]
%   sensor.vrw_ug     accelerometer velocity random walk, body x y z [ug/sqrt(Hz)]
%                     Both are rounded values chosen within the spread of the white-noise
%                     estimates of RawData.csv (Allan deviation at 1 s and low-frequency PSD
%                     floor, 45 s span and full record; nfIdentifyNoise, paper Table 3).
%   sensor.bg_dph     gyro bias 1-sigma [deg/h]; IFOS-500: 0.1 at fixed temperature, 0.3 over
%                     -40..+60 C. It limits the heading accuracy: 5.04 deg per deg/h at this site.
%   sensor.ba_mg      accelerometer bias 1-sigma [mg] (IFOS-500: 0.5)
%   align.window      [t0 t1) leveling + gyrocompassing span that initialises the filter [s]
%                     (t1 = Inf: end of the record); the first 15 s hold the accelerometer warm-up
%   align.gcWindows   K x 2 [t0 t1) spans of the reported static gyrocompass headings [s]
%   filter.attitude   'quaternion' (q_b^n) or 'dcm' (C_b^n): representation of the nominal attitude;
%                     the error state (rotation vector dpsi) and the filter equations are the same
%   filter.t_start, t_end        filter span [s] (t_end = Inf: end of the record)
%   filter.p0_tilt_deg, p0_head_deg  initial 1-sigma tilt and heading errors [deg]; the heading
%                     prior is wide because the coarse heading comes from the same data
%   filter.p0_vel, p0_pos         initial 1-sigma velocity [m/s] and displacement [m]
%   filter.q_ba_rw    accelerometer bias random walk [m/s^2/sqrt(s)]
%   filter.q_bg_step  gyro bias process noise per step [(rad/s)^2]
%   filter.gm_tau, gm_sigma      horizontal specific-force disturbance (floor vibration), first-
%                     order Gauss-Markov: correlation time [s], 1-sigma [m/s^2]
%   filter.zupt_rate, zupt_sigma ZUPT rate [Hz] and 1-sigma [m/s]; 1e-4 m/s is about 3x the real
%                     residual motion, to cover its ~1 s correlation (mean NIS ~0.44 instead of 3)

cfg.dataFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'RawData.csv');
cfg.site = struct('lat_deg', 40.923061376, 'lon_deg', 29.31453938, 'h_m', 9.9355);
cfg.sensor = struct( ...
    'fs', 200, ...
    'gyroAxes', [0 -1 0; 0 0 -1; 1 0 0], ...
    'accelAxes', [0 1 0; 0 0 1; -1 0 0], ...
    'arw_dph', [0.0140 0.0150 0.0069], ...
    'vrw_ug', [6.9 5.7 4.3], ...
    'bg_dph', 0.1, ...
    'ba_mg', 0.5);
cfg.align = struct('window', [15 Inf], 'gcWindows', [15 Inf; 0 Inf]);
cfg.filter = struct( ...
    'attitude', 'quaternion', ...
    't_start', 15, 't_end', Inf, ...
    'p0_tilt_deg', 0.05, 'p0_head_deg', 30, 'p0_vel', 1e-2, 'p0_pos', 1e-3, ...
    'q_ba_rw', 1e-6, 'q_bg_step', 1e-19, ...
    'gm_tau', 2, 'gm_sigma', 6e-5, ...
    'zupt_rate', 10, 'zupt_sigma', 1e-4);

args = varargin;
if ~isempty(args) && isstruct(args{1})
    cfg = merge(cfg, args{1}, '');
    args = args(2:end);
end
if mod(numel(args), 2)
    error('nf:config:args', 'Overrides must be name-value pairs (after an optional struct).');
end
for i = 1:2:numel(args)
    name = char(args{i});
    parts = strsplit(name, '.');
    if ~hasPath(cfg, parts)
        error('nf:config:field', 'Unknown configuration field ''%s''.', name);
    end
    cfg = setfield(cfg, parts{:}, args{i + 1});
end
validate(cfg);
end

% =====================================================================
function cfg = merge(cfg, S, prefix)
for f = fieldnames(S)'
    name = [prefix f{1}];
    if ~isfield(cfg, f{1})
        error('nf:config:field', 'Unknown configuration field ''%s''.', name);
    end
    if isstruct(cfg.(f{1})) && isstruct(S.(f{1}))
        cfg.(f{1}) = merge(cfg.(f{1}), S.(f{1}), [name '.']);
    else
        cfg.(f{1}) = S.(f{1});
    end
end
end

function ok = hasPath(s, parts)
for i = 1:numel(parts)
    if ~isstruct(s) || ~isfield(s, parts{i})
        ok = false;
        return
    end
    s = s.(parts{i});
end
ok = ~isstruct(s);
end

function validate(c)
num = @(v, n) isnumeric(v) && isreal(v) && any(numel(v) == n) && all(isfinite(v(:)));
pos = @(v) num(v, 1) && v > 0;
nneg = @(v) num(v, 1) && v >= 0;
window = @(v) isnumeric(v) && isreal(v) && size(v, 2) == 2 && ~any(isnan(v(:))) && ...
    all(isfinite(v(:, 1))) && all(v(:, 1) < v(:, 2));
axisMap = @(M) isnumeric(M) && isequal(size(M), [3 3]) && all(ismember(M(:), [-1 0 1])) && ...
    all(sum(abs(M), 1) == 1) && all(sum(abs(M), 2) == 1);

check(ischar(c.dataFile) || (isstring(c.dataFile) && isscalar(c.dataFile)), 'dataFile', 'a file name');
check(num(c.site.lat_deg, 1) && abs(c.site.lat_deg) < 89, 'site.lat_deg', 'a latitude in (-89, 89) deg');
check(num(c.site.lon_deg, 1), 'site.lon_deg', 'a longitude [deg]');
check(num(c.site.h_m, 1), 'site.h_m', 'a height [m]');
check(pos(c.sensor.fs), 'sensor.fs', 'a positive rate [Hz]');
check(axisMap(c.sensor.gyroAxes), 'sensor.gyroAxes', 'a 3x3 signed permutation matrix');
check(axisMap(c.sensor.accelAxes) && (isequal(c.sensor.accelAxes, c.sensor.gyroAxes) || ...
    isequal(c.sensor.accelAxes, -c.sensor.gyroAxes)), 'sensor.accelAxes', ...
    '+gyroAxes or -gyroAxes (both triads map to the same body axes; the sign is the output convention)');
check(num(c.sensor.arw_dph, [1 3]) && all(c.sensor.arw_dph >= 0), 'sensor.arw_dph', '1 or 3 non-negative values');
check(num(c.sensor.vrw_ug, [1 3]) && all(c.sensor.vrw_ug >= 0), 'sensor.vrw_ug', '1 or 3 non-negative values');
check(nneg(c.sensor.bg_dph), 'sensor.bg_dph', 'a non-negative value [deg/h]');
check(nneg(c.sensor.ba_mg), 'sensor.ba_mg', 'a non-negative value [mg]');
check(window(c.align.window) && size(c.align.window, 1) == 1, 'align.window', 'a [t0 t1) window [s]');
check(isempty(c.align.gcWindows) || window(c.align.gcWindows), 'align.gcWindows', 'K x 2 [t0 t1) windows [s]');
check(any(strcmp(c.filter.attitude, {'quaternion', 'dcm'})), 'filter.attitude', '''quaternion'' or ''dcm''');
check(num(c.filter.t_start, 1), 'filter.t_start', 'a time [s]');
check(isnumeric(c.filter.t_end) && isscalar(c.filter.t_end) && c.filter.t_end > c.filter.t_start, ...
    'filter.t_end', 'a time after t_start [s] (Inf: end of the record)');
for f = {'p0_tilt_deg', 'p0_head_deg', 'p0_vel', 'p0_pos', 'q_ba_rw', 'q_bg_step', 'gm_sigma'}
    check(nneg(c.filter.(f{1})), ['filter.' f{1}], 'a non-negative value');
end
for f = {'gm_tau', 'zupt_rate', 'zupt_sigma'}
    check(pos(c.filter.(f{1})), ['filter.' f{1}], 'a positive value');
end
end

function check(ok, name, what)
if ~ok
    error('nf:config:value', 'Invalid configuration field ''%s'': expected %s.', name, what);
end
end
