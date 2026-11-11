function [imu, truth] = nfStaticImu(cfg, varargin)
%NFSTATICIMU Synthetic static IMU record at the cfg site, in the nfReadImu format (sensor axes).
%   [imu, truth] = nfStaticImu(cfg, Name, Value, ...)
%   Body model (C = C_b^n of the true attitude, z down):
%       w_b = C' w_ie + b_g + white(ARW) + tone
%       f_b = C' [0; 0; -g] + b_a + white(VRW) + tone + C' [gm_N; gm_E; 0]
%   and the sensor outputs are gyroAxes' * w_b and accelAxes' * f_b, so nfProcess with the same
%   cfg must recover the truth.
%   Options (defaults):
%     'HeadingDeg' (343.1), 'RollDeg' (0), 'PitchDeg' (0)   true attitude of the body [deg]
%     'Index' (2:12006)       sample indices, t = index/fs (the time base of RawData.csv)
%     'ArwDph' (cfg.sensor.arw_dph), 'VrwUg' (cfg.sensor.vrw_ug)   white noise per body axis
%                             [deg/sqrt(h)], [ug/sqrt(Hz)]
%     'BiasGDph' ([0 0 0]), 'BiasAmg' ([0 0 0])   constant biases, body axes [deg/h], [mg]
%     'GmSigma' (0), 'GmTau' (2)   Gauss-Markov horizontal specific force (N, E) [m/s^2], [s]
%     'ToneHz' (99.948), 'ToneGyroDph' ([0 0 0]), 'ToneAccelUg' ([0 0 0])   rms of a sinusoidal
%                             line per body axis (RawData.csv: [6.65 7.30 0.37], [169 152 130])
%     'Seed' (0)              random stream seed
%   truth: heading_deg, roll_deg, pitch_deg, C, b_g [rad/s], b_a [m/s^2] (body).

p = inputParser;
p.addParameter('HeadingDeg', 343.1);
p.addParameter('RollDeg', 0);
p.addParameter('PitchDeg', 0);
p.addParameter('Index', 2:12006);
p.addParameter('ArwDph', cfg.sensor.arw_dph);
p.addParameter('VrwUg', cfg.sensor.vrw_ug);
p.addParameter('BiasGDph', [0 0 0]);
p.addParameter('BiasAmg', [0 0 0]);
p.addParameter('GmSigma', 0);
p.addParameter('GmTau', 2);
p.addParameter('ToneHz', 99.948);
p.addParameter('ToneGyroDph', [0 0 0]);
p.addParameter('ToneAccelUg', [0 0 0]);
p.addParameter('Seed', 0);
p.parse(varargin{:});
o = p.Results;

D2R = pi / 180;
k = nfEarthConstants();
fs = cfg.sensor.fs;
Ts = 1 / fs;
idx = o.Index(:)';
t = idx / fs;
N = numel(t);
L = cfg.site.lat_deg * D2R;
rs = RandStream('mt19937ar', 'Seed', o.Seed);

C = nfEuler2Dcm(o.RollDeg * D2R, o.PitchDeg * D2R, o.HeadingDeg * D2R);
b_g = o.BiasGDph(:) * k.DPH;
b_a = o.BiasAmg(:) * 1e-3 * k.g_std;
Ng = (o.ArwDph(:) .* ones(3, 1)) * D2R / 60;          % rad/sqrt(s)
Na = (o.VrwUg(:) .* ones(3, 1)) * k.UG;               % m/s^2/sqrt(Hz)
W = C' * nfEarthRate(L) + b_g + (Ng / sqrt(Ts)) .* randn(rs, 3, N);
F = C' * [0; 0; -nfGravity(L, cfg.site.h_m)] + b_a + (Na / sqrt(Ts)) .* randn(rs, 3, N);

ph = 2 * pi * rand(rs, 6, 1);
arg = 2 * pi * o.ToneHz * t;
W = W + sqrt(2) * (o.ToneGyroDph(:) * k.DPH) .* cos(arg + ph(1:3));
F = F + sqrt(2) * (o.ToneAccelUg(:) * k.UG) .* cos(arg + ph(4:6));

if o.GmSigma > 0
    a = exp(-Ts / o.GmTau);
    e = zeros(2, N);
    e(:, 1) = o.GmSigma * randn(rs, 2, 1);                     % stationary start
    wn = o.GmSigma * sqrt(1 - a^2) * randn(rs, 2, N);
    for i = 1:2
        e(i, 2:N) = filter(1, [1 -a], wn(i, 2:N), a * e(i, 1));
    end
    F = F + C' * [e; zeros(1, N)];
end

imu.t = t;
imu.gyro = cfg.sensor.gyroAxes' * W;
imu.accel = cfg.sensor.accelAxes' * F;
imu.index = idx;
imu.fs = fs;
imu.report = struct('source', 'nfStaticImu', 'rows_dropped', zeros(1, 0), 'replaced', zeros(0, 5), ...
    'gaps', zeros(0, 2));
truth = struct('heading_deg', o.HeadingDeg, 'roll_deg', o.RollDeg, 'pitch_deg', o.PitchDeg, 'C', C, ...
    'b_g', b_g, 'b_a', b_a);
end
