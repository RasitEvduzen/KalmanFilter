function n = nfIdentifyNoise(t, W_b, F_b, window, fs)
%NFIDENTIFYNOISE Noise characteristics of a static record, per body axis.
%   n = nfIdentifyNoise(t, W_b, F_b, window, fs) analyses the samples with window(1) <= t <
%   window(2) [s] of the body rates W_b [rad/s] and specific forces F_b [m/s^2] (3xN):
%     n.gyro.tau, n.gyro.adev       Allan deviation curves (3 x K) [deg/h]
%     n.gyro.arw_adev               ADEV at tau = 1 s [deg/sqrt(h)]
%     n.gyro.arw_psd                low-frequency PSD floor (0.05-1 Hz, two-sided) [deg/sqrt(h)]
%     n.gyro.bias_instability       minimum ADEV / 0.664 [deg/h]
%     n.accel.(same fields)         ADEV [ug], vrw_adev and vrw_psd [ug/sqrt(Hz)], bias
%                                   instability [ug]
%   The low-frequency PSD floor avoids the resonances and the near-Nyquist line of the record;
%   the configuration values (nfConfig sensor.arw_dph, vrw_ug) come from this analysis.

k = nfEarthConstants();
m = (t >= window(1) - 1e-9) & (t < window(2) - 1e-9);
toDsh = (180 / pi) * 60;                   % rad/s/sqrt(Hz) -> deg/sqrt(h)
toDph = (180 / pi) * 3600;                 % rad/s -> deg/h
g = analyse(W_b(:, m), fs, toDph, toDsh);
a = analyse(F_b(:, m), fs, 1 / k.UG, 1 / k.UG);
n.gyro = struct('tau', g.tau, 'adev', g.adev, 'arw_adev', g.white_adev, 'arw_psd', g.white_psd, ...
    'bias_instability', g.bias_instability);
n.accel = struct('tau', a.tau, 'adev', a.adev, 'vrw_adev', a.white_adev, 'vrw_psd', a.white_psd, ...
    'bias_instability', a.bias_instability);
n.window = window;
end

% =====================================================================
function a = analyse(X, fs, toAdev, toWhite)
for i = 3:-1:1
    x = X(i, :) - mean(X(i, :));
    [tau, ad] = nfAllan(x, fs);
    [~, ad1] = nfAllan(x, fs, round(fs));
    [f, P] = nfPsd(x, fs);
    lf = f > 0.05 & f < 1.0;
    a.adev(i, :) = ad * toAdev;
    a.white_adev(i) = ad1 * toWhite;                        % ADEV(1 s) * sqrt(1 s)
    a.white_psd(i) = sqrt(mean(P(lf)) / 2) * toWhite;      % two-sided white level
    a.bias_instability(i) = min(ad) / 0.664 * toAdev;
end
a.tau = tau;
end
