function res = nfEskf(t, W_b, F_b, init, cfg)
%NFESKF Closed-loop error-state Kalman filter for quasi-stationary north finding (NED, Groves).
%   res = nfEskf(t, W_b, F_b, init, cfg)
%   t (1xN) record time [s]; W_b (3xN) angular rate [rad/s] and F_b (3xN) specific force
%   [m/s^2] in body axes; init initial attitude with fields roll, pitch, heading [rad] (e.g.
%   nfGyrocompass); cfg configuration (nfConfig). The filter runs on the samples with
%   cfg.filter.t_start <= t < cfg.filter.t_end.
%
%   Nominal state: attitude q_b^n or C_b^n (cfg.filter.attitude), NED velocity v and displacement
%   r, accumulated bias estimates b_a, b_g (body) and the disturbance estimate gm (N, E).
%   Error state (17): x = [dpsi(1:3) dv(4:6) dr(7:9) dba(10:12) dbg(13:15) gm(16:17)]
%     dpsi, dv, dr   estimate - truth (attitude error angles, NED)
%     dba, dbg       residual biases, truth - estimate (body)
%     gm             horizontal specific-force disturbance (floor vibration), Gauss-Markov
%   Continuous model of a static IMU:
%     dpsi' = -[w_ie x] dpsi + C dbg
%     dv'   = -[f x] dpsi + C0 dba + [gm; 0],    f = [0; 0; -g]
%     dr'   = dv,    gm' = -gm/tau + noise
%   C0 is the initial attitude: the accelerometer bias is fixed to the body, and a heading
%   estimate that is still moving would make the unobservable horizontal bias look observable.
%   The static specific force -g keeps vibration out of the attitude-to-velocity coupling.
%   Phi = I + A + A^2/2 + A^3/6 (A = F Ts); Q from the noise densities of cfg.sensor.
%   ZUPT at cfg.filter.zupt_rate: z = -v, H = [0 -I 0 ...], Joseph form. Closed-loop reset:
%     quaternion: q <- q(-dpsi) (x) q, normalised        DCM: C <- orth((I - [dpsi x]) C)
%     v -= dv,  r -= dr,  b_a += dba,  b_g += dbg,  gm += dgm
%   The error state, F, Q and H are the same for both attitude representations. The quaternion
%   reset is the exact rotation; the first-order DCM reset turns by atan|dpsi| instead of |dpsi|,
%   which matters only for the large heading corrections of the first seconds (final heading on
%   RawData.csv: 4e-5 deg apart; with an exact DCM reset the two agree to 1e-8 deg).
%
%   res: time (1x(N+1); state j+1 follows the sample at time(j+1) - Ts), heading_deg (heading
%   of body x, unwrapped, last value in [0, 360)), heading_sigma_deg, roll_deg, pitch_deg,
%   roll_sigma_deg, pitch_sigma_deg (tilt errors rotated into the heading frame), v, r, v_sigma
%   (3x(N+1), NED), b_a [m/s^2], b_g [rad/s], b_a_sigma, b_g_sigma (3x(N+1), body), nis and
%   nis_time (one per ZUPT), P0, P_end, qg, qa (per-step white-noise variances), k0 (first
%   filter sample), N.

D2R = pi / 180;
R2D = 180 / pi;
k = nfEarthConstants();
fp = cfg.filter;
Ts = 1 / cfg.sensor.fs;
L = cfg.site.lat_deg * D2R;
g_D = nfGravity(L, cfg.site.h_m);
w_ie = nfEarthRate(L);

% ---- filter span and input checks
k0 = find(t >= fp.t_start - 1e-9, 1);
if isinf(fp.t_end)
    k1 = numel(t);
else
    k1 = find(t < fp.t_end - 1e-9, 1, 'last');
end
if isempty(k0) || isempty(k1) || k1 <= k0
    error('nf:kf:span', 'Empty filter span [%g, %g) s; the record covers [%g, %g] s.', ...
        fp.t_start, fp.t_end, t(1), t(end));
end
W = W_b(:, k0:k1);
F = F_b(:, k0:k1);
tv = t(k0:k1);
N = numel(tv);
if ~all(isfinite(W(:))) || ~all(isfinite(F(:)))
    error('nf:kf:data', 'Non-finite IMU samples in the filter span.');
end
if any(abs(diff(tv) - Ts) > 1e-6)
    warning('nf:kf:gaps', 'Time gaps in the filter span; the filter assumes Ts = %g s.', Ts);
end
if abs(init.pitch) > 80 * D2R
    error('nf:kf:attitude', 'Body x is %.1f deg from the horizontal: its heading is undefined (check the axis maps).', ...
        init.pitch * R2D);
end

% ---- state layout, initial state and covariance
iA = 1:3; iV = 4:6; iP = 7:9; iBa = 10:12; iBg = 13:15; iGm = 16:17; nx = 17;
quat = strcmp(fp.attitude, 'quaternion');
if quat
    att = nfEuler2Quat(init.roll, init.pitch, init.heading);
    C = nfQuat2Dcm(att);
else
    C = nfEuler2Dcm(init.roll, init.pitch, init.heading);
    att = C;
end
C0 = C;
v = zeros(3, 1); r = zeros(3, 1);
b_a = zeros(3, 1); b_g = zeros(3, 1); gm = zeros(2, 1);
sba = cfg.sensor.ba_mg * 1e-3 * k.g_std .* ones(1, 3);
sbg = cfg.sensor.bg_dph * k.DPH .* ones(1, 3);
P = diag([(fp.p0_tilt_deg * D2R)^2 * [1 1], (fp.p0_head_deg * D2R)^2, fp.p0_vel^2 * [1 1 1], ...
    fp.p0_pos^2 * [1 1 1], sba.^2, sbg.^2, fp.gm_sigma^2 * [1 1]]);
P0 = P;

% ---- process model: constant parts of F and Q
qg = (cfg.sensor.arw_dph .* ones(1, 3) * (D2R / 60)).^2 * Ts;     % gyro white noise per step, body
qa = (cfg.sensor.vrw_ug .* ones(1, 3) * k.UG).^2 * Ts;            % accel white noise per step, body
phi_gm = exp(-Ts / fp.gm_tau);
Fm = zeros(nx);
Fm(iA, iA) = -nfSkew(w_ie);
Fm(iV, iA) = -nfSkew([0; 0; -g_D]);
Fm(iV, iBa) = C0;
Fm(iP, iV) = eye(3);
Fm(4:5, iGm) = eye(2);
Fm(iGm, iGm) = -eye(2) / fp.gm_tau;
Qd = zeros(nx);
Qd(iBa, iBa) = fp.q_ba_rw^2 * Ts * eye(3);
Qd(iBg, iBg) = fp.q_bg_step * eye(3);
Qd(iGm, iGm) = fp.gm_sigma^2 * (1 - phi_gm^2) * eye(2);

% ---- ZUPT
mEvery = max(1, round(cfg.sensor.fs / fp.zupt_rate));
H = zeros(3, nx);
H(:, iV) = -eye(3);
Rm = fp.zupt_sigma^2 * eye(3);

% ---- logs
nUpd = floor(N / mEvery);
time = [tv(1), tv + Ts];
yaw = zeros(1, N + 1); hsig = yaw; roll = yaw; pitch = yaw; rollSig = yaw; pitchSig = yaw;
vLog = zeros(3, N + 1); rLog = vLog; vSig = vLog; baLog = vLog; bgLog = vLog; baSig = vLog; bgSig = vLog;
nis = nan(1, nUpd); nisT = nis;
I = eye(nx); I3 = eye(3);

logState(1);
u = 0;
for j = 1:N
    f = F(:, j) - b_a;
    w = W(:, j) - b_g;
    [att, v, r, C] = nfMechanize(att, v, r, f, w, Ts, w_ie, g_D, [gm; 0]);
    gm = gm * phi_gm;

    Fm(iA, iBg) = C;
    A = Fm * Ts;
    A2 = A * A;
    Phi = I + A + 0.5 * A2 + (A2 * A) / 6;
    Qd(iA, iA) = C * diag(qg) * C';
    Qd(iV, iV) = C * diag(qa) * C';
    P = Phi * P * Phi' + Qd;

    if mod(j, mEvery) == 0
        u = u + 1;
        z = -v;
        S = H * P * H' + Rm;
        nis(u) = z' * (S \ z);
        nisT(u) = time(j + 1);
        K = (P * H') / S;
        dx = K * z;
        IKH = I - K * H;
        P = IKH * P * IKH' + K * Rm * K';
        if quat
            att = nfQuatMultiply(nfRotVec2Quat(-dx(iA)), att);
            att = att / norm(att);
            C = nfQuat2Dcm(att);
        else
            C = nfOrthonormalize((I3 - nfSkew(dx(iA))) * C);
            att = C;
        end
        v = v - dx(iV);
        r = r - dx(iP);
        b_a = b_a + dx(iBa);
        b_g = b_g + dx(iBg);
        gm = gm + dx(iGm);
    end
    logState(j + 1);
end

hdeg = unwrap(yaw) * R2D;
res.time = time;
res.heading_deg = hdeg - 360 * floor(hdeg(end) / 360);
res.heading_sigma_deg = hsig * R2D;
res.roll_deg = roll * R2D;
res.pitch_deg = pitch * R2D;
res.roll_sigma_deg = rollSig * R2D;
res.pitch_sigma_deg = pitchSig * R2D;
res.v = vLog;
res.r = rLog;
res.v_sigma = vSig;
res.b_a = baLog;
res.b_g = bgLog;
res.b_a_sigma = baSig;
res.b_g_sigma = bgSig;
res.nis = nis(1:u);
res.nis_time = nisT(1:u);
res.P0 = P0;
res.P_end = P;
res.qg = qg;
res.qa = qa;
res.k0 = k0;
res.N = N;

    function logState(col)
        yaw(col) = atan2(C(2, 1), C(1, 1));
        hsig(col) = sqrt(max(P(3, 3), 0));
        [roll(col), pitch(col)] = nfDcm2Euler(C);
        vLog(:, col) = v;
        rLog(:, col) = r;
        baLog(:, col) = b_a;
        bgLog(:, col) = b_g;
        dP = diag(P);
        vSig(:, col) = sqrt(max(dP(iV), 0));
        baSig(:, col) = sqrt(max(dP(iBa), 0));
        bgSig(:, col) = sqrt(max(dP(iBg), 0));
        cy = cos(yaw(col)); sy = sin(yaw(col));        % NED tilt errors -> roll/pitch errors
        rollSig(col) = sqrt(max(cy^2*P(1, 1) + 2*cy*sy*P(1, 2) + sy^2*P(2, 2), 0));
        pitchSig(col) = sqrt(max(sy^2*P(1, 1) - 2*cy*sy*P(1, 2) + cy^2*P(2, 2), 0));
    end
end
