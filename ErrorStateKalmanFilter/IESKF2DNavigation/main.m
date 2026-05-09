clc, clear, close all;
% Iterated Error State Kalman Filter - 2D Robot Navigation
% IMU + GPS + Magnetometer
% Written By: Rasit
% Date: 09-May-2026
%% Trajectory Generation
Ts       = 1e-2;
TrajTime = 10;
Tstep    = linspace(0, pi, TrajTime/Ts);
xunit    = 3*cos(Tstep);
yunit    = 3*sin(Tstep);
RealTraj = [linspace(3,3,size(xunit,2))   xunit  linspace(-3,-3,size(xunit,2))  -xunit;
            linspace(-4,0,size(yunit,2))   yunit  linspace(0,-4,size(yunit,2))   -yunit-4];
NoD      = size(RealTraj, 2);

vx_true    = gradient(RealTraj(1,:), Ts);
vy_true    = gradient(RealTraj(2,:), Ts);
ax_true    = gradient(vx_true, Ts);
ay_true    = gradient(vy_true, Ts);
theta_true = atan2(vy_true, vx_true);
omega_true = gradient(theta_true, Ts);

%% Noise Parameters
sigma_ax  = 0.05;
sigma_ay  = 0.05;
sigma_w   = 0.02;
sigma_gps = 0.1;
sigma_mag = 0.05;
GPS_rate  = 10;

%% Sensor Measurements
ax_meas  = ax_true    + sigma_ax  * randn(1, NoD);
ay_meas  = ay_true    + sigma_ay  * randn(1, NoD);
w_meas   = omega_true + sigma_w   * randn(1, NoD);
GPS_meas = RealTraj   + sigma_gps * randn(2, NoD);
Mag_meas = theta_true + sigma_mag * randn(1, NoD);

%% IESKF Initialization
NofState = 5;                                                % [px, py, vx, vy, theta]
MaxIter  = 5;
iter_tol = 1e-6;

x_nom  = [RealTraj(1,1); RealTraj(2,1); vx_true(1); vy_true(1); theta_true(1)];
dx_hat = zeros(NofState, 1);

% Low  Kalman Gain: if R >> P => K~=0 (trust model)
% High Kalman Gain: if R << P => K~=1 (trust measurement)
P     = 1e0 * eye(NofState);
Q     = diag([1e-4, 1e-4, sigma_ax^2, sigma_ay^2, sigma_w^2]);
R_gps = sigma_gps^2 * eye(2);
R_mag = sigma_mag^2;
H_gps = [1 0 0 0 0; 0 1 0 0 0];
H_mag = [0 0 0 0 1];

% Innovation Gate (Chi-Squared): GPS 2DOF=5.99, Mag 1DOF=3.84
gate_gps        = 5.99;
gate_mag        = 3.84;
gate_gps_reject = 0;
gate_mag_reject = 0;

%% Storage
x_hist    = zeros(NofState, NoD);
P_hist    = zeros(NofState, NoD);
x_hist(:,1) = x_nom;
P_hist(:,1) = diag(P);

%% IESKF Loop
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for k = 1:NoD-1

    % Prediction
    x_nom(1) = x_nom(1) + x_nom(3)*Ts;
    x_nom(2) = x_nom(2) + x_nom(4)*Ts;
    x_nom(3) = x_nom(3) + ax_meas(k)*Ts;
    x_nom(4) = x_nom(4) + ay_meas(k)*Ts;
    x_nom(5) = x_nom(5) + w_meas(k)*Ts;

    Ft = [1 0 Ts 0  0;
          0 1 0  Ts 0;
          0 0 1  0  0;
          0 0 0  1  0;
          0 0 0  0  1];
    Fw = [0  0  0;
          0  0  0;
          Ts 0  0;
          0  Ts 0;
          0  0  Ts];

    P      = Ft*P*Ft' + Fw*Q(3:5,3:5)*Fw';
    x_pred = x_nom;

    % Iterated Correction: GPS
    if mod(k, GPS_rate) == 0
        z_gps   = GPS_meas(1:2, k);
        P_pred  = P;                                         % Freeze P across iterations
        x_iter  = x_nom;
        inn_gps = z_gps - H_gps*x_iter;
        S_gps   = H_gps*P_pred*H_gps' + R_gps;
        NIS_gps = inn_gps'/S_gps*inn_gps;

        if NIS_gps < gate_gps
            K_gps = P_pred*H_gps' / S_gps;                  % Compute Kalman Gain!
            for j = 1:MaxIter
                inn_j   = z_gps - H_gps*x_iter;
                x_new   = x_iter + K_gps*inn_j;
                x_new(5)= wrapToPi(x_new(5));
                if norm(x_new - x_iter) < iter_tol
                    x_iter = x_new;
                    break
                end
                x_iter = x_new;
            end
            P     = (eye(NofState)-K_gps*H_gps)*P_pred*(eye(NofState)-K_gps*H_gps)' + K_gps*R_gps*K_gps'; % Joseph Form
            x_nom = x_iter;
        else
            gate_gps_reject = gate_gps_reject + 1;
        end
    end

    % Iterated Correction: Magnetometer
    z_mag   = Mag_meas(k);
    P_pred  = P;                                             % Freeze P across iterations
    x_iter  = x_nom;
    inn_mag = wrapToPi(z_mag - H_mag*x_iter);
    S_mag   = H_mag*P_pred*H_mag' + R_mag;
    NIS_mag = inn_mag'/S_mag*inn_mag;

    if NIS_mag < gate_mag
        K_mag = P_pred*H_mag' / S_mag;                      % Compute Kalman Gain!
        for j = 1:MaxIter
            inn_j   = wrapToPi(z_mag - H_mag*x_iter);
            x_new   = x_iter + K_mag*inn_j;
            x_new(5)= wrapToPi(x_new(5));
            if norm(x_new - x_iter) < iter_tol
                x_iter = x_new;
                break
            end
            x_iter = x_new;
        end
        P     = (eye(NofState)-K_mag*H_mag)*P_pred*(eye(NofState)-K_mag*H_mag)' + K_mag*R_mag*K_mag'; % Joseph Form
        x_nom = x_iter;
    else
        gate_mag_reject = gate_mag_reject + 1;
    end

    % Reset
    x_nom(5) = wrapToPi(x_nom(5));
    dx_hat   = zeros(NofState, 1);

    x_hist(:,k+1) = x_nom;
    P_hist(:,k+1) = diag(P);

    if mod(k, 1e2) == 0 || k == NoD-1
        tVec    = (0:k-1) * Ts;
        gps_idx = 1:GPS_rate:k;
        clf

        subplot(2,2,1), grid on, hold on
        plot(RealTraj(1,:), RealTraj(2,:), 'g-', 'LineWidth', 2)
        plot(x_hist(1,1:k), x_hist(2,1:k), 'r-', 'LineWidth', 2)
        plot(GPS_meas(1,gps_idx), GPS_meas(2,gps_idx), 'b.', 'MarkerSize', 8)
        legend('True','IESKF Estimate','GPS','Location','northwest')
        xlabel('X [m]'), ylabel('Y [m]')
        title(sprintf('IESKF  |  MaxIter=%d  |  GPS Rej: %d  |  Mag Rej: %d', MaxIter, gate_gps_reject, gate_mag_reject))
        axis equal

        subplot(2,2,2)
        pos_err = sqrt((x_hist(1,1:k)-RealTraj(1,1:k)).^2 + (x_hist(2,1:k)-RealTraj(2,1:k)).^2);
        plot(tVec, pos_err, 'r', 'LineWidth', 2), grid on
        xlabel('Time [s]'), ylabel('Error [m]'), title('Position Error')

        subplot(2,2,3)
        plot(tVec, rad2deg(theta_true(1:k)),  'g-',  'LineWidth', 2), hold on, grid on
        plot(tVec, rad2deg(Mag_meas(1:k)),    'bs',  'LineWidth', 1, 'MarkerSize', 2)
        plot(tVec, rad2deg(x_hist(5,1:k)),    'r--', 'LineWidth', 2)
        legend('True','Magnetometer','IESKF','Location','northwest')
        xlabel('Time [s]'), ylabel('Heading [deg]'), title('Heading Estimation')

        subplot(2,2,4)
        plot(tVec, P_hist(1,1:k), 'LineWidth', 2), hold on, grid on
        plot(tVec, P_hist(2,1:k), 'LineWidth', 2)
        plot(tVec, P_hist(5,1:k), 'LineWidth', 2)
        legend('\sigma^2_{px}','\sigma^2_{py}','\sigma^2_{\theta}','Location','northwest')
        xlabel('Time [s]'), ylabel('Variance'), title('Covariance P')

        drawnow
    end
end