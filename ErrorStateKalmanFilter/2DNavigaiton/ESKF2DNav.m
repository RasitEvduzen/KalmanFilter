clc, clear, close all;
% Error State Kalman Filter 2D Robot Navigation
% IMU (accel + gyro) + GPS + Magnetometer
% Written By: Rasit
% Date: 27-Mar-2026
%% Trajectory Generation (ConstantAccelModel Reference)
Ts       = 1e-2;
TrajTime = 10;
Tstep    = linspace(0, pi, TrajTime/Ts);
xunit    = 3*cos(Tstep);
yunit    = 3*sin(Tstep);
RealTraj = [linspace(3,3,size(xunit,2))   xunit  linspace(-3,-3,size(xunit,2))  -xunit;
            linspace(-4,0,size(yunit,2))   yunit  linspace(0,-4,size(yunit,2))   -yunit-4];
NoD      = size(RealTraj, 2);

% True Heading, Velocity, Acceleration 
vx_true    = gradient(RealTraj(1,:), Ts);
vy_true    = gradient(RealTraj(2,:), Ts);
ax_true    = gradient(vx_true, Ts);
ay_true    = gradient(vy_true, Ts);
theta_true = atan2(vy_true, vx_true);
omega_true = gradient(theta_true, Ts);

%% Noise Parameters
sigma_ax  = 0.05;    % Accel Noise Std   [m/s^2]
sigma_ay  = 0.05;    % Accel Noise Std   [m/s^2]
sigma_w   = 0.02;    % Gyro  Noise Std   [rad/s]
sigma_gps = 0.1;     % GPS   Noise Std   [m]
sigma_mag = 0.05;    % Mag   Noise Std   [rad]
GPS_rate  = 10;      % GPS every N steps

%% Sensor Measurements 
ax_meas    = ax_true    + sigma_ax  * randn(1, NoD);
ay_meas    = ay_true    + sigma_ay  * randn(1, NoD);
w_meas     = omega_true + sigma_w   * randn(1, NoD);
GPS_meas   = RealTraj   + sigma_gps * randn(2, NoD);
Mag_meas   = theta_true + sigma_mag * randn(1, NoD);

%% ESKF Initialization
NofState   = 5;                                              % [px, py, vx, vy, theta]
x_nom      = [RealTraj(1,1); RealTraj(2,1); vx_true(1); vy_true(1); theta_true(1)];
dx_hat     = zeros(NofState, 1);                             % Error State Always 0

% ---------------------- P~Q~R Matrix ----------------------
% Low  Kalman Gain: if R >> P => K = P/(P+R) | K ~= 0 (Algorithm belief kalman model)
% High Kalman Gain: if R << P => K = P/(P+R) | K ~= 1 (Algorithm belief measurement)
P          = 1e0 * eye(NofState);                            % Error Noise Cov Matrix
Q          = diag([1e-4, 1e-4, sigma_ax^2, sigma_ay^2, sigma_w^2]); % Processes Noise Cov Matrix
R_gps      = sigma_gps^2 * eye(2);                          % GPS Measurement Noise Cov
R_mag      = sigma_mag^2;                                    % Magnetometer Noise Cov
H_gps      = [1 0 0 0 0; 0 1 0 0 0];                        % GPS Observation Matrix
H_mag      = [0 0 0 0 1];                                    % Mag Observation Matrix
% -----------------------------------------------------------
% Innovation Gate (Chi-Squared Consistency Check)
% NIS = innovation' * S^-1 * innovation
% GPS (2 DOF): chi2inv(0.95, 2) = 5.99
% Mag (1 DOF): chi2inv(0.95, 1) = 3.84
gate_gps        = 5.99;                                      % GPS Innovation Gate Threshold
gate_mag        = 3.84;                                      % Mag Innovation Gate Threshold
gate_gps_reject = 0;                                         % GPS Rejection Counter
gate_mag_reject = 0;                                         % Mag Rejection Counter

%% Storage
x_hist     = zeros(NofState, NoD);
P_hist     = zeros(NofState, NoD);
x_hist(:,1)   = x_nom;
P_hist(:,1)   = diag(P);

%% ESKF 
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for k = 1:NoD-1
    % PREDICTION (Nominal State Propagation)
    x_nom(1) = x_nom(1) + x_nom(3)*Ts;                      % px
    x_nom(2) = x_nom(2) + x_nom(4)*Ts;                      % py
    x_nom(3) = x_nom(3) + ax_meas(k)*Ts;                    % vx
    x_nom(4) = x_nom(4) + ay_meas(k)*Ts;                    % vy
    x_nom(5) = x_nom(5) + w_meas(k)*Ts;                     % theta

    % Error State Jacobian 
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

    % Uncertainty Propogation
    P  = Ft*P*Ft' + Fw*Q(3:5,3:5)*Fw';

    % CORRECTION: GPS
    if mod(k, GPS_rate) == 0
        z_gps   = GPS_meas(1:2, k);
        inn_gps = z_gps - H_gps*x_nom;                      % Innovation (residual)
        S_gps   = H_gps*P*H_gps' + R_gps;
        % Innovation Gate Check (Chi-Squared, 2 DOF)
        NIS_gps = inn_gps' / S_gps * inn_gps;
        if NIS_gps < gate_gps
            K_gps   = P*H_gps' / S_gps;                     % Compute Kalman Gain!
            dx_hat  = dx_hat + K_gps * inn_gps;              % Update State with Measurement & Kalman Gain
            P       = (eye(NofState)-K_gps*H_gps)*P*(eye(NofState)-K_gps*H_gps)' + K_gps*R_gps*K_gps'; % Update Estimation Uncertainty
        else
            gate_gps_reject = gate_gps_reject + 1;           % Outlier Rejected!
        end
    end

    % CORRECTION: Magnetometer
    z_mag   = Mag_meas(k);
    inn_mag = wrapToPi(z_mag - H_mag*x_nom);                % Innovation (residual)
    S_mag   = H_mag*P*H_mag' + R_mag;
    % Innovation Gate Check (Chi-Squared, 1 DOF)
    NIS_mag = inn_mag' / S_mag * inn_mag;
    if NIS_mag < gate_mag
        K_mag   = P*H_mag' / S_mag;                          % Compute Kalman Gain!
        dx_hat  = dx_hat + K_mag * inn_mag;                   % Update State with Measurement & Kalman Gain
        P       = (eye(NofState)-K_mag*H_mag)*P*(eye(NofState)-K_mag*H_mag)' + K_mag*R_mag*K_mag'; % Update Estimation Uncertainty
    else
        gate_mag_reject = gate_mag_reject + 1;                % Outlier Rejected!
    end

    % RESET (Inject Error State → Nominal)
    x_nom(5) = wrapToPi(x_nom(5));
    x_nom    = x_nom + dx_hat;                               % True State = Nominal + Error
    x_nom(5) = wrapToPi(x_nom(5));
    dx_hat   = zeros(NofState, 1);                           % Reset Error State to Zero!

    % Store History
    x_hist(:,k+1) = x_nom;
    P_hist(:,k+1) = diag(P);

    % Plot
    if mod(k, 1e2) == 0 || k == NoD-1
        tVec     = (0:k-1) * Ts;                             % Time Vector [s]
        gps_idx  = 1:GPS_rate:k;                             % GPS sample indices
        clf

        subplot(2,2,1), grid on, hold on
        plot(RealTraj(1,:), RealTraj(2,:), 'g-', 'LineWidth', 2)
        plot(x_hist(1,1:k), x_hist(2,1:k), 'r-', 'LineWidth', 2)
        plot(GPS_meas(1,gps_idx), GPS_meas(2,gps_idx), 'b.', 'MarkerSize', 8)
        legend('True Trajectory', 'ESKF Estimate', 'GPS Measurement', 'Location', 'northwest')
        xlabel('X [m]'), ylabel('Y [m]')
        title(sprintf('2D Navigation  |  GPS Rejected: %d  |  Mag Rejected: %d', gate_gps_reject, gate_mag_reject))
        axis equal

        subplot(2,2,2)
        pos_err = sqrt((x_hist(1,1:k)-RealTraj(1,1:k)).^2 + (x_hist(2,1:k)-RealTraj(2,1:k)).^2);
        plot(tVec, pos_err, 'r', 'LineWidth', 2), grid on
        xlabel('Time [s]'), ylabel('Error [m]'), title('Position Error')

        subplot(2,2,3)
        plot(tVec, rad2deg(theta_true(1:k)),  'g-',  'LineWidth', 2), hold on, grid on
        plot(tVec, rad2deg(Mag_meas(1:k)),    'bs',  'LineWidth', 1, 'MarkerSize', 2)
        plot(tVec, rad2deg(x_hist(5,1:k)),    'r--', 'LineWidth', 2)
        legend('True', 'Magnetometer', 'ESKF', 'Location', 'northwest')
        xlabel('Time [s]'), ylabel('Heading [deg]'), title('Heading Estimation')

        subplot(2,2,4)
        plot(tVec, P_hist(1,1:k), 'LineWidth', 2), hold on, grid on
        plot(tVec, P_hist(2,1:k), 'LineWidth', 2)
        plot(tVec, P_hist(5,1:k), 'LineWidth', 2)
        legend('\sigma^2_{px}', '\sigma^2_{py}', '\sigma^2_{\theta}', 'Location', 'northwest')
        xlabel('Time [s]'), ylabel('Variance'), title('Covariance P')

        drawnow
    end
end

