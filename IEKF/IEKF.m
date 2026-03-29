clc, clear all, close all;
% Radar Tracking - Iterated Extended Kalman Filter
% Constant Velocity Model + Nonlinear Radar Measurement (Cartesian to Polar)
% Written By: Rasit Evduzen
% 29-Mar-2026
%%
%% Simulation Parameters
Ts       = 1e-2;                                     % Sampling Period    [s]
TrajTime = 10;                                       % Trajectory Time    [s]
N_seg    = TrajTime/Ts;                              % Samples per segment
Tstep    = linspace(0, pi, N_seg);
omega    = pi / ((N_seg-1)*Ts);                      % Angular rate [rad/s]
xunit    = 3*cos(Tstep);
yunit    = 3*sin(Tstep);

RealTraj = [linspace(3,3,N_seg)    xunit  linspace(-3,-3,N_seg)  -xunit;
            linspace(-4,0,N_seg)   yunit  linspace(0,-4,N_seg)   -yunit-4];
NoD      = size(RealTraj, 2);
tspan    = (0:NoD-1) * Ts;                           % Time Vector        [s]

% Analytic Velocities (per segment)
vx_s1 = zeros(1, N_seg);
vy_s1 = ( 4 / ((N_seg-1)*Ts)) * ones(1, N_seg);     % y: -4 -> 0

vx_s2 = -3*sin(Tstep) * omega;
vy_s2 =  3*cos(Tstep) * omega;

vx_s3 = zeros(1, N_seg);
vy_s3 = (-4 / ((N_seg-1)*Ts)) * ones(1, N_seg);     % y: 0 -> -4

vx_s4 =  3*sin(Tstep) * omega;
vy_s4 = -3*cos(Tstep) * omega;

vx_true   = [vx_s1, vx_s2, vx_s3, vx_s4];
vy_true   = [vy_s1, vy_s2, vy_s3, vy_s4];
TrueState = [RealTraj; vx_true; vy_true];            % [px;py;vx;vy] x NoD


%% CV Motion Model (Linear)
F = [1 0 Ts 0;
    0 1 0  Ts;
    0 0 1  0;
    0 0 0  1];                                      % State Transition Matrix
NofState = 4;

%% Radar Measurement (Nonlinear: Cartesian to Polar)
sigma_r   = 0.1;                                     % Range Noise Std     [m]
sigma_phi = 0.01;                                    % Bearing Noise Std   [rad]
R_meas    = diag([sigma_r^2, sigma_phi^2]);          % Measurement Noise Cov Matrix
h_func    = @(x) [sqrt(x(1)^2 + x(2)^2); atan2(x(2), x(1))];

RadarMeas = zeros(2, NoD);
for k = 1:NoD
    RadarMeas(:,k) = h_func(TrueState(:,k)) + [sigma_r*randn; sigma_phi*randn];
end

%% EKF & IEKF Initialization
% Low  Kalman Gain: if R >> P => K~=0 (trust model)
% High Kalman Gain: if R << P => K~=1 (trust measurement)
sigma_q = 1e-2;                                      % Process Noise Std   [m/s^2]
Q  = sigma_q^2 * diag([Ts^4/4, Ts^4/4, Ts^2, Ts^2]); % Process Noise Cov Matrix
P0 = 1e-3 * diag([sigma_r^2, sigma_r^2, 5, 5]);     % Error Noise Cov Matrix
r0 = RadarMeas(1,1);  phi0 = RadarMeas(2,1);
x0 = [r0*cos(phi0); r0*sin(phi0); vx_true(1); vy_true(1)];

P_ekf  = P0;  P_iekf = P0;
XKalman_ekf  = zeros(NofState, NoD);
XKalman_iekf = zeros(NofState, NoD);
XKalman_ekf(:,1)  = x0;
XKalman_iekf(:,1) = x0;

max_iter = 10;                                       % Max Iterations
epsilon  = 1e-6;                                     % Convergence Threshold

%% EKF & IEKF Loop
for k = 1:NoD-1
    z = RadarMeas(:,k+1);

    % EKF Prediction
    x_pred = F * XKalman_ekf(:,k);
    P_ekf  = F*P_ekf*F' + Q;
    % EKF Correction
    H      = calcH(x_pred);
    nu     = z - h_func(x_pred);  nu(2) = wrapToPi(nu(2));
    K      = P_ekf*H' / (H*P_ekf*H' + R_meas);      % Compute Kalman Gain!
    XKalman_ekf(:,k+1) = x_pred + K*nu;              % Update State with Measurement & Kalman Gain
    P_ekf  = (eye(NofState)-K*H)*P_ekf*(eye(NofState)-K*H)' + K*R_meas*K'; % Joseph Form

    % IEKF Prediction
    x_pred  = F * XKalman_iekf(:,k);
    P_iekf  = F*P_iekf*F' + Q;
    % IEKF Correction (Iterated)
    P_pred_iekf = P_iekf;                            % Freeze P_pred for all iterations
    x_j = x_pred;
    K_j = zeros(NofState,2);  H_j = zeros(2,NofState);
    for j = 1:max_iter
        H_j   = calcH(x_j);
        K_j   = P_pred_iekf*H_j' / (H_j*P_pred_iekf*H_j' + R_meas); % Compute Kalman Gain!
        dnu_j = z - h_func(x_j) - H_j*(x_pred - x_j); % IEKF Innovation (Survey eq. 69)
        dnu_j(2) = wrapToPi(dnu_j(2));
        dx_j  = K_j * dnu_j;
        x_j1  = x_j + dx_j;
        if norm(dx_j) < epsilon, break; end           % Convergence Check
        x_j = x_j1;
    end
    XKalman_iekf(:,k+1) = x_j1;                      % Update State with Measurement & Kalman Gain
    P_iekf = (eye(NofState)-K_j*H_j)*P_pred_iekf*(eye(NofState)-K_j*H_j)' + K_j*R_meas*K_j'; % Joseph Form
end

%% RMSE
rmse_ekf  = sqrt(mean(sum((XKalman_ekf(1:2,:)  - TrueState(1:2,:)).^2)));
rmse_iekf = sqrt(mean(sum((XKalman_iekf(1:2,:) - TrueState(1:2,:)).^2)));

RadarCart = [RadarMeas(1,:).*cos(RadarMeas(2,:));
    RadarMeas(1,:).*sin(RadarMeas(2,:))];

%% Plot
plot_interval = 1e2;
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for i = 1: NoD
    clf
    if mod(i, plot_interval) == 0 || i == NoD

        subplot(3,2,[1 2])
        plot(TrueState(1,:),    TrueState(2,:),    'g',  'LineWidth',2.5), hold on, grid on
        plot(RadarCart(1,1:i),    RadarCart(2,1:i),    'bo', 'MarkerSize',3)
        plot(XKalman_ekf(1,1:i),  XKalman_ekf(2,1:i), 'r',  'LineWidth',2)
        plot(XKalman_iekf(1,1:i), XKalman_iekf(2,1:i),'k',  'LineWidth',2)
        legend('True Trajectory','Radar Measurement','EKF','IEKF','Location','northwest')
        xlabel('X [m]'), ylabel('Y [m]')
        title(sprintf('Radar Tracking  |  EKF RMSE: %.4f m  |  IEKF RMSE: %.4f m', rmse_ekf, rmse_iekf))
        axis equal

        subplot(3,2,3)
        plot(tspan, TrueState(1,:),    'g',  'LineWidth',2), hold on, grid on
        plot(tspan(1:i), XKalman_ekf(1,1:i),  'r',  'LineWidth',2)
        plot(tspan(1:i), XKalman_iekf(1,1:i), 'k',  'LineWidth',2)
        legend('True','EKF','IEKF','Location','northeast')
        xlabel('Time [s]'), ylabel('X [m]'), title('X Position Comparison')

        subplot(3,2,4)
        plot(tspan, TrueState(2,:),    'g',  'LineWidth',2), hold on, grid on
        plot(tspan(1:i), XKalman_ekf(2,1:i),  'r',  'LineWidth',2)
        plot(tspan(1:i), XKalman_iekf(2,1:i), 'k',  'LineWidth',2)
        legend('True','EKF','IEKF','Location','northeast')
        xlabel('Time [s]'), ylabel('Y [m]'), title('Y Position Comparison')

        subplot(3,2,5)
        plot(tspan, vx_true,             'g',  'LineWidth',2), hold on, grid on
        plot(tspan(1:i), XKalman_ekf(3,1:i),   'r',  'LineWidth',2)
        plot(tspan(1:i), XKalman_iekf(3,1:i),  'k',  'LineWidth',2)
        legend('True','EKF','IEKF','Location','northeast')
        xlabel('Time [s]'), ylabel('Vx [m/s]'), title('Vx Comparison')

        subplot(3,2,6)
        plot(tspan, vy_true,             'g',  'LineWidth',2), hold on, grid on
        plot(tspan(1:i), XKalman_ekf(4,1:i),   'r',  'LineWidth',2)
        plot(tspan(1:i), XKalman_iekf(4,1:i),  'k',  'LineWidth',2)
        legend('True','EKF','IEKF','Location','northeast')
        xlabel('Time [s]'), ylabel('Vy [m/s]'), title('Vy Comparison')
        drawnow
    end
end
%% Utility Functions
function H = calcH(x)
% Jacobian of h(x) = [sqrt(px^2+py^2); atan2(py,px)]
r2 = max(x(1)^2 + x(2)^2, 1e-6);                   % Singularity guard
r  = sqrt(r2);
H  = [ x(1)/r,   x(2)/r,  0, 0;
    -x(2)/r2,  x(1)/r2, 0, 0];
end