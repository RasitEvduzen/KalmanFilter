clc, clear all, close all;
% RTS Smoother - Linear Kalman Filter + Backward Pass
% Constant Velocity Model, Linear Measurement
% Written By: Rasit
% Date: 04-Apr-2026
%% Simulation Parameters
Ts       = 5e-2;                                     % Sampling Period    [s]
TrajTime = 10;                                       % Trajectory Time    [s]
SimStep  = 20;                                       % Plot every N steps
Tstep    = linspace(0, pi, TrajTime/Ts);
xunit    = 3*cos(Tstep);
yunit    = 3*sin(Tstep);
RealTraj = [linspace(3,3,size(xunit,2))   xunit  linspace(-3,-3,size(xunit,2))  -xunit;
            linspace(-4,0,size(yunit,2))   yunit  linspace(0,-4,size(yunit,2))   -yunit-4];
NoD      = size(RealTraj, 2);
tspan    = (0:NoD-1) * Ts;

vx_true   = gradient(RealTraj(1,:), Ts);
vy_true   = gradient(RealTraj(2,:), Ts);
TrueState = [RealTraj; vx_true; vy_true];            % [4 x NoD]

%% Noisy Measurement
sigma_meas = 0.3;                                    % Measurement Noise Std [m]
Z          = TrueState(1:2,:) + sigma_meas*randn(2, NoD);

%% System Matrices
NofState = 4;
F = [1 0 Ts 0;
     0 1 0  Ts;
     0 0 1  0;
     0 0 0  1];                                      % State Transition Matrix
H = [1 0 0 0;
     0 1 0 0];                                       % Observation Matrix

% Low  Kalman Gain: if R >> P => K = P/(P+R) | K ~= 0 (Algorithm belief kalman model)
% High Kalman Gain: if R << P => K = P/(P+R) | K ~= 1 (Algorithm belief measurement)
P   = 1e1 * eye(NofState);                           % Error Noise Cov Matrix
Q   = 1e-1 * [Ts^4/4 0      Ts^3/2 0;
              0      Ts^4/4 0      Ts^3/2;
              Ts^3/2 0      Ts^2   0;
              0      Ts^3/2 0      Ts^2];             % Process Noise Cov Matrix
R   = sigma_meas^2 * eye(2);                         % Measurement Noise Cov Matrix

%% Storage
X_kk  = zeros(NofState, NoD);                        % x̂_{k|k}   Forward
P_kk  = zeros(NofState, NofState, NoD);              % P_{k|k}
X_kk1 = zeros(NofState, NoD);                        % x̂_{k|k-1} Predicted
P_kk1 = zeros(NofState, NofState, NoD);              % P_{k|k-1}
X_rts = zeros(NofState, NoD);                        % x̂_{k|N}   Smoothed

%% Live Simulation
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
x_est = [Z(1,1); Z(2,1); 0; 0];
for k = 1:NoD
    %% Forward Pass — KF
    x_pred = F * x_est;
    P_pred = F * P * F' + Q;                         % Uncertainty Propagation
    X_kk1(:,k)   = x_pred;
    P_kk1(:,:,k) = P_pred;

    K     = P_pred * H' / (H * P_pred * H' + R);    % Compute Kalman Gain!
    x_est = x_pred + K * (Z(:,k) - H * x_pred);      % Update State with Measurement & Kalman Gain
    P     = (eye(NofState)-K*H)*P_pred*(eye(NofState)-K*H)' + K*R*K'; % Update Estimation Uncertainty
    X_kk(:,k)   = x_est;
    P_kk(:,:,k) = P;

    %% Backward Pass — RTS (up to current k)
    X_rts(:,k) = X_kk(:,k);
    P_rts_k    = P_kk(:,:,k);
    for j = k-1:-1:1
        G             = P_kk(:,:,j) * F' / P_kk1(:,:,j+1);
        X_rts(:,j)   = X_kk(:,j) + G * (X_rts(:,j+1) - X_kk1(:,j+1));
        P_rts_j       = P_kk(:,:,j) + G * (P_rts_k - P_kk1(:,:,j+1)) * G';
        P_rts_k       = P_rts_j;
    end

    %% Plot
    if mod(k, SimStep) == 0 || k == NoD
        rmse_kf  = sqrt(mean(sum((X_kk(1:2,1:k)  - TrueState(1:2,1:k)).^2)));
        rmse_rts = sqrt(mean(sum((X_rts(1:2,1:k) - TrueState(1:2,1:k)).^2)));
        clf

        subplot(2,2,[1 2])
        plot(TrueState(1,:),    TrueState(2,:),    'g-',  'LineWidth', 2), hold on, grid on
        plot(Z(1,1:k),          Z(2,1:k),          'b.',  'MarkerSize', 6)
        plot(X_kk(1,1:k),       X_kk(2,1:k),      'r--', 'LineWidth', 2)
        plot(X_rts(1,1:k),      X_rts(2,1:k),     'k-',  'LineWidth', 2)
        legend('True Trajectory','Noisy Measurement','KF','RTS Smoother','Location','northwest')
        xlabel('X [m]'), ylabel('Y [m]')
        title(sprintf('RTS Smoother  |  KF RMSE: %.4f m  |  RTS RMSE: %.4f m', rmse_kf, rmse_rts))
        axis equal

        subplot(2,2,3)
        plot(tspan(1:k), TrueState(3,1:k), 'g-',  'LineWidth', 2), hold on, grid on
        plot(tspan(1:k), X_kk(3,1:k),      'r--', 'LineWidth', 2)
        plot(tspan(1:k), X_rts(3,1:k),     'k-',  'LineWidth', 2)
        legend('True','KF','RTS','Location','northeast')
        xlabel('Time [s]'), ylabel('Vx [m/s]'), title('X Velocity')

        subplot(2,2,4)
        plot(tspan(1:k), TrueState(4,1:k), 'g-',  'LineWidth', 2), hold on, grid on
        plot(tspan(1:k), X_kk(4,1:k),      'r--', 'LineWidth', 2)
        plot(tspan(1:k), X_rts(4,1:k),     'k-',  'LineWidth', 2)
        legend('True','KF','RTS','Location','northeast')
        xlabel('Time [s]'), ylabel('Vy [m/s]'), title('Y Velocity')

        drawnow
    end
end