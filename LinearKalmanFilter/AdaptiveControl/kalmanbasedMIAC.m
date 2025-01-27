% Indirect Model Identification Adaptive Control (MIAC) using
% Kalman Filter First-Order SISO System
% Written By: Rasit
% Date: 27-Jan-2025
clc, clear, close all;

%% System Parameters (True and Reference Model)
a = 1; b = 2; theta = 0.2;  % True system parameters
am = -1; bm = 1;            % Reference model parameters

% Basis Functions
Phi = @(x) x * x;
Psi = @(x) [x; Phi(x)];

% Simulation Parameters
Ts = 1e-2;                 % Time step
t = 0:Ts:50;               % Simulation time
ref = zeros(length(t),1);  % Input signal
ref(t > 1) = square(t(t > 1));  % Square input after 1 sec

x_ref = zeros(length(t), 1);   % Reference model state
x = zeros(length(t), 1);       % Actual system state
u = zeros(length(t), 1);       % Control input
xdot = 0;

% Estimated Parameters (Kalman Initialization)
a_est = zeros(length(t), 1);
a_est(1) = 0;  % Initial guess
theta_est = zeros(length(t), 1);
theta_est(1) = 0;

% Kalman Filter Initialization
X_est = [a_est(1); theta_est(1)]; % State estimation [a; theta]
P = 1e2*eye(2);     % Initial covariance matrix    
Q = 5e1*[1 1; 1 1]; % Process noise covariance     
R = 1e4;            % Measurement noise covariance 

g_pdf = @(x,mu,sigma) (1/(sqrt(2*sigma^2*pi)))*exp(-((x-mu).^2)/(2*sigma^2)); % Gauss PDF
x_pdf = -2:Ts:2;

%% Simulation Loop
figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w')
for i = 2:length(t)
    if (i*Ts) > 25
        ref(i) = cos(i*Ts);
    end
    e = x_ref(i-1) - x(i-1);

    % Adaptive Control Law
    kx = (am - X_est(1)) / b;
    kr = bm / b;
    u(i-1) = kx * x(i-1) + kr * ref(i) - X_est(2) * Phi(x(i-1));

    % System Dynamics
    xdot_ref = am * x_ref(i-1) + bm * ref(i);
    xdot = a * x(i-1) + b * u(i-1) + b * theta * Phi(x(i-1));

    % Kalman Filter Prediction
    X_pred = X_est;
    P_pred = P+Q;

    % Kalman Filter Correction
    H = Psi(x(i-1))';                                  % Measurement matrix
    K = (P_pred * H') / (H * P_pred * H' + R);         % Kalman Gain
    y_meas = xdot - (X_pred(1) * x(i-1) + b * u(i-1)); % Measurement residual
    X_est = X_pred + K * (y_meas - H * X_pred);        % Update Kalman State
    P = (eye(2) - K * H) * P_pred * (eye(2) - K * H)' + K*R*K'; % Update covariance via Joseph Formula

    % System Integration (Euler Method)
    x_ref(i) = x_ref(i-1) + xdot_ref * Ts;
    x(i) = x(i-1) + xdot * Ts + .01*randn;
    a_est(i) = X_est(1);
    theta_est(i) = X_est(2);

    % Plot Result
    if mod(i,5e1) == 0 || mod(i,length(t)) == 0
        clf
        subplot(221); % State Space
        plot(t(1:i), x_ref(1:i), 'r', 'LineWidth', 2); hold on;
        plot(t(1:i), x(1:i), 'k', 'LineWidth', 1.5);
        legend('Reference', 'Adaptive Control');
        xlabel('Time (s)'); ylabel('x(t)');

        subplot(222)
        plot(g_pdf(x_pdf,mean(x),std(x)),x_pdf,"b",'LineWidth', 2);
        grid minor

        subplot(2,2,[3,4]); % Estimated Parameters
        plot(t(1:i), a_est(1:i), 'r', 'LineWidth', 2); hold on;
        plot(t(1:i), theta_est(1:i), 'g', 'LineWidth', 2);
        legend('Estimated a', 'Estimated \theta');
        xlabel('Time (s)'); ylabel('Parameter Values'); title('Adaptive Parameter Estimation');

        sgtitle('Kalman Filter Based Indirect Model Identification Adaptive Control');
        drawnow
    end
end
