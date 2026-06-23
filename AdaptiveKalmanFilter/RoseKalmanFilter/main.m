%% Adaptive Kalman Filter Estimation vs Standard Kalman Filter (fixed Q/R) for comparison
% Written By: Rasit
% 30-Jul-2024
clc, clear all, close all;
%% System and Filter Hyperparameters
Gamma   = 8.0;             % Measurement noise scaling factor
AlphaR  = 0.5;             % Smoothing factor for measurement noise variance
AlphaM  = 0.4;             % Smoothing factor for innovation variance M
Ts = 0.1;                  % Sampling time
A  = 1;                    % System matrix
B  = 0;                    % Input matrix
C  = 1;                    % Output matrix
D  = 0;                    % Direct transmission matrix

%% Standard Kalman Filter Parameters (fixed)
Q_std = 0.01;              % Fixed process noise covariance
R_std = 0.05;              % Fixed measurement noise covariance

%% Measurement Data Generator
t      = Ts*(1:1:1200);
x_true = [zeros(1,100) +.5*ones(1,200) zeros(1,150) +1*ones(1,150) ...
          zeros(1,200) -.3*ones(1,200) zeros(1,200)];
v  = [sqrt(1*5e-3)*randn(1,250) sqrt(10*5e-3)*randn(1,350) sqrt(1*5e-3)*randn(1,600)];
y  = x_true + v;
u  = zeros(1, length(y));

%% Init Filter (Adaptive)
x_hat   = y(1);
p_tilde = 0;
E1      = y(1);
EE1     = y(:,1)*y(:,1)';
M       = 0;

%% Init Standard Kalman
x_hat_std   = y(1);        % Initial state estimate (standard)
p_tilde_std = 0.1;         % Initial error covariance (standard)

%% Arrays for Live Plotting
x_hat_plot     = zeros(1, length(y));
x_hat_std_plot = zeros(1, length(y));
Q_plot         = zeros(1, length(y));
R_plot         = zeros(1, length(y));

%% Figure and Graphic Handles
figure('units', 'normalized', 'Position', [0.1 0.1 0.8 0.8], 'color', 'w');
subplot(2,2,[1, 2]);
h_y    = plot(NaN, NaN, 'g'); hold on; grid on;
h_true = plot(t, x_true, 'k--', 'LineWidth', 2);
h_std  = plot(NaN, NaN, 'b', 'LineWidth', 1.2);
h_hat  = plot(NaN, NaN, 'r', 'LineWidth', 1.5);
title('Adaptive Filter vs Standard Kalman Filter')
xlabel('Time (s)'); ylabel('State / Measurement')
legend('Measurement (y)', 'True State', 'Standard KF', 'Adaptive Filter', 'Location', 'best')
xlim([0 t(end)]);
ylim([min(y)-0.2, max(y)+0.2]);

subplot(2,2,3);
h_Q = plot(NaN, NaN, 'b', 'LineWidth', 1.5); grid on;
title('Process Noise Covariance (Q) - Adaptive Kalman')
xlabel('Time (s)'); ylabel('Covariance')
xlim([0 t(end)]);
ylim([-0.01, 0.3]);

subplot(2,2,4);
h_R = plot(NaN, NaN, 'm', 'LineWidth', 1.5); grid on;
title('Measurement Noise Covariance (R) - Adaptive Kalman')
xlabel('Time (s)'); ylabel('Covariance')
xlim([0 t(end)]);
ylim([-0.01, 3.5]);

%% Filtering Loop
for k = 1:length(y)
    % Filter (Adaptive)
    % -- Determine R using a 1st-order IIR filter --
    E1  = AlphaR*y(:,k)         + (1-AlphaR)*E1;
    EE1 = AlphaR*y(:,k)*y(:,k)' + (1-AlphaR)*EE1;
    R   = Gamma*(EE1 - E1*E1');

    % -- Determine M using a 1st-order IIR filter --
    dy = y(:,k) - C*x_hat - D*u(k);
    M  = AlphaM.*dy*dy' + (1-AlphaM).*M;

    % -- Determine Q --
    Q(k) = M - R - p_tilde;
    if Q(k) < 0
        Q(k) = 0;
    end

    % -- Kalman Equations --
    p_hat   = A*p_tilde*A' + Q(k);
    K       = p_hat*C'*pinv(C*p_hat*C' + R);
    x_tilde = x_hat + K*dy;
    p_tilde = (eye(length(B)) - K*C)*p_hat;
    x_hat   = A*x_tilde + B*u(k);

    % Standard Kalman Filter (fixed Q/R)
    % -- Prediction --
    x_pred_std = A*x_hat_std + B*u(k);
    p_hat_std  = A*p_tilde_std*A' + Q_std;

    % -- Innovation --
    dy_std = y(k) - C*x_pred_std - D*u(k);

    % -- Update --
    S_std       = C*p_hat_std*C' + R_std;
    K_std       = p_hat_std*C'/S_std;
    x_hat_std   = x_pred_std + K_std*dy_std;
    p_tilde_std = (1 - K_std*C)*p_hat_std;

    % Save data for plotting
    x_hat_plot(k)     = x_hat;
    x_hat_std_plot(k) = x_hat_std;
    Q_plot(k) = Q(k);
    R_plot(k) = R;

    % -- Plot Result --
    if mod(k, 5) == 0
        set(h_y,    'XData', t(1:k), 'YData', y(1:k));
        set(h_std,  'XData', t(1:k), 'YData', x_hat_std_plot(1:k));
        set(h_hat,  'XData', t(1:k), 'YData', x_hat_plot(1:k));
        set(h_Q,    'XData', t(1:k), 'YData', Q_plot(1:k));
        set(h_R,    'XData', t(1:k), 'YData', R_plot(1:k));
        drawnow
    end
end