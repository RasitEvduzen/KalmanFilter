openExample('control_deeplearning/TrainDDPGAgentToControlFlyingRobotExample')clc,clear all, close;
% Written By: Rasit
% Kalman Filter Based Fourier Series Decomposition
% Date: 06-Jul-2024
%%
fs = 100;             % Sampling Freq (Hz)
t_span = 0:1/fs:10;   % Time Vector (Sn)
f1 = 1;             % First Harmonics  (Hz)
f2 = 2;             % Second Harmonics (Hz)
f3 = 3;             % Third Harmonics  (Hz)
y = 1*sin(2*pi*f1*t_span) + 2*sin(2*pi*f2*t_span) + 3*sin(2*pi*f3*t_span); % Signal
y = y + 0.5*randn(size(y));  % Noise
noh = 3;  % Number of Harmonics

% Kalman filter Parameters
P = 1e3*eye(2*noh);
Q = 1e3*eye(2*noh);   
R = 0.5;  % Measurement Noise Cov
x_est = zeros(2*noh, 1);  % Kalman Initial State Value

F = eye(2*noh);       % State Transition Matrix
H = zeros(1, 2*noh);  % Measurement (State Selection) Matrix

x_estimates = zeros(length(y), 2*noh);

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for k = 1:length(y)
    for n = 1:noh
        H(1, 2*n-1) = cos(2*pi*n*t_span(k));   % Signal Even Component
        H(1, 2*n) = sin(2*pi*n*t_span(k));     % Signal Odd  Component
    end 

    % Prediction Phase
    x_pred = F * x_est;
    P_pred = F * P * F' + Q;
    K = P_pred * H' / (H * P_pred * H' + R);

    % Correction Phase 
    x_est = x_pred + K * (y(k) - H * x_pred);
    P = (eye(2*noh) - K * H) * P_pred;
    x_estimates(k, :) = x_est';
end

% Signal Reconstruction
y_est = zeros(1, length(t_span));
for k = 1:length(t_span)
    for n = 1:noh
        y_est(k) = y_est(k) + x_estimates(k, 2*n-1)*cos(2*pi*n*t_span(k)) + x_estimates(k, 2*n)*sin(2*pi*n*t_span(k));
    end
end

% Plot Phase
subplot(311);
plot(t_span, y, 'b', 'DisplayName', 'Original Signal',LineWidth=2),hold on
plot(t_span, y_est, 'r--', 'DisplayName', 'Kalman Estimated Signal',LineWidth=2);
xlabel('Time (s)');
ylabel('Magnitude');
title("Original Signal & Estimated Signal")
legend;

subplot(312);
plot(t_span, y - y_est, 'k', 'DisplayName', 'Error');
xlabel('Time (s)');
title("Error")
legend;

subplot(313)
hold on
plot(x_estimates(:,2)),plot(x_estimates(:,4)),plot(x_estimates(:,6))
yline(1),yline(2),yline(3)
xlabel('Number of Samples');
title("Harmonics")
