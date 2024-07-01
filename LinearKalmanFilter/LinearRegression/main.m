clc,clear all,close all;
% Written By: Rasit Evduzen
% 24-May-2023
% Kalman Filter Linear Regression
%% Create Data
NoD = 1e2;
X = linspace(0, 10, NoD);  % Bağımsız değişken (x)
true_m = 2.5;  % Gerçek eğim
true_b = 1.0;  % Gerçek kesişim
Y = true_m * X + true_b + normrnd(0, 0.5, size(X));  % Gürültülü çıkış (y)

% Kalman filtresi parametreleri
theta = [0.0; 0.0];  % Initial Parameters [m; b]
P = eye(2);          % Estimation variance Cov
Q = eye(2) * 1e-5;   % Processes Noise Cov
R = 0.25;            % Measurement Noise Cov

% Kalman Estimated Parameters
theta_estimates = zeros(2, NoD);
kalman_gain = [];

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for i = 1:NoD
    clf
    x = X(i);
    y = Y(i);

    % Prediction Phase
    theta = theta;  % F = I
    P = P + Q;      % F = I

    % Correction Phase
    H = [x, 1.0];
    K = P * H' / (H * P * H' + R);
    theta = theta + K * (y - H * theta);
    P = (eye(2) - K * H) * P;
    
    kalman_gain(:,i) = K;
    theta_estimates(:, i) = theta;
    
    subplot(121)
    scatter(X, Y, 'b', 'filled', 'DisplayName', 'Data',LineWidth=10), hold on,grid on
    plot(X, true_m * X + true_b, 'k', 'DisplayName', 'Real Curve',LineWidth=3);
    plot(X, theta_estimates(1, i) * X + theta_estimates(2, i), 'r--', 'DisplayName', 'Predicted Curve',LineWidth=2);
    legend;xlabel('X');ylabel('Y');title('Kalman Filter Based Polynomial Regression');

    subplot(122)
    plot(kalman_gain'),grid
    title("Kalman Gain")
    drawnow
end

