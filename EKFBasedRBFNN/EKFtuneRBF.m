clc; clear; close all;
% RBF NN tune EKF
% Written By: Rasit 
% Date: 10-Nov-2025
%%
% Initial parameters
N = 4;                              % Number of RBF neurons
w = 0.1 * rand(N,1);
c = linspace(-pi, pi, N)' + 0.1 * rand(N,1);
sigma = 0.5 * ones(N,1) + 0.1 * rand(N,1);
theta = [w; c; sigma];              % State vector


P = eye(3*N);                       % Initial covariance matrix
Q = blkdiag(1e-4*eye(N), 1e-6*eye(N), 1e-6*eye(N));  % Process noise covariance
R = 1e-2;                           % Measurement noise covariance

%% Training Data Generation
NoD = 5e2;                          % Number of data 
x_data = linspace(-pi, pi, NoD)';
y_true = sin(x_data);               % Nonlinear target function 
y_meas = y_true + 0.05 * randn(NoD,1);

%% Utility Functions
phi_func = @(x, c, sigma) exp(-((x - c').^2) ./ (2*sigma'.^2));
rbf_func = @(x, ci, si, wi) wi * exp(-((x - ci).^2) ./ (2*si^2));
x_vis = linspace(-pi, pi, 200)';    % Visualization grid

% Parameter History Storage 
w_hist = zeros(NoD, N);
c_hist = zeros(NoD, N);
sigma_hist = zeros(NoD, N);

%% Online EKF Learning and Visualization 
figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w')
for k = 1:NoD
    w = theta(1:N);
    c = theta(N+1:2*N);
    sigma = abs(theta(2*N+1:end));        % sigma > 0 
    phi = phi_func(x_data(k), c, sigma)'; % Nx1 feature vector
    y_est = w' * phi;                     % Estimated output

    % Jacobian (1x3N)
    dYdW = phi';
    dYdC = (w' .* phi') .* ((x_data(k) - c') ./ (sigma'.^2));
    dYdS = (w' .* phi') .* ((x_data(k) - c').^2 ./ (sigma'.^3));
    H = [dYdW, dYdC, dYdS];               % Measurement Jacobian (1x3N)

    % Prediction step
    theta_pred = theta;                   % State prediction (F = I)
    P_pred = P + Q;                       % Covariance prediction

    % Kalman gain
    S = H * P_pred * H' + R;
    K = (P_pred * H') / S;

    % Correction step
    e = y_meas(k) - y_est;                % Innovation (residual)
    theta = theta_pred + K * e;           % Updated state
    P = (eye(3*N) - K*H) * P_pred *(eye(3*N) - K*H)' + K*R*K';        % Updated covariance Joseph Formula

    w_hist(k,:) = theta(1:N);
    c_hist(k,:) = theta(N+1:2*N);
    sigma_hist(k,:) = abs(theta(2*N+1:end));

    y_pred_grid = zeros(size(x_vis));
    for i = 1:N
        y_pred_grid = y_pred_grid + w(i) * exp(-((x_vis - c(i)).^2) / (2*sigma(i)^2));
    end

    % Visualization 
    if mod(k,1e1) == 0 || mod(k,NoD) == 0
        clf
        subplot(2,2,1), hold on
        for i = 1:N
            y_neuron = rbf_func(x_vis, c(i), sigma(i), w(i));
            plot(x_vis, y_neuron, 'Color', [.5 .5 .5]);
            plot(c(i), 0, 'ko', 'MarkerFaceColor', 'r'); % neuron center
        end

        plot(x_vis, y_pred_grid, 'b', 'LineWidth', 2);
        plot(x_data, y_meas, 'r.');
        plot(x_data, y_true, 'k--', 'LineWidth', 1.5);
        title(sprintf('Step %d / %d', k, NoD));
        xlabel('x'); ylabel('y');
        xlim([-pi pi]); ylim([-1.5 1.5]);

        subplot(2,2,2)
        plot(1:k, w_hist(1:k,:), 'LineWidth', 1.2); hold on;
        title('Estimated Parameters (Weights)');
        xlabel('Number of Data Points');
        ylabel('Parameter Values');

        subplot(2,2,3)
        plot(1:k, c_hist(1:k,:), 'LineWidth', 1.2);
        title('Estimated Parameters (Centers)');
        xlabel('Number of Data Points');
        ylabel('Parameter Values');

        subplot(2,2,4)
        plot(1:k, sigma_hist(1:k,:), 'LineWidth', 1.2);
        title('Estimated Parameters (\sigma)');
        xlabel('Number of Data Points');
        ylabel('Parameter Values');

        drawnow;
        
    end
end
