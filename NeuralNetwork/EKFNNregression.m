clc, clear, close all;
% EKF Tuned Neural Network - Nonlinear Regression
% Written By: Rasit
% Date: 27-Mar-2026
%%
%% Training Data Generation
NoD       = 1e3;                                     % Number Of Data
sigma_n   = 0.15;                                     % Measurement Noise Std
x_data    = linspace(-3, 3, NoD)';                   % Input  [NoD x 1]
y_true    = 0.01*x_data.*x_data + 0.1*exp(x_data) + sin(2*x_data); % Target Function
y_meas    = y_true + sigma_n*randn(NoD, 1);          % Noisy Measurement

% Min-Max Normalization → output ∈ [0, 1]
y_min    = min(y_meas);
y_max    = max(y_meas);
y_meas_n = (y_meas - y_min) / (y_max - y_min);      % Normalized Measurement

% Shuffle  uniform sampling across input range
rng(42);
shuf_idx  = randperm(NoD);
x_shuf    = x_data(shuf_idx);
y_shuf    = y_meas_n(shuf_idx);


%% Neural Network Architecture
R          = 1;                                      % Input  Dimension
S          = 75;                                    % Number of Hidden Neurons
act_func   = 'tanh';                                 % Activation: 'tanh' | 'relu' | 'sigmoid' | 'elu'
% -----------------------------------------------------------
% State vector: theta = [Wg(SxR); bh(Sx1); Wc(Sx1); bc(1x1)] -> P_size x 1
% -----------------------------------------------------------
P_size = S*(R+2) + 1;                                % Number Of State

%% NN Weight Initialization
Wg    = 8*rand(S, R)-4;             % Hidden Layer Weights  [S x R]
bh    = 8*rand(S, 1)-4;             % Hidden Layer Biases   [S x 1]
Wc    = 8-rand(1, S)-4;                           % Output Layer Weights  [1 x S]
bc    = 8*rand(1, 1)-4;                           % Output Layer Bias     [1 x 1]
theta = [Wg(:); bh; Wc(:); bc];                      % Full State Vector     [P_size x 1]

% ---------------------- P~Q~R Matrix ----------------------
% Low  Kalman Gain: if R >> P => K = P/(P+R) | K ~= 0 (Algorithm belief kalman model)
% High Kalman Gain: if R << P => K = P/(P+R) | K ~= 1 (Algorithm belief measurement)
P       = 1e-2 * eye(P_size);                         % Error Noise Cov Matrix
Q       = 1e-9 * eye(P_size);                        % Processes Noise Cov Matrix
R_noise = (sigma_n / (y_max - y_min))^2;             % Measurement Noise Cov (min-max normalized space)

%% Visualization Grid (consistent with x_data range)
x_vis  = linspace(min(x_data), max(x_data), NoD)';  % Visualization Grid
y_vis  = zeros(NoD, 1);

%% Parameter History Storage
theta_hist = zeros(NoD, P_size);
plot_interval = 10;
%% EKF Learning
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for k = 1:NoD

    x_k = x_shuf(k);                                 % Current Input  (shuffled)
    z_k = y_shuf(k);                                 % Current Measurement (normalized)

    % Forward Pass + Jacobian
    [y_est, H] = nn_forward(theta, x_k, S, R, act_func);
    
    % Prediction
    theta_pred = theta;                               
    P_pred     = P + Q;                              

    % Correction 
    S_inn = H * P_pred * H' + R_noise;
    K     = (P_pred * H') / S_inn;                   % Compute Kalman Gain!
    e     = z_k - y_est;                             % Innovation (residual)
    theta = theta_pred + K * e;                       % Update State with Measurement & Kalman Gain
    P     = (eye(P_size)-K*H)*P_pred*(eye(P_size)-K*H)' + K*R_noise*K'; % Update Estimation Uncertainty

    % Store History
    theta_hist(k,:) = theta';

    % Plot Regression
    if mod(k, plot_interval) == 0 || k == NoD
        for j = 1:1e3
            y_vis(j) = nn_forward(theta, x_vis(j), S, R, act_func) * (y_max - y_min) + y_min; % Denormalize
        end
        y_vis_full = zeros(NoD, 1);
        for j = 1:NoD
            y_vis_full(j) = nn_forward(theta, x_data(j), S, R, act_func) * (y_max - y_min) + y_min;
        end
        RMSE = sqrt(mean((y_vis_full - y_true).^2));

        clf
        plot(x_data, y_true, 'b', 'LineWidth', 2), hold on, grid on
        scatter(x_data, y_meas, 'k')
        plot(x_vis,  y_vis,  'r-', 'LineWidth', 2)
        xlabel('x'), ylabel('y')
        legend('True Function', 'Noisy Measurement', 'EKF-NN Estimate', 'Location', 'northwest')
        title(sprintf('EKF Tuned Neural Network  (Nonlinear Regression)  |  Step %d / %d  |  RMSE: %.4f', k, NoD, RMSE))
        axis([min(x_data) max(x_data) min(y_meas)-0.5 max(y_meas)+0.5])
        drawnow
    end
end

%% Parameter History
figure('units','normalized','outerposition',[0 0 1 1],'color','w')

subplot(2,2,1)
plot(1:NoD, theta_hist(:, 1:S), 'LineWidth', 1.2), grid on
title('Hidden Layer Weights (Wg)')
xlabel('Number of Data Points'), ylabel('Parameter Values')

subplot(2,2,2)
plot(1:NoD, theta_hist(:, S*R+1:S*R+S), 'LineWidth', 1.2), grid on
title('Hidden Layer Biases (bh)')
xlabel('Number of Data Points'), ylabel('Parameter Values')

subplot(2,2,3)
plot(1:NoD, theta_hist(:, S*R+S+1:S*R+2*S), 'LineWidth', 1.2), grid on
title('Output Layer Weights (Wc)')
xlabel('Number of Data Points'), ylabel('Parameter Values')

subplot(2,2,4)
plot(1:NoD, theta_hist(:, end), 'k', 'LineWidth', 2), grid on
title('Output Bias (bc)')
xlabel('Number of Data Points'), ylabel('Parameter Values')

%% Utility Functions
function [y_est, H] = nn_forward(theta, x_k, S, R, act_func)
% nn_forward - NN Forward Pass + Analytical Jacobian
Wg_cur = reshape(theta(1:S*R), S, R);                % Hidden Weights  [S x R]
bh_cur = theta(S*R+1:S*R+S);                         % Hidden Biases   [S x 1]
Wc_cur = theta(S*R+S+1:S*R+2*S)';                    % Output Weights  [1 x S]
bc_cur = theta(end);                                  % Output Bias     scalar
net_k  = Wg_cur * x_k + bh_cur;                      % Hidden net input  [S x 1]
[act_k, d_act] = activation(net_k, act_func);         % Activation + Derivative [S x 1]
y_est  = Wc_cur * act_k + bc_cur;                     % Estimated Output  scalar
dYdWg  = (Wc_cur' .* d_act * x_k)';                  % [1 x S*R]
dYdbh  = (Wc_cur' .* d_act)';                         % [1 x S]
dYdWc  = act_k';                                      % [1 x S]
dYdbc  = 1;                                           % [1 x 1]
H      = [dYdWg, dYdbh, dYdWc, dYdbc];               % Measurement Jacobian [1 x P_size]
end

function [a, da] = activation(z, act_func)
% activation - Activation Function and Its Derivative
switch act_func
    case 'tanh'
        a  = tanh(z);
        da = 1 - a.^2;                                % tanh'(z) = 1 - tanh²(z)
    case 'relu'
        a  = max(0, z);
        da = double(z > 0);                           % relu'(z) = 1 if z>0 else 0
    case 'sigmoid'
        a  = 1 ./ (1 + exp(-z));
        da = a .* (1 - a);                            % sigmoid'(z) = σ(z)(1-σ(z))
    case 'elu'
        alpha = 1.0;
        a     = z .* (z >= 0) + alpha*(exp(z)-1) .* (z < 0);
        da    = double(z >= 0) + alpha*exp(z) .* double(z < 0);
    otherwise
        error('Unknown activation: %s. Choose tanh | relu | sigmoid | elu', act_func)
end
end