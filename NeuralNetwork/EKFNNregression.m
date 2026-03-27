clc, clear, close all;
% EKF Tuned Neural Network - Nonlinear Regression
% Written By: Rasit
% Date: 27-Mar-2026
%%
%% Training Data Generation
NoD    = 500;                                        % Number Of Data
x_data = linspace(-pi, pi, NoD)';                   % Input  [NoD x 1]
y_true = sin(x_data) + 0.5*sin(3*x_data);           % Target Function
y_meas = y_true + 0.1*randn(NoD, 1);                % Noisy Measurement

%% Neural Network Architecture
R = 1;    % Input  dimension
S = 300;   % Number of hidden neurons
% -----------------------------------------------------------
% State vector: theta = [Wg(S x R); bh(S x 1); Wc(S x 1); bc(1 x 1)]
% Total parameters: P = S*R + S + S + 1 = S*(R+2) + 1
% -----------------------------------------------------------
P_size = S*(R+2) + 1;                                % Number Of State (total parameters)

%% NN Weight Initialization
Wg = 0.1 * randn(S, R);                             % Hidden Layer Weights  [S x R]
bh = 0.1 * randn(S, 1);                             % Hidden Layer Biases   [S x 1]
Wc = 0.1 * randn(1, S);                             % Output Layer Weights  [1 x S]
bc = 0.1 * randn(1, 1);                             % Output Layer Bias     [1 x 1]
theta = [Wg(:); bh; Wc(:); bc];                     % Full State Vector     [P x 1]

% ---------------------- P~Q~R Matrix ----------------------
% Low  Kalman Gain: if R >> P => K = P/(P+R) | K ~= 0 (Algorithm belief kalman model)
% High Kalman Gain: if R << P => K = P/(P+R) | K ~= 1 (Algorithm belief measurement)
P = 1 * eye(P_size);                              % Error Noise Cov Matrix
Q = 1e-5 * eye(P_size);                             % Processes Noise Cov Matrix
R_noise = 1e-1;                                     % Measurement Noise Cov Matrix

%% Parameter History Storage
theta_hist = zeros(NoD, P_size);                    % Full parameter history

%% EKF Learning and Visualization
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for k = 1:NoD
    % Extract Current Parameters from State Vector
    Wg_cur = reshape(theta(1:S*R), S, R);           % Hidden Weights  [S x R]
    bh_cur = theta(S*R+1:S*R+S);                    % Hidden Biases   [S x 1]
    Wc_cur = theta(S*R+S+1:S*R+2*S)';               % Output Weights  [1 x S]
    bc_cur = theta(end);                             % Output Bias     scalar

    x_k = x_data(k);                                % Current Input   (scalar)
    z_k = y_meas(k);                                % Current Label   (Measurement)

    % Forward Pass
    net_k  = Wg_cur * x_k + bh_cur;                 % Hidden net input  [S x 1]
    act_k  = tanh(net_k);                            % Hidden activation [S x 1]
    y_est  = Wc_cur * act_k + bc_cur;               % Estimated Output  (scalar)

    % Jacobian H (1 x P_size) -> H = [dY/dWg | dY/dbh | dY/dWc | dY/dbc]
    d_act  = 1 - act_k.^2;                          % tanh derivative   [S x 1]
    dYdWg  = (Wc_cur' .* d_act * x_k)';             % [1 x S*R]
    dYdbh  = (Wc_cur' .* d_act)';                   % [1 x S]
    dYdWc  = act_k';                                 % [1 x S]
    dYdbc  = 1;                                      % [1 x 1]
    H = [dYdWg, dYdbh, dYdWc, dYdbc];               % Measurement Jacobian (1 x P_size)

    % Time Update (Prediction) Phase
    theta_pred = theta;                             
    P_pred     = P + Q;                             

    % Measurement Update (Correction) Phase
    S_inn = H * P_pred * H' + R_noise;
    K     = (P_pred * H') / S_inn;                  % Compute Kalman Gain!
    e     = z_k - y_est;                             % Innovation (residual)
    theta = theta_pred + K * e;                      % Update State with Measurement & Kalman Gain
    P     = (eye(P_size)-K*H)*P_pred*(eye(P_size)-K*H)' + K*R_noise*K';  % Update Estimation Uncertainty

    % Store History
    theta_hist(k,:) = theta';

    % Visualize Current NN Fit
    if mod(k, 1e1) == 0 || k == NoD
        % Evaluate NN 
        x_vis  = linspace(-pi, pi, 300)';
        y_vis  = zeros(300, 1);
        Wg_vis = reshape(theta(1:S*R), S, R);
        bh_vis = theta(S*R+1:S*R+S);
        Wc_vis = theta(S*R+S+1:S*R+2*S)';
        bc_vis = theta(end);
        for j = 1:300
            y_vis(j) = Wc_vis * tanh(Wg_vis * x_vis(j) + bh_vis) + bc_vis;
        end

        clf
        subplot(2,2,[1 2])
        plot(x_data, y_true, 'b', 'LineWidth', 2), hold on, grid on
        plot(x_data, y_meas, 'k', 'LineWidth', 1)
        plot(x_vis, y_vis, 'r-', 'LineWidth', 2)
        xlabel('x'), ylabel('y')
        legend('True Function', 'Noisy Measurement', 'EKF-NN Estimate', 'Location', 'northwest')
        title(sprintf('EKF Tuned Neural Network - 1D Regression  |  Step %d / %d', k, NoD))
        axis([-4 4 -2 2])

        subplot(2,2,3)
        plot(1:k, theta_hist(1:k, 1:S), 'LineWidth', 1.2), grid on
        title('Hidden Layer Weights (Wg)')
        xlabel('Number of Data Points'), ylabel('Parameter Values')

        subplot(2,2,4)
        plot(1:k, theta_hist(1:k, S*R+S+1:S*R+2*S), 'LineWidth', 1.2), grid on
        title('Output Layer Weights (Wc)')
        xlabel('Number of Data Points'), ylabel('Parameter Values')

        drawnow
    end
end