clc, clear, close all;
% EKF Tuned RBF Neural Network - Nonlinear Classification
% Written By: Rasit
% Date: 27-Mar-2026
%% Spiral Training Data Generation
B = 2;
N = 3e2;
Tall = [];
for i = 1:N/2
    theta_sp = pi/2 + (i-1)*[(2*B-1)/N]*pi;
    Tall = [Tall, [theta_sp*cos(theta_sp); theta_sp*sin(theta_sp)]];
end
Tall     = [Tall, -Tall];
Tmax     = pi/2 + [(N/2-1)*(2*B-1)/N]*pi;
xtrain   = [Tall]'/Tmax;                             % [N x 2] Input Data
ytrain   = [-ones(1, N/2), ones(1, N/2)]';           % [N x 1] Class Labels {-1, +1}

NoD      = length(xtrain);                           % Number Of Data
InputDim = 2;                                        % Input Dimension

%% Test Data
xtest = [];
for t1 = -1:1e-2:1
    for t2 = -1:1e-2:1
        xtest = [xtest; [t1, t2]];
    end
end

%% RBF NN Parameter Initialization
% State vector theta = [w(N x 1); c(2N x 1); sigma(N x 1)] -> 4N x 1
NumberOfNeuron = 160;                                 % Number of RBF Neurons
sigma_init     = 1 / (sqrt(NumberOfNeuron));         

w     = 0.1 * randn(NumberOfNeuron, 1);              % Output Weights
idx_c = randperm(NoD, NumberOfNeuron);
c     = xtrain(idx_c, :);                            % Centers init from training data [N x 2]
sigma = sigma_init * ones(NumberOfNeuron, 1);        % RBF Widths  scaled to data spread

theta    = [w; c(:); sigma];                         % Full State Vector [4N x 1]
NofState = length(theta);                            % Number Of State = 4N

% ---------------------- P~Q~R Matrix ----------------------
% Low  Kalman Gain: if R >> P => K = P/(P+R) | K ~= 0 (Algorithm belief kalman model)
% High Kalman Gain: if R << P => K = P/(P+R) | K ~= 1 (Algorithm belief measurement)
P = 1e2*eye(NofState);                                                     % Error Noise Cov Matrix
Q = blkdiag(1e-3*eye(NumberOfNeuron), ...
            1e-3*eye(2*NumberOfNeuron), ...
            1e-3*eye(NumberOfNeuron));                                      % Processes Noise Cov Matrix
R = 1e-1;                                                                  % Measurement Noise Cov Matrix

%% Parameter Storage
w_hist     = zeros(NoD, NumberOfNeuron);
c_hist     = zeros(NoD, 2*NumberOfNeuron);
sigma_hist = zeros(NoD, NumberOfNeuron);

%% EKF Learning 
figure('units', 'normalized', 'color', 'w')
for k = 1:NoD
    w_cur     = theta(1:NumberOfNeuron);                                   % Output Weights  [N x 1]
    c_cur     = reshape(theta(NumberOfNeuron+1:3*NumberOfNeuron), NumberOfNeuron, InputDim);  % Centers [N x 2]
    sigma_cur = abs(theta(3*NumberOfNeuron+1:end));                        % RBF Widths > 0  [N x 1]

    x_k = xtrain(k,:)';   % Current Input  [2 x 1]
    z_k = ytrain(k);      % Current Label  (Measurement)

    % RBF Activations
    diff_k  = x_k' - c_cur;                                                % [N x 2] diff vectors
    dist2_k = sum(diff_k.^2, 2);                                           % [N x 1] squared distances
    phi_k   = exp(-dist2_k ./ (2*sigma_cur.^2));                           % [N x 1] activations
    y_est   = w_cur' * phi_k;                                              % Estimated Output

    % Jacobian H (1 x 4N)  ->  H = [dY/dw | dY/dc | dY/dsigma]
    dYdW = phi_k';                                                         % [1 x N]
    dYdC = zeros(1, 2*NumberOfNeuron);                                     % [1 x 2N]
    for n = 1:NumberOfNeuron
        dYdC(1, 2*n-1:2*n) = w_cur(n) * phi_k(n) * diff_k(n,:) / sigma_cur(n)^2;
    end
    dYdS = (w_cur .* phi_k .* dist2_k ./ (sigma_cur.^3))';                % [1 x N]
    H    = [dYdW, dYdC, dYdS];                                             % Measurement Jacobian (1 x 4N)

    % Time Update (Prediction) Phase
    theta_pred = theta;              
    P_pred     = P + Q;              

    % Measurement Update (Correction) Phase
    S = H * P_pred * H' + R;
    K = (P_pred * H') / S;          
    e = z_k - y_est;                 % Innovation 
    theta = theta_pred + K * e;      
    P = (eye(NofState)-K*H)*P_pred*(eye(NofState)-K*H)' + K*R*K';        

    w_hist(k,:)     = theta(1:NumberOfNeuron)';
    c_hist(k,:)     = reshape(theta(NumberOfNeuron+1:3*NumberOfNeuron), 1, 2*NumberOfNeuron);
    sigma_hist(k,:) = abs(theta(3*NumberOfNeuron+1:end))';

    % Classify Test Grid
    w_vis     = theta(1:NumberOfNeuron);
    c_vis     = reshape(theta(NumberOfNeuron+1:3*NumberOfNeuron), NumberOfNeuron, InputDim);
    sigma_vis = abs(theta(3*NumberOfNeuron+1:end));
    ytest     = zeros(size(xtest,1), 1);
    for j = 1:size(xtest,1)
        diff_t   = xtest(j,:) - c_vis;
        dist2_t  = sum(diff_t.^2, 2);
        phi_t    = exp(-dist2_t ./ (2*sigma_vis.^2));
        ytest(j) = w_vis' * phi_t;
    end
    ytest = sign(ytest);

    % Plot Result
    if mod(k, 1e1) == 0|| k == NoD
        clf
        plot(xtest(ytest==+1,1), xtest(ytest==+1,2), 'r.'), hold on, grid on
        plot(xtest(ytest==-1,1), xtest(ytest==-1,2), 'b.')
        plot(xtrain(ytrain==+1,1), xtrain(ytrain==+1,2), 'k*', 'LineWidth', 5)
        plot(xtrain(ytrain==-1,1), xtrain(ytrain==-1,2), 'y*', 'LineWidth', 5)
        title(sprintf('EKF Tuned RBF NN - Nonlinear Classification  |  Step %d / %d', k, NoD))
        drawnow
    end
end

% figure('units','normalized','outerposition',[0 0 1 1],'color','w')
% 
% subplot(3,1,1)
% plot(1:NoD, w_hist, 'LineWidth', 1.2), grid on
% title('Estimated Parameters (Weights)')
% xlabel('Number of Data Points'), ylabel('Parameter Values')
% 
% subplot(3,1,2)
% plot(1:NoD, c_hist, 'LineWidth', 1.2), grid on
% title('Estimated Parameters (Centers)')
% xlabel('Number of Data Points'), ylabel('Parameter Values')
% 
% subplot(3,1,3)
% plot(1:NoD, sigma_hist, 'LineWidth', 1.2), grid on
% title('Estimated Parameters (\sigma  -  RBF Widths)')
% xlabel('Number of Data Points'), ylabel('Parameter Values')
