%% Kalman Filter Tuned Extreme Learning Machine
% Written By: Rasit
% 30-Jul-2024
clc, clear all, close all;
%% Nonlinear Classification
% Create Data Spiral Data
B = 4;
N = 200;
Tall = [];
for i = 1:N/2
    theta = pi/2 + (i-1)*[(2*B-1)/N]*pi;
    Tall = [Tall, [theta*cos(theta); theta*sin(theta)]];
end
Tall = [Tall, -Tall];
Tmax = pi/2 + [(N/2-1)*(2*B-1)/N]*pi;
xtrain = [Tall]'/Tmax;
ytrain = [-ones(1, N/2), ones(1, N/2)]';

NoD = length(xtrain);
NumberOfNeuron = 100;  % Number of Neurons

input_w = 2*randn(NumberOfNeuron, size(xtrain, 2));  % generate random input weights
H = tanh(xtrain * input_w');  % TanH Activation Function

%% Kalman Filter Based ELM 
output_w = rand(size(H, 2), 1);  % Initial state vector for Kalman filter (ELM Weight)
P = 1e6*eye(size(H, 2));   % Initial error covariance matrix
Q = 1e-3*eye(size(H, 2));  % Process noise covariance matrix
R = 1e-2;  % Measurement noise covariance

xtest = []; % Test Input
ytest = []; % Test Output
for t1 = -1:1e-2:1
    for t2 = -1:1e-2:1
        xtest = [xtest; [t1, t2]];
    end
end

figure('units', 'normalized', 'color', 'w')
for k = 1:NoD
    % Kalman Filter Prediction Phase
    h_k = H(k, :)';      % Output Selection Matrix
    z_k = ytrain(k, :);  % Measurement

    % Correction Phase
    K = (P*h_k)/(h_k'*P*h_k+R); % Compute Kalman Gain
    output_w = output_w+K*(z_k-h_k'*output_w);  % Update state estimate
    P = (eye(size(P))-K*h_k')*P+Q;  % Update error covariance

    % ELM ~ Prediction
    Ht = tanh(xtest * input_w'); % ELM Hidden Layer
    ytest = sign(Ht * output_w); % ELM Output

    % Plot Result
    clf
    plot(xtest(ytest == +1, 1), xtest(ytest == +1, 2), 'r.'), hold on, grid on
    plot(xtest(ytest == -1, 1), xtest(ytest == -1, 2), 'b.')

    plot(xtrain(ytrain == +1, 1), xtrain(ytrain == +1, 2), 'k*', 'LineWidth', 5)
    plot(xtrain(ytrain == -1, 1), xtrain(ytrain == -1, 2), 'y*', 'LineWidth', 5)
    title('Kalman Filter Tuned ELM - Nonlinear Classification')
    drawnow
end
