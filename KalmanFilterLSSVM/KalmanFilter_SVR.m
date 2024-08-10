% Kalman-tuned Least Squares Support Vector for Nonlinear Regression
% Original Written By: Rasit
% 07-Mar-2024
clc,clear all,close all;
%% Input and Output definition
xtrain = [1:0.1:20]';
NoD = length(xtrain);

% ytrain = sin(xtrain) + 2*sin(2*xtrain)+0.5*randn(NoD,1);
ytrain = 0.01*xtrain.*xtrain + 0.1*exp(-xtrain) + sin(xtrain) + 0.1*randn(NoD,1);

%% Train RLS-SVM
C = 100;     % Over Fitting Param (C=100)
gamma = 5e-1; % RBF param, equal to 1/2sigma^2 (g=1e-2)

kernelSelect = 'rbf';
K = Kernel(kernelSelect,xtrain,gamma);

A = [0, ones(1,NoD);
     ones(NoD,1), (K + 1/C*eye(NoD))];
b = [0; ytrain];

%% Kalman-Tuned RLS-SVR
x_kalman = rand(size(A,1),1);  % Random kalman state vector

% Kalman Filter Parameters
P = 1e1*eye(size(A,1),size(A,1));
Q = 1e-5*eye(size(A,1)); % Process noise covariance
R = 1e-2;                  % Measurement noise covariance

% Prediction Stage
xpred = xtrain; % Prediction Input
ypred = zeros(NoD,1);
tmp = zeros(NoD,1);

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for k=1:NoD
    [x_kalman,K,P] = kalman_filter(A(k,:),b(k,:),x_kalman,P,Q,R);
    % Bias and alpha LaGrange Multipliers
    bias = x_kalman(1);
    alpha = x_kalman(2:end);

    for j=1:NoD
        for i=1:NoD
            tmp(i,1) = alpha(i,1).*exp(-gamma*(xpred(j,1)-xpred(i,1))^2);
        end
        ypred(j,1) = sum(tmp) + bias;
    end
    % Plot Result
    clf
    plot(xtrain, ytrain,'ko-', 'LineWidth', 2), hold on;
    plot(xpred, ypred,'r.-', 'LineWidth', 2);
    axis([0 21 -2 5.5])
    title('Kalman-Tuned SVR Nonlinear Regression');
    legend('Noisy Data', 'Kalman-Tuned SVR ');
    drawnow
end

function [x,K,P] = kalman_filter(a_k,b_k,x,P,Q,R)
a_k = a_k(:);
b_k = b_k(:);
% Predict
x_pred = x;            % Predicted state (No process model here)
P_pred = P + Q;        % Predicted covariance

% Correction
K = (P_pred*a_k)/(a_k'*P_pred*a_k + R); % Compute Kalman Gain
x = x_pred + K*(b_k - a_k'*x_pred);     % State Update
P = (eye(size(K,1)) - K*a_k')*P_pred;   % Covariance Update
end
