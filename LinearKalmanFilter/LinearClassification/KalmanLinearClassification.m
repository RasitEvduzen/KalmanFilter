clc,clear all, close;
% Kalman Filter Based Binary Classification
% Written By: Rasit
% 07-Jul-2024
%% Create Data
load data
NoD = size(Y,1); % Number of Data

% Kalman filter Parameters
Q = 1e-4*eye(3);
R = 4;
P = 1e4*eye(3);
theta = zeros(3,1);  % Kalman Parameters (w1, w2, b)

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for i = 1:NoD
    clf
    % Prediction Phase
    H = [X(i,1), X(i,2), 1]; % Measurement Matrix Update
    theta_pred = theta;  % F = I
    P_pred = P + Q;

    % Correction Phase
    S = H * P_pred * H' + R;
    K = P_pred * H' / S;

    y_pred = sign(H * theta_pred);
    y_diff = Y(i) - y_pred;
    theta = theta_pred + K * y_diff;
    P = (eye(3) - K * H) * P_pred;


    [p, q] = meshgrid(min(X(:,1)):1e-1:max(X(:,1)), min(X(:,2)):1e-1:max(X(:,2)));
    decision_boundary = theta(1) * p + theta(3) + q * theta(2);

    contourf(p, q, sign(decision_boundary), 'LineColor', 'none'),hold on
    scatter(X(Y == 1, 1), X(Y == 1, 2), 'k',"filled")
    scatter(X(Y == -1, 1), X(Y == -1, 2), 'y',"filled");

    x_line = linspace(min(min(X)), max(max(X)), 100);
    y_line = -(theta(1) * x_line + theta(3)) / theta(2);
    plot(x_line, y_line, 'r',LineWidth=4);

    xlabel('X1'),ylabel('X2'),title('Kalman Filter Based Binary Classification');
    axis([min(min(X)), max(max(X)), min(min(X)), max(max(X))]),axis equal

    drawnow
end
