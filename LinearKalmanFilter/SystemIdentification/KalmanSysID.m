clc,clear all,close all
% Kalman Filter Based System Identification
% Written By: Rasit
% Date: 21-Jul-2024
%% ARX Model: y(t)=a(t)*y(t-1)+b(t)*u(t-1)+e(t)
true_a = 0.5;
true_b = 0.7;
e_t = .1*randn;   % Output Noise
NoD = 100;
u = randn(NoD, 1);  % Random Input Signal 
y = zeros(NoD, 1);  % System Response


theta = zeros(3,NoD);  % Kalman States
P = 1e3*eye(3);    % Estimation Noise Cov Matrix
Q = 1e1*eye(3);    % Process Noice Cov Matrix
R = 1e-1;          % Measurement Noise Cov Matrix
theta_hat = [0 0 0]';

figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w')
for t = 2:NoD
    y(t) = true_a*y(t-1)+true_b*u(t-1)+e_t;  % System Response

    if t > 50    % Time Varying Parameters
        true_a = 1;
        true_b = 1.4;
    end
    phi = [y(t-1); u(t-1); 1];            % Regressor Matrix

    % Prediction
    theta_hat = theta_hat;
    P = P+Q;

    % Correction
    K = P*phi/(phi'*P*phi+R);    % Kalman Gain
    S = y(t)-(phi'*theta_hat);   % Innovation
    theta_hat = theta_hat+K*S;
    P = (eye(size(theta,1))-K*phi')*P;
    theta(:, t) = theta_hat;
    
    if mod(t,10) == 0
    clf
    subplot(221)
    plot(1:t,y(1:t))
    title("System Output: y(t)=a(t)*y(t-1)+b(t)*u(t-1)+e(t)")
    xlabel('Sample'),ylabel('y(t)');

    subplot(222)
    plot(1:t, theta(1,1:t),"r",LineWidth=2),hold on
    yline(true_a,"b--",LineWidth=2),xline(50,"b--",LineWidth=2)
    title('Estimated System Parameters');
    xlabel('Sample'),ylabel('a(t)'),axis([-1 101 -.1 1.5])
    legend("Estimated Parameter","Real Parameter")

    subplot(223)
    plot(1:t, theta(2,1:t),"r",LineWidth=2),hold on
    yline(true_b,"b--",LineWidth=2),xline(50,"b--",LineWidth=2)
    title('Estimated System Parameters');
    xlabel('Sample'),ylabel('b(t)'),axis([-1 101 -.1 1.75])
    legend("Estimated Parameter","Real Parameter")

    subplot(224)
    plot(1:t, theta(3,1:t),"r",LineWidth=2),hold on
    yline(e_t,"b--",LineWidth=2),xline(50,"b--",LineWidth=2)
    title('Estimated System Parameters');
    xlabel('Sample'),ylabel('e(t)'),axis([-1 101 -.5 1.5])
    legend("Estimated Parameter","Real Parameter")
    drawnow
    end
end



