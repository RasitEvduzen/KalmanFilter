clc,clear,close all;
% Uncertainty Linear and Nonlinear Propogation
% Written By: Rasit
% 06-Nov-2024
%%
NoD = 1e4;   
mu = [1 pi/2];
sigma = [0.05 0; 0 0.5];
r = sigma(1)*randn(NoD,1)+mu(1);
Theta = sigma(4)*randn(NoD,1)+mu(2);

%% Input Distribution Nonlinear Mapped
x = r.*cos(Theta);
y = r.*sin(Theta);


%% Linearized Covariance 
Jacobian = [0 -1; 1 0];
phat = Jacobian*sigma*Jacobian';
muhat = [0 1];
xlinear = phat(1)*randn(NoD,1)+muhat(1);
ylinear = phat(4)*randn(NoD,1)+muhat(2);

%% Output Covariance via Unscented Transform
% Calculate Sigma Point
NoDim = 2;   % One Dimension Random Distribution
NoSigma = 2*NoDim+1;  % Number of Sigma Points
N_Kappa = 3;   % For Gauss-PDF
Kappa = N_Kappa - NoDim;
P = sqrt(N_Kappa*[std(r)^2 0; 0 std(Theta)^2]);
X_1 = mean([r Theta])';   % Sigma Point-1: Mean Input Distribution
X_2 = mean([r Theta])' + P(:,1); % Sigma Point-2
X_3 = mean([r Theta])' + P(:,2); % Sigma Point-3
X_4 = mean([r Theta])' - P(:,1); % Sigma Point-2
X_5 = mean([r Theta])' - P(:,2); % Sigma Point-3

% Nonlinear Mapped Sigma Points
nonlinear_func = @(x) [x(1)*cos(x(2)); x(1)*sin(x(2))];
X1 = nonlinear_func(X_1);  % Nonlinear Mapped Sigma Point
X2 = nonlinear_func(X_2);  % Nonlinear Mapped Sigma Point
X3 = nonlinear_func(X_3);  % Nonlinear Mapped Sigma Point
X4 = nonlinear_func(X_4);  % Nonlinear Mapped Sigma Point
X5 = nonlinear_func(X_5);  % Nonlinear Mapped Sigma Point
Xsigmas = [X1 X2 X3 X4 X5];
ut_weights = [Kappa/(N_Kappa) ones(1,4)*(1/(2*N_Kappa))];
ut_mean = Xsigmas*ut_weights';
ut_cov  = (Xsigmas-ut_mean)*diag(ut_weights)*(Xsigmas-ut_mean)';
x_ut = sqrt(ut_cov(1))*randn(NoD,1)+ut_mean(1);
y_ut = sqrt(ut_cov(4))*randn(NoD,1)+ut_mean(2);


%%  Plot Result
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
subplot(121)
scatter(r,Theta,[],[0.5 0.5 0.5], 'filled'),hold on
[result1, result2] = plotCovarianceEllipse([r Theta]);
plot(result1, result2, 'r', 'LineWidth', 4);
xlabel("r"),ylabel("theta")
title("Random Sample in Polar Coordinate")
legend("Random Sample","Input Covariance",Location="southwest")
axis tight
% axis off

subplot(122)
scatter(x,y,[],[0.5 0.5 0.5], 'filled'),hold on
[result3, result4] = plotCovarianceEllipse([x y]);  % Nonlinear Mapped data
plot(result3, result4, 'r', 'LineWidth', 4);        % Transformed Data Ellipse
plot(result1.*cos(result2),result1.*sin(result2),'g','LineWidth',4) % Uncertainty Mapped Nonlinear Transform
[result5, result6] = plotCovarianceEllipse([xlinear ylinear]);  % Linearized Covariance
plot(result5, result6, 'k-', 'LineWidth', 4);        % Transformed Data Ellipse
[result7, result8] = plotCovarianceEllipse([x_ut y_ut]);  % Unscented Transform Ellipse
plot(result7, result8, 'b--', 'LineWidth', 4);
scatter(Xsigmas(1,:),Xsigmas(2,:),100,"y","filled")
xlabel("x"),ylabel("y")
title("Random Sample in Cartesian Coordinate")
legend("Random Sample","Output True Covariance","Nonlinear Mapped Input Covariance", ...
    "EKF Linearized Covariance","Unscented Transformed Covariance","Sigma Points",Location="southwest")
axis tight
% axis off
