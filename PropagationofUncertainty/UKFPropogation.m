clc, clear, close all;
% Uncertainty Propagation - Statistical Linearization
% Written By: Rasit
% Date: 27-Mar-2026
%%
NoD   = 1e4;
mu    = [1; pi/2];
sigma = [0.05 0; 0 0.5];
Ng = 100;
r     = sigma(1,1)*randn(NoD,1) + mu(1);
Theta = sigma(2,2)*randn(NoD,1) + mu(2);

%% Nonlinear Transform: Polar → Cartesian
nonlinear_func = @(x) [x(1)*cos(x(2)); x(1)*sin(x(2))];
x_cart = r.*cos(Theta);
y_cart = r.*sin(Theta);

%% UKF Sigma Points
NoDim    = 2;
N_Kappa  = 3;
Kappa    = N_Kappa - NoDim;
P        = sqrt(N_Kappa * [std(r)^2 0; 0 std(Theta)^2]);
X_1      = [mean(r); mean(Theta)];
X_2      = X_1 + P(:,1);
X_3      = X_1 + P(:,2);
X_4      = X_1 - P(:,1);
X_5      = X_1 - P(:,2);

X1 = nonlinear_func(X_1);
X2 = nonlinear_func(X_2);
X3 = nonlinear_func(X_3);
X4 = nonlinear_func(X_4);
X5 = nonlinear_func(X_5);
Xsigmas  = [X1 X2 X3 X4 X5];
ut_weights = [Kappa/N_Kappa, ones(1,4)*(1/(2*N_Kappa))];
ut_mean  = Xsigmas * ut_weights';
ut_cov   = (Xsigmas-ut_mean)*diag(ut_weights)*(Xsigmas-ut_mean)';



%% Plot
figure('units','normalized','outerposition',[0 0 1 1],'color','w')

subplot(1,2,1)
x1v = linspace(min(r)-0.1,    max(r)+0.1,    Ng);
x2v = linspace(min(Theta)-0.1, max(Theta)+0.1, Ng);
Z1  = kde2d(r, Theta, x1v, x2v);
[X1g, X2g] = meshgrid(x1v(1:end-1), x2v(1:end-1));
surf(X1g, X2g, Z1, 'EdgeColor','none')
xlabel('r'), ylabel('theta [rad]'), zlabel('p(r, theta)')
title('Input Distribution in Polar Space')
grid on, view(40,30)


subplot(1,2,2)
x1v = linspace(min(x_cart)-0.1, max(x_cart)+0.1, Ng);
x2v = linspace(min(y_cart)-0.1, max(y_cart)+0.1, Ng);
Z2  = kde2d(x_cart, y_cart, x1v, x2v);
[X1g, X2g] = meshgrid(x1v(1:end-1), x2v(1:end-1));
surf(X1g, X2g, Z2, 'EdgeColor','none', 'FaceAlpha', 0.6), hold on
colormap(gca, cool)

% UKF Gaussian fit
r_gauss = 1;                                       % Gaussian mesh radius
x1v2 = linspace(ut_mean(1)-r_gauss, ut_mean(1)+r_gauss, Ng);
x2v2 = linspace(ut_mean(2)-r_gauss, ut_mean(2)+r_gauss, Ng);
Z_ukf = gauss3d(ut_mean, ut_cov, x1v2, x2v2);
[X1g2, X2g2] = meshgrid(x1v2, x2v2);
surf(X1g2, X2g2, Z_ukf*0.6, 'EdgeColor','k', 'FaceColor',[1 0 0], 'FaceAlpha', 0.85)

% Propagated sigma points
z_sp = max(Z2(:))*ones(1,5);
stem3(Xsigmas(1,:), Xsigmas(2,:), z_sp, 'Color', [0 1 0], ...
      'MarkerFaceColor', [1 1 0], 'MarkerSize', 10, 'LineWidth', 2)

xlabel('x'), ylabel('y'), zlabel('p(x, y)')
title('Output Distribution in Cartesian Space (Polar to Cartesian)')
legend('True PDF', 'UKF Gaussian Fit', 'Sigma Points', 'Location','northwest')
grid on, view(60,5)


function Z = kde2d(data1, data2, x1v, x2v)
    [cnt,~] = histcounts2(data1, data2, x1v, x2v);
    Z = cnt' / (sum(cnt(:)) * (x1v(2)-x1v(1)) * (x2v(2)-x2v(1)));
end

function Z = gauss3d(mu, S, x1v, x2v)
    [X1,X2] = meshgrid(x1v, x2v);
    XG = [X1(:), X2(:)]';
    d  = XG - mu;
    Z  = zeros(size(XG,2),1);
    for i = 1:size(XG,2)
        Z(i) = exp(-0.5*d(:,i)'/S*d(:,i));
    end
    Z = reshape(Z/(2*pi*sqrt(det(S))), length(x2v), length(x1v));
end