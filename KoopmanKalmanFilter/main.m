clc, clear all, close all;
% Inverted Pendulum  Koopman KF
% Offline eDMD lifting + Koopman operator as Jacobian substitute
% Lifting: 'physics' 'poly' 'rbf'
% Written By: Rasit Evduzen
% 28-Apr-2026
%% System Parameters
m  = 1;       % Pendulum Mass      [kg]
M  = 1;       % Cart Mass          [kg]
L  = 2;       % Pendulum Length    [m]
g  = -9.81;   % Gravitational Force [m/s^2]
d  = 0.1;     % Cart Damping       [N·s/m]
Ts = 1e-2;    % Sampling Period    [s]
PendInitialAngle = 75;  % Initial Angle [deg]

%% Lifting Selection
liftType = 'poly';                                % 'physics' | 'poly' | 'rbf'
Nc       = 8;                                        % RBF centers (only for rbf)
rbf_eps  = 2.0;                                      % RBF width

%% True State Simulation (RK4)
xini  = [0; 0; deg2rad(PendInitialAngle); 0];
u     = 0;
tspan = 0:Ts:10;
NoD   = length(tspan);
state = zeros(NoD, 4);
state(1,:) = xini';
for k = 1:NoD-1
    state(k+1,:) = RK4(state(k,:)', m, M, L, g, d, u, Ts)';
end

%% Noisy Measurement
NoisyMeasurement = state' + 1e-1*randn(4, NoD);

%% Offline eDMD Training (PRBS excitation)
N_data   = 3000;
alpha_kp = 1e-6;                                     % Tikhonov regularization
rng(42);
x_buf = xini;
X_data = zeros(4, N_data);
for k = 1:N_data
    X_data(:,k) = x_buf;
    u_prbs      = 0.1*(2*rand-1);
    x_buf       = RK4(x_buf, m, M, L, g, d, u_prbs, Ts);
end

% RBF centers
if strcmp(liftType, 'rbf')
    centers = X_data(:, round(linspace(1,N_data,Nc)));
else
    centers = [];
end

% Lift all snapshots
n_lift = GetNlift(liftType, Nc);
G  = zeros(n_lift, N_data-1);
Gp = zeros(n_lift, N_data-1);
for k = 1:N_data-1
    G(:,k)  = LiftState(X_data(:,k),   liftType, centers, rbf_eps);
    Gp(:,k) = LiftState(X_data(:,k+1), liftType, centers, rbf_eps);
end

% eDMD: Koopman operator K_koop
K_koop = Gp * G' / (G*G' + alpha_kp*eye(n_lift));   % [n_lift x n_lift]

%% Kalman Initialization
NofState = 4;
XKalman  = zeros(NofState, NoD);
XKalman(:,1) = xini;

% Low  Kalman Gain: if R >> P => K = P/(P+R) | K ~= 0 (Algorithm belief kalman model)
% High Kalman Gain: if R << P => K = P/(P+R) | K ~= 1 (Algorithm belief measurement)
P = 1e-5 * eye(NofState);                            % Error Noise Cov Matrix
Q = 1e-5 * eye(NofState);                            % Process Noise Cov Matrix
R = diag([var(NoisyMeasurement(1,:)), ...
          var(NoisyMeasurement(2,:)), ...
          var(NoisyMeasurement(3,:)), ...
          var(NoisyMeasurement(4,:))]);               % Measurement Noise Cov Matrix
C = eye(NofState);                                   % Output Selection Matrix

%% Koopman KF Loop
for i = 1:NoD-1
    % Koopman operator in original state space (top-left block)
    g_cur  = LiftState(XKalman(:,i), liftType, centers, rbf_eps);
    g_next = K_koop * g_cur;
    A_koop = K_koop(1:NofState, 1:NofState);         % Local linearization from Koopman

    % Time Update (Prediction) Phase
    XKalman(:,i) = RK4(XKalman(:,i), m, M, L, g, d, u, Ts); % Nonlinear State Propagation
    P = A_koop*P*A_koop' + Q;                        % Uncertainty Propagation (Koopman)

    % Measurement Update (Correction) Phase
    K     = P*C' / (C*P*C' + R);                     % Compute Kalman Gain!
    XKalman(:,i+1) = XKalman(:,i) + K*(NoisyMeasurement(:,i) - C*XKalman(:,i)); % Update State with Measurement & Kalman Gain
    P = (eye(NofState)-K*C)*P*(eye(NofState)-K*C)' + K*R*K'; % Update Estimation Uncertainty
end

%% Plot Simulation
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
PlotIdx = round(linspace(1, NoD-1, 30));
for i = PlotIdx
    clf
    PlotPend(state, i, L, tspan, NoisyMeasurement, XKalman, liftType)
    drawnow
end

%% Utility Functions
function n = GetNlift(liftType, Nc)
switch liftType
    case 'physics'; n = 6;                           % [x1,x2,x3,x4, sin(x3),cos(x3)]
    case 'poly';    n = 10;                          % [x1..x4, x1²,x2²,x3²,x4², x1x3, 1]
    case 'rbf';     n = 4 + Nc;                     % [x1..x4, rbf_1..rbf_Nc]
end
end

function g = LiftState(x, liftType, centers, rbf_eps)
switch liftType
    case 'physics'
        % Physics-informed: original states + trigonometric terms
        g = [x(1); x(2); x(3); x(4); sin(x(3)); cos(x(3))];
    case 'poly'
        % Polynomial basis up to 2nd order
        g = [x(1); x(2); x(3); x(4); ...
             x(1)^2; x(2)^2; x(3)^2; x(4)^2; ...
             x(1)*x(3); 1];
    case 'rbf'
        % Gaussian RBF
        Nc   = size(centers, 2);
        rbfs = zeros(Nc, 1);
        for k = 1:Nc
            dv     = x - centers(:,k);
            rbfs(k) = exp(-rbf_eps^2 * (dv'*dv));
        end
        g = [x(1); x(2); x(3); x(4); rbfs];
end
end

function x_next = RK4(x, m, M, L, g, d, u, Ts)
k1 = PendDynamics(x,           m, M, L, g, d, u);
k2 = PendDynamics(x + Ts/2*k1, m, M, L, g, d, u);
k3 = PendDynamics(x + Ts/2*k2, m, M, L, g, d, u);
k4 = PendDynamics(x + Ts*k3,   m, M, L, g, d, u);
x_next = x + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function dx = PendDynamics(x, m, M, L, g, d, u)
D       = m*L*L*(M + m*(1-cos(x(3))^2));
dx(1,1) = x(2);
dx(2,1) = (1/D)*(-m^2*L^2*g*cos(x(3))*sin(x(3)) + m*L^2*(m*L*x(4)^2*sin(x(3)) - d*x(2))) + m*L*L*(1/D)*u;
dx(3,1) = x(4);
dx(4,1) = (1/D)*((m+M)*m*g*L*sin(x(3)) - m*L*cos(x(3))*(m*L*x(4)^2*sin(x(3)) - d*x(2))) - m*L*cos(x(3))*(1/D)*u;
end

function PlotPend(state, k, L, t, NoisyMeasurement, XKalman, liftType)
x     = state(k,1);
th    = state(k,3);
W  = 1.5;  H  = 0.5;  wr = 0.5;  mr = 0.67;
y     = wr/2 + H/2;
pendx = x + L*sin(th);
pendy = y - L*cos(th);

subplot(3,2,[1 2])
yline(0,'k--','LineWidth',2), hold on
rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',0.1,'FaceColor',[0.4940 0.1840 0.5560],'LineWidth',1.5)
rectangle('Position',[x-0.9*W/2,0,wr,wr],'Curvature',1,'FaceColor',[1 1 0],'LineWidth',1.5)
rectangle('Position',[x+0.9*W/2-wr,0,wr,wr],'Curvature',1,'FaceColor',[1 1 0],'LineWidth',1.5)
plot([x pendx],[y pendy],'k','LineWidth',2)
rectangle('Position',[pendx-mr/2,pendy-mr/2,mr,mr],'Curvature',1,'FaceColor',[1 0.1 0.1],'LineWidth',1.5)
axis equal, axis([-12 12 -3 3]), grid on
title(sprintf('Inverted Pendulum  Koopman KF [%s]  |  t = %.2f s', liftType, t(k)))

subplot(3,2,3)
plot(t, NoisyMeasurement(1,:), 'Color',[.7 .7 .7], 'LineWidth',2), hold on, grid on
plot(t, XKalman(1,:), 'k', 'LineWidth',2)
plot(t, state(:,1)', 'r--', 'LineWidth',2)
xline(t(k), 'b--', 'LineWidth',1)
legend('Noisy','Koopman KF','Real','Location','northeast')
xlabel('Time [s]'), ylabel('Cart Position [m]')

subplot(3,2,4)
plot(t, NoisyMeasurement(2,:), 'Color',[.7 .7 .7], 'LineWidth',2), hold on, grid on
plot(t, XKalman(2,:), 'k', 'LineWidth',2)
plot(t, state(:,2)', 'r--', 'LineWidth',2)
xline(t(k), 'b--', 'LineWidth',1)
legend('Noisy','Koopman KF','Real','Location','northeast')
xlabel('Time [s]'), ylabel('Cart Velocity [m/s]')

subplot(3,2,5)
plot(t, NoisyMeasurement(3,:), 'Color',[.7 .7 .7], 'LineWidth',2), hold on, grid on
plot(t, XKalman(3,:), 'k', 'LineWidth',2)
plot(t, state(:,3)', 'r--', 'LineWidth',2)
xline(t(k), 'b--', 'LineWidth',1)
legend('Noisy','Koopman KF','Real','Location','northeast')
xlabel('Time [s]'), ylabel('Pendulum Angle [rad]')

subplot(3,2,6)
plot(t, NoisyMeasurement(4,:), 'Color',[.7 .7 .7], 'LineWidth',2), hold on, grid on
plot(t, XKalman(4,:), 'k', 'LineWidth',2)
plot(t, state(:,4)', 'r--', 'LineWidth',2)
xline(t(k), 'b--', 'LineWidth',1)
legend('Noisy','Koopman KF','Real','Location','northeast')
xlabel('Time [s]'), ylabel('Pendulum Angular Velocity [rad/s]')
end