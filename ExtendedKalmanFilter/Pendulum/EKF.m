clc, clear all, close all;
% Inverted Pendulum - Extended Kalman Filter
% Written By: Rasit Evduzen
% 01-Apr-2024
%% System Parameters
m  = 1;       % Pendulum Mass      [kg]
M  = 1;       % Cart Mass          [kg]
L  = 2;       % Pendulum Length    [m]
g  = -9.81;   % Gravitational Force [m/s^2]
d  = 0.1;     % Cart Damping       [N·s/m]
b  = -1;      % Pendulum Down (b=-1)
Ts = 1e-2;    % Sampling Period    [s]
PendInitialAngle = 10;  % Initial Angle [deg]  (180 -> Unstable | 0 -> Stable)

%% True State Simulation (RK4)
xini  = [0; 0; deg2rad(PendInitialAngle); 0];
u     = 0;                                        % Free Response
tspan = 0:Ts:10;
NoD   = length(tspan);
state = zeros(NoD, 4);
state(1,:) = xini';
for k = 1:NoD-1
    state(k+1,:) = RK4(state(k,:)', m, M, L, g, d, u, Ts)';
end

%% Noisy Measurement
NoisyMeasurement = state' + 5e-2*randn(4, NoD);  % [4 x NoD]

%% EKF Initialization
NofState = 4;
XKalman  = zeros(NofState, NoD);

% ---------------------- P~Q~R Matrix ----------------------
% Low  Kalman Gain: if R >> P => K = P/(P+R) | K ~= 0 (Algorithm belief kalman model)
% High Kalman Gain: if R << P => K = P/(P+R) | K ~= 1 (Algorithm belief measurement)
P = 1e-5 * eye(NofState);                         % Error Noise Cov Matrix
Q = 1e-5 * eye(NofState);                         % Processes Noise Cov Matrix
R = diag([var(NoisyMeasurement(1,:)), ...
          var(NoisyMeasurement(2,:)), ...
          var(NoisyMeasurement(3,:)), ...
          var(NoisyMeasurement(4,:))]);             % Measurement Noise Cov Matrix
C = eye(NofState);                                 % Output Selection Matrix

%% EKF Loop
for i = 1:NoD-1
    % Time Update (Prediction) Phase
    XKalman(:,i) = RK4(XKalman(:,i), m, M, L, g, d, u, Ts); % Nonlinear State Propagation
    Phi = calcJacobian(XKalman(:,i), m, M, L, g, d, Ts);     % Jacobian at current estimate
    P   = Phi*P*Phi' + Q;                          % Uncertainty Propagation (Jacobian → Gaussian preserved)

    % Measurement Update (Correction) Phase
    K = P*C' / (C*P*C' + R);                      % Compute Kalman Gain!
    XKalman(:,i+1) = XKalman(:,i) + K*(NoisyMeasurement(:,i) - C*XKalman(:,i)); % Update State with Measurement & Kalman Gain
    P = (eye(NofState)-K*C)*P*(eye(NofState)-K*C)' + K*R*K'; % Update Estimation Uncertainty
end

%% Plot Simulation
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
PlotIdx = round(linspace(1, NoD-1, 30));           % 30 animation frames
for i = PlotIdx
    clf
    PlotPend(state, i, L, tspan, NoisyMeasurement, XKalman)
    drawnow
end

%% Utility Functions
function x_next = RK4(x, m, M, L, g, d, u, Ts)
% RK4 - 4th Order Runge-Kutta Integration
k1 = PendDynamics(x,           m, M, L, g, d, u);
k2 = PendDynamics(x + Ts/2*k1, m, M, L, g, d, u);
k3 = PendDynamics(x + Ts/2*k2, m, M, L, g, d, u);
k4 = PendDynamics(x + Ts*k3,   m, M, L, g, d, u);
x_next = x + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4);
end

function dx = PendDynamics(x, m, M, L, g, d, u)
% PendDynamics - Nonlinear Pendulum State Space
D       = m*L*L*(M + m*(1-cos(x(3))^2));
dx(1,1) = x(2);
dx(2,1) = (1/D)*(-m^2*L^2*g*cos(x(3))*sin(x(3)) + m*L^2*(m*L*x(4)^2*sin(x(3)) - d*x(2))) + m*L*L*(1/D)*u;
dx(3,1) = x(4);
dx(4,1) = (1/D)*((m+M)*m*g*L*sin(x(3)) - m*L*cos(x(3))*(m*L*x(4)^2*sin(x(3)) - d*x(2))) - m*L*cos(x(3))*(1/D)*u;
end

function Phi = calcJacobian(x, m, M, L, g, d, Ts)
% calcJacobian - Analytical Jacobian of pendulum dynamics (linearized at x)
D      = m*L*L*(M + m*(1-cos(x(3))^2));
df2dx3 = (1/D)*(-m^2*L^2*g*(cos(x(3))^2 - sin(x(3))^2) + m*L^2*m*L*x(4)^2*cos(x(3)));
df4dx3 = (1/D)*((m+M)*m*g*L*cos(x(3)) - m*L*(-sin(x(3)))*(m*L*x(4)^2*sin(x(3)) - d*x(2)) ...
         - m*L*cos(x(3))*m*L*x(4)^2*cos(x(3)));
Ac  = [0  1      0       0;
       0 -d/M   df2dx3   0;
       0  0      0       1;
       0  0     df4dx3   0];
Phi = eye(4) + Ac*Ts;                              % Euler discretization
end

function PlotPend(state, k, L, t, NoisyMeasurement, XKalman)
% PlotPend - Pendulum Animation + State Plots
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
title(sprintf('Inverted Pendulum - Extended Kalman Filter  |  t = %.2f s', t(k)))

subplot(3,2,3)
plot(t, NoisyMeasurement(1,:), 'Color',[.7 .7 .7], 'LineWidth',2), hold on, grid on
plot(t, XKalman(1,:),  'k', 'LineWidth',2)
plot(t, state(:,1)', 'r--', 'LineWidth',2)
xline(t(k), 'b--', 'LineWidth',1)
legend('Noisy','EKF','Real','Location','northeast')
xlabel('Time [s]'), ylabel('Cart Position [m]')

subplot(3,2,4)
plot(t, NoisyMeasurement(2,:), 'Color',[.7 .7 .7], 'LineWidth',2), hold on, grid on
plot(t, XKalman(2,:),  'k', 'LineWidth',2)
plot(t, state(:,2)', 'r--', 'LineWidth',2)
xline(t(k), 'b--', 'LineWidth',1)
legend('Noisy','EKF','Real','Location','northeast')
xlabel('Time [s]'), ylabel('Cart Velocity [m/s]')

subplot(3,2,5)
plot(t, NoisyMeasurement(3,:), 'Color',[.7 .7 .7], 'LineWidth',2), hold on, grid on
plot(t, XKalman(3,:),  'k', 'LineWidth',2)
plot(t, state(:,3)', 'r--', 'LineWidth',2)
xline(t(k), 'b--', 'LineWidth',1)
legend('Noisy','EKF','Real','Location','northeast')
xlabel('Time [s]'), ylabel('Pendulum Angle [rad]')

subplot(3,2,6)
plot(t, NoisyMeasurement(4,:), 'Color',[.7 .7 .7], 'LineWidth',2), hold on, grid on
plot(t, XKalman(4,:),  'k', 'LineWidth',2)
plot(t, state(:,4)', 'r--', 'LineWidth',2)
xline(t(k), 'b--', 'LineWidth',1)
legend('Noisy','EKF','Real','Location','northeast')
xlabel('Time [s]'), ylabel('Pendulum Angular Velocity [rad/s]')
end