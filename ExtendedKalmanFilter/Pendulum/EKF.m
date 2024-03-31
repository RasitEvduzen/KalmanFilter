clc,clear all,close all;
% Inverted Pendulum Extended Kalman Filter
% 01-Apr-2024
% Rasit EVDUZEN
%%
% System Parameter
m = 1;   % pendulum mass
M = 1;   % cart mass
L = 2;   % pendulum length
g = -9.81; % gravitational force
d = .1;   % cart damping
b = -1;   % pendulum down (b = -1)
PendInitialAngle = 10;  % 180 -> Unstable Fixed Point | 0 -> Stable Fixed Point
Ts = 5e-3;
tspan = 0:Ts:15;
xini = [0; 0; 0+deg2rad(PendInitialAngle); 0];  % initial condition Small Angle Approach
u = 0;   % Free Response
[t,state] = ode45(@(t,x) PendStateSpace(x,m,M,L,g,d,u),tspan,xini);  % Runge Kutta4 for Forward Dynamics Simulation
NoisyMeasurement = state + 5e-2*randn(size(tspan,2),4);
NoisyMeasurement = NoisyMeasurement';

% Linear State Space Model for EKF
A = [0  1          0               0;
     0 -d/M        b*m*g/M         0;
     0  0          0               1;
     0 -b*d/(M*L) -b*(m+M)*g/(M*L) 0];
NofState = size(A,1);                  % Number Of State
C = eye(NofState);    % Output Selection Matrix
sys = ss(A,zeros(4,1),C,0);
sysd = c2d(sys,Ts);
Phi = sysd.A; % Discrete Time State Transition Matrix

% ---------------------- P~Q~R Matrix Random value --------------------
% Low Kalman Gain:   if R >> P => K = P/(P+R) | K ~= 0 (Algorithm belief kalman model)  XKalman(:,i+1) = XKalman(:,i);
% High Kalman Gain:  if R << P => K = P/(P+R) | K ~= 1 (Algorithm belief measurement)   XKalman(:,i+1) = Measurement(:,i);
P = 1e-5*eye(NofState,NofState);          % Error Noise Cov Matrix
Q = 1e-5*eye(NofState,NofState);          % Processes Noise Cov  Matrix
carposeVar = var(NoisyMeasurement(:,1));
carVelVar  = var(NoisyMeasurement(:,2));
pendAngVar = var(NoisyMeasurement(:,3));
pendVelVar = var(NoisyMeasurement(:,4));
R = eye(NofState);   % Measurement Noise Cov Matrix
R(1,1) = R(1,1)*carposeVar;
R(2,2) = R(2,2)*carVelVar;
R(3,3) = R(3,3)*pendAngVar;
R(4,4) = R(4,4)*pendVelVar;
XKalman(:,1) = zeros(NofState,1);         % Kalman State Initial Value

%% Kalman Filter Loop
for i=1:length(t)-1
    % Time Update (Prediction) Phase
    XKalman(:,i) = Phi*XKalman(:,i);   % State Extrapolation Equation
    P = Phi*P*Phi'+Q;                    % Uncertainty Extrapolation Equation
    % Measurement Update (Correction) Phase
    K = P*C'*inv(C*P*C'+R);           % Compute Kalman Gain!
    XKalman(:,i+1) = XKalman(:,i)+K*(NoisyMeasurement(:,i)-C*XKalman(:,i));  % Update State with Measurement & Kalman Gain
    P = (eye(NofState,NofState)-K*C)*P*(eye(NofState,NofState)-K*C)'+K*R*K';  % Update Estimation Uncertainty
end


%% PLOT Simulation
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for i=1:length(t)-1
    if mod(i,5e1) == 0
        clf
        PlotPend(state,i,L,t,NoisyMeasurement,XKalman)
        drawnow
    end
end

function dx = PendStateSpace(x,m,M,L,g,d,u)
D = m*L*L*(M+m*(1-cos(x(3))^2));
dx(1,1) = x(2);
dx(2,1) = (1/D)*(-m^2*L^2*g*cos(x(3))*sin(x(3)) + m*L^2*(m*L*x(4)^2*sin(x(3)) - d*x(2))) + m*L*L*(1/D)*u;
dx(3,1) = x(4);
dx(4,1) = (1/D)*((m+M)*m*g*L*sin(x(3)) - m*L*cos(x(3))*(m*L*x(4)^2*sin(x(3)) - d*x(2))) - m*L*cos(x(3))*(1/D)*u;
end