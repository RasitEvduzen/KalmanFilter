clc,clear all,close all;
% Multi Dimension Linear Kalman Filter
% 2D Vehicle Pose Estimation with constant acceleration motion model
% Written By: Rasit Evduzen
% 20-Apr-2023
%% Create Data 
Ts = 1e-1;
TrajTime = 10;
NoD = 50;
Mean = 0;
Posestd  = 0.1;
Tstep = linspace(0,pi,TrajTime/Ts);
xunit = 3*cos(Tstep);
yunit = 3*sin(Tstep);
Real = [linspace(3,3,size(xunit,2))  xunit linspace(-3,-3,size(xunit,2)) -xunit;
    linspace(-4,0,size(yunit,2)) yunit linspace(0,-4,size(yunit,2))  -yunit-4] ;

Meas = [Real(1,:)+random('Normal',Mean,Posestd,size(Real,2),1)';
    Real(2,:)+random('Normal',Mean,Posestd,size(Real,2),1)'];
Measurement = Meas;  % Noisy Trajectory
RealTraj = Real;     % Real Trajectory
%% Algorithm Model
F = [1 Ts .5*Ts^2 0 0 0;
    0 1 Ts 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 Ts .5*Ts^2;
    0 0 0 0 1 Ts;
    0 0 0 0 0 1]; % System Matrix (State Space)

H = [1 0 0 0 0 0;
    0 0 0 1 0 0]; % Measurement Matrix  (Output Selection)

% X -> pose, vel, acc  |  Y -> pose, vel, acc [M, M/Sn, M/Sn^2];
NofState = size(F,1);                  % Number Of State
% ---------------------- P~Q~R Matrix Random value --------------------
% High Kalman Gain: if R >> P => P/(P+R) ~= K ~= 0 (Algorithm belief measurement)
% Low Kalman Gain:  if R << P => P/(P+R) ~= K ~= 1 (Algorithm belief kalman model)
P = 1e-3*eye(NofState,NofState);          % High Estimate Uncertainty
Q = 1e-3*eye(NofState,NofState);          % Processes Noise Cov
R = Posestd^2*eye(size(H,1),size(H,1));   % Measurement Noise Cov
XKalman(:,1) = [RealTraj(1) 0 0 RealTraj(2) 0 0]';         
%% Kalman Filter
for i=1:size(Measurement,2)
    % Time Update (Prediction) Phase
    XKalman(:,i) = F * XKalman(:,i);   % State Extrapolation Equation
    P = F*P*F' + Q;                    % Uncertainty Extrapolation Equation
    % Measurement Update (Correction) Phase
    K = P*H'*inv(H*P*H'+R);           % Compute Kalman Gain!
    XKalman(:,i+1) = XKalman(:,i) + K * (Measurement(:,i) - H*XKalman(:,i));  % Update State with Measurement & Kalman Gain
    P = (eye(NofState,NofState)-K*H)*P*(eye(NofState,NofState)-K*H)' + K*R*K';  % Update Estimation Uncertainty
end

%% Plot Data
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
TimeVec = linspace(0,TrajTime,size(XKalman,2));
for i=1:5:size(Measurement,2)
    clf
    subplot(2,3,[1 4])
    plot(RealTraj(1,:),RealTraj(2,:),'g','LineWidth',3),grid on,hold on
    plot(Measurement(1,1:i),Measurement(2,1:i),'r','LineWidth',1)
    plot(XKalman(1,1:i),XKalman(4,1:i),'k','LineWidth',2),axis equal,axis([-5 5 -8 5])
    legend("Real Trajectory", "Measurement Trajectory","Kalman Estimation")
    title("Vehicle Pose Estimation")
    xlabel("X [M]"), ylabel("Y [M]")

    subplot(2,3,2)
    plot(TimeVec(1:i),XKalman(1,1:i),'k','LineWidth',2),grid on, hold on
    plot(TimeVec(1:i),Measurement(1,1:i),'r','LineWidth',1)
    xlabel("Time [Sn]"),ylabel("Px [Mp]")
    legend("Kalman Trajectory","Measurement Trajectory",'Location','northeast')

    subplot(2,3,5)
    plot(TimeVec(1:i),XKalman(4,1:i),'k','LineWidth',2),grid on, hold on
    plot(TimeVec(1:i),Measurement(2,1:i),'r','LineWidth',1)
    xlabel("Time [Sn]"),ylabel("Py [M]")
    legend("Kalman Trajectory","Measurement Trajectory",'Location','northeast')

    subplot(2,3,3)
    plot(TimeVec(1:i),XKalman(2,1:i),'k','LineWidth',2),grid on, hold on
    plot(TimeVec(1:i),XKalman(3,1:i),'r','LineWidth',2)
    xlabel("Time [Sn]"),ylabel("Vx & Ax [M/Sn & M/Sn^2]")
    legend("Estimated Velocity","Estimated Acceleration",'Location','northeast')

    subplot(2,3,6)
    plot(TimeVec(1:i),XKalman(5,1:i),'k','LineWidth',2),grid on, hold on
    plot(TimeVec(1:i),XKalman(6,1:i),'r','LineWidth',2)
    xlabel("Time [Sn]"),ylabel("Vy & Ay [M/Sn & M/Sn^2]")
    legend("Estimated Velocity","Estimated Acceleration",'Location','northeast')

    drawnow
end
