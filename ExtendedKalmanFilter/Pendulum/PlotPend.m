function [] = PlotPend(state,k,L,t,NoisyMeasurement,XKalman)
x = state(k,1);
th = state(k,3);
% dimensions
W  = 1.5;  % cart width
H  = .5;   % cart height
wr = .5;             % wheel radius
mr = 0.6708;     % mass radius
% positions
y = wr/2+H/2; % cart vertical position
pendx = x + L*sin(th);
pendy = y - L*cos(th);

% Plot Pendulum
subplot(3,2,[1,2])
yline(0,'k--','LineWidth',2), hold on
rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[0.4940 0.1840 0.5560],'LineWidth',1.5); % Draw cart
rectangle('Position',[x-.9*W/2,0,wr,wr],'Curvature',1,'FaceColor',[1 1 0],'LineWidth',1.5);    % Draw wheel
rectangle('Position',[x+.9*W/2-wr,0,wr,wr],'Curvature',1,'FaceColor',[1 1 0],'LineWidth',1.5); % Draw wheel
plot([x pendx],[y pendy],'k','LineWidth',2); % Draw pendulum
rectangle('Position',[pendx-mr/2,pendy-mr/2,mr,mr],'Curvature',1,'FaceColor',[1 0.1 .1],'LineWidth',1.5);
axis equal,axis([-12 12 -3 3]),grid on
title('Inverted Pendulum Extended Kalman Filter')

subplot(323)
plot(t,NoisyMeasurement(1,:),Color=[.7 .7 .7],LineWidth=2); hold on
plot(t,XKalman(1,:)','k',LineWidth=2);
plot(t,state(:,1),'r--',LineWidth=2)
legend('Noisy Data','EKF Data','Real Data','Location','NorthEast');
xlabel('Time'),ylabel('Pendulum Linear Pose'),grid on

subplot(324)
plot(t,NoisyMeasurement(2,:),Color=[.7 .7 .7],LineWidth=2); hold on
plot(t,XKalman(2,:)','k',LineWidth=2);
plot(t,state(:,2),'r--',LineWidth=2)
legend('Noisy Data','EKF Data','Real Data','Location','NorthEast');
xlabel('Time'),ylabel('Pendulum Linear Velocity'),grid on

subplot(325)
plot(t,NoisyMeasurement(3,:),Color=[.7 .7 .7],LineWidth=2); hold on
plot(t,XKalman(3,:)','k',LineWidth=2);
plot(t,state(:,3),'r--',LineWidth=2)
legend('Noisy Data','EKF Data','Real Data','Location','NorthEast');
xlabel('Time'),ylabel('Pendulum Angle'),grid on

subplot(326)
plot(t,NoisyMeasurement(4,:),Color=[.7 .7 .7],LineWidth=2); hold on
plot(t,XKalman(4,:)','k',LineWidth=2);
plot(t,state(:,4),'r--',LineWidth=2)
legend('Noisy Data','EKF Data','Real Data','Location','NorthEast');
xlabel('Time'),ylabel('Pendulum Angular Vel'),grid on
end