clc, clear, close all;
% Ball Video Creator Script
% Written By: Rasit
% Date: 11-Aug-2024
%% 
videoFile = 'LinearTraj.avi'; % Video file name
% videoFile = 'InfTraj.avi'; % Video file name
writerObj = VideoWriter(videoFile, 'Uncompressed AVI'); % VideoWriter object
writerObj.FrameRate = 30; % Frame rate
open(writerObj);

%% Infinity Trajectory
% t = linspace(0,2*pi,200)'; 
% NoD = size(t,2);       % Number Of Data
% x = 5e2*sin(t)+500;   
% y = 5e2*sin(2*t)+500;
% dt = t(2)-t(1);
% xdot = diff(x)/dt; % dx/dt
% ydot = diff(y)/dt; % dy/dt
% xdot = [xdot; xdot(end)];
% ydot = [ydot; ydot(end)];
% 
% BallTraj.x = x;
% BallTraj.xnoise = x+randn(NoD,1);
% BallTraj.xd = xdot;
% BallTraj.y = y;
% BallTraj.ynoise = y+randn(NoD,1);
% BallTraj.yd = ydot;
% save BallTrajInf.mat BallTraj

%% Linear Trajectory
Ts = 1/30;               % Sampling Time  Ms
SimTime = 10;        % Simulation Time sec
t = 0:Ts:SimTime;   % Time Vector
NoD = size(t,2);       % Number Of Data
xini = 100;      % Initial Position (X)
xfinal =  900;  % Final Position (X)
vx = ((xfinal-xini)/SimTime)*1.5;    % Velocity (X) Meter/sec
[x,xd,xdd] = SCurveTrajectory(xini,xfinal,SimTime,vx,Ts);

yini = 100;      % Initial Position (Y)
yfinal = 900;    % Final Position (Y)
vy = ((yfinal-yini)/SimTime)*1.5;    % Velocity (Y) Meter/sec (Manual Velocity)
[y,yd,ydd] = SCurveTrajectory(yini,yfinal,SimTime,vy,Ts);

BallTraj.x = x;
BallTraj.xnoise = x+randn(NoD,1);
BallTraj.xd = xd;
BallTraj.xdd = xdd;
BallTraj.y = y;
BallTraj.ynoise = y+randn(NoD,1);
BallTraj.yd = yd;
BallTraj.ydd = ydd;
save BallTrajLinear.mat BallTraj


%% Video Record
ballDiameter = 30;
for k = 1:length(t)
    fig = figure('Visible', 'off', 'Position', [0, 0, 1920, 1080]); % High resolution figure
    set(gca, 'Color', 'w');
    plot(BallTraj.xnoise(k), BallTraj.ynoise(k), 'ro', 'MarkerSize', ballDiameter, 'MarkerFaceColor', 'r'); % Positioning the ball in the center of the screen
    xlim([0, 1920]);
    ylim([0, 1080]);
    axis off;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    close(fig);
end

close(writerObj);
