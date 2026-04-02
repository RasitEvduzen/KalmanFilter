clc, clear, close all;
% Ball Video Creator Script
% Written By: Rasit
% Date: 11-Aug-2024
%%
videoFile = 'LinearTraj.avi';
% videoFile = 'InfTraj.avi';
writerObj = VideoWriter(videoFile, 'Uncompressed AVI');
writerObj.FrameRate = 30;
open(writerObj);

%% Linear Trajectory
Ts      = 1/30;
SimTime = 10;
t       = 0:Ts:SimTime;
NoD     = size(t, 2);

xini = 100;  xfinal = 900;
vx   = ((xfinal-xini)/SimTime) * 1.5;
[x, xd, xdd] = SCurveTrajectory(xini, xfinal, SimTime, vx, Ts);

yini = 100;  yfinal = 900;
vy   = ((yfinal-yini)/SimTime) * 1.5;
[y, yd, ydd] = SCurveTrajectory(yini, yfinal, SimTime, vy, Ts);

noise = 5;                                           % Measurement Noise Std [pixel]
clc, clear, close all;
% Ball Video Creator Script
% Written By: Rasit
% Date: 11-Aug-2024
%%
videoFile = 'LinearTraj.avi';
% videoFile = 'InfTraj.avi';
writerObj = VideoWriter(videoFile, 'Uncompressed AVI');
writerObj.FrameRate = 30;
open(writerObj);

%% Linear Trajectory
Ts      = 1/30;
SimTime = 10;
t       = 0:Ts:SimTime;
NoD     = size(t, 2);

xini = 100;  xfinal = 900;
vx   = ((xfinal-xini)/SimTime) * 1.5;
[x, xd, xdd] = SCurveTrajectory(xini, xfinal, SimTime, vx, Ts);

yini = 100;  yfinal = 900;
vy   = ((yfinal-yini)/SimTime) * 1.5;
[y, yd, ydd] = SCurveTrajectory(yini, yfinal, SimTime, vy, Ts);

BallTraj.noise = 1;                                           % Measurement Noise Std [pixel]
BallTraj.x      = x;
BallTraj.xnoise = x + BallTraj.noise*randn(NoD,1);
BallTraj.xd     = xd;  BallTraj.xdd = xdd;
BallTraj.y      = y;
BallTraj.ynoise = y + BallTraj.noise*randn(NoD,1);
BallTraj.yd     = yd;  BallTraj.ydd = ydd;
save BallTrajLinear.mat BallTraj

%% Video Record
fig = figure('Units','pixels','Position',[0 0 1920 1080],'Color','w', ...
             'MenuBar','none','ToolBar','none');   
ax            = gca;
ax.Units      = 'pixels';
ax.Position   = [0 0 1920 1080];                 % Exact pixel fill, no margins
ballDiameter = 30;

for k = 1:length(t)
    cla
    plot(BallTraj.xnoise(k), BallTraj.ynoise(k), 'ro', ...
         'MarkerSize', ballDiameter, 'MarkerFaceColor', 'r');
    axis off
    xlim([0 1920]), ylim([0 1080])
    writeVideo(writerObj, getframe(gca));
end
close(writerObj);x      = x;
BallTraj.xnoise = x + noise*randn(NoD,1);
BallTraj.xd     = xd;  BallTraj.xdd = xdd;
BallTraj.y      = y;
BallTraj.ynoise = y + noise*randn(NoD,1);
BallTraj.yd     = yd;  BallTraj.ydd = ydd;
save BallTrajLinear.mat BallTraj

%% Video Record
fig = figure('Units','pixels','Position',[0 0 1920 1080],'Color','w', ...
             'MenuBar','none','ToolBar','none');   % Remove toolbar → exact canvas size
ax            = gca;
ax.Units      = 'pixels';
ax.Position   = [0 0 1920 1080];                 % Exact pixel fill, no margins
ballDiameter = 30;

for k = 1:length(t)
    cla
    plot(BallTraj.xnoise(k), BallTraj.ynoise(k), 'ro', ...
         'MarkerSize', ballDiameter, 'MarkerFaceColor', 'r');
    axis off
    xlim([0 1920]), ylim([0 1080])
    writeVideo(writerObj, getframe(gca));
end
close(writerObj);