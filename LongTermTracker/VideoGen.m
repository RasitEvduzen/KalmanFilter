clc, clear, close all;
% Ball Video Creator Script
% Written By: Rasit
% Date: 02-Apr-2026
%%
videoFile = 'LinearTraj.avi';
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

BallTraj.noise = 1;
BallTraj.x      = x;   BallTraj.xnoise = x + BallTraj.noise*randn(NoD,1);
BallTraj.xd     = xd;  BallTraj.xdd    = xdd;
BallTraj.y      = y;   BallTraj.ynoise = y + BallTraj.noise*randn(NoD,1);
BallTraj.yd     = yd;  BallTraj.ydd    = ydd;

% Occlusion Box centered on trajectory, wide on x-axis
occ_box = [300 350 700 650];                 % [x1 y1 x2 y2] world coords
BallTraj.occ_box = occ_box;
save BallTrajLinear.mat BallTraj

%% Video Record
fig          = figure('Units','pixels','Position',[0 0 1920 1080],'Color','w', ...
                      'MenuBar','none','ToolBar','none');
ax           = gca;
ax.Units     = 'pixels';
ax.Position  = [0 0 1920 1080];
ballDiameter = 30;

for k = 1:length(t)
    cla, hold on

    % Draw occlusion box FIRST (background)
    fill([occ_box(1) occ_box(3) occ_box(3) occ_box(1)], ...
         [occ_box(2) occ_box(2) occ_box(4) occ_box(4)], ...
         [0.4 0.4 0.4], 'EdgeColor', 'none');

    % Check if ball is inside occlusion box
    in_occ = BallTraj.xnoise(k) >= occ_box(1) && BallTraj.xnoise(k) <= occ_box(3) && ...
             BallTraj.ynoise(k) >= occ_box(2) && BallTraj.ynoise(k) <= occ_box(4);

    % Draw ball ON TOP of box (hidden only when inside)
    if ~in_occ
        plot(BallTraj.xnoise(k), BallTraj.ynoise(k), 'ro', ...
             'MarkerSize', ballDiameter, 'MarkerFaceColor', 'r');
    end

    axis off
    xlim([0 1920]), ylim([0 1080])
    writeVideo(writerObj, getframe(gca));
end
close(writerObj);
