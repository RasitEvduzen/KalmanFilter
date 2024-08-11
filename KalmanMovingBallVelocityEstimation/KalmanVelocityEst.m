clc, clear, close all;
% Kalman Filter Based Tracked Object Velocity Estimation
% Written By: Rasit
% Date: 11-Aug-2024
%% Load the trajectory data

% Linear Trajectory
videoFile = 'LinearTraj.avi';  
load BallTrajLinear.mat

% Infinity Trajectory
% videoFile = 'InfTraj.avi';     
% load BallTrajInf.mat


videoReader = VideoReader(videoFile);
trackedCoordinates = [];
velocityX = [];
velocityY = [];
time = [];
Ts = 1 / videoReader.FrameRate; % Sampling Periode
idx = 0;  % Data Index

A = [1 0 Ts 0 0.5*Ts^2 0;
    0 1 0 Ts 0 0.5*Ts^2;
    0 0 1 0 Ts 0;
    0 0 0 1 0 Ts;
    0 0 0 0 1 0;
    0 0 0 0 0 1]; % State transition Matrix

H = [1 0 0 0 0 0;
    0 1 0 0 0 0]; % Measurement matrix (position measurement only)

P = 1e6*eye(6); % Initial covariance matrix
Q = 1e3*eye(6); % Process noise covariance matrix
R = 1e0*eye(2); % Measurement noise covariance matrix
x_est = [0; 0; 0; 0; 0; 0]; % Initial estimate (position, velocity, acceleration)

figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w')
sim = "on";
while hasFrame(videoReader)
    frame = readFrame(videoReader);
    idx = idx + 1;
    % Detect the red ball (color thresholding)
    redChannel = frame(:,:,1);
    greenChannel = frame(:,:,2);
    blueChannel = frame(:,:,3);
    redMask = redChannel > 150 & greenChannel < 100 & blueChannel < 100;
    stats = regionprops(redMask, 'BoundingBox', 'Centroid');
    if ~isempty(stats)
        boundingBox = stats(1).BoundingBox;
        centroid = stats(1).Centroid;
        trackedCoordinates = [trackedCoordinates; centroid];

        % Kalman filter 
        z = centroid'; % Measurement (current position)
        % Prediction Phase
        x_pred = A * x_est;
        P_pred = A * P * A' + Q;

        % Correction Phase
        K = P_pred * H' / (H * P_pred * H' + R); % Kalman gain
        x_est = x_pred + K * (z - H * x_pred); 
        P = (eye(6) - K * H) * P_pred; 
        estimatedVelocity = x_est(3:4);  % Velocity State

        velocityX = [velocityX; estimatedVelocity(1)];
        velocityY = [velocityY; estimatedVelocity(2)];
        time = [time; Ts * length(time)];

        if sim == "on"
            frame = insertShape(frame, 'Rectangle', boundingBox, 'Color', 'blue', 'LineWidth', 3);
            frame = insertMarker(frame, centroid, 'x', 'Color', 'yellow', 'Size', 20);

            position = [boundingBox(1)+40, boundingBox(2)+40]; 
            textStr = sprintf('x: %.1f, y: %.1f\nVx: %.2f, Vy: %.2f', centroid(1), centroid(2), estimatedVelocity(1), estimatedVelocity(2));
            frame = insertText(frame, position, textStr, 'FontSize', 18, 'TextColor', 'black', 'BoxColor', 'white', 'BoxOpacity', 0.6);
            
            imshow(frame), hold on
            title('Kalman Filter Based Ball Velocity Estimation');
            xlabel('X Coordinate'),ylabel('Y Coordinate');            
            plot(centroid(1), centroid(2), 'k.', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
            plot(trackedCoordinates(1:end,1),trackedCoordinates(1:end,2), 'k',LineWidth=3)
            quiver(centroid(1), centroid(2), estimatedVelocity(1), estimatedVelocity(2), 'g', 'LineWidth', .5, 'MaxHeadSize', 1);
            drawnow;
        end
    end
end

% Plot Velocity
figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w')
subplot(211);
plot(time, velocityX, 'r', 'DisplayName', 'Vx'),hold on
plot(linspace(0,max(time),size(BallTraj.xd,1)), BallTraj.xd, 'k', 'DisplayName', 'Real Traj Vx')
title('Velocity Over Time');xlabel('Time (s)');ylabel('Velocity (pixels/frame)');grid on;
legend("Estimated Traj","Real Vel Traj")

subplot(212);
plot(time, -velocityY, 'r', 'DisplayName', 'Vy '),hold on
plot(linspace(0,max(time),size(BallTraj.yd,1)),BallTraj.yd, 'k', 'DisplayName', 'Real Traj Vy')
title('Velocity Over Time');xlabel('Time (s)');ylabel('Velocity (pixels/frame)');grid on;
legend("Estimated Traj","Real Vel Traj")
