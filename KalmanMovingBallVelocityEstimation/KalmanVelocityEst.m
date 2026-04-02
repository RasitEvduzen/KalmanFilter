clc, clear, close all;
% Kalman Filter Based Tracked Object Velocity Estimation
% Written By: Rasit
% Date: 11-Aug-2024
%%
%% Load Data
videoFile = 'LinearTraj.avi';
load BallTrajLinear.mat

videoReader = VideoReader(videoFile);
Ts          = 1 / videoReader.FrameRate;
frameH      = videoReader.Height;            % 1080 px

%% Kalman Filter Initialization
NoS = 4;                                     % [px, py, vx, vy]
A   = [1 0 Ts 0;
    0 1 0  Ts;
    0 0 1  0;
    0 0 0  1];                            % State Transition Matrix
H   = [1 0 0 0;
    0 1 0 0];                             % Observation Matrix

% Low  Kalman Gain: if R >> P => K~=0 (trust model)
% High Kalman Gain: if R << P => K~=1 (trust measurement)
P   = 1e0 * eye(NoS);                       % Error Noise Cov Matrix
Q   = 1e4 * [Ts^4/4 0      Ts^3/2 0;
    0      Ts^4/4 0      Ts^3/2;
    Ts^3/2 0      Ts^2   0;
    0      Ts^3/2 0      Ts^2];    % Process Noise Cov Matrix
R   = BallTraj.noise* eye(2);                         % Measurement Noise Cov Matrix

x_est       = [0; 0; 0; 0];                 % Initial State
initialized = false;                          % First detection flag

%% Storage
poseX = [];  poseY = [];
velX  = [];  velY  = [];
time  = [];
trackedCoordinates = [];

%% Kalman Loop
figure('Units','pixels','Position',[0 0 1920 1080],'Color','w','MenuBar','none','ToolBar','none');
idx = 0;
SimPlot = 1;
while hasFrame(videoReader)
    frame = readFrame(videoReader);
    idx   = idx + 1;

    % Detect red ball
    redMask = frame(:,:,1) > 150 & frame(:,:,2) < 100 & frame(:,:,3) < 100;
    stats   = regionprops(redMask, 'BoundingBox', 'Centroid');

    if ~isempty(stats)
        centroid    = stats(1).Centroid;
        boundingBox = stats(1).BoundingBox;
        trackedCoordinates = [trackedCoordinates; centroid];

        % Coordinate transform: pixel → world (y-flip)
        z = [centroid(1); frameH - centroid(2)];

        % Initialize from first detection
        if ~initialized
            x_est       = [z(1); z(2); 0; 0];
            initialized = true;
        end

        % Time Update (Prediction) Phase
        x_pred = A * x_est;
        P_pred = A * P * A' + Q;                         % Uncertainty Propagation

        % Measurement Update (Correction) Phase
        K     = P_pred * H' / (H * P_pred * H' + R);    % Compute Kalman Gain!
        x_est = x_pred + K * (z - H * x_pred);           % Update State with Measurement & Kalman Gain
        P     = (eye(NoS)-K*H)*P_pred*(eye(NoS)-K*H)' + K*R*K'; % Joseph Form

        % Store (world coordinates)
        poseX(end+1) = x_est(1);
        poseY(end+1) = x_est(2);
        velX(end+1)  = x_est(3);
        velY(end+1)  = x_est(4);
        time(end+1)  = Ts * idx;
        if SimPlot == 1
            % Visualization (image coordinates)
            frame = insertShape(frame, 'Rectangle', boundingBox, 'Color','blue','LineWidth',3);
            frame = insertMarker(frame, centroid, 'x', 'Color','yellow','Size',20);
            textStr = sprintf('x: %.1f  y: %.1f\nVx: %.2f  Vy: %.2f', ...
                x_est(1), x_est(2), x_est(3), x_est(4));
            frame = insertText(frame, [boundingBox(1)+40, boundingBox(2)+40], textStr, ...
                'FontSize',18,'TextColor','black','BoxColor','white','BoxOpacity',0.6);
            imshow(frame), hold on
            title('Kalman Filter Based Ball Velocity Estimation')
            plot(centroid(1), centroid(2), 'k.', 'MarkerSize',5)
            plot(trackedCoordinates(:,1), trackedCoordinates(:,2), 'k', 'LineWidth',2)
            quiver(centroid(1), centroid(2), x_est(3), -x_est(4), 'g', 'LineWidth',1.5, 'MaxHeadSize',1)
            drawnow
        end
    end
end

%% Plot Results
t_ref = linspace(0, max(time), length(BallTraj.x));

figure('units','normalized','outerposition',[0 0 1 1],'color','w')

subplot(2,2,1)
plot(time, poseX, 'r', 'LineWidth',2), hold on, grid on
plot(t_ref, BallTraj.x, 'k', 'LineWidth',2)
legend('KF Estimate','True','Location','northwest')
xlabel('Time [s]'), ylabel('X [pixel]'), title('X Position')

subplot(2,2,2)
plot(time, poseY, 'r', 'LineWidth',2), hold on, grid on
plot(t_ref, BallTraj.y, 'k', 'LineWidth',2)
legend('KF Estimate','True','Location','northwest')
xlabel('Time [s]'), ylabel('Y [pixel]'), title('Y Position')

subplot(2,2,3)
plot(time, velX, 'r', 'LineWidth',2), hold on, grid on
plot(t_ref, BallTraj.xd, 'k', 'LineWidth',2)
legend('KF Estimate','True','Location','northwest')
xlabel('Time [s]'), ylabel('Vx [pixel/s]'), title('X Velocity')

subplot(2,2,4)
plot(time, velY, 'r', 'LineWidth',2), hold on, grid on
plot(t_ref, BallTraj.yd, 'k', 'LineWidth',2)
legend('KF Estimate','True','Location','northwest')
xlabel('Time [s]'), ylabel('Vy [pixel/s]'), title('Y Velocity')