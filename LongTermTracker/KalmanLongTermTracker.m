clc, clear, close all;
% Kalman Filter Based Tracked Object Velocity Estimation
% Written By: Rasit
% Date: 02-Apr-2026

%% Load Data
videoFile = 'LinearTraj.avi';
load BallTrajLinear.mat

videoReader = VideoReader(videoFile);
Ts          = 1 / videoReader.FrameRate;
frameH      = videoReader.Height;
frameW      = videoReader.Width;
xlim_gen    = 1920;
ylim_gen    = 1080;
occ_box = BallTraj.occ_box;                 % [300 350 700 650]

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
R   = BallTraj.noise^2 * eye(2);            % Measurement Noise Cov Matrix
P_cap = 1e4 * eye(NoS);                     % P cap during occlusion (prevents spike on exit)

x_est       = [0; 0; 0; 0];                 % Initial State
initialized = false;                         % First detection flag

%% Storage
poseX = [];  poseY = [];
velX  = [];  velY  = [];
time  = [];
trackedCoordinates = [];

%% Kalman Loop
figure('Units','pixels','Position',[0 0 1920 1080],'Color','w','MenuBar','none','ToolBar','none');
ax          = gca;
ax.Units    = 'pixels';
ax.Position = [0 0 1920 1080];
SimPlot = 1;
idx = 0;
while hasFrame(videoReader)
    frame = readFrame(videoReader);
    idx   = idx + 1;

    % Detect red ball
    redMask  = frame(:,:,1) > 150 & frame(:,:,2) < 100 & frame(:,:,3) < 100;
    stats    = regionprops(redMask, 'BoundingBox', 'Centroid');
    detected = ~isempty(stats);

    % Time Update (Prediction)
    if initialized
        x_pred = A * x_est;                              % State Propagation
        P_pred = A * P * A' + Q;                         % Uncertainty Propagation
    end

    if detected
        centroid    = stats(1).Centroid;
        boundingBox = stats(1).BoundingBox;
        trackedCoordinates = [trackedCoordinates; centroid];

        z = [centroid(1); frameH - centroid(2)];

        % Initialize from first detection
        if ~initialized
            x_est       = [z(1); z(2); 0; 0];
            initialized = true;
            x_pred      = x_est;
            P_pred      = P;
        end

        %  Collision check 
        in_occ = z(1) >= occ_box(1) && z(1) <= occ_box(3) && ...
                 z(2) >= occ_box(2) && z(2) <= occ_box(4);

        if in_occ
            % Inside box: prediction only, cap P to prevent exit spike
            x_est = x_pred;
            P     = min(P_pred, P_cap);
        else
            % Outside box: normal Kalman correction
            K     = P_pred * H' / (H * P_pred * H' + R); % Compute Kalman Gain!
            x_est = x_pred + K * (z - H * x_pred);        % Update State with Measurement & Kalman Gain
            P     = (eye(NoS)-K*H)*P_pred*(eye(NoS)-K*H)' + K*R*K'; % Joseph Form
        end

    elseif initialized
        x_est = x_pred;
        P     = min(P_pred, P_cap);
    end

    if initialized
        poseX(end+1) = x_est(1);
        poseY(end+1) = x_est(2);
        velX(end+1)  = x_est(3);
        velY(end+1)  = x_est(4);
        time(end+1)  = Ts * idx;

        if SimPlot == 1
            pred_img_x = x_est(1) / xlim_gen * frameW;
            pred_img_y = (1 - x_est(2)/ylim_gen) * frameH;

            if detected
                frame = insertShape(frame, 'Rectangle', boundingBox, 'Color','blue','LineWidth',3);
                frame = insertMarker(frame, centroid, 'x', 'Color','blue','Size',20);
                arrow_x = centroid(1);              
                arrow_y = centroid(2);
            else
                arrow_x = pred_img_x;               
                arrow_y = pred_img_y;
            end

            textStr = sprintf('x: %.1f px  y: %.1f px\nVx: %.2f px/s  Vy: %.2f px/s', ...
                              x_est(1), x_est(2), x_est(3), x_est(4));
            frame = insertText(frame, [30 30], textStr, ...
                              'FontSize',18,'TextColor','black','BoxColor','white','BoxOpacity',0.6);

            imshow(frame), hold on
            if ~isempty(trackedCoordinates)
                plot(trackedCoordinates(:,1), trackedCoordinates(:,2), 'k', 'LineWidth',2)
            end
            quiver(arrow_x, arrow_y, x_est(3), -x_est(4), 'g', 'LineWidth',2, 'MaxHeadSize',2, 'ShowArrowHead','on', 'AutoScale','off')
            drawnow
        end
    end
end

%% Plot Results
t_ref = linspace(0, max(time), length(BallTraj.x));

figure('units','normalized','outerposition',[0 0 1 1],'color','w')

subplot(2,2,1)
plot(time, poseX, 'r', 'LineWidth',3), hold on, grid on
plot(t_ref, BallTraj.x, 'k--', 'LineWidth',2)
legend('KF Estimate','True','Location','northwest')
xlabel('Time [s]'), ylabel('X [px]'), title('X Position')

subplot(2,2,2)
plot(time, poseY, 'r', 'LineWidth',3), hold on, grid on
plot(t_ref, BallTraj.y, 'k--', 'LineWidth',2)
legend('KF Estimate','True','Location','northwest')
xlabel('Time [s]'), ylabel('Y [px]'), title('Y Position')

subplot(2,2,3)
plot(time, velX, 'r', 'LineWidth',3), hold on, grid on
plot(t_ref, BallTraj.xd, 'k--', 'LineWidth',2)
legend('KF Estimate','True','Location','northwest')
xlabel('Time [s]'), ylabel('Vx [px/s]'), title('X Velocity')

subplot(2,2,4)
plot(time, velY, 'r', 'LineWidth',3), hold on, grid on
plot(t_ref, BallTraj.yd, 'k--', 'LineWidth',2)
legend('KF Estimate','True','Location','northwest')
xlabel('Time [s]'), ylabel('Vy [px/s]'), title('Y Velocity')