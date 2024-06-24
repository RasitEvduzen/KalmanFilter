clc,clear all,close all;
% Written By: Rasit Evduzen
% 24-May-2023
% Quasi Stationary North Finding 
%% Load Data
global W_ie Eq_rad Pol_rad Flat Eccen mu Deg_to_rad Rad_to_deg
Deg_to_rad = pi/180;   % Degree to Radian convert constant
Rad_to_deg = 180/pi;   % Radian to Degree convert constant
% --------------------------- Load Data ----------------------------------------
addpath("Functions");
% ----------------------- Load Data ---------------------------
load RawData.csv
AlgorithmFreq = 200; % Sampling Frequency [Hz]
time_vec = RawData(:,1)/AlgorithmFreq;
GX = RawData(:,2)*Deg_to_rad;
GY = RawData(:,3)*Deg_to_rad;
GZ = RawData(:,4)*Deg_to_rad;
AX = RawData(:,5);
AY = RawData(:,6);
AZ = RawData(:,7);
clear RawData
% -------------------------------- Anomaly Delete --------------------------------
thmin = 1e-1;
thmax = 2e1;
[AX] = DataClean(AX,thmin,thmax);
[AY] = DataClean(AY,thmin,thmax);
[AZ] = DataClean(AZ,thmin,thmax);
[GX] = DataClean(GX,thmin,thmax);
[GY] = DataClean(GY,thmin,thmax);
[GZ] = DataClean(GZ,thmin,thmax);
% --------------------------- Algorithm Init Params ---------------------------
W_ie =  (360 / ((((23*60)+56)*60)+4)) * Deg_to_rad;  % Earth Angular Velocity [Rad/Sn]
Eq_rad = 6378137;        % Equatorial radius, R0  [M]    Standards -> WGS84
Pol_rad = 6356752.31425; % Polar      radius, Rp  [M]    Standards -> WGS84
Flat = 1/298.257223563;  % Flattening   No Unit          Standards -> WGS84
Eccen = 0.0818191908425; % Eccentricity No Unit          Standards -> WGS84
mu = 3.986004418E14;     % Earth grav const(m^3 s^-2)    Standards -> WGS84
L_b = 40.923061376 * Deg_to_rad;       % latitude  [rad] Teknopark
Lambda_b = 29.31453938 * Deg_to_rad;   % longitude [rad] Teknopark
H_b = 9.9355;           % altitude [m]
avarage_time = 10;      % Avarage Time Stop Value [Sn]
Ts = 1 / AlgorithmFreq; % Sampling Periode [Sn]
% --------------------------- Kalman Init Param ---------------------------
NofState = 17;     % Number Of State "Attitude(3x1), Pose(3x1), Velocity(3x1),  WSin(1x1), WCos(1x1), Ba(3x1), Bg(3x1)
model_order = 3;   % State Transition Taylor Series Truncation Order
H = zeros(NofState,NofState);   % Measurement Matrix
H(7:9,7:9) = -eye(3);            % Pose State Selection
V_eb_w = zeros(3,1);            % NED Frame Velocity Initial Value
R_eb_w = zeros(3,1);            % Pose Initial Value
W_ie_w = zeros(3,1);            % Earth Angular Velocity, wander frame resolution
F = zeros(NofState,NofState);   % Error State Space Matrix
Phi = zeros(NofState,NofState); % Error State Transition Matrix
wander_ang = 0;                 % Wander Angle Initial Value
wander_check = 1;               % Wander Correction Check
b_a = zeros(3,1);               % Accel Dynamic Bias
b_g = zeros(3,1);               % Gyro Dynamic Bias
Nominal_State_Vec = zeros(NofState,1);   % Nominal State Vector (Measurement Value and Mechanization)
Error_State_Vec = zeros(NofState,1);     % Error State Vector   (ESKF)
True_State_Vec = zeros(NofState,1);      % True State Vector = Nominal State Vector - Error State Vector
K = zeros(NofState,NofState);            % Kalman Gain Initialization

% ----------------------------- Algorithm Init ---------------------------------
avarage_sample = avarage_time*AlgorithmFreq;
Gx_mean = mean(GX(1:avarage_sample));
Gy_mean = mean(GY(1:avarage_sample));
Gz_mean = mean(GZ(1:avarage_sample));
Ax_mean = mean(AX(1:avarage_sample));
Ay_mean = mean(AY(1:avarage_sample));
Az_mean = mean(AZ(1:avarage_sample));

% ------------- Leveling and direct gyrocompassing calculations -------------
% [roll_ini pitch_ini] = accel2attitude(-[Ay_mean Az_mean Ax_mean]'); % roll, pitch angles
% [yaw_ini]  = stationaryGyrocompassing([Gy_mean Gz_mean Gx_mean], [roll_ini pitch_ini]);


% ---------------------- Triad Settings ----------------------
F_val = [-AY(avarage_sample:end) AZ(avarage_sample:end) -AX(avarage_sample:end)]';
W_val = [-GY(avarage_sample:end) GZ(avarage_sample:end) -GX(avarage_sample:end)]';
[roll_ini pitch_ini] = accel2attitude([-Ay_mean Az_mean -Ax_mean]'); % roll, pitch angles
[yaw_ini]  = stationaryGyrocompassing([-Gy_mean Gz_mean -Gx_mean], [roll_ini pitch_ini]);
Euler_Ang = [roll_ini pitch_ini yaw_ini]';          % Attitude Initial value
C_b_w = euler_to_ctm(roll_ini,pitch_ini,yaw_ini);   % Base to wander frame transformation initial val


NofData_Offline = size(F_val,2);
% ---------------------- P~Q~R Matrix value from Hakan --------------------
Patt = 1e-5; Pvel = 1e-4; Ppos = 1e-6; Pwander = 1e-2; Pba = 1e-4; Pbg = 1e-12;
P(1:3,1:3) = eye(3,3)*Patt;
P(4:6,4:6) = eye(3,3)*Pvel;
P(7:9,7:9) = eye(3,3)*Ppos;
P(10:11,10:11) = eye(2,2)*Pwander;
P(12:14,12:14) = eye(3,3)*Pba;
P(15:17,15:17) = eye(3,3)*Pbg; % Estimation Noise

Qatt = 1e-7; Qvel = 1e-7; Qpos = 1e-18; Qwander = 1e-8; Qba = 1e-9; Qbg = 1e-19;
Q(1:3,1:3) = eye(3,3)*Qatt;
Q(4:6,4:6) = eye(3,3)*Qvel;
Q(7:9,7:9) = eye(3,3)*Qpos;
Q(10:11,10:11) = eye(2,2)*Qwander;
Q(12:14,12:14) = eye(3,3)*Qba;
Q(15:17,15:17) = eye(3,3)*Qbg;  % Processs Noise

R = eye(NofState)*1e-9;        % Measurement Noise

% ---------------------- P~Q~R Matrix Random value --------------------
% High Kalman Gain: if R >> P => P/(P+R) ~= K ~= 0 (Algorithm belief measurement)
% Low Kalman Gain:  if R << P => P/(P+R) ~= K ~= 1 (Algorithm belief kalman model)
% P = eye(NofState)*1e-1;        % Estimation variance
% Q = eye(NofState)*1e-15;       % Processes Noise
% R = eye(NofState)*1e-9;        % Measurement Noise
yaw_offset = 360;  % Yaw angle offset for sign
True_State_Vec(1:3) = [roll_ini pitch_ini yaw_ini]';
% ---------------------- ESKF Loop Code --------------------
for j=1:NofData_Offline
        
            % ESKF Loop Code
    F_ib_b = F_val(:,j) - b_a(:,j);   % Accelerometer Data [M/Sn^2]  (Fx Fy Fz)
    W_ib_b = W_val(:,j) - b_g(:,j);   % Gyroscope Data [Rad/Sn]      (Wx Wy Wz)

%     if (j*Ts) == 20
%         W_ib_b(3) = (W_ib_b(3) + (1));  % Small Perturbation 
%     end

    % Navigation Mechanization Wander Frame
    [C_b_w(:,:,j+1),V_eb_w(:,j+1),R_eb_w(:,j+1),W_ie_w(:,j+1),wander_ang(j+1)] = Navigation_Update_Wander(Ts,L_b,Lambda_b,H_b,C_b_w(:,:,j),V_eb_w(:,j),R_eb_w(:,j),F_ib_b,W_ib_b,wander_ang(j));
    [roll,pitch,yaw] = ctm_to_euler(C_b_w(:,:,j+1));

%     if (j*Ts) == 25
%         wander_ang(j) = wander_ang(j) + (10 * Deg_to_rad);   % Small Perturbation 
%         % yaw = yaw + (1 * Deg_to_rad); 
%     end

    [Euler_Ang(:,j+1)] = [roll,pitch,yaw];
    Nominal_State_Vec(:,j+1) = [Euler_Ang(:,j+1); V_eb_w(:,j+1); R_eb_w(:,j+1); sin(wander_ang(j+1)); cos(wander_ang(j+1)); b_a(:,j); b_g(:,j)]; % Update!
    % ------------------------- Time Update (Prediction) Stage -------------------------
    F(:,:,j+1) = [-[skew_symmetric(W_ie_w(:,j+1))] zeros(3) zeros(3) [0 W_ie*cos(L_b) 0]' -[W_ie*cos(L_b) 0 0]' zeros(3) [C_b_w(:,:,j+1)];
                  -[skew_symmetric(C_b_w(:,:,j+1)*F_ib_b)]  zeros(3) zeros(3) zeros(3,1) zeros(3,1) [C_b_w(:,:,j+1)] zeros(3);
                    zeros(3) eye(3) zeros(3) zeros(3,1) zeros(3,1) zeros(3) zeros(3);
                    zeros(1,3) zeros(1,3) zeros(1,3) zeros(1) zeros(1) zeros(1,3) zeros(1,3);
                    zeros(1,3) zeros(1,3) zeros(1,3) zeros(1) zeros(1) zeros(1,3) zeros(1,3);
                    zeros(3) zeros(3) zeros(3) zeros(3,1) zeros(3,1) zeros(3) zeros(3);
                    zeros(3) zeros(3) zeros(3) zeros(3,1) zeros(3,1) zeros(3) zeros(3)];

    Phi(:,:,j+1) = state_transition(F(:,:,j+1),Ts,model_order);             % F --> Phi (State Transition Matrix) Power Series Expansion
%     Error_State_Vec(:,j) = Phi(:,:,j+1) * Error_State_Vec(:,j);             % State Extrapolation Equation
%     Q = (Q*Ts) + (0.5*F(:,:,j+1)*Q) + (0.5*Q*F(:,:,j+1)');                  % Measurement noise discrization
    P(:,:,j+1) = Phi(:,:,j+1)*P(:,:,j)*Phi(:,:,j+1)'+Q;                                     % Error State Covariance

    % ------------------------- Measurement Update (Correction) Stage -------------------------
    K(:,:,j+1) = P(:,:,j+1)*H'*inv(H*P(:,:,j+1)*H'+R);                                        % Compute Kalman Gain!
    Innovation = -[zeros(6,1); Nominal_State_Vec(7:9,j+1); zeros(8,1)];
    Error_State_Vec(:,j+1) = Error_State_Vec(:,j)+K(:,:,j+1)*(Innovation);    % Update State with Measurement & Kalman Gain
    P(:,:,j+1) = (eye(NofState,NofState)-K(:,:,j+1)*H)*P(:,:,j+1)*(eye(NofState,NofState)-K(:,:,j+1)*H)'+K(:,:,j+1)*R*K(:,:,j+1)';         % Update Error State Covariance

    % ---------- True State = Nominal State - Error State ----------
    True_State_Vec(:,j+1) = Nominal_State_Vec(:,j+1) - Error_State_Vec(:,j+1);

    % ------------------------- Close Loop Correction Phase -------------------------
%     b_a(:,j+1) = b_a(:,j) + Error_State_Vec(12:14,j+1);  % Close loop Accelerometer Dynamic Bias correction
%     b_g(:,j+1) = b_g(:,j) + Error_State_Vec(15:end,j+1); % Close loop Gyro Dynamic Bias correction    
    b_a(:,j+1) = Error_State_Vec(12:14,j+1);  % Close loop Accelerometer Dynamic Bias correction
    b_g(:,j+1) = Error_State_Vec(15:end,j+1); % Close loop Gyro Dynamic Bias correction 

    % ------------------------- True State Mechanization --------------------
    C_b_w_temp = euler_to_ctm(True_State_Vec(1,j+1),True_State_Vec(2,j+1),True_State_Vec(3,j+1));
    C_b_w_temp = ortho_normalisation(C_b_w_temp);              % Orthonormalization decomposition
    V_eb_w_temp = True_State_Vec(4:6,j+1);
    R_eb_w_temp = True_State_Vec(7:9,j+1);
    wander_ang(j+1) = atan2(True_State_Vec(10),True_State_Vec(11));    % Wander Angle [Rad/Sn] Update
    wander_check(j+1) = (True_State_Vec(10)^2+True_State_Vec(11)^2);   % sin^2 + cos^2 = 1

    [C_b_w(:,:,j+1),V_eb_w(:,j+1),R_eb_w(:,j+1),W_ie_w(:,j+1),wander_ang(j+1)] = Navigation_Update_Wander(Ts,L_b,Lambda_b,H_b,C_b_w_temp,V_eb_w_temp,R_eb_w_temp,F_ib_b,W_ib_b,wander_ang(j+1));
    Error_State_Vec(:,j+1) = zeros(NofState,1); % Set the Error State Zero for close loop correction


end

%% -------------------------------- Anomaly Delete --------------------------------

g_pdf = @(x,mu,sigma) (1/(sqrt(2*sigma^2*pi)))*exp(-((x-mu).^2)/(2*sigma^2)); % Gauss PDF

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
subplot(231)
std_Gx = std(GX);
mean_Gx = mean(GX);
plot(time_vec,GX),grid minor,xlabel("Time [Sn]"),ylabel("[Deg/Sn]"),hold on
title(["Gyro Gx"," Mean: "+string(mean_Gx), " Std: "+string(std_Gx)])
yline(mean_Gx,"r")


subplot(232)
std_Gy = std(GY);
mean_Gy = mean(GY);
plot(time_vec,GY),grid minor,xlabel("Time [Sn]"),ylabel("[Deg/Sn]"),hold on
title(["Gyro Gy"," Mean: "+string(mean_Gy), " Std: "+string(std_Gy)])
yline(mean_Gy,"r")

subplot(233)
std_Gz = std(GZ);
mean_Gz = mean(GZ);
plot(time_vec,GZ),grid minor,xlabel("Time [Sn]"),ylabel("[Deg/Sn]"),hold on
title(["Gyro Gz"," Mean: "+string(mean_Gz), " Std: "+string(std_Gz)])
yline(mean_Gz,"r")

subplot(234)
std_Ax = std(AX);
mean_Ax = mean(AX);
plot(time_vec,AX),grid minor,xlabel("Time [Sn]"),ylabel("[M/Sn^2]"),hold on
title(["Accel Ax"," Mean: "+string(mean_Ax), " Std: "+string(std_Ax)])
yline(mean_Ax,"r")

subplot(235)
std_Ay = std(AY);
mean_Ay = mean(AY);
plot(time_vec,AY),grid minor,xlabel("Time [Sn]"),ylabel("[M/Sn^2]"),hold on
title(["Accel Ay"," Mean: "+string(mean_Ay), " Std: "+string(std_Ay)])
yline(mean_Ay,"r")

subplot(236)
std_Az = std(AZ);
mean_Az = mean(AZ);
plot(time_vec,AZ),grid minor,xlabel("Time [Sn]"),ylabel("[M/Sn^2]"),hold on
title(["Accel Az"," Mean: "+string(mean_Az), " Std: "+string(std_Az)])
yline(mean_Az,"r")

% % ---------------------- Plot Histogram & PDF -----------------------------

% figure('units','normalized','outerposition',[0 0 1 1],'color','w')
% Ax_hist = histogram(AX,BinWidth=5e-4,Normalization='pdf');grid on,hold on
% t_pdf_ax = linspace(Ax_hist.BinLimits(1),Ax_hist.BinLimits(2),size(AX,1))';
% plot(t_pdf_ax,g_pdf(t_pdf_ax,mean_Ax,std_Ax),LineWidth=2,color="r")
% title(["Accel Ax"," Mean: "+string(mean_Ax)," Std: "+string(std_Ax)])
% ylabel("Gauss PDF")
% 
% figure('units','normalized','outerposition',[0 0 1 1],'color','w')
% Ay_hist = histogram(AY,BinWidth=5e-4,Normalization='pdf');grid on,hold on
% t_pdf_ay = linspace(Ay_hist.BinLimits(1),Ay_hist.BinLimits(2),size(AY,1))';
% plot(t_pdf_ay,g_pdf(t_pdf_ay,mean_Ay,std_Ay),LineWidth=2,color="r")
% title(["Accel Ay"," Mean: "+string(mean_Ay)," Std: "+string(std_Ay)])
% ylabel("Gauss PDF")
% 
% figure('units','normalized','outerposition',[0 0 1 1],'color','w')
% Az_hist = histogram(AZ,BinWidth=5e-4,Normalization='pdf');grid on,hold on
% t_pdf_az = linspace(Az_hist.BinLimits(1),Az_hist.BinLimits(2),size(AZ,1))';
% plot(t_pdf_az,g_pdf(t_pdf_az,mean_Az,std_Az),LineWidth=2,color="r")
% title(["Accel Az"," Mean: "+string(mean_Az)," Std: "+string(std_Az)])
% ylabel("Gauss PDF")
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1],'color','w')
% Gx_hist = histogram(GX,BinWidth=5e-4,Normalization='pdf');grid on,hold on
% t_pdf_gx = linspace(Gx_hist.BinLimits(1),Gx_hist.BinLimits(2),size(GX,1))';
% plot(t_pdf_gx,g_pdf(t_pdf_gx,mean_Gx,std_Gx),LineWidth=2,color="r")
% title(["Gyro Gx"," Mean: "+string(mean_Gx)," Std: "+string(std_Gx)])
% ylabel("Gauss PDF")
% 
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1],'color','w')
% Gy_hist = histogram(GY,BinWidth=5e-4,Normalization='pdf');grid on,hold on
% t_pdf_gy = linspace(Gy_hist.BinLimits(1),Gy_hist.BinLimits(2),size(GY,1))';
% plot(t_pdf_gy,g_pdf(t_pdf_gy,mean_Gy,std_Gy),LineWidth=2,color="r")
% title(["Gyro Gy"," Mean: "+string(mean_Gy)," Std: "+string(std_Gy)])
% ylabel("Gauss PDF")
% 
% figure('units','normalized','outerposition',[0 0 1 1],'color','w')
% Gz_hist = histogram(GZ,BinWidth=5e-4,Normalization='pdf');grid on,hold on
% t_pdf_gz = linspace(Gz_hist.BinLimits(1),Gz_hist.BinLimits(2),size(GZ,1))';
% plot(t_pdf_gz,g_pdf(t_pdf_gz,mean_Gz,std_Gz),LineWidth=2,color="r")
% title(["Gyro Gz"," Mean: "+string(mean_Gz)," Std: "+string(std_Gz)])
% ylabel("Gauss PDF")

%% ESKF Output Plot
eskf_time = (1:size(True_State_Vec,2))/AlgorithmFreq;
% --------------------- Attitude ---------------------------
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
subplot(331)
plot(eskf_time,True_State_Vec(1,:)*Rad_to_deg,LineWidth=2),hold on,grid minor
scatter(eskf_time(1),True_State_Vec(1,1)*Rad_to_deg,'red','filled')
scatter(eskf_time(end),True_State_Vec(1,end)*Rad_to_deg,'green','filled')
ylabel("Roll Angle [Deg]")
xlabel("Time [Sn]")
title("Euler Angle")

subplot(332)
plot(eskf_time,True_State_Vec(2,:)*Rad_to_deg,LineWidth=2),hold on,grid minor
scatter(eskf_time(1),True_State_Vec(2,1)*Rad_to_deg,'red','filled')
scatter(eskf_time(end),True_State_Vec(2,end)*Rad_to_deg,'green','filled')
ylabel("Pitch Angle [Deg]")
xlabel("Time [Sn]")
title("Euler Angle")

subplot(333)
plot(eskf_time,mod(((True_State_Vec(3,:)*Rad_to_deg)+yaw_offset),360),LineWidth=2),hold on,grid minor
scatter(eskf_time(1),mod(((True_State_Vec(3,1)*Rad_to_deg)+yaw_offset),360),'red','filled')
scatter(eskf_time(end),mod(((True_State_Vec(3,end)*Rad_to_deg)+yaw_offset),360),'green','filled')
ylabel("Yaw Angle [Deg]")
xlabel("Time [Sn]")
title("Euler Angle")
% --------------------- Velocity ---------------------------
subplot(334)
plot(eskf_time,True_State_Vec(4,:),LineWidth=2)
grid minor
ylabel("North Velocity (X) [M/Sn]")
xlabel("Time [Sn]")
title("NED Velocity")

subplot(335)
plot(eskf_time,True_State_Vec(5,:),LineWidth=2)
grid minor
ylabel("East Velocity (Y) [M/Sn]")
xlabel("Time [Sn]")
title("NED Velocity")

subplot(336)
plot(eskf_time,True_State_Vec(6,:),LineWidth=2)
grid minor
ylabel("Down Velocity (Z) [M/Sn]")
xlabel("Time [Sn]")
title("NED Velocity")
% --------------------- Pose ---------------------------
subplot(337)
plot(eskf_time,True_State_Vec(7,:),LineWidth=2)
grid minor
ylabel("North Pose (X) [M]")
xlabel("Time [Sn]")
title("NED Pose")

subplot(338)
plot(eskf_time,True_State_Vec(8,:),LineWidth=2)
grid minor
ylabel("East Pose (Y) [M]")
xlabel("Time [Sn]")
title("NED Pose")

subplot(339)
plot(eskf_time,True_State_Vec(9,:),LineWidth=2)
grid minor
ylabel("Down Pose (Z) [M]")
xlabel("Time [Sn]")
title("NED Pose")
% --------------------- wander Plot ---------------------------

figure('units','normalized','outerposition',[0 0 1 1],'color','w')
subplot(211)
plot(eskf_time,True_State_Vec(10,:),'k',LineWidth=2)
grid minor
ylabel("Sin(Wander) [Rad]")
xlabel("Time [Sn]")
title("Wander Angle")

subplot(212)
plot(eskf_time,True_State_Vec(11,:),'k',LineWidth=2)
grid minor
ylabel("cos(Wander) [Rad]")
xlabel("Time [Sn]")
title("Wander Angle")
% --------------------- Bias Plot ---------------------------
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
subplot(231)
plot(eskf_time,True_State_Vec(12,:),'r',LineWidth=2)
grid minor
ylabel("Ba X")
xlabel("Time [Sn]")
title("Accel Dynamic Bias")

subplot(232)
plot(eskf_time,True_State_Vec(13,:),'r',LineWidth=2)
grid minor
ylabel("Ba Y")
xlabel("Time [Sn]")
title("Accel Dynamic Bias")

subplot(233)
plot(eskf_time,True_State_Vec(14,:),'r',LineWidth=2)
grid minor
ylabel("Ba Z")
xlabel("Time [Sn]")
title("Accel Dynamic Bias")

subplot(234)
plot(eskf_time,True_State_Vec(15,:),'r',LineWidth=2)
grid minor
ylabel("Bg X")
xlabel("Time [Sn]")
title("Gyro Dynamic Bias")

subplot(235)
plot(eskf_time,True_State_Vec(16,:),'r',LineWidth=2)
grid minor
ylabel("Bg Y")
xlabel("Time [Sn]")
title("Gyro Dynamic Bias")

subplot(236)
plot(eskf_time,True_State_Vec(17,:),'r',LineWidth=2)
grid minor
ylabel("Bg Z")
xlabel("Time [Sn]")
title("Gyro Dynamic Bias")

% -------------------- Heading Angle --------------------
figure('units','normalized','outerposition',[0 0 1 1],'color','w')

subplot(311)
plot(eskf_time,mod(((True_State_Vec(3,:)*Rad_to_deg)+yaw_offset),360),LineWidth=2),hold on,grid minor
scatter(eskf_time(1),mod(((True_State_Vec(3,1)*Rad_to_deg)+yaw_offset),360),'red','filled')
scatter(eskf_time(end),mod(((True_State_Vec(3,end)*Rad_to_deg)+yaw_offset),360),'green','filled')
ylabel("Yaw Angle [Deg]")
xlabel("Time [Sn]")
title("Euler Angle")

subplot(312)
plot(eskf_time,(wander_ang*Rad_to_deg),LineWidth=2),hold on,grid minor
ylabel("Wander Angle [Deg]")
xlabel("Time [Sn]")
title("Wander Angle")

subplot(313)
plot(eskf_time,(wander_ang*Rad_to_deg)+(mod(((True_State_Vec(3,:)*Rad_to_deg)+yaw_offset),360)),LineWidth=2),hold on,grid minor
ylabel("Heading Angle [Deg]")
xlabel("Time [Sn]")
title("Heading Angle")