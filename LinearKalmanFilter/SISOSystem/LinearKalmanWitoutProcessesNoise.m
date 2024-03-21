clc,clear all,close all
% Linear Kalman Without Processes Noise
% Written By:Rasit Evduzen
% 11-Apr-2023
%% Estimating the height of building
% Initialization
MeasurementVal = [49.03 48.44 55.21 49.98 50.6 52.61 45.87 42.64 48.26 55.84]';  % Measurement Value
NoD = size(MeasurementVal,1);  % Number Of Iteration
TrueVal = 50*ones(NoD,1);   % True High of building
Xini = 60;     % Amprical Estimation of building hight
Pini = 225;    % Estimation of building hight Variance
Rn = 25;       % Sensor Variance
XNext = [Xini]; % Predict Next State Based On system Dynamics
PNext = [Pini]; % Extrapolate Estimate Uncertainty
% These two parameters depend on system dynamic model (State Space Model)
K = zeros(NoD,1);  % Kalman Gain Vector
%% Calculation Kalman Filter
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for i=1:NoD
    % Prediction Calculation
    XNext(i+1) = XNext(i);
    PNext(i+1) = PNext(i);

    % Correction Stage Calculation
    K(i) = PNext(i) / (PNext(i) + Rn); % Update Kalman Gain!
    XNext(i+1) = XNext(i) + K(i) * (MeasurementVal(i) - XNext(i)); % Current State Estimate!
    PNext(i+1) = (1 - K(i)) * PNext(i);   % Update Current Estimate Uncertainty!
    % Plot Result
    clf
    tSpan = (1:NoD)';
    subplot(211)
    plot(K(1:i),'k','LineWidth',2),grid on
    xlabel("Measurement Number"),ylabel("Kalman Gain"),title("Kalman Gain")

    subplot(212)
    plot(tSpan,TrueVal,'gd-','LineWidth',2),grid on,hold on
    plot(tSpan,MeasurementVal,'bs-','LineWidth',2)
    plot(XNext(1:i),'ro-','LineWidth',2)
    xlabel("Measurement Number"),ylabel("Height [M]"),title("Building height")
    legend("True Value","Measurement Value","Kalman Estimate",'Location','northwest')
    drawnow
end


