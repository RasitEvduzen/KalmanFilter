clc,clear all,close all
% Linear Kalman Without Processes Noise
% Written By:Rasit Evduzen
% 11-Apr-2023
%% Estimating the temperature of the liquid in a tank "Const Dynamics ~ Static System"
TrueVal = [50.005  49.994  49.993  50.001  50.006  49.998  50.021  50.005  50  49.997]';   % True value of Liquid temp in a tank
MeasurementVal = [50.486  50.963  51.597  52.001  52.518  53.05  53.438  53.858  54.465  55.114]';  % Measurement Value
Xini = 50;     % Amprical Estimation of liquid temp initial value
Pini = 1e-6;    % Estimation variance temp of liquid
Rn = 1e-2;     % Measurement Noise Variance
Qn = 1e-4;      % Processes Noise Variance
% -----------------------------------------------------------
NoD = size(MeasurementVal,1);  % Number Of Iteration
Ts = 5;        % Sampling Periode [Sn]
XNext = [Xini]; % Predict Next State Based On system Dynamics
PNext = [Pini + Qn]; % Extrapolate Estimate Uncertainty
% These two parameters depend on system dynamic model (State Space Model)
K = zeros(NoD,1);  % Kalman Gain Vector
tSpan = (1:NoD)';
%% Calculation Kalman Filter
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
for i=1:NoD
    % Prediction Stage Calculation
    XNext(i+1) = XNext(i);
    PNext(i+1) = PNext(i) + Qn;
    % Correction Stage Calculation
    K(i) = PNext(i) / (PNext(i) + Rn); % Update Kalman Gain!
    XNext(i+1) = XNext(i) + K(i) * (MeasurementVal(i) - XNext(i)); % Current State Estimate!
    PNext(i+1) = (1 - K(i)) * PNext(i);   % Update Current Estimate Uncertainty!
    
    % Plot Result
    clf
    subplot(211)
    plot(K(1:i),'k','LineWidth',2),grid on
    xlabel("Measurement Number"),ylabel("Kalman Gain"),title("Kalman Gain")

    subplot(212)
    plot(tSpan,TrueVal,'gd-','LineWidth',2),hold on,grid on
    plot(tSpan,MeasurementVal,'bs-','LineWidth',2)
    plot(XNext(1:i),'ro-','LineWidth',2)
    xlabel("Measurement Number"),ylabel("Temperature [^oC]"),title("Estimating the temperature")
    legend("True Value","Measurement Value","Kalman Estimate",'Location','northwest')
    drawnow
end


