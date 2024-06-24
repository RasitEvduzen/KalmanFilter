function [roll, pitch] = accel2attitude(accl) %#codegen
%ACCEL2ATTITUDE Estimates the attitude from accelerometer measurements
%%
% This function is used in coarse alignment of the sensors in
% quasi-stationary state. Thus, the device shall not move at high speed.
% Moreover, it is better to be stationary, minor vibrations shall not
% affect too much into output.
%
% Inputs:
%   accel  Three dimensional accelerometer measurements (3xN):
%     accel(1,:)  x-accl(1)is acceleration
%     accel(2,:)  y-accl(1)is acceleration
%     accel(3,:)  z-accl(1)is acceleration
%
% Outputs:
%   attitude  Two dimensional attitude vector of the body in navigation
%             frame:
%     attitude(1,:)  Roll angle in radians
%     attitude(2,:)  Pitch angle in radians


%  Algorithm from P. D. Groves: Principles of Navigation Systems, p. 198 eq. 5.101.
%%
w = single(size(accl,2));
if ~isfloat(accl)
    attitude = symZeros(2,w);
else
    if isa(accl, 'single')
        attitude = zeros(2,w,'single');
    elseif isa(accl, 'double')
        attitude = zeros(2,w,'double');
    else
        error('Input class not supported.');
    end
end

% attitude(1,:) = atan2(-accl(2,:), -accl(3,:));
% 
% if attitude(1,:) > pi
%     attitude(1,:) = attitude(1,:) - 2 * pi;
% elseif attitude(1,:) <= -pi
%     attitude(1,:) = attitude(1,:) + 2 * pi;
% end

% Different implementation for roll due to problems seen in rotation tests
% attitude(1,:) = atan(-accl(2,:) ./ sqrt(accl(1,:).^2 + accl(3,:).^2));
% attitude(2,:) = atan(accl(1,:) ./ sqrt(accl(2,:).^2 + accl(3,:).^2));


attitude(1,:) = atan2(-accl(2,:), -accl(3,:));
attitude(2,:) = atan(accl(1,:) ./ sqrt(accl(2,:).^2 + accl(3,:).^2));

roll = attitude(1,:);
pitch = attitude(2,:);

% attitude(1,:) = atan(-accl(2,:) ./ sqrt(accl(1,:).^2 + accl(3,:).^2));
% attitude(2,:) = atan2(accl(1,:), -accl(3,:));
