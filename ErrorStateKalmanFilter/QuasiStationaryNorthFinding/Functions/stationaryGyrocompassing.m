function [heading] = stationaryGyrocompassing(bodyRate, attitude)
%STATIONARYGYROCOMPASSING Estimates the heading after leveling the body
%   Input:
%   bodyRate,   the angular rate of the body in body frame
%   attude,     roll and pitch attitude of the body
%   Output:
%   heading,    estimated heading
%   headingError,   estimated heading error

roll = attitude(1);
pitch = attitude(2);

w_x = bodyRate(1);
w_y = bodyRate(2);
w_z = bodyRate(3);



sinPsi = -w_y * cos(roll) + w_z * sin(roll);
cosPsi = w_x * cos(pitch) + w_y * sin(roll) * sin(pitch) + w_z * cos(roll) * sin(pitch);

heading = atan2(sinPsi, cosPsi);

% PG 198, eq [5.105]
end

