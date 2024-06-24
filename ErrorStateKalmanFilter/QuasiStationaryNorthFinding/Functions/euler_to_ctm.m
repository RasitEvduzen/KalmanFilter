function [ctm] = euler_to_ctm(roll,pitch,yaw)

% rtx = [1  0           0; 
%        0  cos(roll)  sin(roll);
%        0 -sin(roll)  cos(roll)]; 
% 
% rty = [cos(pitch)  0  -sin(pitch); 
%        0            1   0;
%        sin(pitch)  0  cos(pitch)];  
% 
% rtz = [ cos(yaw)  sin(yaw)  0;
%        -sin(yaw)  cos(yaw)  0;
%         0          0          1]; 
% ctm = rtx*rty*rtz; 

ctm = eul2rotm([roll,pitch,yaw],'XYZ'); % The default order for Euler angle rotations is "ZYX"
end