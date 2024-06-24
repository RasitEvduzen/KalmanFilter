function C = euler2RotMat(angles)
% EULER2ROTMAT Converts Euler angles to rotation matrix (aka direction cosine matrix, DCM).
% The Euler angles rotate the frame n (navigation) to the frame b (body) according to 'zyx' sequence.
% INPUT:
%   * angles,           Euler angles (angles = [roll; pitch; yaw])              (3 x 1) vector      [rad]
%
% OUTPUT:
%   * C,                Coordinate transformation matrix from n to b            (3 x 3) matrix      []
%                       NB: That is v_b  = C * v_n.
%                           ('_b' or '_n' mean the vector 'v' is expressed in the frame b or n)


    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(angles), [3 1]))
        error('angles must be a (3 x 1) vector.');
    end
    
    r1 = angles(1);      % roll
    r2 = angles(2);      % pitch
    r3 = angles(3);      % yaw
    
    C = rotX(r1) * rotY(r2) * rotZ(r3);
    
end
