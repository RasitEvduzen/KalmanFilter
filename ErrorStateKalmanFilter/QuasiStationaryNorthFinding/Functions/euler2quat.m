function q = euler2quat(angles)
% INPUT:
%   * angles,           Euler angles (angles = [roll; pitch; yaw])              (3 x 1) vector      [rad]
%
% OUTPUT:
%   * q,                Quaternion (q = [q0; q1; q2; q3], q0 is the scalar)     (4 x 1) vector      []


    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(angles), [3 1]))
        error('angles must be a (3 x 1) vector.');
    end
    
%     angles = wrapToPi(angles);
    
    c1 = cos(angles(1) / 2);
    c2 = cos(angles(2) / 2);
    c3 = cos(angles(3) / 2);
    s1 = sin(angles(1) / 2);
    s2 = sin(angles(2) / 2);
    s3 = sin(angles(3) / 2);
    
    
    q = [c1 * c2 * c3 + s1 * s2 * s3;
         s1 * c2 * c3 - c1 * s2 * s3;
         c1 * s2 * c3 + s1 * c2 * s3;
         c1 * c2 * s3 - s1 * s2 * c3];
        
        
    % Normalize the quaternion
    q = q / norm(q);
    
    % Ensure quaternion scalar part is non-negative
    if (q(1) < 0)
        q = -q;
    end
            
end