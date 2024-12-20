function r = RotMat2euler(C)
% ROTMAT2EULER Converts direction cosine matrix to Euler angles ('zyx' sequence).
%  "We’ll follow the notational conventions of Shoemake’s “Euler Angle Conversion”
% INPUT:
%   * C,                Coordinate transformation matrix                        (3 x 3) matrix      []
%
% OUTPUT:
%   * r,                Euler angles (r = [roll; pitch; yaw])                   (3 x 1) vector      [rad]

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(C), [3 3]))
        error('C must be a (3 x 3) matrix.');
    end

    C11 = C(1, 1);
    C12 = C(1, 2);
    C13 = C(1, 3);
    C21 = C(2, 1);
    C22 = C(2, 2);
    C23 = C(2, 3);
    C31 = C(3, 1);
    C32 = C(3, 2);
    C33 = C(3, 3);
    
    
    r1 = atan2(C23, C33);                                   % roll
    
    c2 = sqrt(C11^2 + C12^2);
    r2 = atan2(-C13, c2);                                   % pitch
    
    s1 = sin(r1);
    c1 = cos(r1);
    r3 = atan2(s1 * C31 - c1 * C21, c1 * C22 - s1 * C32);   % yaw
%     r3 = atan2(C12, C11);               
    
    r = [r1; r2; r3];
        
end