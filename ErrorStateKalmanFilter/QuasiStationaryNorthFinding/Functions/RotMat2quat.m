function q = RotMat2quat(C)
% [Eq. 98a - 98b - 99a - 99b Trawny]
% ROTMAT2QUAT Converts direction cosine matrix to quaternion.
% INPUT:
%   * C,                Coordinate transformation matrix                        (3 x 3) matrix      
% OUTPUT:
%   * q,                Quaternion (q = [q0; q1; q2; q3], q0 is the scalar)     (4 x 1) vector      

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(C), [3 3]))
        error('C must be a (3 x 3) matrix.');
    end

    % Norm check
    if (abs(1 - norm(C)) > 1e-15)
        warning(['C should have norm 1:    norm(C) = ', num2str(norm(C)), ',    norm(C) - 1 = ', num2str(norm(C)-1)])
    end
        
    for i=1:3
        for j=1:3
            if (C(i,j) < -1)
                warning('C[%i,%i]=%f is smaller than -1', i, j, C(i,j));
            elseif (C(i,j) > 1)
                warning('C[%i,%i]=%f is greater than 1', i, j, C(i,j));
            end
        end
    end
    
    c11 = C(1,1);
    c12 = C(1,2);
    c13 = C(1,3);
    c21 = C(2,1);
    c22 = C(2,2);
    c23 = C(2,3);
    c31 = C(3,1);
    c32 = C(3,2);
    c33 = C(3,3);
    
    % [Eq. 96 Trawny]
    T = trace(C);
    [~, max_index] = max([c11, c22, c33, T]);
    
    switch(max_index)
        
        case 1
            % [Eq. 98a Trawny]
            q1 = sqrt(1 + 2 * c11 - T) / 2;
            q = [q1;
                 (c12 + c21) / (4 * q1);
                 (c13 + c31) / (4 * q1);
                 (c23 - c32) / (4 * q1)];
            
        case 2
            % [Eq. 98b Trawny]
            q2 = sqrt(1 + 2 * c22 - T) / 2;
            q = [(c12 + c21) / (4 * q2);
                 q2;
                 (c23 + c32) / (4 * q2);
                 (c31 - c13) / (4 * q2)];
             
        case 3
            % [Eq. 99a Trawny]
            q3 = sqrt(1 + 2 * c33 - T) / 2;
            q = [(c13 + c31) / (4 * q3);
                 (c23 + c32) / (4 * q3);
                 q3;
                 (c12 - c21) / (4 * q3)];
             
        case 4
            % [Eq. 99b Trawny]
            q4 = sqrt(1 + T) / 2;
            q = [(c23 - c32) / (4 * q4);
                 (c31 - c13) / (4 * q4);
                 (c12 - c21) / (4 * q4);
                 q4];
             
        otherwise
            error('Unexpected behavior.')
    end
    
    % quaternion normalization
    q = q / norm(q);

    % Ensure quaternion scalar part is non-negative
    if (q(4) < 0)
        q = -q;
    end
    
    q = [q(4); q(1); q(2); q(3)]; 
    
end