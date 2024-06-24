function q_dot = quatDerivative(q, omega)
% [Eq. 106 Trawny] OR [Eq. 2 Suh] (and [Eq. 4 Suh])
% QUATDERIVATIVE Computes the quaternion derivative.
% INPUT:
%   * q,                Quaternion (q = [q0; q1; q2; q3], q0 is the scalar)     (4 x 1) vector      []
%   * omega,            Angular velocity omega(k)                               (3 x 1) vector      [rad / s]
% OUTPUT:
%   * q_dot,            Quaternion derivative of q                              (4 x 1) vector      [1 / s]

    % Check number of arguments
    narginchk(2,2);
    if (~isequal(size(q), [4 1]))
        error('q must be a (4 x 1) vector.');
    end
    if (~isequal(size(omega), [3 1]))
        error('omega must be a (3 x 1) vector.');
    end

    % [Eq. 85 Trawny]
    q_dot = 0.5 * Omega(omega) * q;

    % [Eq. 2 Suh] (and [Eq. 4 Suh])
    q_dot_suh = 0.5 * quatMultiplication(q, omega2quat(omega));  %        

    q_dot - q_dot_suh                                                                 
    

%     q_dot = quaternProd(0.5 * q, omega2quat(omega));     
%     % Update orientation
%     % Compute the quaternion derivative
%     Qdot=quaternProd(0.5*Q(i-1,:),Sw);
% 
%     % Update the estimated position
%     Q(i,:)=Q(i-1,:)+Qdot*dt;
% 
%     % Normalize quaternion
%     Q(i,:)=Q(i,:)/norm(Q(i,:));

end