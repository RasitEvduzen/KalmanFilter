function q_times_p = quatMultiplicationTrawny(q, p)
% INPUT:
%   * q,Quaternion q(k) (q = [q0; q1; q2; q3], q0 is the scalar)        (4 x 1) vector      []
%   * p,Quaternion p(k) (p = [p0; p1; p2; p3], p0 is the scalar)        (4 x 1) vector      []
    narginchk(2,2);
    if (~isequal(size(q), [4 1]) && ~isequal(size(p), [4 1]))
        error('q and p must be (4 x 1) vectors.');
    end
    if (~isequal(size(p), [4 1]))
        error('p must be a (4 x 1) vector.');
    end
    if (~isequal(size(q), [4 1]))
        error('q must be a (4 x 1) vector.');
    end
    % Trawny Eq 8 L(q)
    q_times_p = [q(4)*eye(3)-skew(q(1:3)) q(1:3); -q(1:3)' q(4)] * p;
  
    % Trawny Eq 10 L(p)
%     q_times_p = [p(4)*eye(3)-skew(p(1:3)) p(1:3); -p(1:3)' p(4)] * q;
end