function r = nfQuatMultiply(p, q)
%NFQUATMULTIPLY Hamilton product r = p (x) q of scalar-first quaternions [q0; q1; q2; q3].
%   With C(q) from nfQuat2Dcm, C(p (x) q) = C(p) C(q): the right factor is applied first.

p = p(:); q = q(:);
r = [p(1)*q(1) - p(2:4)'*q(2:4);
     p(1)*q(2:4) + q(1)*p(2:4) + cross(p(2:4), q(2:4))];
end
