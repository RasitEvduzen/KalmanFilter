function q = nfRotVec2Quat(rho)
%NFROTVEC2QUAT Quaternion of a rotation vector rho [rad]: angle |rho| about rho/|rho|.
%   q = [cos(|rho|/2); sin(|rho|/2) rho/|rho|]; the Taylor series is used for |rho| < 1e-8.

a = norm(rho);
if a < 1e-8
    q = [1 - a^2/8; (0.5 - a^2/48) * rho(:)];
else
    q = [cos(a/2); sin(a/2) / a * rho(:)];
end
end
