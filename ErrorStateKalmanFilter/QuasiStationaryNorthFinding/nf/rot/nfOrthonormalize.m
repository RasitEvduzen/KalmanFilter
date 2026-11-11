function C = nfOrthonormalize(C)
%NFORTHONORMALIZE Nearest rotation matrix (Frobenius norm) of a near-orthonormal 3x3 matrix.
%   C = nfOrthonormalize(C) returns U*V' from the SVD C = U*S*V'. A rotation matrix is
%   returned unchanged (to round-off). The Q factor of qr(C), used by the 2023 code, is not a
%   valid replacement: LAPACK does not force diag(R) > 0, so two axes can flip (a 180 deg yaw jump).

[U, ~, V] = svd(C);
C = U*V';
if det(C) < 0                   % C was far from a rotation: flip the weakest singular direction
    U(:, 3) = -U(:, 3);
    C = U*V';
end
end
