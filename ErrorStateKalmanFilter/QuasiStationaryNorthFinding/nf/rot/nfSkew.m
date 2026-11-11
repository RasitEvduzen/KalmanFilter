function S = nfSkew(a)
%NFSKEW Skew-symmetric (cross-product) matrix: nfSkew(a)*b = cross(a, b).

S = [    0, -a(3),  a(2);
      a(3),     0, -a(1);
     -a(2),  a(1),     0];
end
