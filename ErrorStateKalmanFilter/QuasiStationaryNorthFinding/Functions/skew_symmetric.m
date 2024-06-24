function A = skew_symmetric(a)
%Skew_symmetric - Calculates skew-symmetric matrix

A = [    0, -a(3),  a(2);
      a(3),     0, -a(1);
     -a(2),  a(1),     0];
end