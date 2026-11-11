function w = nfEarthRate(L)
%NFEARTHRATE Earth rotation rate resolved in NED at latitude L [rad]: W_ie*[cos(L); 0; -sin(L)] [rad/s].

k = nfEarthConstants();
w = k.W_ie * [cos(L); 0; -sin(L)];
end
