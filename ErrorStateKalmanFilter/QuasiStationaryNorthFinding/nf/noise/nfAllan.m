function [tau, adev] = nfAllan(x, fs, m)
%NFALLAN Overlapping Allan deviation of a rate signal (result in the units of x).
%   [tau, adev] = nfAllan(x, fs)      cluster sizes log-spaced from 1 to (N-1)/3 samples
%   [tau, adev] = nfAllan(x, fs, m)   given cluster sizes m [samples]
%   With theta = [0, cumsum(x)]/fs and tau = m/fs:
%       AVAR(tau) = sum((theta(k+2m) - 2 theta(k+m) + theta(k))^2) / (2 tau^2 (N - 2m + 1))
%   White rate noise of density N gives ADEV = N/sqrt(tau): ADEV(1 s) is the angle (velocity)
%   random walk.

x = x(:)';
n = numel(x);
if nargin < 3
    m = unique(round(logspace(0, log10(floor((n - 1) / 3)), 60)));
end
th = [0, cumsum(x)] / fs;
tau = m / fs;
adev = zeros(size(m));
for i = 1:numel(m)
    d = th(2*m(i) + 1:end) - 2*th(m(i) + 1:end - m(i)) + th(1:end - 2*m(i));
    adev(i) = sqrt(sum(d.^2) / (2 * tau(i)^2 * numel(d)));
end
end
