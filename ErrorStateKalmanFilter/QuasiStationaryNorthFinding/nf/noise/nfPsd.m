function [f, P] = nfPsd(x, fs, nseg)
%NFPSD One-sided power spectral density by Welch's method [units^2/Hz].
%   [f, P] = nfPsd(x, fs)         segments of min(4096, largest power of 2 <= N) samples
%   [f, P] = nfPsd(x, fs, nseg)   segment length (even)
%   Periodic Hann window, 50 % overlap, mean removed per segment, density scaling; all bins
%   except DC and Nyquist are doubled. White noise of two-sided density N gives P = 2 N^2.

x = x(:)';
n = numel(x);
if nargin < 3
    nseg = min(4096, 2^floor(log2(n)));
end
w = 0.5 - 0.5*cos(2*pi*(0:nseg - 1)/nseg);
starts = 1:nseg/2:n - nseg + 1;
P = zeros(1, nseg/2 + 1);
for s = starts
    seg = x(s:s + nseg - 1);
    X = fft((seg - mean(seg)) .* w);
    P = P + abs(X(1:nseg/2 + 1)).^2;
end
P = P / (numel(starts) * fs * sum(w.^2));
P(2:end - 1) = 2 * P(2:end - 1);
f = (0:nseg/2) * fs / nseg;
end
