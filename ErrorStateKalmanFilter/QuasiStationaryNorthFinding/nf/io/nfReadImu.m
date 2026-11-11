function imu = nfReadImu(src, fs)
%NFREADIMU Read and clean an IMU record in the RawData.csv format (sensor axes, SI units).
%   imu = nfReadImu(file, fs)   CSV rows [index, GX GY GZ (deg/s), AX AY AZ (m/s^2)]
%   imu = nfReadImu(M, fs)      the same data as an N x 7 matrix
%   fs is the sampling rate [Hz] (default 200); the record time is t = index/fs.
%
%   Cleaning: rows without data are dropped (all six sensor values exactly 0, e.g. the first row
%   of RawData.csv, a non-finite index, or no finite sensor value). Then, per channel, samples
%   that are non-finite, outside the sensor range (400 deg/s, 200 m/s^2) or more than 30 robust
%   standard deviations (1.4826 MAD) from the median are replaced by the next valid sample of
%   the channel (on RawData.csv only the AY = -11236 m/s^2 spike of data row 2).
%
%   imu.t (1xN) [s], imu.gyro (3xN) [rad/s], imu.accel (3xN) [m/s^2], imu.index (1xN), imu.fs,
%   imu.report: source, rows_dropped (data rows of the file), replaced (one row per replaced
%   sample: [data row, sample index, channel 1..6 = GX..AZ, old value, new value]) and gaps
%   (index gaps of the kept rows, [index before, index after]).

if nargin < 2
    fs = 200;
end
if ischar(src) || isstring(src)
    source = char(src);
    raw = readCsv(source);
else
    source = '<matrix>';
    raw = double(src);
end
if size(raw, 2) ~= 7
    error('nf:io:format', 'Expected N x 7 data [index GX GY GZ AX AY AZ], got %d columns.', size(raw, 2));
end

% ---- rows without data
rowNum = (1:size(raw, 1))';
drop = all(raw(:, 2:7) == 0, 2) | ~isfinite(raw(:, 1)) | all(~isfinite(raw(:, 2:7)), 2);
R = raw(~drop, :);
rows = rowNum(~drop);
if isempty(R)
    error('nf:io:empty', 'No valid rows in %s.', source);
end
index = R(:, 1);
X = R(:, 2:7);                                   % deg/s and m/s^2

% ---- spikes, non-finite and out-of-range samples
range = [400 400 400 200 200 200];
replaced = zeros(0, 5);
for c = 1:6
    x = X(:, c);
    fin = isfinite(x);
    if ~any(fin)
        error('nf:io:channel', 'Channel %d of %s has no finite samples.', c, source);
    end
    med = median(x(fin));
    madv = median(abs(x(fin) - med));
    bad = ~fin | abs(x - med) > 30 * max(1.4826*madv, 1e-4) | abs(x) > range(c);
    good = find(~bad);
    if isempty(good)
        error('nf:io:channel', 'All samples of channel %d of %s are invalid.', c, source);
    end
    for i = find(bad)'
        nxt = good(find(good > i, 1));
        if isempty(nxt)
            nxt = good(find(good < i, 1, 'last'));
        end
        replaced(end + 1, :) = [rows(i), index(i), c, x(i), x(nxt)]; %#ok<AGROW>
        X(i, c) = x(nxt);
    end
end

% ---- time base
t = index / fs;
d = diff(index);
if any(d <= 0)
    warning('nf:io:index', 'The sample index of %s is not strictly increasing: t = (index(1) + (0:N-1))/fs.', source);
    t = (index(1) + (0:numel(index) - 1)') / fs;
end

imu.t = t';
imu.gyro = X(:, 1:3)' * (pi / 180);
imu.accel = X(:, 4:6)';
imu.index = index';
imu.fs = fs;
imu.report = struct('source', source, 'rows_dropped', rowNum(drop)', 'replaced', replaced, ...
    'gaps', [index(d ~= 1), index(find(d ~= 1) + 1)]);
end

% =====================================================================
function raw = readCsv(file)
% fast path: plain numeric lines; anything else (header, empty fields) goes to READMATRIX
fid = fopen(file, 'r');
if fid < 0
    error('nf:io:file', 'Cannot open the IMU record ''%s''.', file);
end
c = textscan(fid, '%f%f%f%f%f%f%f', 'Delimiter', ',');
complete = feof(fid) && all(cellfun(@numel, c) == numel(c{1}));
fclose(fid);
if complete
    raw = [c{:}];
else
    raw = readmatrix(file, 'FileType', 'text', 'NumHeaderLines', 0);
end
end
