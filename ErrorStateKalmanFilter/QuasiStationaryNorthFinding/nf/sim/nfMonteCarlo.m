function mc = nfMonteCarlo(cfg, n, varargin)
%NFMONTECARLO Heading error statistics of nfProcess on n synthetic static records.
%   mc = nfMonteCarlo(cfg, n, Name, Value, ...) simulates records with nfStaticImu
%   (seeds 1..n; Name/Value options are passed on) and runs the full pipeline on each.
%   Extra options:
%     'RandomHeading' (false)  true heading uniform in [0, 360) deg per run
%     'RandomBias' (false)     constant biases drawn from the filter priors (cfg.sensor.bg_dph,
%                              ba_mg) per run, so that the heading sigma can be checked
%   mc: heading_deg (truth), err_deg (estimate - truth, wrapped to +-180), sigma_deg, z =
%   err/sigma, gc_err_deg (gyrocompass on cfg.align.window), gc_sigma_deg, and the statistics
%   mean_err_deg, std_err_deg, rms_z (consistent filter: ~1).

p = inputParser;
p.KeepUnmatched = true;
p.addParameter('RandomHeading', false);
p.addParameter('RandomBias', false);
p.parse(varargin{:});
o = p.Results;
pass = namedargs(p.Unmatched);

mc.heading_deg = zeros(1, n);
mc.err_deg = zeros(1, n); mc.sigma_deg = zeros(1, n);
mc.gc_err_deg = zeros(1, n); mc.gc_sigma_deg = zeros(1, n);
wrap = @(x) mod(x + 180, 360) - 180;
for i = 1:n
    rs = RandStream('mt19937ar', 'Seed', 100000 + i);
    extra = {};
    if o.RandomHeading
        extra = [extra, {'HeadingDeg', 360 * rand(rs)}]; %#ok<AGROW>
    end
    if o.RandomBias
        extra = [extra, {'BiasGDph', cfg.sensor.bg_dph * randn(rs, 1, 3), ...
            'BiasAmg', cfg.sensor.ba_mg * randn(rs, 1, 3)}]; %#ok<AGROW>
    end
    [imu, truth] = nfStaticImu(cfg, pass{:}, extra{:}, 'Seed', i);
    res = nfProcess(cfg, imu);
    mc.heading_deg(i) = truth.heading_deg;
    mc.err_deg(i) = wrap(res.summary.heading_deg - truth.heading_deg);
    mc.sigma_deg(i) = res.summary.sigma_deg;
    mc.gc_err_deg(i) = wrap(res.summary.coarse_heading_deg - truth.heading_deg);
    mc.gc_sigma_deg(i) = res.summary.coarse_sigma_deg;
end
mc.z = mc.err_deg ./ mc.sigma_deg;
mc.mean_err_deg = mean(mc.err_deg);
mc.std_err_deg = std(mc.err_deg);
mc.rms_z = sqrt(mean(mc.z.^2));
end

% =====================================================================
function c = namedargs(s)
c = reshape([fieldnames(s)'; struct2cell(s)'], 1, []);
end
