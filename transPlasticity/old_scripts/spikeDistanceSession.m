function R = spikeDistanceSession(spikeFile, varargin)
%SPIKEDISTANCESESSION  SPIKE-distance profile over a whole session plus a set
%   of interval-jitter nulls.
%
%   R = spikeDistanceSession(spikeFile, <options>)
%
%   OBSERVED. Multivariate SPIKE-distance S(t) (spikeDistanceProfile) of the
%   selected units over [0 dur], dur = last spike of any unit. Stored exactly
%   (event times E, interval start/end values y1, y2) and as exact bin means
%   on a fixed grid of width 'BinWidth' starting at t = 0.
%
%   NULL. INTERVAL JITTER on a fixed grid anchored at the start of the
%   recording (t = 0), not at each spike: time is cut into bins [k*J,(k+1)*J),
%   J = 'Jitter', and every spike is redrawn uniformly inside its own bin (the
%   last bin cut at dur). Each unit's spike count in every J-bin is preserved;
%   timing finer than J is destroyed. Surrogate s uses
%   RandStream('mt19937ar','Seed', seed + s): the set reproduces exactly.
%   Each surrogate's profile is reduced to exact bin means on the SAME grid as
%   the observed one (its event times differ, so the exact form is not kept).
%
%   OPTIONS
%     'Jitter'    1        s
%     'NSurr'     100
%     'BinWidth'  0.1      s, grid for the stored bin means (observed and null)
%     'MatchN'    10       units, seeded draw as in syncBatch ([] = all in band)
%     'RateRange' [0.1 10] Hz over the whole session ([] = no rate exclusion)
%     'RI', 'MRTS'         passed to spikeDistanceProfile (defaults: standard)
%     'Seed'      []       default: the syncBatch rule on the file path
%
%   OUTPUT R
%     .E .y1 .y2     exact observed profile (see spikeDistanceProfile)
%     .edges         bin edges (s from recording start); last bin may be short
%     .S             observed bin means;  .Snull  nBins x NSurr (single)
%     .D             observed time-averaged SPIKE-distance; .Dnull NSurr x 1
%     .N .keepIdx .detIdx .uid .dur .params .ok .why .file
%   S is a DISSIMILARITY: lower = more synchronous.

p = inputParser;
p.addParameter('Jitter', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('NSurr', 100, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('BinWidth', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('MatchN', 10);
p.addParameter('RateRange', [0.1 10]);
p.addParameter('RI', false);
p.addParameter('MRTS', 0);
p.addParameter('Seed', []);
p.addParameter('KeepExact', true);   % true | false | 'auto' (exact only if pairs x spikes <= ExactLimit)
p.addParameter('ExactLimit', 2e9);   % pairs x pooled spikes; ~1e9 took ~2 min
p.addParameter('MaxPairs', Inf);     % > this many unit pairs: average over a seeded random subset
p.addParameter('Verbose', true);
p.parse(varargin{:});
o = p.Results;

R = struct('ok',false,'why','','file',spikeFile,'E',[],'y1',[],'y2',[],'edges',[], ...
    'S',[],'Snull',[],'D',NaN,'Dnull',[],'N',0,'keepIdx',[],'detIdx',[],'uid',[], ...
    'dur',NaN,'params',o);

S = load(spikeFile, 'spikes');
if ~isfield(S,'spikes') || ~isfield(S.spikes,'times'), R.why = 'no spikes.times'; return; end
ts = S.spikes.times(:);
ts = cellfun(@(x) unique(reshape(x(~isnan(x)), [], 1)), ts, 'UniformOutput', false);
lastSpk = cellfun(@(x) max([x; -Inf]), ts);
dur = max(lastSpk(isfinite(lastSpk)));
if isempty(dur) || ~isfinite(dur) || dur <= 0, R.why = 'bad duration'; return; end
R.dur = dur;

% ---- units ---------------------------------------------------------------
M = cellfun(@numel, ts);
rr = o.RateRange; if isempty(rr), rr = [0 Inf]; end   % [] = no rate exclusion (any unit with spikes)
keep = M / dur > rr(1) & M / dur < rr(2) & M > 0;
R.keepIdx = find(keep);
if numel(R.keepIdx) < 2, R.why = sprintf('%d units in rate band', numel(R.keepIdx)); return; end
seed = o.Seed;
if isempty(seed), seed = mod(sum(double(spikeFile) .* (1:numel(spikeFile))), 2^31 - 1); end
sel = 1:numel(R.keepIdx);
if ~isempty(o.MatchN)
    if numel(R.keepIdx) < o.MatchN
        R.why = sprintf('%d units in rate band < MatchN %d', numel(R.keepIdx), o.MatchN); return
    end
    rs = RandStream('mt19937ar', 'Seed', seed);
    sel = randperm(rs, numel(R.keepIdx), o.MatchN);
end
R.detIdx = R.keepIdx(sel);
if isfield(S.spikes, 'UID') && numel(S.spikes.UID) >= max(R.detIdx), R.uid = S.spikes.UID(R.detIdx); end
trains = ts(R.detIdx);
R.N = numel(trains);
if R.N < 2, R.why = 'fewer than 2 units'; return; end
iv = [0 dur];
% unit pairs: all, or a seeded random subset (same for data and every surrogate)
PL = nchoosek(1:R.N, 2); nAll = size(PL, 1);
if nAll > o.MaxPairs
    rsP = RandStream('mt19937ar', 'Seed', seed + 1e6);
    PL = sortrows(PL(randperm(rsP, nAll, o.MaxPairs), :));
end
R.pairs = R.detIdx(PL);                      % pairs as indices into spikes.times
R.params.nPairsAll = nAll; R.params.nPairsUsed = size(PL, 1);
pa = {'Interval', iv, 'RI', o.RI, 'MRTS', o.MRTS, 'Pairs', PL};

g = (0:o.BinWidth:dur)';
if g(end) < dur, g(end+1) = dur; end
R.edges = g;

nSpk = sum(cellfun(@numel, trains));
if o.Verbose
    fprintf('    %d units, %d spikes, %.1f h, %d pairs (of %d), %d bins\n', R.N, nSpk, dur/3600, size(PL,1), nAll, numel(g)-1);
end

% ---- observed --------------------------------------------------------------
tObs = tic;
doExact = o.KeepExact;
if ischar(doExact) || isstring(doExact)      % 'auto'
    doExact = size(PL, 1) * nSpk <= o.ExactLimit;
    if o.Verbose && ~doExact, fprintf('    exact observed profile skipped (pairs x spikes = %.2g)\n', size(PL,1)*nSpk); end
end
if doExact
    % exact profile on the pooled grid: cost ~ pairs x total spikes; done once
    [R.E, R.y1, R.y2, R.D] = spikeDistanceProfile(trains, pa{:});
    R.S = spikeDistanceBinMean(R.E, R.y1, R.y2, g);
else
    [R.S, ~, ~, R.D] = spikeDistanceProfile(trains, pa{:}, 'Bins', g);
end
if o.Verbose, fprintf('    observed: %.0f s\n', toc(tObs)); end

% ---- interval-jitter surrogates ---------------------------------------------
% Fast binned path: each pair evaluated on its own spikes only, integrated to
% the bins (identical numbers to binning the exact profile; see
% spikeDistanceProfile 'Bins'). Surrogate 1 runs first to time it.
nS = o.NSurr; J = o.Jitter; N = R.N;
cols = cell(1, nS); Dn = nan(nS, 1);
if nS > 0
    t1s = tic;
    [m1, ~, ~, D1] = spikeDistanceProfile(local_jitter(trains, J, dur, seed + 1), pa{:}, 'Bins', g);
    tS = toc(t1s);
    if o.Verbose
        nw = 1; pp = gcp('nocreate'); if ~isempty(pp), nw = pp.NumWorkers; end
        fprintf('    surrogate: %.0f s each -> ~%.1f h for %d on %d worker(s)\n', tS, tS*(nS-1)/nw/3600, nS, nw);
    end
end
parfor s = 2:nS
    [ms, ~, ~, Ds] = spikeDistanceProfile(local_jitter(trains, J, dur, seed + s), pa{:}, 'Bins', g);
    cols{s} = single(ms);
    Dn(s) = Ds;
end
if nS > 0, cols{1} = single(m1); Dn(1) = D1; end   % set after parfor: sliced outputs
if nS > 0, R.Snull = [cols{:}]; else, R.Snull = zeros(numel(g) - 1, 0, 'single'); end
R.Dnull = Dn;
R.params.Seed = seed;
R.ok = true;
end


function tj = local_jitter(trains, J, dur, seed)
% interval jitter on the fixed grid [k*J, (k+1)*J) anchored at t = 0
rs = RandStream('mt19937ar', 'Seed', seed);
tj = cell(numel(trains), 1);
for n = 1:numel(trains)
    x = trains{n};
    lo = floor(x / J) * J;
    hi = min(lo + J, dur);
    tj{n} = sort(lo + rand(rs, numel(x), 1) .* (hi - lo));
end
end