function R = syncProfileSession(spikeFile, varargin)
%SYNCPROFILESESSION  Spike-resolved SPIKE-synchronization over a whole session,
%   plus a set of interval-jitter nulls on the same spikes.
%
%   R = syncProfileSession(spikeFile, <options>)
%
%   OBSERVED. spikeSync's pooled profile C(t_k), one value per spike of the
%   selected units, over [0 dur] (dur = last spike of any unit). Stored as the
%   integer count of the other N-1 units with a coincident spike,
%   count = C * (N-1), so C = count/(N-1).
%
%   NULL. INTERVAL JITTER on a fixed grid anchored at the START OF THE
%   RECORDING (t = 0), not at each spike: time is cut into bins
%   [k*J, (k+1)*J), J = 'Jitter' (default 1 s), and every spike is redrawn
%   uniformly inside ITS OWN bin (the last bin is cut at dur). This preserves
%   each unit's spike count in every bin exactly -- the slow rate profile at
%   resolution J -- and destroys timing finer than J. A spike-centred jitter
%   (t + U(-J/2, J/2)) is a different surrogate: it lets the jittered spike
%   leave its bin, so it does not preserve the binned counts, and the amount
%   of structure it can mistake for fine timing depends on where the spike
%   sat. Surrogate k is drawn from RandStream('mt19937ar','Seed', seed + k),
%   so the whole set reproduces from the file name and the options.
%
%   ALIGNMENT. Bin counts are preserved, so the i-th spike of a unit stays in
%   the same bin. Surrogate profiles are therefore stored on the OBSERVED
%   spike axis: R.null(k, s) is the coincidence count of the spike that
%   started as observed spike k, after jitter, in surrogate s. Its true time in
%   the surrogate differs from R.t(k) by less than J. Any smoothing of the
%   profile (spike-count windows, kernels, time windows) can be applied
%   identically to the observed column and to every null column.
%
%   OPTIONS
%     'Jitter'     1       interval width (s)
%     'NSurr'      100     number of jittered surrogates
%     'MatchN'     10      units per session, drawn once with the same seeded
%                          randperm as syncBatch ([] = all units in the rate band)
%     'RateRange'  [0.1 10] Hz, over the whole session [0 dur]
%     'Seed'       []      default: sum(double(file).*(1:numel(file))) mod (2^31-1),
%                          the syncBatch rule, so the subset draw is the same
%                          function of the path
%
%   OUTPUT R
%     .ok .why .file
%     .t        M x 1 spike time (s from recording start), observed order
%     .unit     M x 1 index into .detIdx of the unit that fired
%     .idx      M x 1 rank of the spike within its unit's train
%     .count    M x 1 observed coincidence count (0..N-1), uint8 (uint16 if N > 256)
%     .null     M x NSurr  jittered coincidence counts, same class
%     .N        units used;  .nPartners = N-1
%     .C        scalar observed SPIKE-synchronization (mean of C over spikes)
%     .CnullMean, .CnullSd   scalar C of the surrogates
%     .keepIdx  units in the rate band;  .detIdx the N used (indices into spikes.times)
%     .uid      spikes.UID(detIdx) when present
%     .dur, .params
%
%   See also SPIKESYNC, RUNSYNCPROFILES.

p = inputParser;
p.addParameter('Jitter', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('NSurr', 100, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('MatchN', 10);
p.addParameter('RateRange', [0.1 10]);
p.addParameter('Seed', []);
p.parse(varargin{:});
o = p.Results;

R = struct('ok',false,'why','','file',spikeFile,'t',[],'unit',[],'idx',[],'count',[], ...
    'null',[],'N',0,'nPartners',NaN,'C',NaN,'CnullMean',NaN,'CnullSd',NaN, ...
    'keepIdx',[],'detIdx',[],'uid',[],'dur',NaN,'params',o);

S = load(spikeFile, 'spikes');
if ~isfield(S,'spikes') || ~isfield(S.spikes,'times'), R.why = 'no spikes.times'; return; end
ts = S.spikes.times(:);
ts = cellfun(@(x) unique(x(~isnan(x(:)))), ts, 'UniformOutput', false);   % sorted, no ties
lastSpk = cellfun(@(x) max([x; -Inf]), ts);
dur = max(lastSpk(isfinite(lastSpk)));
if isempty(dur) || ~isfinite(dur) || dur <= 0, R.why = 'bad duration'; return; end
R.dur = dur;

% ---- units: rate band over the whole session, then a seeded MatchN draw ----
M = cellfun(@numel, ts);
FR = M / dur;
keep = FR > o.RateRange(1) & FR < o.RateRange(2) & M > 0;
R.keepIdx = find(keep);
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
N = numel(trains); R.N = N; R.nPartners = N - 1;
if N < 2, R.why = 'fewer than 2 units'; return; end
cls = 'uint8'; if N - 1 > 255, cls = 'uint16'; end
iv = [0 dur];

% ---- observed profile ------------------------------------------------------
[C, ~, prof] = spikeSync(trains, 'Interval', iv);
Mtot = numel(prof.t);
R.C = C;
R.t = prof.t(:);
R.unit = uint8(prof.train(:));
R.idx = uint32(prof.idx(:));
R.count = cast(round(prof.C(:) * (N - 1)), cls);
pos = cell(N, 1);                       % pos{n}(i) = row of unit n's i-th spike
for n = 1:N
    k = find(prof.train == n);
    pos{n} = zeros(numel(trains{n}), 1);
    pos{n}(prof.idx(k)) = k;
end

% ---- interval-jitter surrogates --------------------------------------------
nS = o.NSurr; J = o.Jitter;
cols = cell(1, nS); Cs = nan(nS, 1);
parfor s = 1:nS
    rsS = RandStream('mt19937ar', 'Seed', seed + s);
    tj = cell(N, 1);
    for n = 1:N
        x = trains{n};
        lo = floor(x / J) * J;
        hi = min(lo + J, dur);
        y = sort(lo + rand(rsS, numel(x), 1) .* (hi - lo));
        d = find(diff(y) <= 0);                 % vanishingly rare tie: keep order defined
        while ~isempty(d)
            y(d + 1) = y(d) + 1e-9; y = sort(y); d = find(diff(y) <= 0); %#ok<FXUP>
        end
        tj{n} = y;
    end
    [Cj, ~, pj] = spikeSync(tj, 'Interval', iv);
    col = zeros(Mtot, 1, cls);
    cj = round(pj.C(:) * (N - 1));
    for n = 1:N
        k = find(pj.train == n);
        col(pos{n}(pj.idx(k))) = cj(k);
    end
    cols{s} = col; Cs(s) = Cj;
end
if nS > 0, R.null = [cols{:}]; else, R.null = zeros(Mtot, 0, cls); end
R.CnullMean = mean(Cs); R.CnullSd = std(Cs);
R.params.Seed = seed; R.params.cls = cls;
R.ok = true;
end
