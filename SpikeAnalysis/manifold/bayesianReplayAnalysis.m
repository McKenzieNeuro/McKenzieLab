function out = bayesianReplayAnalysis(spikes, behavior, ripples, opts)
% bayesianReplayAnalysis  Bayesian replay decoding + shuffle significance over ripples.
%
%   out = bayesianReplayAnalysis(spikes, behavior, ripples, opts)
%
% Pipeline (Davidson, Kloosterman & Wilson 2009; Karlsson & Frank 2009;
% Grosmark & Buzsaki 2016):
%   1. Build 1D linearized place fields from RUN periods.
%   2. For each ripple, bin spikes (tau), Bayesian-decode position per bin.
%   3. Score the posterior for sequential structure (weighted corr + line/radon).
%   4. Compare to shuffles -> per-event p-value and z-score.
%
% INPUTS
%   spikes   : Buzcode spikes struct (uses spikes.times) OR a {nCells x 1} cell of
%              spike-time vectors.
%   behavior : Buzcode behavior struct (uses .timestamps and .position.lin) OR []
%              if you instead pass opts.posTime / opts.posLin.
%   ripples  : Buzcode ripples struct (uses .timestamps) OR an [N x 2] matrix of
%              [start stop] times (s).  These are the candidate events.
%   opts     : struct (all optional). Key fields:
%       .cellIdx       indices of cells to use (e.g. pyramidal); default all.
%       .runEpoch      [start stop] used to BUILD place fields; default full.
%       .rateMaps,.binCenters  precomputed fields; if given, step 1 is skipped.
%       .posTime,.posLin       override behavior fields directly.
%       --- decoding / fields ---
%       .tau (0.02 s) .binSize (2 cm) .smoothSD (4 cm) .speedThresh (5 cm/s)
%       .rateFloor (1e-2 Hz)
%       --- event window ---
%       .peakWindow  fixed analysis window around each ripple peak, in seconds:
%                    scalar w -> peak +/- w/2 ; [pre post] -> peak-pre..peak+post ;
%                    [] (default) -> use detected ripple [start stop].
%                    (Buzcode ripples.peaks is used as the anchor; if you pass an
%                     [N x 2] matrix, the window midpoint is used instead.)
%       --- event inclusion ---
%       .minActiveCells (5) .minBins (5) .minSpkBins (5)
%       --- scoring / significance ---
%       .scoreForP    'wcorr' (default) | 'line'   (which score drives the p-value)
%       .shuffleMethod 'columnCycle' (default) | 'permute'
%       .nShuffles    (500)
%       .lineBand (15 cm) .nSlope (120) .computeLine (true)
%       .savePosteriors (true)
%       .verbose (true)
%
% OUTPUT struct out:
%   .results    : table, one row per candidate event (see column names below).
%   .posteriors : {N x 1} cell of [nPos x nT] posteriors (empty if excluded /
%                 if savePosteriors=false).
%   .rateMaps, .binCenters : the place fields used.
%   .opts       : resolved options.
%
% Table columns: eventID peakTime startTime stopTime duration nActiveCells
%   nTimeBins nSpikeBins included weightedCorr lineScore lineSlope_cmps
%   pValue zScore scoreType
%   (startTime/stopTime are the WINDOW actually decoded -- equal to the padded
%    peakWindow when set, otherwise the detected ripple bounds.)
%
% NOTE on units & direction:
%   * posLin / binSize / speedThresh must share units (cm assumed). If your
%     behavior.position.lin is normalized or in trial units, rescale first.
%   * The sign of weightedCorr / lineSlope encodes decoded direction along the
%     LINEARIZED axis; forward vs reverse replay depends on how lin is oriented
%     relative to the animal's travel during RUN. Interpret with that in mind.

if nargin < 4, opts = struct(); end
opts = setdefault(opts, 'tau',            0.02);
opts = setdefault(opts, 'binSize',        2);
opts = setdefault(opts, 'smoothSD',       4);
opts = setdefault(opts, 'speedThresh',    5);
opts = setdefault(opts, 'rateFloor',      1e-2);
opts = setdefault(opts, 'peakWindow',     []);    % [] = use ripple bounds; see below
opts = setdefault(opts, 'minActiveCells', 5);
opts = setdefault(opts, 'minBins',        5);
opts = setdefault(opts, 'minSpkBins',     5);
opts = setdefault(opts, 'scoreForP',      'wcorr');
opts = setdefault(opts, 'shuffleMethod',  'columnCycle');
opts = setdefault(opts, 'nShuffles',      500);
opts = setdefault(opts, 'lineBand',       15);
opts = setdefault(opts, 'nSlope',         120);
opts = setdefault(opts, 'computeLine',    true);
opts = setdefault(opts, 'savePosteriors', true);
opts = setdefault(opts, 'verbose',        true);

% ---- resolve spikes ----
if isstruct(spikes) && isfield(spikes, 'times')
    spikeTimes = spikes.times(:);
elseif iscell(spikes)
    spikeTimes = spikes(:);
else
    error('spikes must be a Buzcode spikes struct (.times) or a cell array.');
end
if isfield(opts, 'cellIdx') && ~isempty(opts.cellIdx)
    spikeTimes = spikeTimes(opts.cellIdx);
end
nCells = numel(spikeTimes);

% ---- resolve behavior ----
if isfield(opts,'posTime') && ~isempty(opts.posTime)
    posTime = opts.posTime(:); posLin = opts.posLin(:);
elseif isstruct(behavior)
    posTime = behavior.timestamps(:);
    if isfield(behavior,'position') && isfield(behavior.position,'lin')
        posLin = behavior.position.lin(:);
    elseif isfield(behavior,'position') && isfield(behavior.position,'linearized')
        posLin = behavior.position.linearized(:);
    else
        error('behavior.position.lin not found; pass opts.posLin instead.');
    end
else
    posTime = []; posLin = [];   % allowed only if rateMaps supplied
end

% ---- resolve ripples (timestamps + optional peak anchor) ----
peaks = [];
if isstruct(ripples) && isfield(ripples, 'timestamps')
    evt = ripples.timestamps;
    if isfield(ripples, 'peaks') && ~isempty(ripples.peaks)
        peaks = ripples.peaks(:);
    end
elseif isnumeric(ripples) && size(ripples,2) == 2
    evt = ripples;
else
    error('ripples must be a Buzcode ripples struct (.timestamps) or an [N x 2] matrix.');
end
N = size(evt, 1);
if isempty(peaks), peaks = mean(evt, 2); end       % fall back to window midpoint

% ---- optional fixed analysis window centered on each ripple ----
% opts.peakWindow = w        -> peak +/- w/2   (total width w, seconds)
%                 = [pre post]-> peak-pre .. peak+post  (asymmetric, seconds)
%                 = []        -> use the detected ripple [start stop] (default)
% Use this when raw ripple bounds are too short for >= minBins decoding bins
% (e.g. a 60-100 ms window at tau=20 ms gives 3-5 bins). The window pulls in
% spikes just outside the detected ripple and may overlap neighbors if events
% are <peakWindow apart; that is expected for candidate-event decoding.
if ~isempty(opts.peakWindow)
    w = opts.peakWindow(:).';
    if isscalar(w), pre = w/2; post = w/2; else, pre = w(1); post = w(2); end
    winStart = peaks - pre;
    winStop  = peaks + post;
else
    winStart = evt(:,1);
    winStop  = evt(:,2);
end

% ---- place fields ----
if isfield(opts,'rateMaps') && ~isempty(opts.rateMaps)
    rateMaps = opts.rateMaps; binCenters = opts.binCenters;
    if isfield(opts,'cellIdx') && ~isempty(opts.cellIdx) && size(rateMaps,1) ~= nCells
        rateMaps = rateMaps(opts.cellIdx, :);
    end
else
    if isempty(posTime), error('No behavior provided and no precomputed rateMaps.'); end
    fopts = struct('binSize',opts.binSize,'smoothSD',opts.smoothSD, ...
                   'speedThresh',opts.speedThresh,'rateFloor',opts.rateFloor);
    if isfield(opts,'runEpoch') && ~isempty(opts.runEpoch), fopts.epoch = opts.runEpoch; end
    [rateMaps, binCenters] = computeRateMaps1D(spikeTimes, posTime, posLin, fopts);
end
nPos = size(rateMaps, 2);

sopts = struct('computeLine',opts.computeLine,'lineBand',opts.lineBand,'nSlope',opts.nSlope);

% ---- preallocate ----
eventID      = (1:N).';
peakTime     = peaks;
startTime    = winStart;
stopTime     = winStop;
duration     = stopTime - startTime;
nActiveCells = zeros(N,1);
nTimeBins    = zeros(N,1);
nSpikeBins   = zeros(N,1);
included     = false(N,1);
weightedCorr = nan(N,1);
lineScore    = nan(N,1);
lineSlope    = nan(N,1);
pValue       = nan(N,1);
zScore       = nan(N,1);
posteriors   = cell(N,1);

for e = 1:N
    edges = startTime(e) : opts.tau : stopTime(e);
    if numel(edges) < 2, continue; end
    nT = numel(edges) - 1;

    countMat = zeros(nCells, nT);
    for c = 1:nCells
        spk = spikeTimes{c};
        countMat(c,:) = histcounts(spk, edges);
    end

    nActiveCells(e) = sum(any(countMat > 0, 2));
    nSpikeBins(e)   = sum(sum(countMat,1) > 0);
    nTimeBins(e)    = nT;
    included(e) = nT >= opts.minBins && ...
                  nActiveCells(e) >= opts.minActiveCells && ...
                  nSpikeBins(e) >= opts.minSpkBins;
    if ~included(e), continue; end

    posterior = bayesianDecode(countMat, rateMaps, opts.tau);
    S = scoreTrajectory(posterior, binCenters, opts.tau, sopts);
    weightedCorr(e) = S.weightedCorr;
    lineScore(e)    = S.lineScore;
    lineSlope(e)    = S.lineSlope;
    if opts.savePosteriors, posteriors{e} = posterior; end

    % observed score that drives significance
    obs = pickScore(S, opts.scoreForP);
    shvals = zeros(opts.nShuffles, 1);
    needLine = strcmpi(opts.scoreForP, 'line');
    shOpts = sopts; shOpts.computeLine = needLine;   % skip line score in shuffles unless needed
    for s = 1:opts.nShuffles
        Psh = shufflePosterior(posterior, opts.shuffleMethod);
        Ssh = scoreTrajectory(Psh, binCenters, opts.tau, shOpts);
        shvals(s) = pickScore(Ssh, opts.scoreForP);
    end
    pValue(e) = (1 + sum(abs(shvals) >= abs(obs))) / (1 + opts.nShuffles);
    zScore(e) = (abs(obs) - mean(abs(shvals))) / (std(abs(shvals)) + eps);

    if opts.verbose && mod(e, 25) == 0
        fprintf('  event %d/%d\n', e, N);
    end
end

scoreType = repmat({opts.scoreForP}, N, 1);
results = table(eventID, peakTime, startTime, stopTime, duration, nActiveCells, ...
    nTimeBins, nSpikeBins, included, weightedCorr, lineScore, lineSlope, ...
    pValue, zScore, scoreType, ...
    'VariableNames', {'eventID','peakTime','startTime','stopTime','duration', ...
    'nActiveCells','nTimeBins','nSpikeBins','included','weightedCorr', ...
    'lineScore','lineSlope_cmps','pValue','zScore','scoreType'});

out = struct('results', results, 'posteriors', {posteriors}, ...
             'rateMaps', rateMaps, 'binCenters', binCenters, 'opts', opts);

if opts.verbose
    nInc = sum(included);
    nSig = sum(included & pValue < 0.05);
    fprintf('Decoded %d/%d candidate events; %d significant at p<0.05 (%s, %s shuffle).\n', ...
        nInc, N, nSig, opts.scoreForP, opts.shuffleMethod);
end
end

% ---------------------------------------------------------------------------
function s = setdefault(s, f, val)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = val; end
end

function v = pickScore(S, which)
switch lower(which)
    case 'wcorr', v = S.weightedCorr;
    case 'line',  v = S.lineScore;
    otherwise, error('scoreForP must be ''wcorr'' or ''line''.');
end
end

function Q = shufflePosterior(P, method)
[nPos, nT] = size(P);
switch lower(method)
    case 'columncycle'                       % per-bin random circular shift
        Q = P;
        for j = 1:nT
            Q(:,j) = circshift(P(:,j), randi(nPos) - 1);
        end
    case 'permute'                           % shuffle temporal order of bins
        Q = P(:, randperm(nT));
    otherwise
        error('shuffleMethod must be ''columnCycle'' or ''permute''.');
end
end
