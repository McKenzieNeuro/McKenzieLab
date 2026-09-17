function [rateMaps, binCenters, occupancy, info] = computeRateMaps1D(spikeTimes, posTime, posLin, opts)
% computeRateMaps1D  Linearized 1D place fields (rate maps) from running periods.
%
%   [rateMaps, binCenters, occupancy, info] = computeRateMaps1D(spikeTimes, posTime, posLin, opts)
%
% INPUTS
%   spikeTimes : {nCells x 1} cell array, each a vector of spike times (s).
%                (Buzcode: pass spikes.times.)
%   posTime    : [nSamp x 1] behavior timestamps (s).  (Buzcode: behavior.timestamps)
%   posLin     : [nSamp x 1] linearized position, SAME UNITS as opts.binSize / speedThresh
%                (assumed cm here).            (Buzcode: behavior.position.lin)
%   opts       : struct, fields (all optional):
%       .binSize     spatial bin width            (default 2  cm)
%       .smoothSD    Gaussian SD for smoothing    (default 4  cm)
%       .speedThresh run inclusion speed          (default 5  cm/s)
%       .epoch       [start stop] to build fields (default full posTime range)
%       .rateFloor   minimum rate, avoids log(0)  (default 1e-2 Hz)
%       .speed       precomputed speed aligned to posTime (default [], computed here)
%
% OUTPUTS
%   rateMaps   : [nCells x nPosBins] firing rate (Hz), floored at rateFloor.
%   binCenters : [1 x nPosBins] position bin centers.
%   occupancy  : [1 x nPosBins] dwell time during running (s).
%   info       : struct with speed, runMask, edges, peakRate, nSpikesUsed.
%
% Method: smooth spike counts and occupancy with the same kernel, then divide
% (standard occupancy-normalized rate map; edge attenuation cancels in ratio).

if nargin < 4, opts = struct(); end
opts = setdefault(opts, 'binSize',     2);
opts = setdefault(opts, 'smoothSD',    4);
opts = setdefault(opts, 'speedThresh', 5);
opts = setdefault(opts, 'rateFloor',   1e-2);
opts = setdefault(opts, 'speed',       []);

posTime = posTime(:); posLin = posLin(:);
if ~isfield(opts,'epoch') || isempty(opts.epoch)
    opts.epoch = [posTime(1) posTime(end)];
end

% defensively sort/unique time for interpolation
[posTime, ord] = unique(posTime, 'stable');
posLin = posLin(ord);
[posTime, si] = sort(posTime); posLin = posLin(si);

% restrict to epoch
ein = posTime >= opts.epoch(1) & posTime <= opts.epoch(2);
posTime = posTime(ein); posLin = posLin(ein);
dt = median(diff(posTime));

% speed
if isempty(opts.speed)
    v = abs(diff(posLin)) ./ diff(posTime);
    v = [v(1); v];                                  % pad to length
    v = gaussSmooth1D(v(:).', max(1, round(0.1/dt))).';   % ~100 ms smoothing
else
    v = opts.speed(:); v = v(ein);
end
runMask = v >= opts.speedThresh;

% position bins
lo = floor(min(posLin)); hi = ceil(max(posLin));
edges = lo : opts.binSize : (hi + opts.binSize);
binCenters = edges(1:end-1) + opts.binSize/2;
nPos = numel(binCenters);
nCells = numel(spikeTimes);

% occupancy during running
occRaw = histcounts(posLin(runMask), edges) * dt;    % 1 x nPos (s)

% spike counts during running, mapped to position
spikeCount = zeros(nCells, nPos);
nSpikesUsed = zeros(nCells, 1);
for c = 1:nCells
    spk = spikeTimes{c}(:);
    spk = spk(spk >= opts.epoch(1) & spk <= opts.epoch(2));
    if isempty(spk), continue; end
    spkPos = interp1(posTime, posLin, spk, 'linear', NaN);
    spkVel = interp1(posTime, v,      spk, 'linear', NaN);
    keep = ~isnan(spkPos) & spkVel >= opts.speedThresh;
    spikeCount(c,:) = histcounts(spkPos(keep), edges);
    nSpikesUsed(c)  = sum(keep);
end

% smooth-then-divide
sdBins = opts.smoothSD / opts.binSize;
smOcc = gaussSmooth1D(occRaw, sdBins);
rateMaps = zeros(nCells, nPos);
for c = 1:nCells
    smSC = gaussSmooth1D(spikeCount(c,:), sdBins);
    r = smSC ./ smOcc;
    r(~isfinite(r) | smOcc <= 0) = 0;
    rateMaps(c,:) = r;
end
rateMaps = max(rateMaps, opts.rateFloor);

occupancy = occRaw;
info = struct('speed', v, 'runMask', runMask, 'edges', edges, ...
              'dt', dt, 'peakRate', max(rateMaps, [], 2), ...
              'nSpikesUsed', nSpikesUsed, 'opts', opts);
end

% ---------------------------------------------------------------------------
function s = setdefault(s, f, val)
if ~isfield(s, f) || isempty(s.(f)), s.(f) = val; end
end

function y = gaussSmooth1D(x, sdBins)
if sdBins <= 0, y = x; return; end
r = ceil(3*sdBins);
k = exp(-((-r:r).^2) / (2*sdBins^2));
k = k / sum(k);
y = conv(x, k, 'same');
end
