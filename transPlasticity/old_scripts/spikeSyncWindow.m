function [Ct, tc, nSpk, Cnull] = spikeSyncWindow(prof, varargin)
%SPIKESYNCWINDOW  Time-resolved SPIKE-Synchronization from the pooled profile.
%
%   [Ct, tc, nSpk, Cnull] = spikeSyncWindow(prof, ...)
%
%   Computes C over sliding windows EXACTLY, by re-summing the numerator and
%   denominator of Eq. 19 within each window:
%
%       C(window) = sum_{k: t_k in window} C(t_k) / #{k: t_k in window}
%
%   This is not a smoothed version of prof.C. Convolving prof.C with a kernel
%   would weight every spike equally regardless of how many spikes fall in the
%   window; re-summing keeps the spike-count weighting that defines C, so the
%   windowed values are the same quantity as the scalar, just restricted in
%   time. Windows are half-open (a,b].
%
%   INPUT
%     prof : the third output of spikeSync.
%
%   NAME/VALUE OPTIONS
%     'Width'     : window width in the time units of the spike data. Required.
%     'Step'      : window step (default Width/4).
%     'Range'     : [t0 t1] over which windows are placed (default prof.t range).
%     'MinSpikes' : windows with fewer pooled spikes return NaN (default 10).
%                   Do not set this to 0. C is defined as 1 for an empty set
%                   (common silence = perfect synchrony), which is defensible
%                   for a whole recording and actively misleading in a trace:
%                   quiescent windows become synchrony peaks. Low-count windows
%                   are also high-variance, since the denominator IS the spike
%                   count.
%
%   OUTPUT
%     Ct     : nWin x 1 windowed SPIKE-Synchronization (NaN where under MinSpikes)
%     tc     : nWin x 1 window centres
%     nSpk   : nWin x 1 pooled spike count per window
%     Cnull  : [] unless surrogates are requested via spikeSyncSurrogate.
%
%   CAVEAT worth stating in any figure legend that uses this: the window width
%   is a free time scale, and the whole point of this measure class is that the
%   coincidence detection itself has none. Windowing does not corrupt the
%   coincidence criterion (that is still rate-adaptive per spike pair), it only
%   sets the resolution at which you read the result out. Choose Width from the
%   process you are asking about, not from the data.
%
%   See also SPIKESYNC, SPIKESYNCSURROGATE.

p = inputParser;
p.addParameter('Width', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Step', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Range', [], @(x) isempty(x) || numel(x) == 2);
p.addParameter('MinSpikes', 10, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.parse(varargin{:});
W = p.Results.Width;
if isempty(W), error('spikeSyncWindow:width','''Width'' is required.'); end
step = p.Results.Step; if isempty(step), step = W/4; end
minSpk = p.Results.MinSpikes;

t = prof.t(:); c = prof.C(:);
if isempty(t), Ct = []; tc = []; nSpk = []; Cnull = []; return; end

rng_ = p.Results.Range;
if isempty(rng_), rng_ = [min(t) max(t)]; else, rng_ = sort(rng_(:))'; end

a = (rng_(1) : step : rng_(2) - W)';
if isempty(a), a = rng_(1); end
b = a + W;
tc = a + W/2;

% cumulative sums -> O(1) per window
cs = [0; cumsum(c)];
kA = local_countLE(t, a);
kB = local_countLE(t, b);
nSpk = kB - kA;
Ct = (cs(kB+1) - cs(kA+1)) ./ nSpk;
Ct(nSpk < minSpk) = NaN;
Cnull = [];

end


function k = local_countLE(t, q)
% #{t <= q} for sorted t, vectorised over q, safe against duplicate times
% (pooled spike times tie whenever two trains fire at the same instant).
nt = numel(t); nq = numel(q);
v    = [t(:); q(:)];
flag = [zeros(nt,1); ones(nq,1)];
[~, ord] = sortrows([v flag]);        % ties: t before q, so "<=" is respected
isq = flag(ord) == 1;
run = cumsum(~isq);                   % running count of t elements
k = zeros(nq,1);
k(ord(isq) - nt) = run(isq);
end
