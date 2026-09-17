function [rOverall, rInBout, ratio, boutSec] = boutRateBatch(T, width, step)
%BOUTRATEBATCH  boutRate applied across every session in a syncBatch output.
%
%   [rOverall, rInBout, ratio, boutSec] = boutRateBatch(T, width, step)
%
%   T          : the .ok-filtered struct array from syncBatch (e.g. `a`, `b`)
%   width,step : same values used for detection
%
%   Reloads raw spikes per session (boutRate does), so this is slower than
%   the other batch metrics, which only read syncBatch's already-computed
%   output. Not expensive on the scale of the original surrogate computation
%   -- no spikeSync calls here, just spike-time loading and counting -- but
%   not free either across 100+ sessions on a network drive.
%
%   See also BOUTRATE.

n = numel(T);
rOverall = nan(n,1); rInBout = nan(n,1); ratio = nan(n,1); boutSec = nan(n,1);
for i = 1:n
    try
        o = boutRate(T(i), width, step);
        rOverall(i) = o.rateOverall;
        rInBout(i)  = o.rateInBout;
        ratio(i)    = o.ratio;
        boutSec(i)  = o.boutSeconds;
    catch ME
        warning('boutRateBatch:session','session %d (%s): %s', i, T(i).file, ME.message);
    end
end
end
