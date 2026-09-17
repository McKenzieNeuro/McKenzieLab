function m = meanPeakMagnitude(trans)
%MEANPEAKMAGNITUDE  Mean peak dC (height above surrogate median) over
%   significant transients in one session's syncBatch output.
%
%   m = meanPeakMagnitude(trans)   where trans = <session>.transients
%
%   Returns NaN if the session has no significant transients -- there is no
%   magnitude to average, same censoring as transientStats' meanDur. Any
%   group comparison on this quantity is conditional on a transient having
%   occurred, exactly like duration, and should be read with the same caveat.
%
%   peakDC is always >= 0 by construction (it's Ct minus the surrogate
%   median, evaluated only inside a window that exceeded the surrogate's
%   UPPER band, so it's guaranteed positive relative to the median too).
%
%   This lives in its own file, not as a local function inside
%   runSyncAnalysis.m, for the same reason optoPulseIntervals was moved:
%   local functions in a script are only in scope while that script runs.
%
%   See also TRANSIENTSTATS, NESTINGCHECK, LMECOMPARE.

if isempty(trans), m = NaN; return; end
sig = [trans.significant];
if ~any(sig), m = NaN; return; end
m = mean([trans(sig).peakDC]);
end
