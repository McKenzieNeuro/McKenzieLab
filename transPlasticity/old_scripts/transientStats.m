function [rate, duty, meanDur, allDur] = transientStats(r, step, width)
%TRANSIENTSTATS  Rate, duty cycle, and duration from one syncBatch session.
%
%   [rate, duty, meanDur, allDur] = transientStats(r, step, width)
%
%   r      : one element of a syncBatch output struct array (must have
%            .transients and .dur)
%   step   : the 'Step' value used for detection (seconds)
%   width  : the 'Width' value used for detection (seconds)
%
%   OUTPUT (all counting SIGNIFICANT clusters only)
%     rate    : transients per hour
%     duty    : fraction of the session spent in a significant transient
%     meanDur : mean transient duration in seconds (NaN if none)
%     allDur  : per-transient durations in seconds (empty if none)
%
%   Duration of an nWin-window cluster is (nWin-1)*step + width, not
%   nWin*step: consecutive windows overlap by (width-step), so a cluster's
%   real elapsed time is width (for the first window) plus step for each
%   additional window in the run. Using nWin*step alone overstates duration
%   by (width-step) per cluster; using nWin*width overstates it far more.
%
%   Pulling this into one function exists because the same formula was
%   written inline four times across an analysis and got the wrong variable
%   (Width instead of Step) once. One place to get it right.

if isempty(r.transients)
    rate = 0; duty = 0; meanDur = NaN; allDur = [];
    return
end

sig = [r.transients.significant];
nWin = [r.transients(sig).nWin];

allDur = (nWin - 1) * step + width;
rate = numel(allDur) / (r.dur / 3600);
duty = sum(allDur) / r.dur;
meanDur = mean(allDur);   % NaN for empty allDur, as intended
end
