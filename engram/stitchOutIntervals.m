function [tsOut, durOut] = stitchOutIntervals(ts, dur, excludeIntervals, pad)
%STITCHOUTINTERVALS  Remove intervals from spike trains, closing the gaps.
%
%   [tsOut, durOut] = stitchOutIntervals(ts, dur, excludeIntervals, pad)
%
%   ts               : {N x 1} cell of spike times (s)
%   dur              : recording duration (s)
%   excludeIntervals : K x 2 [start stop] to remove (e.g. light pulses)
%   pad              : seconds to additionally exclude on each side of every
%                       interval (light artifact, opsin-driven afterdischarge)
%
%   Excluded windows are not merely blanked -- they are CLOSED, i.e. spikes
%   after an excluded window are shifted earlier by that window's duration,
%   so remaining spikes form a continuous stitched timeline. This matters
%   specifically for SPIKE-Synchronization: its coincidence window derives
%   from local inter-spike intervals (Eq. 15 in Kreuz et al.), and simply
%   deleting spikes in the excluded window while leaving the gap in place
%   would let that ISI-based window bridge across the gap, silently
%   widening the effective coincidence criterion right at every pulse
%   boundary. Stitching removes the gap itself, so the ISI structure the
%   surrounding spikes see is the same as if the excluded period had never
%   existed.
%
%   Returns durOut = dur minus total excluded time (after merging overlaps
%   and clipping to [0, dur]).
%
%   See also OPTOPULSEINTERVALS, SYNCBATCH.

if isempty(excludeIntervals)
    tsOut = ts; durOut = dur; return
end

ex = sortrows(excludeIntervals(:,1:2));
ex(:,1) = max(ex(:,1) - pad, 0);
ex(:,2) = min(ex(:,2) + pad, dur);

% merge overlapping/touching padded intervals
ms = ex(1,1); me = ex(1,2); merged = zeros(0,2);
for k = 2:size(ex,1)
    if ex(k,1) <= me
        me = max(me, ex(k,2));
    else
        merged(end+1,:) = [ms me]; %#ok<AGROW>
        ms = ex(k,1); me = ex(k,2);
    end
end
merged(end+1,:) = [ms me];
merged = merged(merged(:,2) > merged(:,1), :);

if isempty(merged), tsOut = ts; durOut = dur; return; end

starts = merged(:,1); stops = merged(:,2);

tsOut = cell(size(ts));
for n = 1:numel(ts)
    x = ts{n}(:);
    inExcl = false(size(x));
    for k = 1:size(merged,1)
        inExcl = inExcl | (x >= starts(k) & x < stops(k));
    end
    x = x(~inExcl);
    tsOut{n} = local_shift(x, starts, stops);
end

durOut = dur - sum(stops - starts);
end


function y = local_shift(x, starts, stops)
% Shift each remaining spike time earlier by the total excluded duration
% strictly before it (spikes inside an excluded window are assumed already
% removed by the caller).
if isempty(x), y = x; return; end
excl_before = zeros(size(x));
for k = 1:numel(starts)
    excl_before = excl_before + (stops(k) - starts(k)) .* (x >= stops(k));
end
y = x - excl_before;
end
