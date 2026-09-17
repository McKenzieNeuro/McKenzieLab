function out = boutRate(r, width, step)
%BOUTRATE  Mean per-neuron firing rate during significant transients vs overall.
%
%   out = boutRate(r, width, step)
%
%   r          : ONE session struct from syncBatch (run with DetectTransients=true)
%   width,step : the SAME 'Width'/'Step' values used to produce r.transients
%
%   syncBatch never stores raw spike times in its output (only derived
%   summaries), so this reloads spikes.times from r.file and reconstructs
%   the EXACT unit subset actually used for detection via r.detIdx (indices
%   into the raw spikes.times cell array -- see syncBatch's local_session:
%   keepIdx = units passing the rate filter, detIdx = the seeded MatchN draw
%   from keepIdx). Using r.rate here would be wrong: r.rate covers all
%   Nfull rate-qualifying units, not the specific MatchN subset that
%   generated r.transients.
%
%   If the session was duration-truncated (r.durTruncated), spikes past
%   r.dur are dropped here too -- otherwise "overall" rate would be computed
%   over a longer window than the one transients were actually detected in.
%
%   BOUT WINDOWS: a significant transient's TRUE elapsed interval is
%   [tStart - width/2, tEnd + width/2], not [tStart, tEnd] -- tStart/tEnd in
%   r.transients are WINDOW-CENTER times of the first/last window in the
%   cluster (see local_transientsFromSurrogates in syncBatch.m), the same
%   convention transientStats.m uses for duration. Adjacent clusters can have
%   overlapping or touching true intervals after this padding; they are
%   merged before counting so no bout time or spike is double-counted.
%
%   OUTPUT (all rates in Hz, per neuron, pooled across the N detection units)
%     rateOverall : total detection-unit spikes / (N * r.dur)
%     rateInBout  : spikes inside a merged significant-bout interval, divided
%                   by (N * total bout duration). NaN if there are no
%                   significant transients -- same censoring as
%                   transientStats' duration and meanPeakMagnitude's peakDC.
%     ratio       : rateInBout / rateOverall. Self-normalizing: controls for
%                   any baseline rate difference between datasets, answers
%                   "how much do these neurons speed up during a bout" rather
%                   than "how fast are these neurons overall."
%     boutSeconds : total merged bout duration -- sanity-check how much data
%                   rateInBout rests on (one short bout = a noisy estimate).
%     N           : detection units actually found in the file (should equal
%                   numel(r.detIdx); a mismatch is flagged, not silently
%                   ignored -- it means the file changed since syncBatch ran).
%
%   See also SYNCBATCH, TRANSIENTSTATS, MEANPEAKMAGNITUDE.

S = load(r.file, 'spikes');
allTs = S.spikes.times(:);

if isempty(r.detIdx)
    error('boutRate:noDetIdx', ...
        'r.detIdx is empty -- was this session run with MatchN set?');
end
if max(r.detIdx) > numel(allTs)
    error('boutRate:mismatch', ...
        '%s: detIdx references unit %d but the file only has %d units -- ' , ...
        'file has changed since syncBatch ran on it.', r.file, max(r.detIdx), numel(allTs));
end

ts = allTs(r.detIdx);
N = numel(ts);
if N ~= numel(r.detIdx)
    warning('boutRate:countMismatch','%s: expected %d detection units, got %d', ...
        r.file, numel(r.detIdx), N);
end

% Align raw spikes to the SAME time base as r.transients / r.tc.
% When syncBatch ran with 'Window', it re-zeroed spike times to the window
% start, so transient tStart/tEnd are window-relative while the spikes
% reloaded here are in absolute recording time. Without this shift every
% bout-membership test below would compare two different clocks and return
% near-zero in-bout spike counts -- wrong, and quietly so.
if ~isempty(r.windowAbs)
    w0 = r.windowAbs(1);
    ts = cellfun(@(x) x(x >= r.windowAbs(1) & x < r.windowAbs(2)) - w0, ts, ...
        'UniformOutput', false);
elseif r.durTruncated
    ts = cellfun(@(x) x(x <= r.dur), ts, 'UniformOutput', false);
end

totalSpikes = sum(cellfun(@numel, ts));
out.rateOverall = totalSpikes / (N * r.dur);
out.N = N;

sig = [r.transients.significant];
if ~any(sig)
    out.rateInBout = NaN;
    out.ratio = NaN;
    out.boutSeconds = 0;
    return
end

trans = r.transients(sig);
starts = [trans.tStart] - width/2;
ends   = [trans.tEnd]   + width/2;
starts = max(starts, 0);
ends   = min(ends, r.dur);

% merge overlapping/touching padded intervals
[starts, ord] = sort(starts);
ends = ends(ord);
ms = starts(1); me = ends(1);
mergedS = []; mergedE = [];
for k = 2:numel(starts)
    if starts(k) <= me
        me = max(me, ends(k));
    else
        mergedS(end+1) = ms; mergedE(end+1) = me; %#ok<AGROW>
        ms = starts(k); me = ends(k);
    end
end
mergedS(end+1) = ms; mergedE(end+1) = me; %#ok<AGROW>

boutSeconds = sum(mergedE - mergedS);
inBoutSpikes = 0;
for n = 1:N
    x = ts{n};
    for k = 1:numel(mergedS)
        inBoutSpikes = inBoutSpikes + sum(x >= mergedS(k) & x < mergedE(k));
    end
end

out.rateInBout = inBoutSpikes / (N * boutSeconds);
out.ratio = out.rateInBout / out.rateOverall;
out.boutSeconds = boutSeconds;
end