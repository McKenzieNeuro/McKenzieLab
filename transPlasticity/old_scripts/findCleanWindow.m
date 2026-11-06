function w = findCleanWindow(excludeIntervals, dur, L, prefer)
%FINDCLEANWINDOW  A length-L window containing none of the excluded intervals.
%
%   w = findCleanWindow(excludeIntervals, dur, L, prefer)
%
%   excludeIntervals : K x 2 [start stop] to avoid (e.g. opto pulses, padded)
%   dur   : recording duration (s)
%   L     : required window length (s)
%   prefer: 'first' (default), 'last', or 'longest' -- which clean gap to
%           take the window from when several qualify
%
%   Returns [t0 t1], or [] if no gap of length L exists.
%
%   WHY NOT JUST MASK PULSE WINDOWS. syncBatch's 'Exclude' NaNs any analysis
%   window overlapping a pulse, which is right for a whole-session trace but
%   wrong when comparing epochs ACROSS datasets: the masked epoch ends up
%   with fewer valid windows than an unmasked one, and transient-cluster
%   significance depends on window count. Selecting a genuinely clean gap
%   keeps the window count identical on both sides.

if nargin < 4, prefer = 'first'; end
w = [];
if ~isfinite(dur) || dur < L, return; end

if isempty(excludeIntervals)
    gaps = [0 dur];
else
    e = sortrows(reshape(excludeIntervals, [], 2));
    % merge overlapping exclusions
    ms = e(1,1); me = e(1,2); merged = zeros(0,2);
    for k = 2:size(e,1)
        if e(k,1) <= me
            me = max(me, e(k,2));
        else
            merged(end+1,:) = [ms me]; %#ok<AGROW>
            ms = e(k,1); me = e(k,2);
        end
    end
    merged(end+1,:) = [ms me];

    edges = [0; merged(:,2)];
    stops = [merged(:,1); dur];
    gaps = [edges stops];
    gaps = gaps(gaps(:,2) - gaps(:,1) >= L, :);
end

gaps = gaps(gaps(:,2) - gaps(:,1) >= L, :);
if isempty(gaps), return; end

switch lower(prefer)
    case 'last'
        g = gaps(end,:);  w = [g(2)-L, g(2)];
    case 'longest'
        [~,i] = max(gaps(:,2)-gaps(:,1)); g = gaps(i,:); w = [g(1), g(1)+L];
    otherwise
        g = gaps(1,:);    w = [g(1), g(1)+L];
end
w(1) = max(w(1), 0); w(2) = min(w(2), dur);
if w(2) - w(1) < L, w = []; end
end
