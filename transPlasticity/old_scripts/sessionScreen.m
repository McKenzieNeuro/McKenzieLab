function [dur, ok, nKeep] = sessionScreen(fils, rateRange, minUnits, cacheDir)
%SESSIONSCREEN  Recording duration and usability per session, without
%   computing any synchrony.
%
%   [dur, ok, nKeep] = sessionScreen(fils, rateRange, minUnits, cacheDir)
%
%   Replaces the "screen" syncBatch call (NullMethod 'none', no MatchN) the
%   analysis scripts used only for .dur and .ok. That call still ran a full
%   spikeSync over EVERY unit of the whole recording, which runs out of memory
%   on sessions with hundreds of units. The rule here is the same as
%   syncBatch's:
%     dur   = last spike of any unit (s)
%     keep  = rate in (rateRange(1), rateRange(2)) over [0 dur] and >= 1 spike
%     ok    = sum(keep) >= max(2, minUnits)
%   Results are cached per file (keyed on path, file date and the rule), so
%   reruns only load spikes for new or changed sessions. cacheDir = [] disables.

if nargin < 4, cacheDir = []; end
if isempty(rateRange), rateRange = [0 Inf]; end   % [] = no rate exclusion
if ~isempty(cacheDir) && ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
n = numel(fils); dur = nan(n,1); ok = false(n,1); nKeep = zeros(n,1);
for i = 1:n
    f = fils{i};
    d = dir(f);
    cf = '';
    if ~isempty(cacheDir) && ~isempty(d)
        key = sprintf('%s|%.6f|%g|%g|%d', f, d.datenum, rateRange(1), rateRange(2), minUnits);
        h = typecast(int32(java.lang.String(key).hashCode()), 'uint32');
        cf = fullfile(cacheDir, sprintf('screen_%010u.mat', h));
        if exist(cf, 'file')
            C = load(cf, 'v', 'key');
            if strcmp(C.key, key), dur(i) = C.v(1); ok(i) = C.v(2) > 0; nKeep(i) = C.v(3); continue; end
        end
    end
    try
        S = load(f, 'spikes');
        ts = S.spikes.times(:);
        M = cellfun(@(x) sum(isfinite(x)), ts);
        last = cellfun(@(x) max([x(isfinite(x)); -Inf]), ts);
        dur(i) = max(last(isfinite(last)));
        FR = M / dur(i);
        keep = FR > rateRange(1) & FR < rateRange(2) & M > 0;
        nKeep(i) = sum(keep);
        ok(i) = isfinite(dur(i)) && dur(i) > 0 && nKeep(i) >= max(2, minUnits);
    catch ME
        warning('sessionScreen:fail', '%s: %s', f, ME.message);
    end
    if ~isempty(cf)
        v = [dur(i), double(ok(i)), nKeep(i)]; %#ok<NASGU>
        save(cf, 'v', 'key');
    end
end
end
