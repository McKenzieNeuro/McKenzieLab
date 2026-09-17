function ex = optoPulseIntervals(spikeFile)
%OPTOPULSEINTERVALS  Stimulation intervals to exclude, for syncBatch's 'Exclude'.
%
%   ex = optoPulseIntervals(spikeFile)  ->  K x 2 [start stop] in SECONDS
%
%   Reads *.pulseInfo.mat from the same folder as the given spikes file.
%
%   This lives in its own file, not as a local function at the bottom of
%   runSyncAnalysis.m, because local functions in a SCRIPT are only in scope
%   while that script is running. Passing @optoPulseIntervals to syncBatch
%   from the command line -- e.g. for a MatchN sweep or any one-off
%   diagnostic -- would throw "Undefined function 'optoPulseIntervals'" for
%   every session, which syncBatch's per-session try/catch then swallowed as
%   an ordinary session failure. That produced 66 silently "failed" OptoRAM
%   sessions that looked like a data problem and were not.

ex = zeros(0,2);
d = fileparts(spikeFile);
c = dir(fullfile(d, '*.pulseInfo.mat'));
if isempty(c), c = dir(fullfile(d, '*pulseInfo*.mat')); end
if isempty(c)
    warning('optoPulseIntervals:missing','no pulseInfo file in %s', d);
    return
end

S = load(fullfile(c(1).folder, c(1).name));
f = fieldnames(S);
P = S.(f{1});
if ~isstruct(P)
    warning('optoPulseIntervals:shape','%s: top-level variable is not a struct', c(1).name);
    return
end
if numel(P) > 1
    % informational: pooling across elements is intended, but if the elements
    % are different stimulation protocols you may want only some of them.
    fprintf('optoPulseIntervals: %s has %d struct elements; pooling all.\n', ...
        c(1).name, numel(P));
end

% P may be a STRUCT ARRAY (e.g. one element per stimulation channel or per
% pulse train). Indexing P.time on a struct array yields a comma-separated
% list, which is why size(P.time,2) errored. Loop over elements and pool.
ex = zeros(0,2);
for e = 1:numel(P)
    Pe = P(e);
    if isfield(Pe,'time') && ~isempty(Pe.time)
        v = Pe.time;
        if size(v,2) >= 2, ex = [ex; v(:,1:2)]; else, ex = [ex; v(:) v(:)]; end %#ok<AGROW>
    elseif isfield(Pe,'timestamps') && ~isempty(Pe.timestamps)
        v = Pe.timestamps;
        if size(v,2) >= 2, ex = [ex; v(:,1:2)]; else, ex = [ex; v(:) v(:)]; end %#ok<AGROW>
    elseif isfield(Pe,'ints') && ~isempty(Pe.ints)
        v = Pe.ints;
        if size(v,2) >= 2, ex = [ex; v(:,1:2)]; end %#ok<AGROW>
    end
end
if isempty(ex)
    warning('optoPulseIntervals:fields', ...
        '%s: no usable time/timestamps/ints across %d struct element(s) (fields: %s)', ...
        c(1).name, numel(P), strjoin(fieldnames(P(1))', ', '));
    return
end
ex = sortrows(ex);

% Sanity check: pulse times must be SECONDS and inside the recording. If they
% are sample indices (values in the millions) every window overlap test fails
% silently and you get 0 excluded windows.
if ~isempty(ex) && max(ex(:)) > 1e5
    warning('optoPulseIntervals:units', ...
        '%s: max pulse time %.4g looks like samples, not seconds. Divide by Fs.', ...
        c(1).name, max(ex(:)));
end

% pad by the stimulation transient you want to exclude, not by zero
ex = [ex(:,1) - 0.5, ex(:,2) + 2.0];
end