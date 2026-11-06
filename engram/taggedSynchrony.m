function out = taggedSynchrony(spikeFile, varargin)
%TAGGEDSYNCHRONY  Are light-modulated units more synchronous with each other,
%   during periods with the light off, than a random same-size group of
%   units from the same session?
%
%   out = taggedSynchrony(spikeFile, ...)
%
%   DESIGN. Light-on epochs (padded) are removed by STITCHING (not masking)
%   before anything is computed -- see stitchOutIntervals.m for why masking
%   would bias the coincidence criterion at every pulse boundary. Units are
%   tagged from the light-ON data but ALL synchrony is then computed only on
%   the light-OFF, stitched timeline, so tag identity and the synchrony
%   estimate never share data.
%
%   REBOUND CAVEAT (ArchT). Units silenced by ArchT often fire a rebound
%   burst at light offset. Any rebound that falls outside the padded
%   exclusion stays in the "light-off" data, and because the suppressed
%   units rebound together after EVERY pulse, it would raise their C for
%   reasons that have nothing to do with light-off coordination. With a
%   small 'Pad' this is a live confound: check the post-offset PSTH of the
%   tagged units and rerun with 'Pad' covering the rebound.
%
%   The comparison is size-matched BY CONSTRUCTION: observed C for the
%   tagged group (size k) is compared against a null built from NDrawNull
%   random draws of k units (pool set by 'NullPool'), same session, same
%   light-off interval, same spikeSync call.
%
%   NAME/VALUE
%     'LatWin','MinProb','NPerm','Alpha','MaxLatency','MaxJitter' -- passed
%         to optoTagUnits (used only if 'TagMethod' is 'latency').
%     'TagMethod'  'rate' (default) or 'latency'. 'rate' uses
%                  optoTagUnitsRate (paired Poisson GLMM on binned counts,
%                  sm_getOptoMod); 'latency' uses optoTagUnits.
%     'ModulationSign' 'decrease' (default), 'increase', or 'either' --
%                  which significantly modulated units form the tested
%                  group: tag.isSuppressed, tag.isTagged (beta>0), or their
%                  union. 'decrease' and 'either' require TagMethod 'rate'
%                  (optoTagUnits has no suppression field).
%     'Dt'         bin duration (s) for 'rate' tagging (default 0.090).
%     'PreOffset','PostOffset' passed to optoTagUnitsRate (defaults 1, 0.005).
%     'Pad'        s excluded around each pulse beyond its own on/off,
%                  default 0.002. See REBOUND CAVEAT.
%     'RateRange'  [lo hi] Hz inclusion band on the LIGHT-OFF data, default
%                  [0.1 10].
%     'MinSpikesLightOff' minimum spikes on the stitched light-off timeline
%                  to keep a unit, default 20.
%     'NDrawNull'  random same-size draws for the null, default 500.
%     'NullPool'   'all' (default) or 'untagged'. 'all' draws the null from
%                  tagged+untagged combined, size k -- avoids degeneracy and
%                  does not require trusting the untagged pool as an
%                  unbiased reference. 'untagged' draws ONLY from non-tagged
%                  qualifying units and becomes DEGENERATE when nTagged
%                  approaches nUntaggedPool (draws nearly identical, Cnull
%                  variance collapses, |z| inflates). A warning is printed
%                  when nTagged/nUntaggedPool > 0.5 and NullPool is 'untagged'.
%     'MinTagged'  minimum tagged units to attempt the comparison, default 2.
%     'MinUntagged' minimum untagged pool size, default = tagged-group size.
%
%   OUTPUT struct (same fields whether or not the session succeeds):
%     .ok, .why, .file
%     .nTagged, .nUntaggedPool, .lightOffDur
%     .Cobs        observed multivariate C among tagged units, light-off
%     .Cnull       C for the random same-size draws
%     .z, .p       (Cobs - mean(Cnull)) / std(Cnull); permutation p
%     .tag         tagging output for every unit in the session
%     .taggedIdx, .untaggedIdx  indices into spikes.times
%     .rawPulses   K x 2 [on off] pulse times used for tagging AND stitching
%     .pulseFile   the pulseInfo file they came from
%     .params      every parsed option (incl. Pad, ModulationSign, NullPool),
%                  so downstream code can rebuild the exact light-off timeline
%
%   See also OPTOTAGUNITSRATE, OPTOTAGUNITS, STITCHOUTINTERVALS, SPIKESYNC.

p = inputParser;
p.addParameter('LatWin', 0.010);
p.addParameter('MinProb', 0.1);
p.addParameter('NPerm', 1000);
p.addParameter('Alpha', 0.01);
p.addParameter('MaxLatency', []);
p.addParameter('MaxJitter', []);
p.addParameter('TagMethod', 'rate', @(s) any(strcmpi(s,{'rate','latency'})));
p.addParameter('ModulationSign', 'decrease', @(s) any(strcmpi(s,{'increase','decrease','either'})));
p.addParameter('Dt', 0.090);
p.addParameter('PreOffset', 1);
p.addParameter('PostOffset', .005);
p.addParameter('Pad', 0.002);
p.addParameter('RateRange', [0.1 10]);
p.addParameter('MinSpikesLightOff', 20);
p.addParameter('NDrawNull', 500);
p.addParameter('NullPool', 'all', @(s) any(strcmpi(s,{'untagged','all'})));
p.addParameter('MinTagged', 2);
p.addParameter('MinUntagged', []);
p.parse(varargin{:});
o = p.Results;

out = taggedSynchronyBlank();
out.file = spikeFile;
out.params = o;

S = load(spikeFile, 'spikes');
if ~isfield(S,'spikes') || ~isfield(S.spikes,'times')
    out.why = 'no spikes.times'; return
end
tsAll = cellfun(@(x) sort(x(:)), S.spikes.times(:), 'UniformOutput', false);

% raw pulse on/off, for tagging (ms-scale latency logic) -- separate from
% optoPulseIntervals, which pads for second-scale EXCLUSION and is the
% wrong pad to reason about spike latency with.
d = fileparts(spikeFile);
c = dir(fullfile(d,'pulseInfo.mat'));
if isempty(c), out.why = 'no pulseInfo file'; return; end
out.pulseFile = fullfile(c(1).folder, c(1).name);
Praw = load(out.pulseFile);
fn = fieldnames(Praw); Pstruct = Praw.(fn{1});
rawPulses = zeros(0,2);
for e = 1:numel(Pstruct)
    if isfield(Pstruct(e),'time') && ~isempty(Pstruct(e).time)
        v = Pstruct(e).time;
        if size(v,2) >= 2, rawPulses = [rawPulses; v(:,1:2)]; %#ok<AGROW>
        else, rawPulses = [rawPulses; v(:) v(:)]; end %#ok<AGROW>
    end
end
if isempty(rawPulses), out.why = 'no usable raw pulse times'; return; end
rawPulses = sortrows(rawPulses);
out.rawPulses = rawPulses;

%% ---- tag from full (light-containing) data -------------------------------
if ~strcmpi(o.TagMethod,'rate') && ~strcmpi(o.ModulationSign,'increase')
    out.why = '''ModulationSign'' other than ''increase'' requires TagMethod ''rate''';
    return
end
if strcmpi(o.TagMethod, 'rate')
    tag = optoTagUnitsRate(tsAll, rawPulses, o.Dt, ...
        'PreOffset',o.PreOffset, 'PostOffset',o.PostOffset, 'Alpha',o.Alpha);
else
    tag = optoTagUnits(tsAll, rawPulses, 'LatWin',o.LatWin, 'MinProb',o.MinProb, ...
        'NPerm',o.NPerm, 'Alpha',o.Alpha, 'MaxLatency',o.MaxLatency, 'MaxJitter',o.MaxJitter);
end
out.tag = tag;
switch lower(o.ModulationSign)
    case 'decrease', isTagged = [tag.isSuppressed]';
    case 'increase', isTagged = [tag.isTagged]';
    case 'either',   isTagged = [tag.isTagged]' | [tag.isSuppressed]';
end

%% ---- stitch out light-on (+pad) for the synchrony computation -----------
lastSpk = cellfun(@(x) max([x; -Inf]), tsAll);
dur = max(lastSpk(isfinite(lastSpk)));
if isempty(dur) || ~isfinite(dur) || dur <= 0
    out.why = 'bad duration'; return
end

[tsOff, durOff] = stitchOutIntervals(tsAll, dur, rawPulses, o.Pad);
out.lightOffDur = durOff;
if durOff < 60
    out.why = sprintf('only %.0f s light-off after stitching', durOff);
    return
end

M = cellfun(@numel, tsOff);
FR = M / durOff;
qualifies = FR > o.RateRange(1) & FR < o.RateRange(2) & M >= o.MinSpikesLightOff;

taggedIdx = find(isTagged & qualifies);
untaggedIdx = find(~isTagged & qualifies);
out.taggedIdx = taggedIdx; out.untaggedIdx = untaggedIdx;
out.nTagged = numel(taggedIdx); out.nUntaggedPool = numel(untaggedIdx);

minUnt = o.MinUntagged; if isempty(minUnt), minUnt = out.nTagged; end
if out.nTagged < o.MinTagged
    out.why = sprintf('%d qualifying tagged units < MinTagged=%d', out.nTagged, o.MinTagged);
    return
end
if out.nUntaggedPool < minUnt
    out.why = sprintf('%d qualifying untagged units < required %d', out.nUntaggedPool, minUnt);
    return
end

%% ---- observed C among tagged units ---------------------------------------
try
    Cobs = spikeSync(tsOff(taggedIdx), 'Interval', [0 durOff]);
catch ME
    out.why = ['spikeSync (tagged) failed: ' ME.message]; return
end
out.Cobs = Cobs;

%% ---- null: random same-size draws --------------------------------------
k = out.nTagged;

if strcmpi(o.NullPool, 'all')
    drawPool = [taggedIdx; untaggedIdx];
else
    drawPool = untaggedIdx;
    ratio = k / max(out.nUntaggedPool, 1);
    if ratio > 0.5
        warning('taggedSynchrony:degenerateNull', ...
            ['nTagged/nUntaggedPool = %.2f: drawing %d of %d units leaves only ', ...
             '%d out per draw. The null draws will be nearly identical sets ', ...
             'and Cnull''s variance will collapse, inflating |z| regardless of ', ...
             'any true effect. Consider ''NullPool'',''all''.'], ...
            ratio, k, out.nUntaggedPool, out.nUntaggedPool - k);
    end
end
nPool = numel(drawPool);

Cnull = nan(o.NDrawNull,1);
for r = 1:o.NDrawNull
    sel = drawPool(randperm(nPool, k));
    try
        Cnull(r) = spikeSync(tsOff(sel), 'Interval', [0 durOff]);
    catch
        % skip a draw that happens to hit a degenerate pair; rare
    end
end
Cnull = Cnull(isfinite(Cnull));
out.Cnull = Cnull;
if numel(Cnull) < 10
    out.why = 'too few valid null draws'; return
end

out.z = (Cobs - mean(Cnull)) / max(std(Cnull), eps);
out.p = (1 + sum(Cnull >= Cobs)) / (1 + numel(Cnull));
out.ok = true;
out.why = '';
end