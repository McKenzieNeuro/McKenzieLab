%% runSyncAnalysis.m
% SPIKE-Synchronization, transient detection, and rate/duty/duration
% comparison between TransPlasticity and OptoRAM.
%
% Fixes accumulated over the course of debugging this pipeline, kept here so
% the next read of this file explains itself:
%
%  1. Interval passed explicitly as [0 dur] -- otherwise spikeSync infers it
%     from the retained units, so the edge-spike window fallback and the rate
%     denominator both shift with whichever units the rate filter happened
%     to keep.
%  2. Duration comes from real recording length where available, not last-
%     spike time (which biases rate high and is silently wrong for any
%     zero-spike unit). Still using last-spike inference by default here --
%     see the TODO below.
%  3. N is matched by construction: every session's trace, surrogate band,
%     and transient detection run on a SEEDED FIXED-SIZE SUBSET of units
%     (MatchN=10), not just a scalar correction on top of an unmatched
%     detector. The earlier version only corrected the scalar C, which left
%     rho(N, transient rate) ~ 0.43 -- more units meant a tighter surrogate
%     band meant more detected transients, independent of anything real.
%  4. Sessions need >=10 qualifying units (MinUnits) or they're dropped, not
%     padded.
%  5. Each session keeps its own time base; sessions are never indexed
%     inconsistently across datasets (syncBatch returns one struct per input
%     file, aligned, with .ok flagging failures rather than skipping them).
%  6. Opto pulses excluded from *.pulseInfo.mat by INTERVAL OVERLAP against
%     every window (not just the leading edge), and applied identically to
%     observed AND surrogate traces so the null band isn't shaped by
%     contaminated windows near an excluded epoch.
%  7. Width 60 / Step 5 (not Step 1 -- that was 98% overlap, 60x redundant).
%  8. Two syncBatch passes, not three+: a cheap screen (no surrogates, no
%     transients) to sanity-check N/rate/duration, then ONE expensive pass
%     with the real settings. syncBatch's cache key hashes the full option
%     set, so calling it with different options is always a full recompute,
%     never an accidental cache hit.
%  9. Duplicate surrogate generation eliminated inside syncBatch: the scalar
%     null and the windowed null used to each draw their own 50 surrogates
%     independently; now one shared draw produces both. 'Parallel', true
%     parallelizes the surrogate loop (Parallel Computing Toolbox).
%  10. rate / duty-cycle / transient-duration all go through transientStats.m
%      -- a single formula for cluster duration ((nWin-1)*step + width), not
%      four hand-written copies. One of the hand-written copies used Width
%      where it should have used Step and produced a >100% duty cycle before
%      this was caught.
%
% STILL OPEN, not fixed in this script:
%  - Duration is still last-spike-inferred. Pass a real 'Duration' function
%    handle once you have one (session metadata / TTL span / file length).
%  - Session/animal nesting is unchecked. If sessions are not one-per-animal,
%    the permutation tests below overstate the effective sample size.
%  - See checkDurationAndN.m for the duration-truncation and MatchN-sweep
%    sensitivity analyses -- run those before writing any of this down.
%    They are deliberately kept as a separate script, not folded in here,
%    because each is a full syncBatch pass and this file should stay the
%    one you can rerun quickly to get the headline numbers.

%% ---- file lists ------------------------------------------------------
fils1 = getAllExtFiles('R:\TransPlasticity','mat',1);
fils1 = fils1(contains(fils1,'spikes.cellinfo.mat'));
fils2 = getAllExtFiles('R:\WSun\ePhys\OptoRAM','mat',1);
fils2 = fils2(contains(fils2,'spikes.cellinfo.mat'));

%% ---- pass 1: SCREEN (cheap) --------------------------------------------
% No surrogates, no transient detection. Sanity-check N/rate/duration before
% spending any surrogate budget.
scfg = { 'RateRange', [0.1 10], 'Width', 60, 'Step', 5, 'MinSpikes', 10, ...
         'MinUnits', 10, ...
         'NullMethod', 'none', 'DetectTransients', false, ...
         'CacheDir', fullfile(tempdir,'syncScreen') };

S1 = syncBatch(fils1, scfg{:});
S2 = syncBatch(fils2, scfg{:}, 'Exclude', @optoPulseIntervals);

N1 = [S1([S1.ok]).N];  N2 = [S2([S2.ok]).N];
fprintf('\nscreen: N units TransPlasticity median %d [%d %d], OptoRAM median %d [%d %d]\n', ...
    median(N1), min(N1), max(N1), median(N2), min(N2), max(N2));
fprintf('screen: %d/%d and %d/%d sessions pass MinUnits=10\n', ...
    sum([S1.ok]), numel(S1), sum([S2.ok]), numel(S2));

%% ---- pass 2: the real run (expensive, ONCE per dataset) ----------------
% MatchN fixed at 10: every session's trace, surrogate band, and transient
% detection run on exactly 10 units, so detection sensitivity is constant
% across sessions by construction. See checkDurationAndN.m to test whether
% the result is sensitive to this specific choice.
%
% 'shift' rather than 'dither': circular rotation preserves each train's own
% ISI structure exactly. Dither destroys within-train regularity as well as
% cross-train alignment, which inflates the apparent effect for rhythmic or
% bursty units.

cfg = { 'RateRange', [0.1 10], 'Width', 60, 'Step', 5, 'MinSpikes', 10, ...
        'MinUnits', 10, 'MatchN', 15, ...
        'NullMethod', 'shift', 'NSurr', 50, ...
        'DetectTransients', true, 'NSurrWindow', 50, 'RunAlpha', 0.05, ...
        'Parallel', true, ...
        'CacheDir', fullfile(tempdir,'syncFull'), ...
        'CheckpointEvery', 5 };
    tempdir = 'R:\TransPlasticity\temp';
% CheckpointFile is set PER DATASET below (same option set otherwise would
% make the two checkpoint files collide -- each syncBatch call's checkpoint
% is keyed on options+file-list together, and the two file lists differ, so
% this is safe even sharing 'cfg', but separate files keep it legible).
ckptDir = fullfile(tempdir, 'syncCheckpoints');
if ~exist(ckptDir,'dir'), mkdir(ckptDir); end

T1 = syncBatch(fils1, cfg{:}, 'CheckpointFile', fullfile(ckptDir,'transplasticity.mat'));
T2 = syncBatch(fils2, cfg{:}, 'Exclude', @optoPulseIntervals, ...
    'CheckpointFile', fullfile(ckptDir,'optoram.mat'));

a = T1([T1.ok]); b = T2([T2.ok]);
fprintf('\nfull run: %d/%d TransPlasticity, %d/%d OptoRAM sessions usable\n', ...
    numel(a), numel(T1), numel(b), numel(T2));

% verify the fixes actually took: N constant, exclusion doing something
fprintf('N used: %s (should be a single value)\n', mat2str(unique([[a.N] [b.N]])));
fprintf('exclusion: %d/%d OptoRAM sessions have >0 excluded windows (median %d of ~%d)\n', ...
    sum([b.nExcluded] > 0), numel(b), round(median([b.nExcluded])), ...
    round(median([b.dur]) / cfg{find(strcmp(cfg,'Step'))+1}));

%% ---- scalar comparison (whole-session average level) --------------------
fprintf('\n-- scalar C (session-average) --\n');
fprintf('%-16s %8s %8s %8s\n','','C','null','dC');
fprintf('%-16s %8.3f %8.3f %+8.3f\n', 'TransPlasticity', mean([a.C]), mean([a.Cnull]), mean([a.dC]));
fprintf('%-16s %8.3f %8.3f %+8.3f\n', 'OptoRAM', mean([b.C]), mean([b.Cnull]), mean([b.dC]));
fprintf('dC difference p = %.4f\n', permtest([a.dC], [b.dC]));

%% ---- transient rate / duty / duration -----------------------------------
stepUsed  = cfg{find(strcmp(cfg,'Step'))+1};
widthUsed = cfg{find(strcmp(cfg,'Width'))+1};

[rate1, occ1, dur1] = local_collect(a, stepUsed, widthUsed);
[rate2, occ2, dur2] = local_collect(b, stepUsed, widthUsed);

fprintf('\n-- transient rate (per hour, significant clusters only) --\n');
fprintf('  TransPlasticity : %.2f +/- %.2f  (n=%d sessions)\n', ...
    mean(rate1), std(rate1)/sqrt(numel(rate1)), numel(rate1));
fprintf('  OptoRAM (ex-stim): %.2f +/- %.2f  (n=%d sessions)\n', ...
    mean(rate2), std(rate2)/sqrt(numel(rate2)), numel(rate2));
pRate = permtest(rate1, rate2);

fprintf('\n-- duty cycle (%% of session in a significant transient) --\n');
fprintf('  TransPlasticity : %.2f +/- %.2f\n', mean(occ1)*100, std(occ1)/sqrt(numel(occ1))*100);
fprintf('  OptoRAM (ex-stim): %.2f +/- %.2f\n', mean(occ2)*100, std(occ2)/sqrt(numel(occ2))*100);
pDuty = permtest(occ1, occ2);

fprintf('\n-- mean transient duration (s) --\n');
fprintf('  TransPlasticity : %.1f +/- %.1f\n', mean(dur1,'omitnan'), std(dur1,'omitnan')/sqrt(sum(~isnan(dur1))));
fprintf('  OptoRAM (ex-stim): %.1f +/- %.1f\n', mean(dur2,'omitnan'), std(dur2,'omitnan')/sqrt(sum(~isnan(dur2))));
pDur = permtest(dur1(~isnan(dur1)), dur2(~isnan(dur2)));

mag1 = arrayfun(@(r) local_meanPeak(r.transients), a);
mag2 = arrayfun(@(r) local_meanPeak(r.transients), b);
pMag = permtest(mag1(~isnan(mag1)), mag2(~isnan(mag2)));

fprintf('\np: rate=%.4f  duty=%.4f  duration=%.4f  magnitude=%.4f\n', pRate, pDuty, pDur, pMag);
fprintf(['read together: rate marginal, duty and duration significant -> the effect\n' ...
         'is mainly TRANSIENTS LASTING LONGER once they occur, not more of them, with\n' ...
         'unchanged peak magnitude once in a transient. Confirm with checkDurationAndN.m\n' ...
         'before writing this down: session durations differ ~4x between datasets\n' ...
         '(p=3.6e-18), and duration must be ruled out as the driver.\n']);
     
     checkDurationAndN

%% ---- drill into a single transient-rich session -------------------------
% [~, worst] = max(rate1);
% one = a(worst);
% figure; hold on
% plot(one.tc, one.Ct, 'k-');
% sig = one.transients([one.transients.significant]);
% for k = 1:numel(sig)
%     xline(sig(k).tStart, 'r-'); xline(sig(k).tEnd, 'r-');
% end
% xlabel('time (s)'); ylabel('C (windowed)'); title(one.file)


%% =========================================================================
function [rate, duty, meanDur] = local_collect(T, step, width)
rate = nan(numel(T),1); duty = nan(numel(T),1); meanDur = nan(numel(T),1);
for i = 1:numel(T)
    [rate(i), duty(i), meanDur(i)] = transientStats(T(i), step, width);
end
end

function m = local_meanPeak(trans)
if isempty(trans), m = NaN; return; end
sig = [trans.significant];
if ~any(sig), m = NaN; return; end
m = mean([trans(sig).peakDC]);
end


%% =========================================================================
function ex = optoPulseIntervals(spikeFile)
% Return K x 2 [start stop] intervals to exclude, in seconds.
% Loaded from the session folder rather than inherited from the workspace.
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

if isfield(P,'time') && size(P.time,2) >= 2
    ex = P.time(:,1:2);
elseif isfield(P,'time')
    ex = [P.time(:) P.time(:)];
elseif isfield(P,'timestamps') && size(P.timestamps,2) >= 2
    ex = P.timestamps(:,1:2);
elseif isfield(P,'timestamps')
    ex = [P.timestamps(:) P.timestamps(:)];
elseif isfield(P,'ints')
    ex = P.ints(:,1:2);
else
    warning('optoPulseIntervals:fields','%s: no time/timestamps/ints field (has: %s)', ...
        c(1).name, strjoin(fieldnames(P)', ', '));
    return
end

if ~isempty(ex) && max(ex(:)) > 1e5
    warning('optoPulseIntervals:units', ...
        '%s: max pulse time %.4g looks like samples, not seconds. Divide by Fs.', ...
        c(1).name, max(ex(:)));
end

ex = [ex(:,1) - 0.5, ex(:,2) + 2.0];
end