%% checkPositionConfound.m
% Does epoch position-in-session explain the null in runCrossDatasetNaive,
% or is it genuinely about drug history? Two checks, run after
% runCrossDatasetNaive.m (uses its saved A, B, posA, posB).
%
% CHECK 1: within TransPlasticity, does the DRUG-NAIVE group's epoch
% position correlate with duration? If position alone drives duration,
% naive epochs recorded LATE in a session should look like the "typical"
% (late, post-hoc) sessions from the original whole-session comparison even
% before any exposure -- which would mean position, not drug history, was
% doing the work in the original effect.
%
% CHECK 2: reproduce the original whole-session TransPlasticity-vs-OptoRAM
% duration comparison, but using ONLY OptoRAM sessions restricted to their
% first LFIX seconds (matching this analysis's OptoRAM epoch) vs ONLY
% TransPlasticity sessions restricted to their first LFIX seconds (an
% EARLY, not pre-infusion, TransPlasticity epoch, so it removes the position
% mismatch while still being drug-naive-or-not depending on the session's
% condition). If the ASV-exposed group still differs from OptoRAM under
% matched EARLY position, drug history survives as the explanation. If the
% original 2x TransPlasticity-vs-OptoRAM effect vanishes once position is
% matched regardless of exposure, position was doing more of the work than
% previously credited.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
WIDTH = 60; STEP = 5; MATCHN = 15; LFIX = 3600;
load(fullfile(CACHE_ROOT,'crossNaiveResults.mat'), 'A','B','aAnim','bAnim','posA','posB');

dA = nan(numel(A),1);
for i = 1:numel(A), [~,~,dA(i)] = transientStats(A(i), STEP, WIDTH); end

fprintf('CHECK 1: does epoch position predict duration WITHIN drug-naive TransPlasticity epochs?\n');
g = isfinite(dA);
[r,p] = corr(posA(g), log(dA(g)), 'type','Spearman');
fprintf('  rho(position, log duration) = %+.3f, p = %.4f (n=%d)\n', r, p, sum(g));
fprintf('  A strong positive rho here means later epochs run longer regardless\n');
fprintf('  of exposure, i.e. position alone could produce the original effect.\n');

%% ---- CHECK 2: match position instead of exposure -------------------------
fprintf('\nCHECK 2: TransPlasticity vs OptoRAM, BOTH restricted to first %.0f min,\n', LFIX/60);
fprintf('regardless of ASV exposure (removes position mismatch, keeps exposure history mixed)\n');

fils1 = getAllExtFiles('R:\TransPlasticity','mat',1);
fils1 = fils1(contains(fils1,'spikes.cellinfo.mat'));
cnd1 = cellfun(@(f) conditionFromPath(f,'TransPlasticity'), fils1, 'UniformOutput',false);
keep1 = ismember(cnd1, {'Baseline','ASV','ASV24','SAL','SAL24'});
fils1 = fils1(keep1);

scfg = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',10,'NullMethod','none','DetectTransients',false, ...
         'CacheDir', fullfile(CACHE_ROOT,'screenPosCheck') };
S1 = syncBatch(fils1, scfg{:});
ok1 = [S1.ok]' & [S1.dur]' >= LFIX;
uf1 = fils1(ok1);
mapEarly = containers.Map(uf1, repmat({[0 LFIX]}, sum(ok1), 1));

cfgB = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',10,'MatchN',MATCHN,'NullMethod','shift','NSurr',50, ...
         'DetectTransients',true,'NSurrWindow',50,'RunAlpha',0.05, ...
         'Parallel', ~isempty(ver('parallel')), 'CheckpointEvery',5 };
T1e = syncBatch(uf1, cfgB{:}, 'Window', @(f) mapEarly(f), ...
    'CacheDir', fullfile(CACHE_ROOT,'earlyTP'), ...
    'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_earlyTP.mat'));
Ae = T1e([T1e.ok]);
aAnimE = arrayfun(@(r) animalFromPath(r.file,'TransPlasticity'), Ae, 'UniformOutput',false);

dAe = nan(numel(Ae),1);
for i = 1:numel(Ae), [~,~,dAe(i)] = transientStats(Ae(i), STEP, WIDTH); end

fprintf('  TransPlasticity (first %.0f min, any exposure): %d epochs, %d animals\n', ...
    LFIX/60, numel(Ae), numel(unique(aAnimE)));

dBvec = nan(numel(B),1);
for i = 1:numel(B), [~,~,dBvec(i)] = transientStats(B(i), STEP, WIDTH); end

Rn = nestingCheck(Ae, B, log(dAe), log(dBvec), ...
    'TransPlasticity','OptoRAM','log duration, position-matched');
fprintf(['\nCompare this p to the ORIGINAL whole-session result (p_animal=0.0073).\n' ...
         'If position-matched-but-exposure-mixed ALSO comes back null, the\n' ...
         'original effect depended on the late timing of typical TransPlasticity\n' ...
         'sessions as much as on drug history, and the two are not yet separable\n' ...
         'with this design.\n']);
