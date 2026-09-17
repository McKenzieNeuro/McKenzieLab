%% runASVAnalysis.m
% Master orchestrator for the ASV synchrony analysis, TransPlasticity only.
%
% REPLACES runAllTests.m (renamed to _deprecated_runAllTests_crossDataset.m,
% kept for reference, not meant to be run). That script was built around
% comparing TransPlasticity against OptoRAM; the cross-dataset comparison was
% retired after the position-in-session confound (checkPositionConfound.m)
% showed the original dataset-level duration difference did not survive
% matching epoch position, and every question this project now asks is
% within TransPlasticity, within animal.
%
% All five steps below share ONE CACHE_ROOT. Each step is also an
% independently runnable script (unchanged) -- this file just sequences them
% and confirms each one's expected output exists before moving on, so a
% failure in step 2 doesn't waste an hour redoing step 1's surrogate compute,
% and doesn't silently let step 3 run against a stale step-2 result.
%
% Re-run this whenever sessions are added. syncBatch's cache is
% content+options keyed, so unchanged sessions are not recomputed; only new
% or changed ones cost time.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';   % LOCAL disk -- not R:\, not tempdir
if ~exist(CACHE_ROOT,'dir'), mkdir(CACHE_ROOT); end

fprintf('\n=========================================================\n');
fprintf(' STEP 1/5: full syncBatch run (TransPlasticity)\n');
fprintf('=========================================================\n');

fils1 = getAllExtFiles('R:\TransPlasticity','mat',1);
fils1 = fils1(contains(fils1,'spikes.cellinfo.mat'));
fprintf('%d candidate session files\n', numel(fils1));

cfg = { 'RateRange', [0.1 10], 'Width', 60, 'Step', 5, 'MinSpikes', 10, ...
        'MinUnits', 10, 'MatchN', 15, ...
        'NullMethod', 'shift', 'NSurr', 50, ...
        'DetectTransients', true, 'NSurrWindow', 50, 'RunAlpha', 0.05, ...
        'Parallel', ~isempty(ver('parallel')), ...
        'CacheDir', fullfile(CACHE_ROOT,'full'), 'CheckpointEvery', 5 };

T1 = syncBatch(fils1, cfg{:}, 'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_tp.mat'));
a = T1([T1.ok]);
fprintf('%d/%d sessions usable, N used = %s\n', numel(a), numel(T1), mat2str(unique([a.N])));
assert(numel(unique([a.N])) == 1 && unique([a.N]) == 15, ...
    'syncBatch:runASVAnalysis','N did not come out matched at 15 -- stop and check before continuing.');

% quick design check before spending time on the rest: how are sessions
% distributed across animals now that more have been added?
anims = arrayfun(@(r) animalFromPath(r.file,'TransPlasticity'), a, 'UniformOutput', false);
[ua, ~, gi] = unique(anims);
counts = accumarray(gi, 1);
fprintf('\nsessions per animal (%d animals):\n', numel(ua));
[cs, ord] = sort(counts, 'descend');
for i = 1:numel(ua), fprintf('  %-10s %d\n', ua{ord(i)}, cs(i)); end
if max(counts) / sum(counts) > 0.25
    fprintf('\n>> %s alone is >25%% of usable sessions -- run a leave-one-out on\n', ua{ord(1)});
    fprintf('   any result that depends on it before trusting the group estimate.\n');
end

fprintf('\n=========================================================\n');
fprintf(' STEP 2/5: acute within-session infusion analysis\n');
fprintf('=========================================================\n');
runInfusionAnalysis
assert(logical(exist(fullfile(CACHE_ROOT,'infusionResults.mat'),'file')), ...
    'runASVAnalysis:step2','infusionResults.mat was not produced -- check the output above.');

fprintf('\n=========================================================\n');
fprintf(' STEP 3/5: persistence (drug-free epoch) analysis\n');
fprintf('=========================================================\n');
runPersistenceAnalysis
assert(logical(exist(fullfile(CACHE_ROOT,'persistenceResults.mat'),'file')), ...
    'runASVAnalysis:step3','persistenceResults.mat was not produced -- check the output above.');

fprintf('\n=========================================================\n');
fprintf(' STEP 4/5: step-vs-trend shape test\n');
fprintf('=========================================================\n');
fprintf('NOTE: add ''overall rate'' to stepAnalysis''s MET list if not already\n');
fprintf('present -- this is the outstanding check on whether the firing-rate\n');
fprintf('decline in post-exposure epochs is a step (like duty cycle) or a\n');
fprintf('smooth trend (consistent with recording drift instead of drug).\n');
stepAnalysis

fprintf('\n=========================================================\n');
fprintf(' STEP 5/5: firing rate within synchrony bouts\n');
fprintf('=========================================================\n');
runBoutRateAnalysis

fprintf('\n=========================================================\n');
fprintf(' DONE. Cross-reference against methods_results_draft before\n');
fprintf(' updating any number in the manuscript.\n');
fprintf('=========================================================\n');
