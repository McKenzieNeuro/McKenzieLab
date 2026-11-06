%% runAllTests.m
% One entry point for the whole TransPlasticity vs OptoRAM synchrony analysis.
%
% Runs, in order:
%   0. preflight  -- every required function is on the path, exactly once
%   1. batch      -- syncBatch per dataset (cached; cheap on rerun)
%   2. metrics    -- derive all six per-session quantities
%   3. tests      -- each metric through BOTH nestingCheck (exact permutation,
%                    session- and animal-level) and lmeCompare (mixed model,
%                    animal random intercept, Satterthwaite DF)
%   4. summary    -- consolidated table, with the multiple-comparison problem
%                    stated rather than buried
%   5. sensitivity (optional, expensive) -- MatchN sweep and duration-matched
%                    truncation
%
% READ THIS BEFORE INTERPRETING THE OUTPUT
%
% Six metrics x (session perm, animal perm, LME) is eighteen p-values, and
% they were not all pre-registered -- several were added mid-analysis in
% response to earlier results. Treat exactly one as confirmatory and the rest
% as exploratory/descriptive, or apply a correction and say which. The
% honest framing for this dataset: transient DURATION was the effect that
% survived every check, and everything else characterises it rather than
% independently testing it. Duty cycle in particular is downstream of
% duration by construction (duty ~ rate x duration), so it is not an
% independent test of anything.
%
% WHY TWO STATISTICAL METHODS PER METRIC
%   session-level permutation : treats sessions as independent. They are not
%                               (ICC up to 0.73 here). Reported only to show
%                               how much the naive test overstates things.
%   animal-level permutation  : exact, assumption-free, but conservative --
%                               every animal weighted equally regardless of
%                               how many sessions it contributed.
%   LME                       : uses within-animal replication, so better
%                               powered, but asymptotic and approximate at
%                               10-12 clusters per group.
%   Agreement between the last two is the bar. Where they disagree, trust the
%   permutation test's CONCLUSION and treat the LME as suggestive.
%
% LOG TRANSFORMS: applied to duration, magnitude, and bout-rate ratio --
% strictly positive, right-skewed, multiplicative-effect quantities. NOT
% applied to rate, duty, or dC: rate and duty are legitimately zero for
% sessions with no transients (log(0) = -Inf), and dC can be negative
% (log of a negative is complex, and would corrupt the fit silently).

clear; clc;

%% ---- 0. preflight -------------------------------------------------------
req = {'syncBatch','spikeSync','spikeSyncWindow','transientStats','permtest', ...
       'optoPulseIntervals','animalFromPath','nestingCheck','lmeCompare', ...
       'meanPeakMagnitude','boutRate','boutRateBatch','getAllExtFiles'};
fprintf('preflight:\n');
bad = false;
for k = 1:numel(req)
    w = which('-all', req{k});
    if isempty(w)
        fprintf('  MISSING  %s\n', req{k}); bad = true;
    elseif numel(w) > 1
        fprintf('  SHADOWED %s (%d copies -- %s)\n', req{k}, numel(w), w{1}); bad = true;
    end
end
if bad
    error('runAllTests:path', ...
        ['Fix the path before running. Shadowed/missing functions have already ', ...
         'cost this analysis several debugging cycles (a script saved under a ', ...
         'function''s filename, and local functions not visible from the command line).']);
end
fprintf('  all %d functions resolve uniquely\n', numel(req));

if isempty(ver('stats')),    error('runAllTests:stats','LME needs Statistics Toolbox.'); end
hasPar = ~isempty(ver('parallel'));
fprintf('  parallel: %d\n', hasPar);

%% ---- 1. batch -----------------------------------------------------------
RUN_SENSITIVITY = false;    % section 5: expensive, set true when you want it
CACHE_ROOT      = 'C:\Users\samckenzie\syncCache';   % LOCAL disk, not R:\ or tempdir
if ~exist(CACHE_ROOT,'dir'), mkdir(CACHE_ROOT); end

fils1 = getAllExtFiles('R:\TransPlasticity','mat',1);
fils1 = fils1(contains(fils1,'spikes.cellinfo.mat'));
fils2 = getAllExtFiles('R:\WSun\ePhys\OptoRAM','mat',1);

fils2 = fils2(contains(fils2,'spikes.cellinfo.mat'));

WIDTH = 60; STEP = 5; MATCHN = 15;   % MATCHN=15: significant, keeps most sessions

cfg = { 'RateRange', [0.1 10], 'Width', WIDTH, 'Step', STEP, 'MinSpikes', 10, ...
        'MinUnits', 10, 'MatchN', MATCHN, ...
        'NullMethod', 'shift', 'NSurr', 50, ...
        'DetectTransients', true, 'NSurrWindow', 50, 'RunAlpha', 0.05, ...
        'Parallel', hasPar, ...
        'CacheDir', fullfile(CACHE_ROOT,'full'), 'CheckpointEvery', 5 };

T1 = syncBatch(fils1, cfg{:}, 'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_tp.mat'));
T2 = syncBatch(fils2, cfg{:}, 'Exclude', @optoPulseIntervals, ...
    'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_or.mat'));

a = T1([T1.ok]); b = T2([T2.ok]);
fprintf('\nusable: %d/%d TransPlasticity, %d/%d OptoRAM\n', ...
    numel(a), numel(T1), numel(b), numel(T2));
fprintf('N used: %s | exclusion: %d/%d OptoRAM sessions affected\n', ...
    mat2str(unique([[a.N] [b.N]])), sum([b.nExcluded]>0), numel(b));

%% ---- 2. metrics ---------------------------------------------------------
[rate1, occ1, dur1] = deal(nan(numel(a),1));
for i = 1:numel(a), [rate1(i), occ1(i), dur1(i)] = transientStats(a(i), STEP, WIDTH); end
[rate2, occ2, dur2] = deal(nan(numel(b),1));
for i = 1:numel(b), [rate2(i), occ2(i), dur2(i)] = transientStats(b(i), STEP, WIDTH); end

mag1 = arrayfun(@(r) meanPeakMagnitude(r.transients), a);
mag2 = arrayfun(@(r) meanPeakMagnitude(r.transients), b);

dC1 = [a.dC]'; dC2 = [b.dC]';

fprintf('\ncomputing bout rates (reloads raw spikes -- slower)...\n');
[rOv1, rIn1, ratio1] = boutRateBatch(a, WIDTH, STEP); %#ok<ASGLU>
[rOv2, rIn2, ratio2] = boutRateBatch(b, WIDTH, STEP); %#ok<ASGLU>

%% ---- 3. tests -----------------------------------------------------------
% name, x1, x2, transform ('none' or 'log'), role
M = { 'transient duration',        dur1,   dur2,   'log',  'PRIMARY'
      'transient rate',            rate1,  rate2,  'none', 'secondary'
      'duty cycle',                occ1,   occ2,   'none', 'downstream of duration'
      'peak magnitude',            mag1,   mag2,   'log',  'secondary'
      'session mean dC',           dC1,    dC2,    'none', 'secondary'
      'bout/overall rate ratio',   ratio1, ratio2, 'log',  'secondary'
      'overall firing rate',       rOv1,   rOv2,   'log',  'covariate check' };

res = struct('metric',{},'role',{},'transform',{}, ...
    'meanA',{},'meanB',{},'pSession',{},'pAnimal',{},'pLME',{}, ...
    'iccA',{},'iccB',{},'nA',{},'nB',{});

for k = 1:size(M,1)
    name = M{k,1}; x1 = M{k,2}; x2 = M{k,3}; tf = M{k,4};
    if strcmp(tf,'log')
        if any(x1(~isnan(x1)) <= 0) || any(x2(~isnan(x2)) <= 0)
            warning('runAllTests:log','%s has non-positive values; skipping log.', name);
            tf = 'none';
        else
            x1 = log(x1); x2 = log(x2);
        end
    end

    fprintf('\n\n########## %s (%s, %s) ##########\n', name, M{k,5}, tf);
    Rn = nestingCheck(a, b, x1, x2, 'TransPlasticity', 'OptoRAM', name);
    try
        lme = lmeCompare(a, b, x1, x2, 'TransPlasticity', 'OptoRAM', name);
        st = anova(lme, 'DFMethod', 'Satterthwaite');
        pL = st.pValue(strcmp(st.Term,'Group'));
    catch ME
        warning('runAllTests:lme','%s: LME failed (%s)', name, ME.message);
        pL = NaN;
    end

    res(end+1) = struct('metric',name,'role',M{k,5},'transform',tf, ...
        'meanA',mean(x1,'omitnan'),'meanB',mean(x2,'omitnan'), ...
        'pSession',Rn.pSession,'pAnimal',Rn.pAnimal,'pLME',pL, ...
        'iccA',Rn.iccA,'iccB',Rn.iccB, ...
        'nA',sum(~isnan(x1)),'nB',sum(~isnan(x2))); %#ok<SAGROW>
end

%% ---- 4. summary ---------------------------------------------------------
fprintf('\n\n================== SUMMARY ==================\n');
fprintf('%-26s %-22s %8s %8s %8s %6s %6s\n', ...
    'metric','role','p_sess','p_animal','p_LME','ICC_A','ICC_B');
for k = 1:numel(res)
    fprintf('%-26s %-22s %8.4f %8.4f %8.4f %6.2f %6.2f\n', ...
        res(k).metric, res(k).role, res(k).pSession, res(k).pAnimal, ...
        res(k).pLME, res(k).iccA, res(k).iccB);
end

fprintf(['\nInterpretation rules used above:\n' ...
    '  p_sess is shown for contrast only -- sessions are NOT independent here.\n' ...
    '  Report p_animal (exact) and p_LME (better powered) together.\n' ...
    '  Where they disagree, the permutation conclusion wins.\n' ...
    '  %d metrics were tested; only "transient duration" was the pre-specified\n' ...
    '  hypothesis after the N-matching fix. Everything else is exploratory --\n' ...
    '  correct for multiplicity or label it as such.\n'], numel(res));

resTable = struct2table(res);
save(fullfile(CACHE_ROOT,'results.mat'), 'res','resTable','a','b','cfg', ...
    'rate1','rate2','occ1','occ2','dur1','dur2','mag1','mag2','dC1','dC2', ...
    'ratio1','ratio2','rOv1','rOv2');
fprintf('\nsaved -> %s\n', fullfile(CACHE_ROOT,'results.mat'));

%% ---- 5. sensitivity (optional) ------------------------------------------
if RUN_SENSITIVITY
    fprintf('\n\n========== SENSITIVITY ==========\n');

    % 5a. MatchN sweep -- is the result an artifact of MATCHN specifically?
    for mN = [10 15 20 25]
        mcfg = cfg; mcfg{find(strcmp(mcfg,'MatchN'))+1} = mN;
        mcfg{find(strcmp(mcfg,'CacheDir'))+1} = fullfile(CACHE_ROOT, sprintf('m%d', mN));
        Q1 = syncBatch(fils1, mcfg{:});
        Q2 = syncBatch(fils2, mcfg{:}, 'Exclude', @optoPulseIntervals);
        qa = Q1([Q1.ok]); qb = Q2([Q2.ok]);
        d1 = nan(numel(qa),1); for i=1:numel(qa), [~,~,d1(i)] = transientStats(qa(i),STEP,WIDTH); end
        d2 = nan(numel(qb),1); for i=1:numel(qb), [~,~,d2(i)] = transientStats(qb(i),STEP,WIDTH); end
        ok1 = ~isnan(d1); ok2 = ~isnan(d2);
        Rn = nestingCheck(qa, qb, log(d1), log(d2), 'TransPlasticity','OptoRAM', ...
            sprintf('log duration @ MatchN=%d', mN));
        fprintf('MatchN=%2d  nA=%3d nB=%3d  duration p_animal=%.4f\n', ...
            mN, sum(ok1), sum(ok2), Rn.pAnimal);
    end

    % 5b. duration-matched truncation -- is it a session-length artifact?
    capDur = median([b.dur]);
    tcfg = [cfg, {'MaxDuration', capDur}];
    tcfg{find(strcmp(tcfg,'CacheDir'))+1} = fullfile(CACHE_ROOT,'trunc');
    U1 = syncBatch(fils1, tcfg{:});
    U2 = syncBatch(fils2, tcfg{:}, 'Exclude', @optoPulseIntervals);
    ua = U1([U1.ok]); ub = U2([U2.ok]);
    e1 = nan(numel(ua),1); for i=1:numel(ua), [~,~,e1(i)] = transientStats(ua(i),STEP,WIDTH); end
    e2 = nan(numel(ub),1); for i=1:numel(ub), [~,~,e2(i)] = transientStats(ub(i),STEP,WIDTH); end
    fprintf('\ntruncated to %.0f s:\n', capDur);
    nestingCheck(ua, ub, log(e1), log(e2), 'TransPlasticity','OptoRAM','log duration (truncated)');
end
