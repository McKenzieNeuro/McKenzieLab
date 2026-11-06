%% responderPersistenceSplit.m
% Restrict to animals with an acute dosing-day effect, then test persistence
% -- using SPLIT baselines so selection and test do not share noise.
%
% THE PROBLEM THIS SOLVES
% Selecting on acute change (post-infusion minus pre-infusion baseline) and
% then measuring persistence against that same baseline means regression to
% the mean produces a positive persistence step with nothing having
% persisted. Dropping animals to avoid it (the earlier fix) cost too many
% subjects.
%
% THE SPLIT
% The pre-infusion baseline epoch is partitioned into alternating time
% blocks: A blocks feed the acute contrast, B blocks serve as the drug-naive
% reference for the persistence step. Disjoint minutes, so the sampling
% noise that drives regression to the mean is not shared. Every animal is
% retained.
%
% WHAT THIS DOES NOT FIX
% A and B come from the same epoch, so anything slowly varying across that
% hour -- state, drift, a settling electrode -- is still common to both.
% The split removes shared SAMPLING noise, not shared slow structure. It
% weakens the circularity, it does not eliminate it.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
WIDTH = 60; STEP = 5; BLOCK = 600;

% POST INCLUSION. Broad by default: every drug-free epoch at least
% MIN_DAYS_POST days after a dose, INCLUDING baselines recorded on later
% dosing days. No drug is on board at >=24 h, so these test persistence
% whichever day they fall on, and keeping them is what preserves enough
% animals to look at a responder subgroup at all.
%
% The trade is known and should be reported both ways: on the full group the
% duty step was +0.139 (6/8 animals) with dosing-day baselines included and
% +0.071 (4/7) with them excluded. Flip EXCLUDE_ASVDAY_FROM_POST to see the
% conservative version of everything below.
MIN_DAYS_POST            = 1;
EXCLUDE_ASVDAY_FROM_POST = false;

A = load(fullfile(CACHE_ROOT,'infusionResults.mat'));      % within-session acute
P = load(fullfile(CACHE_ROOT,'persistenceResults.mat'));   % drug-free epochs

%% ---- split the ACUTE pre-infusion baselines ----------------------------
% A.Tpre holds the pre-infusion epoch sessions; A.oP is its duty over the
% whole epoch. Recompute as two half-estimates.
% A.ok is the PAIRWISE mask (usable in both pre and post epochs) that
% runInfusionAnalysis used when it saved aa/cc/oP/oQ -- those arrays are
% already filtered, while Tpre/Tpost were saved whole. Select Tpre with the
% same mask so every array below has matching length and ordering.
Ppre = A.Tpre(A.ok);
nS = numel(Ppre);
aPreA = nan(nS,1); aPreB = nan(nS,1);
for i = 1:nS
    [aPreA(i), aPreB(i)] = splitEpochDuty(Ppre(i), WIDTH, STEP, BLOCK);
end
fprintf('acute baselines split: %d sessions, %d with both halves estimable\n', ...
    nS, sum(isfinite(aPreA) & isfinite(aPreB)));
fprintf('half-to-half agreement (should be high if the split is stable): rho = %.3f\n', ...
    corr(aPreA(isfinite(aPreA)&isfinite(aPreB)), aPreB(isfinite(aPreA)&isfinite(aPreB)), ...
    'type','Spearman'));

acuteAnim = A.aa(:); acuteCond = A.cc(:);   % already filtered by A.ok
acutePost = A.oQ(:);
if numel(acuteAnim) ~= nS
    error('responderPersistenceSplit:align', ...
        'acute arrays (%d) do not match selected Tpre sessions (%d).', ...
        numel(acuteAnim), nS);
end
isASVsess = ~cellfun(@isempty, regexp(acuteCond,'^ASV','once'));

%% ---- per-animal acute change, using the A half only --------------------
ua = unique(acuteAnim(isASVsess));
acuteChangeA = nan(numel(ua),1);
for i = 1:numel(ua)
    m = isASVsess & strcmp(acuteAnim, ua{i});
    acuteChangeA(i) = mean(acutePost(m) - aPreA(m), 'omitnan');
end

%% ---- split the PERSISTENCE drug-naive epochs ---------------------------
E = P.T(P.ok);            % same mask runPersistenceAnalysis used for ea/ec/oE
persAnim = P.ea; persPhase = P.phase; persCond = P.ec;
nE = numel(E);
if numel(P.ea) ~= nE
    error('responderPersistenceSplit:align', ...
        'persistence arrays (%d) do not match selected epochs (%d).', numel(P.ea), nE);
end
preA = nan(nE,1); preB = nan(nE,1);
for i = 1:nE
    [preA(i), preB(i)] = splitEpochDuty(E(i), WIDTH, STEP, BLOCK);
end

isASVday = ~cellfun(@isempty, regexp(persCond,'^ASV$','once'));
persLast = P.dLast(:);
selPre  = strcmp(persPhase,'pre');
selPost = strcmp(persPhase,'post') & (persLast >= MIN_DAYS_POST);
if EXCLUDE_ASVDAY_FROM_POST
    selPost = selPost & ~isASVday;
end
fprintf('\npost epochs: %d (>= %d day since dose%s); conditions: %s\n', ...
    sum(selPost), MIN_DAYS_POST, ...
    local_tern(EXCLUDE_ASVDAY_FROM_POST,', dosing days excluded',''), ...
    strjoin(unique(persCond(selPost))', ', '));

up = unique(persAnim);
persStepB = nan(numel(up),1);      % reference = B half of the naive epochs
for i = 1:numel(up)
    yp = preB(selPre  & strcmp(persAnim, up{i}));      % B half only
    yq = P.oE(selPost & strcmp(persAnim, up{i}));      % post epochs, full
    yp = yp(isfinite(yp)); yq = yq(isfinite(yq));
    if ~isempty(yp) && ~isempty(yq), persStepB(i) = mean(yq) - mean(yp); end
end

%% ---- table + tests ------------------------------------------------------
fprintf('\n%-8s %14s %14s\n','animal','acuteChange(A)','persistStep(B)');
xa = []; yy = []; who = {};
for i = 1:numel(up)
    j = find(strcmp(ua, up{i}), 1);
    if isempty(j), continue; end
    fprintf('%-8s %14s %14s\n', up{i}, local_f(acuteChangeA(j)), local_f(persStepB(i)));
    if isfinite(acuteChangeA(j)) && isfinite(persStepB(i))
        xa(end+1,1) = acuteChangeA(j); yy(end+1,1) = persStepB(i); who{end+1} = up{i}; %#ok<AGROW>
    end
end

resp = xa > 0;
fprintf('\nresponders (acute change > 0, measured on A half): %d/%d -- %s\n', ...
    sum(resp), numel(resp), strjoin(who(resp), ', '));
if sum(resp) >= 2
    d = yy(resp);
    p = local_signflip(d, 1e4);
    fprintf('persistence step in responders: mean %+.3f, %d/%d positive, p = %.4f\n', ...
        mean(d), sum(d>0), numel(d), p);
    if numel(d) < 6
        fprintf('  (exact floor at n=%d is %.3f)\n', numel(d), 2/2^numel(d));
    end
end
if sum(resp) >= 2 && sum(~resp) >= 2
    fprintf('responders vs non: %+.3f vs %+.3f, p = %.4f\n', ...
        mean(yy(resp)), mean(yy(~resp)), permtest(yy(resp), yy(~resp)));
end
if numel(yy) >= 4
    [rho,pr] = corr(xa, yy, 'type','Spearman');
    fprintf('\nacuteChange(A) vs persistStep(B): rho %+.3f, p %.4f (n=%d animals)\n', ...
        rho, pr, numel(yy));
    fprintf(['This is the key number: selection statistic and test statistic now\n' ...
             'come from disjoint minutes, so a positive rho is not regression to\n' ...
             'the mean. Compare it to the unsplit version in responderPersistence.m\n' ...
             '-- if the unsplit rho was much larger, that gap was the bias.\n']);
end

%% =========================================================================
function s = local_tern(c,a,b)
if c, s = a; else, s = b; end
end

function s = local_f(v)
if isnan(v), s = '-'; else, s = sprintf('%+.3f', v); end
end

function p = local_signflip(d, nperm)
d = d(:); d = d(isfinite(d)); n = numel(d); obs = mean(d);
null = zeros(nperm,1);
for k = 1:nperm, null(k) = mean((2*(rand(n,1)>0.5)-1).*d); end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end