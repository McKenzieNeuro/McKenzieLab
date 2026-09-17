%% stepAnalysis.m
% Is the post-exposure change a STEP (shift to a new level that stays) or a
% TREND (progressive change that keeps going)? Run after
% runPersistenceAnalysis.m.
%
% WHY THE EARLIER TEST WAS THE WRONG SHAPE
% Section 4 of runPersistenceAnalysis fits  M ~ Days + Phase + (1|Animal)
% and reads the Phase p-value. When the truth is a step and pre/post barely
% overlap in time, a linear Days term absorbs much of that step, so the
% Phase term is biased toward null. That guard is right when the worry is
% drift, but it makes the model conservative against precisely the
% hypothesis of interest. This script tests the SHAPE instead of assuming it.
%
% THE DISCRIMINATING PREDICTIONS
%   step  : slope within pre  ~ 0
%           slope within post ~ 0
%           level(post) - level(pre) /= 0
%   trend : a single slope across everything, with no special status for the
%           exposure date; level difference is just the slope integrated
%
% So three things get reported per metric:
%   1. per-animal pre vs post means and the within-animal step, with an
%      exact sign-flip test across animals (no time covariate -- the step is
%      the hypothesis, not a nuisance)
%   2. slopes fitted SEPARATELY within pre and within post epochs. A step
%      requires both to be flat.
%   3. AIC comparison of step-only, trend-only, and step+trend.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
load(fullfile(CACHE_ROOT,'persistenceResults.mat'), ...
     'ea','ec','ed','eHasInf','rE','oE','dE','mE','phase','dFirst','dLast','daysInStudy');

MET = {'duty cycle', oE; 'transient rate', rE; 'log duration', log(dE); 'log magnitude', log(mE)};

%% ---- which epochs count as POST ----------------------------------------
% An ASV-day baseline that follows an EARLIER exposure is ambiguous: it is
% post-exposure with respect to the previous dose, but it is also the epoch
% immediately preceding that day's infusion, on a day the animal is being
% handled and dosed. Excluding those from the post group leaves a post group
% made only of days on which no drug was given at all -- ASV24, SAL, SAL24,
% Baseline -- which is the cleaner probe of a persistent state.
%
% PRE epochs are NOT filtered: a first-exposure ASV-day baseline is genuine
% drug-naive data and is exactly what the pre group needs.
% Two inclusion rules, both reported. They are NOT equivalent and the choice
% must be made on principle, not on which gives the larger effect:
%
%   MIN_DAYS_POST = 1, EXCLUDE_ASVDAY_FROM_POST = false
%       Every drug-free epoch at least a day after a dose. No drug on board,
%       so it tests persistence regardless of which day it sits on. Keeps the
%       most data. Includes baselines recorded ON later dosing days, where
%       handling and anticipation are alternative explanations.
%
%   EXCLUDE_ASVDAY_FROM_POST = true
%       Post restricted to non-dosing days (ASV24, SAL, SAL24). Cleaner
%       interpretation, fewer epochs and fewer animals.
%
% For the record, on this dataset the duty-cycle step was +0.139 (6/8 animals
% positive) with ASV-day epochs included and +0.071 (4/7) with them excluded.
% Report both; presenting only the larger one after having seen both is the
% thing a reviewer will ask about.
MIN_DAYS_POST            = 1;      % post epochs must be >= this many days after a dose
EXCLUDE_ASVDAY_FROM_POST = false;

isASVday = ~cellfun(@isempty, regexp(ec, '^ASV$', 'once'));
tooSoon  = strcmp(phase,'post') & ~(dLast >= MIN_DAYS_POST);
dropPost = strcmp(phase,'post') & ( tooSoon | (EXCLUDE_ASVDAY_FROM_POST & isASVday) );

fprintf('post inclusion: >= %d day(s) since last dose', MIN_DAYS_POST);
if EXCLUDE_ASVDAY_FROM_POST, fprintf(', ASV-day epochs excluded'); end
fprintf('\n');
fprintf('  dropped %d post epochs (%d too soon, %d ASV-day); %d remain\n', ...
    sum(dropPost), sum(tooSoon), sum(EXCLUDE_ASVDAY_FROM_POST & strcmp(phase,'post') & isASVday), ...
    sum(strcmp(phase,'post') & ~dropPost));
fprintf('  post group conditions: %s\n', ...
    strjoin(unique(ec(strcmp(phase,'post') & ~dropPost))', ', '));
fprintf('  days since last dose in post group: %s\n', ...
    mat2str(unique(dLast(strcmp(phase,'post') & ~dropPost))'));

sel0 = ismember(phase,{'pre','post'}) & ~dropPost;

for k = 1:size(MET,1)
    y = MET{k,2};
    g = sel0 & isfinite(y);
    isPost = strcmp(phase,'post');

    fprintf('\n================ %s ================\n', MET{k,1});

    %% 1. per-animal step
    ua = unique(ea(g));
    rows = {};
    d = [];
    for i = 1:numel(ua)
        m = g & strcmp(ea, ua{i});
        yp = y(m & ~isPost); yq = y(m & isPost);
        rows(end+1,:) = {ua{i}, numel(yp), numel(yq), ...
            local_m(yp), local_m(yq), local_m(yq)-local_m(yp)}; %#ok<AGROW>
        if ~isempty(yp) && ~isempty(yq), d(end+1,1) = mean(yq)-mean(yp); end %#ok<AGROW>
    end
    fprintf('%-8s %4s %4s %9s %9s %9s\n','animal','nPre','nPost','pre','post','step');
    for i = 1:size(rows,1)
        fprintf('%-8s %4d %4d %9.3f %9.3f %+9.3f\n', rows{i,:});
    end
    if numel(d) >= 2
        p = local_signflip(d, 1e4);
        fprintf('within-animal step: mean %+.3f, %d/%d animals positive, p = %.4f (exact)\n', ...
            mean(d), sum(d>0), numel(d), p);
        if numel(d) < 6
            fprintf('  (n=%d -> exact p floor is %.3f; cannot reach 0.05)\n', ...
                numel(d), 2/2^numel(d));
        end
    else
        fprintf('too few animals with both phases for a paired step test\n');
    end

    %% 2. slopes within each phase -- a step needs both flat
    for ph = {'pre','post'}
        m = g & strcmp(phase, ph{1});
        if sum(m) < 6 || numel(unique(ea(m))) < 2, continue; end
        tb = table(y(m), categorical(ea(m)), daysInStudy(m), ...
            'VariableNames',{'M','Animal','Days'});
        mm = fitlme(tb,'M ~ Days + (1|Animal)','FitMethod','REML');
        aa = anova(mm,'DFMethod','Satterthwaite');
        b = mm.Coefficients.Estimate(strcmp(mm.CoefficientNames,'Days'));
        fprintf('slope within %-4s : %+8.5f /day  p = %.4f  (n=%d epochs, %d animals)\n', ...
            ph{1}, b, aa.pValue(strcmp(aa.Term,'Days')), sum(m), numel(unique(ea(m))));
    end

    %% 3. which shape fits better
    tb = table(y(g), categorical(ea(g)), categorical(phase(g),{'pre','post'}), daysInStudy(g), ...
        'VariableNames',{'M','Animal','Phase','Days'});
    if numel(unique(tb.Phase)) < 2, continue; end
    mStep  = fitlme(tb,'M ~ Phase + (1|Animal)','FitMethod','ML');
    mTrend = fitlme(tb,'M ~ Days + (1|Animal)','FitMethod','ML');
    mBoth  = fitlme(tb,'M ~ Days + Phase + (1|Animal)','FitMethod','ML');
    aS = anova(mStep,'DFMethod','Satterthwaite');
    fprintf('AIC  step %.1f | trend %.1f | both %.1f   -> %s fits best\n', ...
        mStep.ModelCriterion.AIC, mTrend.ModelCriterion.AIC, mBoth.ModelCriterion.AIC, ...
        local_best(mStep.ModelCriterion.AIC, mTrend.ModelCriterion.AIC, mBoth.ModelCriterion.AIC));
    fprintf('step p (no time covariate) = %.4f\n', aS.pValue(strcmp(aS.Term,'Phase')));
end

%% ---- per-animal trajectories: look before trusting any p ---------------
% One panel per animal, drug-free epoch metric against date, exposure dates
% marked. A step should be visible by eye; if it is not, no model will
% rescue it.
figure('Color','w','Name','drug-free epoch trajectories');
y = MET{1,2};                      % duty cycle by default; change index to switch
ua = unique(ea(sel0));
nA = numel(ua); nc = ceil(sqrt(nA)); nr = ceil(nA/nc);
for i = 1:nA
    subplot(nr,nc,i); hold on
    m = sel0 & strcmp(ea, ua{i}) & isfinite(y);
    if ~any(m), title(ua{i}); continue; end
    t = ed(m) - min(ed(m));
    pre = ~strcmp(phase(m),'post');
    yy = y(m);
    plot(t(pre),  yy(pre),  'ko', 'MarkerFaceColor','w');
    plot(t(~pre), yy(~pre), 'ro', 'MarkerFaceColor','r');
    expd = unique(ed(strcmp(ea,ua{i}) & ~cellfun(@isempty, regexp(ec,'^ASV$','once'))));
    for e = expd', xline(e - min(ed(m)), 'b--'); end
    title(sprintf('%s', ua{i})); xlabel('day'); ylabel(MET{1,1}); box off
end
sgtitle('open = drug-naive, filled red = post-exposure, blue = ASV exposure days');

%% =========================================================================
function v = local_m(x)
if isempty(x), v = NaN; else, v = mean(x); end
end

function s = local_best(a,b,c)
[~,i] = min([a b c]); opts = {'STEP','TREND','both'}; s = opts{i};
end

function p = local_signflip(d, nperm)
d = d(:); d = d(isfinite(d)); n = numel(d); obs = mean(d);
null = zeros(nperm,1);
for k = 1:nperm, null(k) = mean((2*(rand(n,1)>0.5)-1).*d); end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end