%% runBoutRateAnalysis.m
% Firing rate INSIDE significant synchrony bouts: does ASV change it, acutely
% and persistently? Run after runInfusionAnalysis.m and runPersistenceAnalysis.m.
%
% THREE QUANTITIES, AND THEY ANSWER DIFFERENT QUESTIONS
%   rateOverall : mean per-neuron firing rate across the whole epoch. This is
%                 a CONTROL, not a result: if ASV simply raised firing rate,
%                 that would need reporting even though C is rate-adaptive.
%   rateInBout  : mean per-neuron firing rate during significant bouts only.
%   ratio       : rateInBout / rateOverall -- how much neurons speed up (or
%                 do not) when the population synchronises. Self-normalising,
%                 so it is unaffected by any baseline rate difference, and it
%                 is the quantity that speaks to mechanism.
%
% Report the ratio as primary and rateOverall alongside it. Reporting
% rateInBout alone would confound "neurons fire faster in bouts" with
% "neurons fire faster overall".
%
% CENSORING. rateInBout and ratio exist only for epochs containing at least
% one significant bout, so both comparisons are conditional on a bout having
% occurred -- the same asymmetry that applies to duration and magnitude.
% rateOverall is defined for every epoch.
%
% A NOTE ON WHAT THE RATIO CANNOT SHOW. Bouts are detected as intervals of
% elevated C, and C is computed from the same spikes. If a rate increase
% itself made coincidence detection more likely, bouts would preferentially
% land on high-rate stretches and the ratio would exceed 1 with no change in
% coordination. The adaptive coincidence window makes C largely rate-
% invariant, which is why the measure was chosen, but the ratio should still
% be interpreted as descriptive of what bouts look like rather than as
% independent evidence that synchrony drives rate.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
WIDTH = 60; STEP = 5;

A = load(fullfile(CACHE_ROOT,'infusionResults.mat'));
P = load(fullfile(CACHE_ROOT,'persistenceResults.mat'));

%% ======================= ACUTE (within session) =========================
Ppre  = A.Tpre(A.ok);
Ppost = A.Tpost(A.ok);
aa    = A.aa(:); cc = A.cc(:);

n = numel(Ppre);
[ovPre,inPre,raPre]    = deal(nan(n,1));
[ovPost,inPost,raPost] = deal(nan(n,1));
for i = 1:n
    try
        o = boutRate(Ppre(i),  WIDTH, STEP);
        ovPre(i)=o.rateOverall; inPre(i)=o.rateInBout; raPre(i)=o.ratio;
    catch ME, warning('pre %d: %s', i, ME.message); end
    try
        o = boutRate(Ppost(i), WIDTH, STEP);
        ovPost(i)=o.rateOverall; inPost(i)=o.rateInBout; raPost(i)=o.ratio;
    catch ME, warning('post %d: %s', i, ME.message); end
end

isASV = ~cellfun(@isempty, regexp(cc,'^ASV','once'));
isSAL = ~cellfun(@isempty, regexp(cc,'^SAL','once'));

fprintf('\n================= ACUTE: firing rate in bouts =================\n');
fprintf('%-14s %-6s %9s %9s %10s %8s\n','quantity','arm','pre','post','change','p_anim');
MET = { 'overall rate (Hz)', ovPre,      ovPost
        'in-bout rate (Hz)', inPre,      inPost
        'log ratio',         log(raPre), log(raPost) };
for k = 1:size(MET,1)
    xp = MET{k,2}; xq = MET{k,3}; d = xq - xp;
    dAA = []; dAS = [];
    for arm = {'ASV','SAL'}
        sel = (strcmp(arm{1},'ASV') & isASV) | (strcmp(arm{1},'SAL') & isSAL);
        g = sel & isfinite(d);
        if sum(g) < 3, continue; end
        au = unique(aa(g)); dA = arrayfun(@(j) mean(d(g & strcmp(aa,au{j}))), 1:numel(au))';
        fprintf('%-14s %-6s %9.3f %9.3f %+10.3f %8.4f\n', MET{k,1}, arm{1}, ...
            mean(xp(g)), mean(xq(g)), mean(d(g)), local_signflip(dA,1e4));
        if strcmp(arm{1},'ASV'), dAA = dA; else, dAS = dA; end
    end
    if numel(dAA) >= 3 && numel(dAS) >= 3
        fprintf('%-14s %-6s %9s %9s %+10.3f %8.4f  <- DiD\n', MET{k,1}, 'DiD', '', '', ...
            mean(dAA)-mean(dAS), permtest(dAA, dAS));
    end
end

%% ==================== PERSISTENCE (drug-free epochs) ====================
E = P.T(P.ok); ea = P.ea(:); ec = P.ec(:); phase = P.phase(:); dLast = P.dLast(:);
m = numel(E);
[ovE,inE,raE] = deal(nan(m,1));
for i = 1:m
    try
        o = boutRate(E(i), WIDTH, STEP);
        ovE(i)=o.rateOverall; inE(i)=o.rateInBout; raE(i)=o.ratio;
    catch ME, warning('epoch %d: %s', i, ME.message); end
end

selPre  = strcmp(phase,'pre');
selPost = strcmp(phase,'post') & dLast >= 1;

fprintf('\n=============== PERSISTENCE: firing rate in bouts ===============\n');
fprintf('%-18s %9s %9s %10s %8s %8s %7s\n', ...
    'quantity','naive','post','step','p_exact','p_LME','nAnim');
MEP = { 'overall rate (Hz)', ovE
        'in-bout rate (Hz)', inE
        'log ratio',         log(raE) };
for k = 1:size(MEP,1)
    y = MEP{k,2};
    ua = unique(ea);
    d = []; 
    for i = 1:numel(ua)
        yp = y(selPre  & strcmp(ea,ua{i}) & isfinite(y));
        yq = y(selPost & strcmp(ea,ua{i}) & isfinite(y));
        if ~isempty(yp) && ~isempty(yq), d(end+1,1) = mean(yq)-mean(yp); end %#ok<AGROW>
    end
    g = (selPre | selPost) & isfinite(y);
    pL = NaN;
    if numel(unique(phase(g))) == 2
        tb = table(y(g), categorical(ea(g)), categorical(phase(g),{'pre','post'}), ...
            'VariableNames',{'M','Animal','Phase'});
        try
            lme = fitlme(tb,'M ~ Phase + (1|Animal)','FitMethod','REML');
            av = anova(lme,'DFMethod','Satterthwaite');
            pL = av.pValue(strcmp(av.Term,'Phase'));
        catch, end
    end
    fprintf('%-18s %9.3f %9.3f %+10.3f %8.4f %8.4f %7d\n', MEP{k,1}, ...
        mean(y(selPre & isfinite(y))), mean(y(selPost & isfinite(y))), ...
        mean(d), local_signflip(d,1e4), pL, numel(d));
end

fprintf(['\nRead the ratio row against the overall-rate row. A ratio change with\n' ...
         'a flat overall rate means bouts became more (or less) distinct from the\n' ...
         'surrounding activity. Both moving together means the whole epoch shifted\n' ...
         'and the bouts are not special.\n']);

save(fullfile(CACHE_ROOT,'boutRateResults.mat'), ...
    'ovPre','inPre','raPre','ovPost','inPost','raPost','aa','cc', ...
    'ovE','inE','raE','ea','ec','phase','dLast');

%% =========================================================================
function p = local_signflip(d, nperm)
d = d(:); d = d(isfinite(d));
if numel(d) < 2, p = NaN; return; end
nn = numel(d); obs = mean(d); null = zeros(nperm,1);
for k = 1:nperm, null(k) = mean((2*(rand(nn,1)>0.5)-1).*d); end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end
