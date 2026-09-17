%% runPersistenceAnalysis.m
% Does the ASV effect persist? Compares DRUG-FREE epochs before an animal's
% first ASV exposure with drug-free epochs after it, then stratifies the
% post-exposure epochs by time since the most recent exposure.
%
% THE KEY MOVE: EVERY SESSION CONTRIBUTES A DRUG-FREE EPOCH
%   infusion session -> the baseline epoch BEFORE infusion
%   no-infusion session -> an epoch from the start of the recording
% So the baseline of ASV session #3 is post-exposure but drug-free, which is
% precisely the state the persistence question is about. It also gives
% animals whose recordings begin at ASV (CP4, CP5, CP7 had zero pre-phase
% sessions in the earlier analysis) a genuine drug-naive data point: the
% baseline of their FIRST ASV session, recorded before any drug was given.
%
% FIXED EPOCH LENGTH ACROSS ALL SESSIONS -- different from the within-session
% infusion analysis, where L could vary because every test was paired inside
% one recording. Here epochs are compared ACROSS sessions, so detection
% sensitivity must be identical everywhere: a single LFIX for all of them.
% Sessions that cannot supply LFIX drug-free seconds are dropped, not
% shortened.
%
% EXPOSURE EVENTS ARE 'ASV' SESSIONS ONLY. ASV24 is the follow-up recording
% of the SAME exposure, not a new dose, so it does not reset "time since
% last exposure". Change EXPOSURE_PAT if that is wrong.
%
% THE CONFOUND, AGAIN. Post-exposure is always later in the study than
% pre-exposure, so a persistent drug effect and slow drift predict the same
% sign. Two things separate them here, neither available before:
%   1. DaysInStudy enters the model alongside Phase (step vs trend), and
%   2. time since LAST exposure is not monotonic with time in study -- an
%      animal can be 1 day post one exposure and 30 days post another later
%      on. A dose-recency relationship that is not explained by DaysInStudy
%      is hard for drift to mimic.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
WIDTH = 60; STEP = 5; MATCHN = 20;
LFIX  = 3600;                 % fixed drug-free epoch length (s)
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
EXPOSURE_PAT = '^ASV$';       % exact: ASV24 is not a new exposure

fils = getAllExtFiles('R:\TransPlasticity','mat',1);
fils = fils(contains(fils,'spikes.cellinfo.mat'));

cnd = cellfun(@(f) conditionFromPath(f,'TransPlasticity'), fils, 'UniformOutput',false);
ani = cellfun(@(f) animalFromPath(f,'TransPlasticity'),    fils, 'UniformOutput',false);
dnm = cellfun(@(f) sessionDateFromPath(f),                 fils);

keep = ismember(cnd, CONDS) & ~cellfun(@isempty, ani) & isfinite(dnm);
fils = fils(keep); cnd = cnd(keep); ani = ani(keep); dnm = dnm(keep);
fprintf('%d sessions in the %d requested conditions\n', numel(fils), numel(CONDS));

%% ---- 1. drug-free epoch per session -------------------------------------
inf_t = nan(numel(fils),1);
for i = 1:numel(fils), inf_t(i) = infusionTime(fils{i}); end
hasInf = ~isnan(inf_t);
fprintf('%d have an infusion marker; %d do not (epoch taken from recording start)\n', ...
    sum(hasInf), sum(~hasInf));

scfg = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',10,'NullMethod','none','DetectTransients',false, ...
         'CacheDir', fullfile(CACHE_ROOT,'screenPersist') };
Sc = syncBatch(fils, scfg{:});
durs = [Sc.dur]';

% available drug-free seconds: up to infusion, or the whole recording
avail = durs;
avail(hasInf) = min(inf_t(hasInf), durs(hasInf));

usable = [Sc.ok]' & avail >= LFIX;
fprintf('\n%d/%d sessions can supply a %.0f-min drug-free epoch\n', ...
    sum(usable), numel(usable), LFIX/60);
for c = unique(cnd(usable))'
    fprintf('   %-10s %2d sessions, %d animals\n', c{1}, sum(usable & strcmp(cnd,c{1})), ...
        numel(unique(ani(usable & strcmp(cnd,c{1})))));
end

uf = fils(usable); ua = ani(usable); uc = cnd(usable);
ud = dnm(usable);  uav = avail(usable);

% epoch = the LFIX seconds immediately BEFORE infusion (or from the start)
w0 = uav - LFIX;
mapW = containers.Map(uf, num2cell([w0, w0 + LFIX], 2));

%% ---- 2. detection on the drug-free epochs -------------------------------
cfgB = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',10,'MatchN',MATCHN,'NullMethod','shift','NSurr',50, ...
         'DetectTransients',true,'NSurrWindow',50,'RunAlpha',0.05, ...
         'Parallel', ~isempty(ver('parallel')), 'CheckpointEvery',5 };

T = syncBatch(uf, cfgB{:}, 'Window', @(f) mapW(f), ...
    'CacheDir', fullfile(CACHE_ROOT,'persistEpoch'), ...
    'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_persist.mat'));

ok = [T.ok]';
E = T(ok); ea = ua(ok); ec = uc(ok); ed = ud(ok);
eHasInf = hasInf(usable); eHasInf = eHasInf(ok);   % epoch precedes an infusion?
fprintf('\n%d/%d drug-free epochs usable after detection\n', sum(ok), numel(ok));

[rE,oE,dE] = deal(nan(numel(E),1));
for i = 1:numel(E), [rE(i),oE(i),dE(i)] = transientStats(E(i), STEP, WIDTH); end
mE = arrayfun(@(r) meanPeakMagnitude(r.transients), E);

%% ---- 3. exposure timing --------------------------------------------------
% For each epoch: days since FIRST exposure (negative = drug-naive) and days
% since the MOST RECENT exposure at or before that date.
uan = unique(ea);
dFirst = nan(numel(E),1); dLast = nan(numel(E),1);
for i = 1:numel(uan)
    m = strcmp(ea, uan{i});
    expDates = unique(ed(m & ~cellfun(@isempty, regexp(ec, EXPOSURE_PAT,'once'))));
    if isempty(expDates)
        continue                      % never exposed: both stay NaN
    end
    idx = find(m);
    for j = idx'
        % Whether an exposure on the SAME date counts depends on where the
        % epoch sits in that recording:
        %   infusion session -> the epoch is taken BEFORE the infusion, so
        %      that day's dose has not been delivered yet as far as this
        %      epoch is concerned. Use strictly-before. Without this, the
        %      pre-infusion baseline of an animal's very first ASV session --
        %      genuinely drug-naive data -- lands in the post-exposure group,
        %      simultaneously starving the naive group and diluting the post
        %      group with pre-dose recordings.
        %   no-infusion session (ASV24, SAL24, Baseline) -> the epoch is from
        %      the start of the recording and any exposure dated that day or
        %      earlier has already happened. Use on-or-before.
        if eHasInf(j)
            prior = expDates(expDates <  ed(j));
        else
            prior = expDates(expDates <= ed(j));
        end
        if ~isempty(prior)
            dFirst(j) = ed(j) - min(prior);
            dLast(j)  = ed(j) - max(prior);
        else
            dFirst(j) = -1;           % drug-naive epoch (any negative value)
            dLast(j)  = NaN;
        end
    end
end
phase = repmat({'naive'}, numel(E), 1);
phase(dFirst >= 0)  = {'post'};
phase(dFirst == -1) = {'pre'};

fprintf('\nepochs: %d pre-exposure (drug-naive), %d post-exposure, %d from never-exposed animals\n', ...
    sum(strcmp(phase,'pre')), sum(strcmp(phase,'post')), sum(strcmp(phase,'naive')));
fprintf('  of the drug-naive epochs, %d are pre-infusion baselines from an\n', ...
    sum(strcmp(phase,'pre') & eHasInf));
fprintf('  animal''s FIRST exposure session -- data that an on-or-before rule\n');
fprintf('  would have misassigned to the post group.\n');
nBoth = sum(arrayfun(@(i) any(strcmp(ea,uan{i}) & strcmp(phase,'pre')) && ...
                          any(strcmp(ea,uan{i}) & strcmp(phase,'post')), 1:numel(uan)));
fprintf('%d/%d animals contribute BOTH pre and post drug-free epochs\n', nBoth, numel(uan));

fprintf('\ndays since last exposure (post epochs): ');
fprintf('%s\n', mat2str(sort(dLast(strcmp(phase,'post') & isfinite(dLast)))'));

%% ---- 4. pre vs post, with a drift control -------------------------------
sel = ismember(phase, {'pre','post'});
daysInStudy = nan(numel(E),1);
for i = 1:numel(uan)
    m = strcmp(ea, uan{i});
    daysInStudy(m) = ed(m) - min(ed(m));
end

MET = {'duty cycle', oE; 'transient rate', rE; 'log duration', log(dE); 'log magnitude', log(mE)};
fprintf('\n%-15s %9s %9s %9s %9s\n','metric','pre','post','p_step','p_step|time');
for k = 1:size(MET,1)
    y = MET{k,2};
    g = sel & isfinite(y);
    if numel(unique(phase(g))) < 2, continue; end
    tb = table(y(g), categorical(ea(g)), categorical(phase(g),{'pre','post'}), daysInStudy(g), ...
        'VariableNames',{'M','Animal','Phase','Days'});
    m1 = fitlme(tb,'M ~ Phase + (1|Animal)','FitMethod','ML');
    m3 = fitlme(tb,'M ~ Days + Phase + (1|Animal)','FitMethod','ML');
    a1 = anova(m1,'DFMethod','Satterthwaite'); a3 = anova(m3,'DFMethod','Satterthwaite');
    fprintf('%-15s %9.3f %9.3f %9.4f %9.4f\n', MET{k,1}, ...
        mean(y(g & strcmp(phase,'pre'))), mean(y(g & strcmp(phase,'post'))), ...
        a1.pValue(strcmp(a1.Term,'Phase')), a3.pValue(strcmp(a3.Term,'Phase')));
end
fprintf(['p_step alone is not evidence -- post is always later. p_step|time is\n' ...
         'the step AFTER adjusting for a smooth trend; that is the one to read.\n']);

%% ---- 5. stratify by time since last exposure ----------------------------
BINS = [0 1; 2 7; 8 30; 31 inf];
fprintf('\nstratified by days since last exposure (drug-free epochs only)\n');
fprintf('%-15s %10s', 'metric', 'pre-expo');
for b = 1:size(BINS,1), fprintf('%10s', sprintf('%g-%g d', BINS(b,1), BINS(b,2))); end
fprintf('\n');
for k = 1:size(MET,1)
    y = MET{k,2};
    fprintf('%-15s %10.3f', MET{k,1}, mean(y(strcmp(phase,'pre') & isfinite(y)),'omitnan'));
    for b = 1:size(BINS,1)
        g = strcmp(phase,'post') & dLast >= BINS(b,1) & dLast <= BINS(b,2) & isfinite(y);
        if sum(g) == 0, fprintf('%10s','-'); else, fprintf('%10.3f', mean(y(g))); end
    end
    fprintf('\n');
end
fprintf('%-15s %10d', 'n epochs', sum(strcmp(phase,'pre')));
for b = 1:size(BINS,1)
    fprintf('%10d', sum(strcmp(phase,'post') & dLast >= BINS(b,1) & dLast <= BINS(b,2)));
end
fprintf('\n');
fprintf('%-15s %10d', 'n animals', numel(unique(ea(strcmp(phase,'pre')))));
for b = 1:size(BINS,1)
    g = strcmp(phase,'post') & dLast >= BINS(b,1) & dLast <= BINS(b,2);
    fprintf('%10d', numel(unique(ea(g))));
end
fprintf('\n');

% continuous recency model within post epochs: does the effect decay?
fprintf('\nrecency (post epochs only): M ~ log(1+daysSinceLast) + Days + (1|Animal)\n');
for k = 1:size(MET,1)
    y = MET{k,2};
    g = strcmp(phase,'post') & isfinite(y) & isfinite(dLast);
    if sum(g) < 8, continue; end
    tb = table(y(g), categorical(ea(g)), log(1+dLast(g)), daysInStudy(g), ...
        'VariableNames',{'M','Animal','LogRecency','Days'});
    mm = fitlme(tb,'M ~ LogRecency + Days + (1|Animal)','FitMethod','REML');
    aa2 = anova(mm,'DFMethod','Satterthwaite');
    ii = find(strcmp(mm.CoefficientNames,'LogRecency'));
    fprintf('  %-15s slope %+7.4f  p = %.4f  (n=%d epochs, %d animals)\n', MET{k,1}, ...
        mm.Coefficients.Estimate(ii), aa2.pValue(strcmp(aa2.Term,'LogRecency')), ...
        sum(g), numel(unique(ea(g))));
end
fprintf(['  A negative slope means the metric returns toward baseline as time\n' ...
         '  since the last dose grows -- i.e. the change is NOT permanent.\n' ...
         '  A flat slope with an elevated post level is what persistence looks like.\n']);

%% ---- 6. focused 24 h comparison -----------------------------------------
% The recency distribution is dominated by short intervals, so the one
% timepoint with real n is 24 h post exposure (dLast == 1: the drug-free
% baseline epochs of ASV24 sessions). Compare those against drug-naive
% epochs from the SAME animals -- paired, so between-animal baseline
% differences cancel rather than being assumed away.
is24  = strcmp(phase,'post') & dLast == 1;
isNv  = strcmp(phase,'pre');
fprintf('\n=== 24 h post exposure vs drug-naive ===\n');
fprintf('%d epochs at 24 h (%d animals); %d drug-naive epochs (%d animals)\n', ...
    sum(is24), numel(unique(ea(is24))), sum(isNv), numel(unique(ea(isNv))));

pairAn = intersect(unique(ea(is24)), unique(ea(isNv)));
fprintf('%d animals have BOTH -> paired comparison uses these\n', numel(pairAn));

if numel(pairAn) >= 3
    fprintf('\n%-15s %9s %9s %9s %8s %8s\n', ...
        'metric','naive','24h','diff','p_pair','n_anim');
    for k = 1:size(MET,1)
        y = MET{k,2};
        dd = nan(numel(pairAn),1);
        for j = 1:numel(pairAn)
            v0 = y(isNv  & strcmp(ea,pairAn{j}) & isfinite(y));
            v1 = y(is24  & strcmp(ea,pairAn{j}) & isfinite(y));
            if ~isempty(v0) && ~isempty(v1), dd(j) = mean(v1) - mean(v0); end
        end
        g = isfinite(dd);
        if sum(g) < 3, continue; end
        % exact sign-flip on within-animal differences
        nn = sum(g); obs = mean(dd(g)); nullv = zeros(1e4,1);
        for q = 1:1e4, nullv(q) = mean((2*(rand(nn,1)>0.5)-1).*dd(g)); end
        pp = (1 + sum(abs(nullv) >= abs(obs)))/(1+1e4);
        fprintf('%-15s %9.3f %9.3f %+9.3f %8.4f %8d\n', MET{k,1}, ...
            mean(y(isNv & isfinite(y))), mean(y(is24 & isfinite(y))), obs, pp, nn);
    end
    fprintf(['\nWith n=%d animals the smallest reachable paired p is %.3f --\n' ...
             'check that before reading a null as absence of an effect.\n'], ...
        numel(pairAn), 2/2^numel(pairAn));
else
    fprintf('too few animals with both -- paired 24h comparison not possible.\n');
end

save(fullfile(CACHE_ROOT,'persistenceResults.mat'), ...
    'T','ok','ea','ec','ed','eHasInf','rE','oE','dE','mE','phase','dFirst','dLast', ...
    'daysInStudy','LFIX','BINS');

%% ---- 6. FOCUSED: drug-naive vs 24 h post-exposure ------------------------
% The single pre-specified contrast the data can actually support. The
% recency stratification above is too imbalanced to read (most post epochs
% sit within a day of a dose and the far bins hold 2-3 animals), but the 24 h
% point has real numbers behind it: these are ASV24-day epochs, recorded a
% full day after the dose with no drug on board.
%
% Both a paired within-animal test (animals contributing BOTH a drug-naive
% epoch and a 24 h epoch) and an LME across all such epochs with the
% time-in-study covariate. The paired test is exact and immune to
% between-animal baseline differences -- which you have independently
% observed -- so it is the one to report.

isPre  = strcmp(phase,'pre');
is24   = strcmp(phase,'post') & dLast == 1;

fprintf('\n=== drug-naive vs 24 h post-exposure ===\n');
fprintf('%d drug-naive epochs (%d animals), %d 24-h epochs (%d animals)\n', ...
    sum(isPre), numel(unique(ea(isPre))), sum(is24), numel(unique(ea(is24))));

bothA = intersect(unique(ea(isPre)), unique(ea(is24)));
fprintf('%d animals contribute BOTH -> paired test uses these\n', numel(bothA));
if numel(bothA) >= 2
    fprintf('   %s\n', strjoin(bothA', ', '));
end

fprintf('\n%-15s %9s %9s %9s %8s %8s %6s\n', ...
    'metric','naive','24h','diff','p_paired','p_LME','n_pair');
for k = 1:size(MET,1)
    y = MET{k,2};
    gP = isPre & isfinite(y); g24 = is24 & isfinite(y);
    if sum(gP) < 2 || sum(g24) < 2, continue; end

    % paired: one value per animal per phase, then sign-flip on the difference
    d = nan(numel(bothA),1);
    for j = 1:numel(bothA)
        d(j) = mean(y(g24 & strcmp(ea,bothA{j}))) - mean(y(gP & strcmp(ea,bothA{j})));
    end
    d = d(isfinite(d));
    if numel(d) >= 2
        pPair = local_signflip2(d, 1e4);
    else
        pPair = NaN;
    end

    % LME across all epochs of both phases, adjusting for time in study
    g = (gP | g24);
    tb = table(y(g), categorical(ea(g)), ...
        categorical(phase(g),{'pre','post'}), daysInStudy(g), ...
        'VariableNames',{'M','Animal','Phase','Days'});
    pL = NaN;
    try
        mm = fitlme(tb,'M ~ Days + Phase + (1|Animal)','FitMethod','ML');
        aa3 = anova(mm,'DFMethod','Satterthwaite');
        pL = aa3.pValue(strcmp(aa3.Term,'Phase'));
    catch
    end

    fprintf('%-15s %9.3f %9.3f %+9.3f %8.4f %8.4f %6d\n', MET{k,1}, ...
        mean(y(gP)), mean(y(g24)), mean(y(g24))-mean(y(gP)), pPair, pL, numel(d));
end

fprintf(['\nn_pair is animals with BOTH a drug-naive and a 24-h epoch. A paired\n' ...
         'test on fewer than ~6 animals cannot reach p<0.05 by sign-flip alone\n' ...
         '(2/2^n floor), so read the effect size and the LME there, not the p.\n']);

function p = local_signflip2(d, nperm)
d = d(:); d = d(isfinite(d));
n = numel(d); obs = mean(d);
null = zeros(nperm,1);
for k = 1:nperm, null(k) = mean((2*(rand(n,1)>0.5)-1) .* d); end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end