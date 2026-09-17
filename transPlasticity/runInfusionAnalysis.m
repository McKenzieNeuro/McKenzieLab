%% runInfusionAnalysis.m
% Pre-infusion baseline vs the 4-6 h post-infusion window, WITHIN session.
%
% DESIGN
% Each infusion session contains an hours-long baseline before infusion and,
% if long enough, the 4-6 h post-infusion window where the drug is expected
% to act. Comparing those two epochs within one recording holds constant
% everything that confounded the earlier analyses: electrodes, sorted units,
% animal, day, implant age, recording order. No normalisation is applied --
% the metrics are absolute and comparable by construction.
%
% EPOCH LENGTHS ARE MATCHED WITHIN SESSION.
% Transient significance comes from a per-epoch surrogate max-run
% distribution, which depends on how many windows the epoch holds. A long
% baseline compared against a short post window would differ for that reason
% alone. So per session:
%     L    = min(baseline available, post-window available, LMAX)
%     pre  = [inf - L, inf]
%     post = [inf + 4h, inf + 4h + L]
% Both epochs then contain the same number of windows and face the same
% detection threshold. L varies across sessions, which is fine because every
% test below is paired WITHIN session.
%
% THE CONTROL THAT MAKES THIS INTERPRETABLE
% An epoch 4-6 h later differs from baseline for reasons unrelated to drug --
% sleep pressure, satiety, electrode settling, time of day. Comparing the
% pre->post CHANGE in ASV sessions against the pre->post change in SAL
% sessions (difference-in-differences) subtracts that common time course.
% The DiD is the estimate to report; the within-ASV pre/post change on its
% own conflates drug with time-in-session.
warning off
CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
WIDTH = 60; STEP = 5; MATCHN = 15;
LAT   = 4*3600;      % start of the treatment window, s after infusion
LMAX  = 2*3600;      % cap on epoch length
LMIN  = 1800;        % need at least this much on both sides (~360 windows)

fils1 = getAllExtFiles('R:\TransPlasticity','mat',1);
fils1 = fils1(contains(fils1,'spikes.cellinfo.mat'));

%% ---- 1. infusion markers + feasibility ----------------------------------
inf_t = nan(numel(fils1),1);
for i = 1:numel(fils1), inf_t(i) = infusionTime(fils1{i}); end
has = ~isnan(inf_t);

cnd = cellfun(@(f) conditionFromPath(f,'TransPlasticity'), fils1, 'UniformOutput', false);
ani = cellfun(@(f) animalFromPath(f,'TransPlasticity'),    fils1, 'UniformOutput', false);

fprintf('%d/%d sessions carry an Infusion_start marker\n', sum(has), numel(fils1));
fprintf('by condition:\n');
for c = unique(cnd(has))'
    fprintf('   %-16s %d\n', c{1}, sum(has & strcmp(cnd,c{1})));
end
fprintf('\ninfusion time (s into recording): p0 %.0f | p50 %.0f | p100 %.0f\n', ...
    min(inf_t(has)), median(inf_t(has)), max(inf_t(has)));
fprintf('   -> baseline epoch available: median %.1f h\n', median(inf_t(has))/3600);
fprintf('CHECK these look like SECONDS. Values in ms or samples would silently\n');
fprintf('mis-split every session.\n');

% session durations, from a cheap screen (also applies the unit criteria)
scfg = { 'RateRange',[0.1 10], 'Width',WIDTH, 'Step',STEP, 'MinSpikes',10, ...
         'MinUnits',10, 'NullMethod','none', 'DetectTransients',false, ...
         'CacheDir', fullfile(CACHE_ROOT,'screenInf') };
Sc = syncBatch(fils1(has), scfg{:});
durs = [Sc.dur]';
it   = inf_t(has); cn = cnd(has); an = ani(has); fl = fils1(has);

Lpre  = it;
Lpost = max(0, min(durs, it + LAT + LMAX) - (it + LAT));
L     = min([Lpre, Lpost, repmat(LMAX, numel(it), 1)], [], 2);

usable = [Sc.ok]' & L >= LMIN;
fprintf('\n%d/%d infusion sessions usable (>=%.0f min matched on both sides)\n', ...
    sum(usable), numel(usable), LMIN/60);
fprintf('matched epoch L: median %.0f min (range %.0f-%.0f)\n', ...
    median(L(usable))/60, min(L(usable))/60, max(L(usable))/60);
fprintf('usable by condition:\n');
for c = unique(cn(usable))'
    fprintf('   %-16s %d sessions, %d animals\n', c{1}, sum(usable & strcmp(cn,c{1})), ...
        numel(unique(an(usable & strcmp(cn,c{1})))));
end
if sum(usable & ~cellfun(@isempty, regexp(cn,'^SAL','once'))) < 3
    warning(['Few/no usable SAL infusion sessions -- without them there is no ' ...
             'difference-in-differences control, and a within-ASV pre/post change ' ...
             'cannot be separated from the normal time course of a long recording.']);
end

uf = fl(usable); uL = L(usable); ut = it(usable);
uc = cn(usable);  ua = an(usable);

%% ---- 2. detection on each epoch -----------------------------------------
mapPre  = containers.Map(uf, num2cell([ut - uL,        ut],             2));
mapPost = containers.Map(uf, num2cell([ut + LAT, ut + LAT + uL],        2));

cfgB = { 'RateRange',[0.1 10], 'Width',WIDTH, 'Step',STEP, 'MinSpikes',10, ...
         'MinUnits',10, 'MatchN',MATCHN, 'NullMethod','shift', 'NSurr',50, ...
         'DetectTransients',true, 'NSurrWindow',50, 'RunAlpha',0.05, ...
         'Parallel', ~isempty(ver('parallel')), 'CheckpointEvery',5 };

Tpre  = syncBatch(uf, cfgB{:}, 'Window', @(f) mapPre(f), ...
    'CacheDir', fullfile(CACHE_ROOT,'infPre'), ...
    'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_infPre.mat'));
Tpost = syncBatch(uf, cfgB{:}, 'Window', @(f) mapPost(f), ...
    'CacheDir', fullfile(CACHE_ROOT,'infPost'), ...
    'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_infPost.mat'));

ok = [Tpre.ok]' & [Tpost.ok]';
fprintf('\n%d/%d sessions usable in BOTH epochs\n', sum(ok), numel(ok));

%% ---- 3. metrics ----------------------------------------------------------
P = Tpre(ok); Q = Tpost(ok); cc = uc(ok); aa = ua(ok);
n = sum(ok);
[rP,oP,dP] = deal(nan(n,1)); [rQ,oQ,dQ] = deal(nan(n,1));
for i = 1:n
    [rP(i),oP(i),dP(i)] = transientStats(P(i), STEP, WIDTH);
    [rQ(i),oQ(i),dQ(i)] = transientStats(Q(i), STEP, WIDTH);
end
mP = arrayfun(@(r) meanPeakMagnitude(r.transients), P);
mQ = arrayfun(@(r) meanPeakMagnitude(r.transients), Q);

isASV = ~cellfun(@isempty, regexp(cc,'^ASV','once'));
isSAL = ~cellfun(@isempty, regexp(cc,'^SAL','once'));

MET = { 'duty cycle',     oP,       oQ
        'transient rate', rP,       rQ
        'log duration',   log(dP),  log(dQ)
        'log magnitude',  log(mP),  log(mQ) };

fprintf('\n%-15s %-6s %8s %8s %9s %8s %8s\n', ...
    'metric','arm','pre','post','change','p_sess','p_anim');
DiD = struct();
for k = 1:size(MET,1)
    xp = MET{k,2}; xq = MET{k,3}; d = xq - xp;
    for arm = {'ASV','SAL'}
        sel = (strcmp(arm{1},'ASV') & isASV) | (strcmp(arm{1},'SAL') & isSAL);
        g = sel & isfinite(d);
        if sum(g) < 3, continue; end
        pS = local_signflip(d(g), 1e4);
        au = unique(aa(g)); dA = nan(numel(au),1);
        for j = 1:numel(au), dA(j) = mean(d(g & strcmp(aa,au{j}))); end
        pA = local_signflip(dA, 1e4);
        fprintf('%-15s %-6s %8.3f %8.3f %+9.3f %8.4f %8.4f\n', ...
            MET{k,1}, arm{1}, mean(xp(g)), mean(xq(g)), mean(d(g)), pS, pA);
    end
    % difference-in-differences, at the animal level
    gA = isASV & isfinite(d); gS = isSAL & isfinite(d);
    if sum(gA) >= 3 && sum(gS) >= 3
        auA = unique(aa(gA)); dAA = arrayfun(@(j) mean(d(gA & strcmp(aa,auA{j}))), 1:numel(auA))';
        auS = unique(aa(gS)); dAS = arrayfun(@(j) mean(d(gS & strcmp(aa,auS{j}))), 1:numel(auS))';
        pDiD = permtest(dAA, dAS);
        fprintf('%-15s %-6s %8s %8s %+9.3f %8s %8.4f  <- DiD (ASV change - SAL change)\n', ...
            MET{k,1}, 'DiD', '', '', mean(dAA)-mean(dAS), '', pDiD);
        DiD.(matlab.lang.makeValidName(MET{k,1})) = struct( ...
            'ASVchange',mean(dAA),'SALchange',mean(dAS),'p',pDiD, ...
            'nASVanimals',numel(auA),'nSALanimals',numel(auS));
    end
end

fprintf(['\nThe DiD row is the estimate to report: it subtracts the pre->post\n' ...
         'change seen in vehicle sessions, so what remains is drug-specific\n' ...
         'rather than the normal time course of a long recording.\n']);

save(fullfile(CACHE_ROOT,'infusionResults.mat'), ...
    'Tpre','Tpost','ok','cc','aa','rP','rQ','oP','oQ','dP','dQ','mP','mQ', ...
    'uf','ut','uL','DiD','LAT','LMAX','LMIN');

%% =========================================================================
function p = local_signflip(d, nperm)
d = d(:); d = d(isfinite(d));
n = numel(d); obs = mean(d);
null = zeros(nperm,1);
for k = 1:nperm, null(k) = mean((2*(rand(n,1)>0.5)-1) .* d); end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end
