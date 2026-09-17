%% checkFirstHourPrePost.m
% First hour PRE-ASV vs OptoRAM, and first hour POST-ASV vs OptoRAM.
% Position held fixed (first hour, both sides) so any difference here is
% about drug/dataset, not about where in the session the epoch sits -- the
% confound checkPositionConfound.m just demonstrated matters a great deal.
%
% "First hour pre-ASV" here means: the first hour of a TransPlasticity
% recording, on a day BEFORE the animal's first ASV exposure (drug-naive).
% "First hour post-ASV" means: the first hour of a TransPlasticity
% recording, on a day AT OR AFTER the first exposure (ASV, ASV24, SAL,
% SAL24 sessions once exposure has started) -- NOT the pre-infusion baseline
% used elsewhere; this is deliberately the start of ANY post-exposure
% session, since the position confound requires "first hour" specifically,
% and an ASV session's own first hour is typically its pre-infusion period
% anyway given multi-hour latency to infusion in these recordings.
%
% Requires checkPositionConfound.m to have been run (reuses B, dBvec, Ae from
% that script's Check 2, which is the first-60-min TransPlasticity/OptoRAM
% comparison already computed -- Ae there mixed pre/post ASV; here it is
% split by phase).

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
WIDTH = 60; STEP = 5; MATCHN = 15; LFIX = 3600;

load(fullfile(CACHE_ROOT,'crossNaiveResults.mat'), 'B','bAnim');
dBvec = nan(numel(B),1);
for i = 1:numel(B), [~,~,dBvec(i)] = transientStats(B(i), STEP, WIDTH); end
[rB,oB,~] = deal(nan(numel(B),1));
for i = 1:numel(B), [rB(i),oB(i)] = transientStats(B(i), STEP, WIDTH); end
mB = arrayfun(@(r) meanPeakMagnitude(r.transients), B);

%% ---- classify every TransPlasticity session by phase, at the SESSION level
fils1 = getAllExtFiles('R:\TransPlasticity','mat',1);
fils1 = fils1(contains(fils1,'spikes.cellinfo.mat'));
cnd = cellfun(@(f) conditionFromPath(f,'TransPlasticity'), fils1, 'UniformOutput',false);
ani = cellfun(@(f) animalFromPath(f,'TransPlasticity'),    fils1, 'UniformOutput',false);
dnm = cellfun(@(f) sessionDateFromPath(f),                 fils1);
keep = ismember(cnd, {'Baseline','ASV','ASV24','SAL','SAL24'}) & ~cellfun(@isempty,ani) & isfinite(dnm);
fils1 = fils1(keep); cnd = cnd(keep); ani = ani(keep); dnm = dnm(keep);

isExp = strcmp(cnd,'ASV');
uan = unique(ani);
sessPhase = repmat({'pre'}, numel(fils1),1);
for i = 1:numel(uan)
    m = strcmp(ani, uan{i});
    expDates = unique(dnm(m & isExp));
    if isempty(expDates), continue; end
    t0 = min(expDates);
    sessPhase(m & dnm >= t0) = {'post'};
end
fprintf('sessions: %d pre-ASV, %d post-ASV (session-level phase, by DATE not by epoch)\n', ...
    sum(strcmp(sessPhase,'pre')), sum(strcmp(sessPhase,'post')));

%% ---- first-hour epoch on each session, screen for length --------------
scfg = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',10,'NullMethod','none','DetectTransients',false, ...
         'CacheDir', fullfile(CACHE_ROOT,'screenFH') };
Sc = syncBatch(fils1, scfg{:});
ok = [Sc.ok]' & [Sc.dur]' >= LFIX;
fprintf('%d/%d sessions have a usable first hour\n', sum(ok), numel(ok));

uf = fils1(ok); ph = sessPhase(ok);
mapEarly = containers.Map(uf, repmat({[0 LFIX]}, numel(uf), 1));

cfgB = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',10,'MatchN',MATCHN,'NullMethod','shift','NSurr',50, ...
         'DetectTransients',true,'NSurrWindow',50,'RunAlpha',0.05, ...
         'Parallel', ~isempty(ver('parallel')), 'CheckpointEvery',5 };
Tfh = syncBatch(uf, cfgB{:}, 'Window', @(f) mapEarly(f), ...
    'CacheDir', fullfile(CACHE_ROOT,'firstHourTP'), ...
    'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_firstHourTP.mat'));

good = [Tfh.ok]';
Fpre  = Tfh(good & strcmp(ph,'pre'));
Fpost = Tfh(good & strcmp(ph,'post'));
aPre  = arrayfun(@(r) animalFromPath(r.file,'TransPlasticity'), Fpre,  'UniformOutput',false);
aPost = arrayfun(@(r) animalFromPath(r.file,'TransPlasticity'), Fpost, 'UniformOutput',false);
fprintf('\nfirst-hour PRE-ASV : %d sessions, %d animals\n', numel(Fpre),  numel(unique(aPre)));
fprintf('first-hour POST-ASV: %d sessions, %d animals\n', numel(Fpost), numel(unique(aPost)));

%% ---- metrics -------------------------------------------------------------
[rPre,oPre,dPre] = deal(nan(numel(Fpre),1));
for i=1:numel(Fpre), [rPre(i),oPre(i),dPre(i)] = transientStats(Fpre(i),STEP,WIDTH); end
mPre = arrayfun(@(r) meanPeakMagnitude(r.transients), Fpre);

[rPost,oPost,dPost] = deal(nan(numel(Fpost),1));
for i=1:numel(Fpost), [rPost(i),oPost(i),dPost(i)] = transientStats(Fpost(i),STEP,WIDTH); end
mPost = arrayfun(@(r) meanPeakMagnitude(r.transients), Fpost);

MET = { 'duty cycle',     oPre, oPost, oB
        'transient rate', rPre, rPost, rB
        'log duration',   log(dPre), log(dPost), log(dBvec)
        'log magnitude',  log(mPre), log(mPost), log(mB) };

fprintf('\n================ PRE-ASV (first hr) vs OptoRAM (first hr) ================\n');
for k = 1:size(MET,1)
    fprintf('\n-- %s --\n', MET{k,1});
    nestingCheck(Fpre, B, MET{k,2}, MET{k,4}, 'TransPlasticity','OptoRAM', [MET{k,1} ' pre']);
    try, lmeCompare(Fpre, B, MET{k,2}, MET{k,4}, 'TransPlasticity','OptoRAM', [MET{k,1} ' pre']); catch ME
        warning('LME failed: %s', ME.message); end
end

fprintf('\n================ POST-ASV (first hr) vs OptoRAM (first hr) ================\n');
for k = 1:size(MET,1)
    fprintf('\n-- %s --\n', MET{k,1});
    nestingCheck(Fpost, B, MET{k,3}, MET{k,4}, 'TransPlasticity','OptoRAM', [MET{k,1} ' post']);
    try, lmeCompare(Fpost, B, MET{k,3}, MET{k,4}, 'TransPlasticity','OptoRAM', [MET{k,1} ' post']); catch ME
        warning('LME failed: %s', ME.message); end
end

fprintf(['\nRead pre-vs-OptoRAM and post-vs-OptoRAM SIDE BY SIDE. Position is now\n' ...
         'matched on both comparisons, so:\n' ...
         '  both null           -> no dataset difference detectable at this position\n' ...
         '  pre null, post not  -> ASV exposure creates a divergence from OptoRAM\n' ...
         '                         that was not there before exposure\n' ...
         '  both differ equally -> a dataset difference unrelated to ASV\n' ...
         'The middle case is the one that supports a drug-specific claim; the\n' ...
         'other two do not.\n']);

save(fullfile(CACHE_ROOT,'firstHourPrePostResults.mat'), ...
    'Fpre','Fpost','B','aPre','aPost','bAnim','rPre','oPre','dPre','mPre', ...
    'rPost','oPost','dPost','mPost','rB','oB','dBvec','mB');
