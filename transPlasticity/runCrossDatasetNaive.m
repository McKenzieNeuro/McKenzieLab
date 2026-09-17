%% runCrossDatasetNaive.m
% TransPlasticity vs OptoRAM using DRUG-NAIVE epochs only.
%
% The earlier cross-dataset comparison pooled TransPlasticity sessions
% regardless of ASV exposure, so any lasting drug effect sat inside the
% TransPlasticity side of a between-dataset contrast. Restricting to
% pre-first-exposure epochs removes that: both sides are then drug-free and
% exposure-naive, and a surviving difference is a property of the
% preparations rather than of ASV.
%
% MATCHING, WHICH IS THE WHOLE DIFFICULTY
%  1. Epoch length. Cluster significance depends on how many windows an
%     epoch holds, so both sides use the same LFIX.
%  2. Stimulation. OptoRAM epochs are taken from a genuinely pulse-FREE gap
%     rather than masking pulse windows after the fact -- masking would leave
%     OptoRAM epochs with fewer valid windows than TransPlasticity's and bias
%     detection between groups.
%  3. Unit count. Same MatchN on both sides, as before.
%
% WHAT STILL DOES NOT MATCH, and should be said in the text: the
% TransPlasticity epoch is the hour immediately BEFORE infusion, so it sits
% some hours into a recording, while the OptoRAM epoch is the first clean
% hour available. Position within the session is therefore not matched. The
% script prints both distributions so the size of that mismatch is visible.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
WIDTH = 60; STEP = 5; MATCHN = 15; LFIX = 3600;

%% ---- TransPlasticity: reuse the drug-naive epochs already computed ------
P = load(fullfile(CACHE_ROOT,'persistenceResults.mat'));
E = P.T(P.ok);
isPre = strcmp(P.phase,'pre');
A = E(isPre);
aAnim = P.ea(isPre);
fprintf('TransPlasticity drug-naive epochs: %d from %d animals\n', ...
    numel(A), numel(unique(aAnim)));

%% ---- OptoRAM: matched clean epochs --------------------------------------
fils2 = getAllExtFiles('R:\WSun\ePhys\OptoRAM','mat',1);
fils2 = fils2(contains(fils2,'spikes.cellinfo.mat'));

scfg = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',10,'NullMethod','none','DetectTransients',false, ...
         'CacheDir', fullfile(CACHE_ROOT,'screenOR') };
S2 = syncBatch(fils2, scfg{:});
dur2 = [S2.dur]';

win2 = cell(numel(fils2),1); okWin = false(numel(fils2),1);
for i = 1:numel(fils2)
    if ~S2(i).ok, continue; end
    ex = optoPulseIntervals(fils2{i});
    w = findCleanWindow(ex, dur2(i), LFIX, 'first');
    if ~isempty(w), win2{i} = w; okWin(i) = true; end
end
fprintf('OptoRAM: %d/%d sessions have a %.0f-min pulse-free window\n', ...
    sum(okWin), numel(fils2), LFIX/60);

uf2 = fils2(okWin);
mapW = containers.Map(uf2, win2(okWin));

cfgB = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',10,'MatchN',MATCHN,'NullMethod','shift','NSurr',50, ...
         'DetectTransients',true,'NSurrWindow',50,'RunAlpha',0.05, ...
         'Parallel', ~isempty(ver('parallel')), 'CheckpointEvery',5 };

T2 = syncBatch(uf2, cfgB{:}, 'Window', @(f) mapW(f), ...
    'CacheDir', fullfile(CACHE_ROOT,'naiveOR'), ...
    'CheckpointFile', fullfile(CACHE_ROOT,'ckpt_naiveOR.mat'));
B = T2([T2.ok]);
bAnim = arrayfun(@(r) animalFromPath(r.file,'OptoRAM'), B, 'UniformOutput', false);
fprintf('OptoRAM usable epochs: %d from %d animals\n', numel(B), numel(unique(bAnim)));

%% ---- position-in-session mismatch ---------------------------------------
posA = arrayfun(@(r) r.windowAbs(1), A);
posB = arrayfun(@(r) r.windowAbs(1), B);
fprintf('\nepoch start (s into recording): TransPlasticity median %.0f [%.0f-%.0f]\n', ...
    median(posA), min(posA), max(posA));
fprintf('                                OptoRAM         median %.0f [%.0f-%.0f]\n', ...
    median(posB), min(posB), max(posB));
fprintf('If these differ a lot, position-in-session is a live alternative\n');
fprintf('explanation for any difference below.\n');

%% ---- metrics + tests ----------------------------------------------------
[rA,oA,dA] = deal(nan(numel(A),1));
for i = 1:numel(A), [rA(i),oA(i),dA(i)] = transientStats(A(i), STEP, WIDTH); end
mA = arrayfun(@(r) meanPeakMagnitude(r.transients), A);
[rB,oB,dB] = deal(nan(numel(B),1));
for i = 1:numel(B), [rB(i),oB(i),dB(i)] = transientStats(B(i), STEP, WIDTH); end
mB = arrayfun(@(r) meanPeakMagnitude(r.transients), B);

MET = { 'duty cycle',      oA,        oB,        'none'
        'transient rate',  rA,        rB,        'none'
        'log duration',    log(dA),   log(dB),   'log'
        'log magnitude',   log(mA),   log(mB),   'log' };

for k = 1:size(MET,1)
    fprintf('\n########## %s ##########\n', MET{k,1});
    nestingCheck(A, B, MET{k,2}, MET{k,3}, 'TransPlasticity','OptoRAM', MET{k,1});
    try
        lme = lmeCompare(A, B, MET{k,2}, MET{k,3}, 'TransPlasticity','OptoRAM', MET{k,1});
    catch ME
        warning('LME failed for %s: %s', MET{k,1}, ME.message);
    end
end

save(fullfile(CACHE_ROOT,'crossNaiveResults.mat'), ...
    'A','B','aAnim','bAnim','rA','oA','dA','mA','rB','oB','dB','mB','posA','posB');
