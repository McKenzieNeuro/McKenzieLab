%% compareInfusionOrder.m
% Is the missing post-infusion rate increase after ASV an ORDER effect?
% Every infused session (SAL, ASV) is scored on its OWN infusion: rate before that session's
% infusion (t < infusion) and 2-8 h after it, so no cross-day reference and no PreASV
% contamination. LR = log10(post / pre) of the session median rate (PYR, INT, PYRdeep).
% Each session is labelled by what the animal had received before it (chronological by date):
%   SAL1 / ASV1          no earlier infusion
%   SALafterSAL          only earlier SAL           ASVafterSAL   only earlier SAL
%   SALafterASV          any earlier ASV            ASVafterASV   any earlier ASV
%   SALrep = SALafterSAL + SALafterASV (saline given as a repeat);  ASVrep = ASVafterSAL + ASVafterASV
% Tests, per sleep state (ALL / WAKE / NREM):
%   within group: sign-flip on per-animal mean LR (exact);
%   between groups: permutation over sessions, LME LR ~ group + (1|Animal), and paired
%       sign-flip on animals with both groups;
%   order-adjusted drug effect: LME LR ~ drug + ordinal + (1|Animal), and ~ drug + repeat + (1|Animal),
%       where ordinal = position of the infusion among that animal's infusions.
% CP17_251015 (SAL1) and CP17_251020 (SALrep) are both analysed.
% Output: compareInfusionOrder.mat (R = per-session table, printed in full).
ROOT        = 'R:\TransPlasticity';
CACHE_ROOT  = 'C:\Users\samckenzie\syncCache';
if exist('SHARED', 'var'), CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root; end
SHORT_WIN = [2 8] * 3600; MIN_EPOCH = 1800; MIN_STRAT = 3; MIN_EPOCH_STATE = 300;
STATES = {'ALL','WAKE','NREM'}; MET = {'PYR','INT','PYRdeep','PYRmatched'}; MCOL = [3 4 5];   % columns of S{q}: 1 hours, 2.. = rate_ALL, PYR, INT, PYRdeep, PYRsup
% PYRmatched (ALL only): median over pyramidal units firing in BOTH windows of log10(post/pre) of
% each unit's rate, i.e. a within-unit change, insensitive to units that are silent in one window.
% A session is FLAGGED when its pre-infusion PYR median < LOW_PRE Hz, fewer than MIN_MATCH units fire in
% both windows, or the pre-vs-post rank correlation of unit rates < MIN_RHO (unit identity / stability).
LOW_PRE = 0.5; MIN_MATCH = 5; MIN_RHO = 0.3;
GRP = {'SAL1','SALrep','SALafterSAL','SALafterASV','ASV1','ASVrep','ASVafterSAL','ASVafterASV','SAL','ASV'};

fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'spikes.cellinfo.mat'));
cnd = cellfun(@(f) conditionFromPath(f, 'TransPlasticity'), fa, 'UniformOutput', false);
ani = cellfun(@(f) animalFromPath(f, 'TransPlasticity'),    fa, 'UniformOutput', false);
dnm = cellfun(@sessionDateFromPath, fa);
inf_t = cellfun(@infusionTime, fa);
k = ismember(cnd, {'SAL','ASV'}) & ~cellfun(@isempty, ani) & isfinite(dnm) & isfinite(inf_t);
fa = fa(k); cnd = cnd(k); ani = ani(k); dnm = dnm(k); inf_t = inf_t(k);
n = numel(fa);

%% ---- order labels --------------------------------------------------------------------------
ord = zeros(n, 1); prior = repmat({''}, n, 1); grp = repmat({''}, n, 1);
for i = 1:n
    e = find(strcmp(ani, ani{i}) & dnm < dnm(i));           % earlier infusion sessions of this animal
    ord(i) = numel(e) + 1;
    if isempty(e), prior{i} = 'none'; elseif any(strcmp(cnd(e), 'ASV')), prior{i} = 'ASV'; else, prior{i} = 'SAL'; end
    switch prior{i}
        case 'none', grp{i} = [cnd{i} '1'];
        otherwise,   grp{i} = [cnd{i} 'after' prior{i}];
    end
end


%% ---- per-session pre / post rates -----------------------------------------------------------
pre = nan(n, numel(STATES), numel(MET)); post = pre; hPre = nan(n, numel(STATES)); hPost = hPre;
nPyr = zeros(n, 1); silPre = zeros(n, 1); silPost = zeros(n, 1); nMatch = zeros(n, 1); rhoPP = nan(n, 1); lrU = nan(n, 1);
for i = 1:n
    ti = inf_t(i);
    [sP, uP] = local_window(fa{i}, @(t, dur) t >= 0 & t < ti, MIN_EPOCH, MIN_STRAT, STATES, MIN_EPOCH_STATE);
    [sQ, uQ] = local_window(fa{i}, @(t, dur) t - ti >= SHORT_WIN(1) & t - ti < SHORT_WIN(2), MIN_EPOCH, MIN_STRAT, STATES, MIN_EPOCH_STATE);
    for q = 1:numel(STATES)
        if ~isempty(sP{q}), hPre(i,q) = sP{q}{1}; for m = 1:numel(MCOL), pre(i,q,m) = sP{q}{MCOL(m)}; end, end
        if ~isempty(sQ{q}), hPost(i,q) = sQ{q}{1}; for m = 1:numel(MCOL), post(i,q,m) = sQ{q}{MCOL(m)}; end, end
    end
    if ~isempty(uP) && ~isempty(uQ) && height(uP) == height(uQ)
        kp = uP.PYR; nPyr(i) = sum(kp); silPre(i) = sum(kp & uP.rate == 0); silPost(i) = sum(kp & uQ.rate == 0);
        km = kp & uP.rate > 0 & uQ.rate > 0; nMatch(i) = sum(km);
        if sum(kp) >= 3, rhoPP(i) = corr(uP.rate(kp), uQ.rate(kp), 'Type', 'Spearman'); end
        if nMatch(i) >= 3
            lrU(i) = median(log10(uQ.rate(km) ./ uP.rate(km)));
            pre(i,1,4) = median(uP.rate(km)); post(i,1,4) = median(uQ.rate(km));
        end
    end
end
LR = log10(post ./ pre); LR(:,1,4) = lrU;
flag = pre(:,1,1) < LOW_PRE | nMatch < MIN_MATCH | rhoPP < MIN_RHO;

[~, sn] = cellfun(@(f) fileparts(fileparts(f)), fa, 'UniformOutput', false);
R = table(ani, cnd, dnm, ord, prior, grp, inf_t / 3600, hPre(:,1), hPost(:,1), pre(:,1,1), post(:,1,1), LR(:,1,1), ...
    nPyr, silPre, silPost, nMatch, rhoPP, lrU, flag, sn, ...
    'VariableNames', {'Animal','Cond','Date','Ord','Prior','Group','InfH','PreH','PostH','PreHz','PostHz','LR', ...
    'nPyr','silentPre','silentPost','nMatched','rhoUnits','LRmatched','Flag','Session'});
[~, o] = sortrows(table(ani, dnm)); R = R(o, :);
fprintf('\n%d infused sessions, %d animals.  PYR rate before the session''s own infusion and 2-8 h after (ALL states)\n', n, numel(unique(ani)));
fprintf('PYR units: n, silent (0 spikes) before / after, firing in both; rho = Spearman of unit rates before vs after;\n');
fprintf('ratio = session median after / before; matched = median within-unit ratio; ** = flagged (pre < %g Hz, < %d matched, or rho < %g)\n', LOW_PRE, MIN_MATCH, MIN_RHO);
fprintf('%-6s %-20s %-4s %-3s %-4s %-12s %5s %5s %6s %6s %5s %5s %5s %5s %6s %6s %6s\n', 'Animal','Session','Cond','Ord','Prev','Group','PreH','PostH','PreHz','PostHz','nPYR','sil0','sil1','both','rho','ratio','match');
for r = 1:height(R)
    fl = '  '; if R.Flag(r), fl = '**'; end
    fprintf('%-6s %-20s %-4s %-3d %-4s %-12s %5.2f %5.2f %6.2f %6.2f %5d %5d %5d %5d %6.2f %6.2f %6.2f %s\n', R.Animal{r}, R.Session{r}, R.Cond{r}, R.Ord(r), R.Prior{r}, R.Group{r}, ...
        R.PreH(r), R.PostH(r), R.PreHz(r), R.PostHz(r), R.nPyr(r), R.silentPre(r), R.silentPost(r), R.nMatched(r), R.rhoUnits(r), 10^R.LR(r), 10^R.LRmatched(r), fl);
end
fprintf('%d of %d sessions flagged\n', sum(R.Flag), height(R));

%% ---- group tables and tests ----------------------------------------------------------------
member = @(g) strcmp(grp, g) | (strcmp(g, 'SALrep') & startsWith(grp, 'SALafter')) | (strcmp(g, 'ASVrep') & startsWith(grp, 'ASVafter')) ...
    | (strcmp(g, 'SAL') & strcmp(cnd, 'SAL')) | (strcmp(g, 'ASV') & strcmp(cnd, 'ASV'));
prs = {'SAL1','SALrep'; 'SAL1','ASV'; 'SALrep','ASV'; 'SAL1','ASVrep'; 'SALrep','ASVrep'; 'SALafterSAL','ASVafterSAL'};
for pass = 1:2
if pass == 2, fprintf('\n\n################ ROBUSTNESS: flagged sessions removed (ALL state) ################\n'); end
for q = 1:numel(STATES)
    if pass == 2 && q > 1, continue; end
    for m = 1:numel(MET)
        x = LR(:,q,m); pr = pre(:,q,m); po = post(:,q,m);
        if pass == 2, x(flag) = NaN; end
        if ~any(isfinite(x)), continue; end
        fprintf('\n================ %s, %s: log10(post/pre), 2-8 h after the session''s own infusion ================\n', STATES{q}, MET{m});
        fprintf('mean +/- SEM over animals (animal = mean of its sessions); log10 ratio = mean of log10(post/pre); ratio = 10^mean\n');
        fprintf('%-12s %9s %15s %15s %16s %7s %12s\n', 'group', 'sess(an)', 'pre Hz', 'post Hz', 'log10 ratio', 'ratio', 'p within');
        for g = GRP
            kk = member(g{1}) & isfinite(x);
            if ~any(kk), fprintf('%-12s %9s\n', g{1}, '0'); continue; end
            [d, ~] = local_animalMeans(x(kk), ani(kk));
            dP = local_animalMeans(pr(kk), ani(kk)); dQ = local_animalMeans(po(kk), ani(kk));
            fprintf('%-12s %4d (%2d) %6.2f+/-%5.2f %6.2f+/-%5.2f %+7.3f+/-%5.3f %7.2f %12.3f\n', g{1}, sum(kk), numel(d), mean(dP), local_sem(dP), mean(dQ), local_sem(dQ), mean(d), local_sem(d), 10^mean(d), local_flip1(d));
        end
        fprintf('between groups (ratio B / ratio A):  perm | LME | paired (animals with both)\n');
        for j = 1:size(prs, 1)
            a = member(prs{j,1}) & isfinite(x); b = member(prs{j,2}) & isfinite(x);
            if sum(a) < 2 || sum(b) < 2, fprintf('  %-12s vs %-12s  too few\n', prs{j,1}, prs{j,2}); continue; end
            pp = permtest(x(a), x(b), 1e4);
            pl = local_lme2(x, ani, a, b);
            [pd, na] = local_paired2(x, ani, a, b);
            dA = local_animalMeans(x(a), ani(a)); dB = local_animalMeans(x(b), ani(b));
            fprintf('  %-12s vs %-12s  log10 ratio %+.3f+/-%.3f / %+.3f+/-%.3f, B-A %+.3f+/-%.3f   perm %.3f | LME %.3f | paired %.3f (%d animals)\n', prs{j,1}, prs{j,2}, mean(dA), local_sem(dA), mean(dB), local_sem(dB), mean(dB) - mean(dA), hypot(local_sem(dA), local_sem(dB)), pp, pl, pd, na);
        end
        kk = isfinite(x);
        T = table(x(kk), categorical(cnd(kk)), ord(kk), double(ord(kk) > 1), categorical(ani(kk)), 'VariableNames', {'LR','Drug','Ord','Rep','Animal'});
        for f = {'Ord', 'Rep'}
            try
                M = fitlme(T, ['LR ~ Drug + ' f{1} ' + (1|Animal)']);
                [~, ~, C] = fixedEffects(M, 'DFMethod', 'satterthwaite'); if isa(C, 'dataset'), C = dataset2table(C); end
                C.Name = cellstr(string(C.Name)); kd = find(strcmp(C.Name, 'Drug_SAL')); ko = find(strcmp(C.Name, f{1}));
                if isempty(kd), kd = find(startsWith(C.Name, 'Drug_')); end
                sg = 1; if strcmp(C.Name{kd}, 'Drug_SAL'), sg = -1; end      % report ASV - SAL
                fprintf('  order-adjusted [%s]: ASV - SAL = %+.3f+/-%.3f log10 (x%.2f), p %.3f | %s coef %+.3f+/-%.3f, p %.3f | n=%d sessions\n', f{1}, sg*C.Estimate(kd), C.SE(kd), 10^(sg*C.Estimate(kd)), C.pValue(kd), f{1}, C.Estimate(ko), C.SE(ko), C.pValue(ko), height(T));
            catch ME
                fprintf('  order-adjusted [%s]: %s\n', f{1}, ME.message);
            end
        end
    end
end
end
save(fullfile(CACHE_ROOT, 'compareInfusionOrder.mat'), 'R', 'pre', 'post', 'LR', 'ani', 'cnd', 'grp', 'ord', 'flag', 'STATES', 'MET');

%% ---- helpers ---------------------------------------------------------------------------------
function [S, u] = local_window(spkFile, inWin, minEpoch, minStrat, states, minEpochState)
useCM = false; %#ok<NASGU>
% rates of every unit inside the window defined by inWin(t, dur) (t = s from recording start), and
% inside (window AND sleep state) for each entry of states. S{q} = one table row for state q, or {}
% when the window (q = 'ALL': minEpoch) or the state inside it (minEpochState) is too short.
% u = per-unit table for the unrestricted window (used by the figure).
% useCM: whole-recording window, so with no SleepState file the state rates may be taken from
% cell_metrics.firingRate_<state>state (CellExplorer); those rows get Hours = NaN.
S = cell(1, numel(states)); u = table(); hasSS = false;
fld = fileparts(spkFile);
dc = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(dc), warning('rw:cm', 'no cell_metrics in %s', fld); return; end
L = load(spkFile, 'spikes'); sp = L.spikes;
C = load(fullfile(dc(1).folder, dc(1).name), 'cell_metrics'); cm = C.cell_metrics;
ts = sp.times(:);
ts = cellfun(@(x) reshape(x(~isnan(x)), [], 1), ts, 'UniformOutput', false);
nU = numel(ts);
dur = max(cellfun(@(x) max([x; -Inf]), ts));
if ~isfinite(dur), return; end
% window on a 1 s grid (inWin may need ts2 evaluation, so evaluate it once, for grid and spikes)
g = (0.5:1:dur)';
inG = inWin(g, dur);
if sum(inG) < minEpoch, return; end
inX = cellfun(@(x) inWin(x, dur), ts, 'UniformOutput', false);
% sleep states
ints = struct();
sf = dir(fullfile(fld, '*.SleepState.states.mat'));
if ~isempty(sf)
    Ls = load(fullfile(sf(1).folder, sf(1).name), 'SleepState');
    if isfield(Ls, 'SleepState') && isfield(Ls.SleepState, 'ints'), ints = Ls.SleepState.ints; hasSS = true; end
end

% labels
nC = numel(cm.putativeCellType);
cmType = cellstr(string(cm.putativeCellType(:)));
cmLay = repmat({''}, nC, 1);
if isfield(cm, 'deepSuperficial') && numel(cm.deepSuperficial) == nC, cmLay = cellstr(string(cm.deepSuperficial(:))); end
clu = []; if isfield(sp, 'UID') && numel(sp.UID) == nU, clu = sp.UID(:); end
if isfield(sp, 'cluID') && numel(sp.cluID) == nU, clu = sp.cluID(:); end
loc = matchCellMetrics(cm, (1:nU)', clu, nU);
if ~any(loc), warning('rw:match', '%s: cannot match %d spike units to %d cell_metrics units', fld, nU, nC); end
ct = repmat({''}, nU, 1); ly = ct; ok = loc > 0;
ct(ok) = cmType(loc(ok)); ly(ok) = cmLay(loc(ok));
isP = contains(ct, 'Pyramidal'); isI = contains(ct, 'Interneuron');
isD = strcmp(ly, 'Deep'); isS = strcmp(ly, 'Superficial');

for q = 1:numel(states)
    if strcmp(states{q}, 'ALL')
        mG = inG; mX = inX; need = minEpoch;
    else
        fn = [states{q} 'state'];
        cmf = ['firingRate_' fn];
        if isfield(ints, fn) && ~isempty(ints.(fn))
            iv = ints.(fn);
            mG = inG & local_inIv(g, iv);
            mX = cellfun(@(x, w) w & local_inIv(x, iv), ts, inX, 'UniformOutput', false);
            need = minEpochState;
        elseif useCM && isfield(cm, cmf) && numel(cm.(cmf)) == nC
            v = double(cm.(cmf)(:)); r = nan(nU, 1); r(ok) = v(loc(ok));
            hasSS = true; mG = []; mX = {}; need = -Inf;
        else
            continue
        end
    end
    if isempty(mX)                           % CellExplorer state rate: no duration available
        winLen = NaN;
    else
        winLen = sum(mG);
        if winLen < need, continue; end
        r = cellfun(@sum, mX) / winLen;
    end
    uq = table(r, true(nU, 1), isP, isI, isP & isD, isP & isS, 'VariableNames', {'rate','ALL','PYR','INT','PYRdeep','PYRsup'});
    st = uq.Properties.VariableNames(2:end);
    rs = nan(1, 5); ns = zeros(1, 5);
    for k = 1:5
        m = uq.(st{k}) & r > 0; ns(k) = sum(m);
        if ns(k) >= minStrat, rs(k) = median(r(m)); end
    end
    fI = NaN; if ns(2) + ns(3) > 0, fI = ns(3) / (ns(2) + ns(3)); end
    fD = NaN; if ns(4) + ns(5) > 0, fD = ns(4) / (ns(4) + ns(5)); end
    S{q} = [{winLen / 3600}, num2cell(rs), num2cell(ns), {fI, fD}];
    if strcmp(states{q}, 'ALL'), u = uq; end
end
end


function k = local_inIv(t, iv)
% t inside any [start end) row of iv
k = false(size(t));
if isempty(iv), return; end
iv = sortrows(iv(all(isfinite(iv), 2) & iv(:,2) > iv(:,1), :));
if isempty(iv), return; end
M = iv(1,:);                                  % merge overlapping / touching intervals
for q = 2:size(iv, 1)
    if iv(q,1) <= M(end,2), M(end,2) = max(M(end,2), iv(q,2)); else, M(end+1,:) = iv(q,:); end %#ok<AGROW>
end
b = discretize(t, reshape(M', [], 1));         % odd bin = inside an interval, even = gap
k = mod(b, 2) == 1;
end


function [d, ids] = local_animalMeans(x, a)
% per-animal mean of x
ids = unique(a); d = arrayfun(@(k) mean(x(strcmp(a, ids{k}))), 1:numel(ids))';
end


function e = local_sem(d)
% standard error of the mean (NaN for fewer than 2 values)
d = d(isfinite(d)); e = NaN; if numel(d) > 1, e = std(d) / sqrt(numel(d)); end
end


function p = local_flip1(d)
% exact sign-flip test of mean(d) = 0 over animals
p = NaN; n = numel(d);
if n < 2, return; end
obs = abs(mean(d));
if n <= 14
    Sg = 2*(dec2bin(0:2^n-1, n) == '1') - 1; p = mean(abs(Sg*d/n) >= obs - 1e-12);
else
    nul = arrayfun(@(k) mean((2*(rand(n,1) > 0.5) - 1) .* d), 1:1e4); p = (1 + sum(abs(nul) >= obs)) / (1 + 1e4);
end
end


function [p, nA] = local_paired2(x, ani, a, b)
% sign-flip on per-animal (mean b - mean a) over animals with both
p = NaN; nA = 0;
ia = unique(ani(a)); ib = unique(ani(b)); ids = intersect(ia, ib);
nA = numel(ids);
if nA < 2, return; end
d = arrayfun(@(k) mean(x(b & strcmp(ani, ids{k}))) - mean(x(a & strcmp(ani, ids{k}))), 1:nA)';
p = local_flip1(d);
end


function p = local_lme2(x, ani, a, b)
p = NaN;
k = a | b;
T = table(x(k), double(b(k)), categorical(ani(k)), 'VariableNames', {'y','g','Animal'});
try
    M = fitlme(T, 'y ~ g + (1|Animal)');
    [~, ~, C] = fixedEffects(M, 'DFMethod', 'satterthwaite'); if isa(C, 'dataset'), C = dataset2table(C); end
    p = C.pValue(strcmp(cellstr(string(C.Name)), 'g'));
catch
end
end
