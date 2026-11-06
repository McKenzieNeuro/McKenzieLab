%% compareRateWindowsByCellType.m
% Firing rate by putative cell type / layer in four groups, rates counted
% from spike times inside each window (cell_metrics supplies the labels):
%   Control    never-infused recordings, whole recording
%   PreASV     any time before the animal's first ASV infusion (ts2 < 0),
%              including SAL days; never-ASV animals added (INCLUDE_NO_ASV)
%   ASVshort   SHORT_WIN after the first ASV infusion (default 2-8 h)
%   ASVlong    > LONG_AFTER after the first ASV infusion (default 24 h)
%   PreSAL     before the first SAL infusion (and before the first ASV)
%   SALshort   SHORT_WIN after the first SAL infusion, before the first ASV
%   SALlong    >= LONG_AFTER after the first SAL infusion, before the first ASV
% SAL groups only for animals whose first SAL came before their first ASV; on
% later infusion days only the pre-infusion part counts (LATER_INFUSIONS).
% PreASV still contains the SAL windows (it is all time before the first ASV).
%   SALaPre    SAL session after the first ASV: its own pre-infusion part
%   SALaShort  the same session, 2-8 h after its infusion (that session's own infusion,
%              not ts3, which counts from the FIRST SAL only). These windows also
%              lie inside ASVlong, which they overlap by construction.
% DiD row 3: (SALaShort - SALaPre) - (SALshort - PreSAL) = does the saline
% response change after ASV exposure (persistent ASV effect, or habituation).
% ts2 (s from first ASV) comes from ts2FromReference: exact when both folder
% names carry HHMMSS ('clock'), otherwise date-based (+-24 h; within the
% first-ASV session it is exact either way).
% LATER_INFUSIONS = 'preonly': on any infusion day after the first ASV, only
% the time before that day's infusion counts toward ASVlong (no re-exposure).
%
% One replicate = one session's portion of a window. Tests as in
% compareRippleWindowsByCellType.m: permutation + LME (1|Animal) for the
% between-group contrasts, paired sign-flip for the within-animal changes,
% and DiD = sign-flip on per-animal (ASV change from PreASV) - (SAL change from PreSAL).

ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
SHORT_WIN  = [2 8] * 3600;
LONG_AFTER = 24 * 3600;
LATER_INFUSIONS = 'preonly';    % 'preonly' | 'all'
INCLUDE_NO_ASV = true;
MIN_EPOCH = 1800;               % s of window a session must contribute
MIN_STRAT = 3; NPERM = 1e4;
ST = {'ALL','PYR','INT','PYRdeep','PYRsup'};
G  = {'Control','PreASV','ASVshort','ASVlong','PreSAL','SALshort','SALlong','SALaPre','SALaShort'};

%% ---- sessions ---------------------------------------------------------------------------
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'spikes.cellinfo.mat'));
cnd = cellfun(@(f) conditionFromPath(f, 'TransPlasticity'), fa, 'UniformOutput', false);
ani = cellfun(@(f) animalFromPath(f, 'TransPlasticity'),    fa, 'UniformOutput', false);
dnm = cellfun(@sessionDateFromPath, fa);
k = ismember(cnd, CONDS) & ~cellfun(@isempty, ani) & isfinite(dnm);
fa = fa(k); cnd = cnd(k); ani = ani(k); dnm = dnm(k);
inf_t = cellfun(@infusionTime, fa);
recStart = cellfun(@recordingStartFromPath, fa);
ref = firstASVReference(fa, ani, cnd, dnm, inf_t, recStart, 'ASV');
refS = firstASVReference(fa, ani, cnd, dnm, inf_t, recStart, 'SAL');

fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, 'spikes.cellinfo.mat'));
fc = fc(~isfinite(cellfun(@infusionTime, fc)));

S = {}; U = {};
for i = 1:numel(fc)
    [s, u] = local_window(fc{i}, @(t, dur) t >= 0 & t <= dur, MIN_EPOCH, MIN_STRAT);
    if isempty(s), continue; end
    S(end+1,:) = [{'Control', controlAnimal(fc{i}), fc{i}}, s]; U{end+1,1} = u; %#ok<AGROW>
end
modes = {};
for i = 1:numel(fa)
    [~, mode, ~, rf] = ts2FromReference(0, ani{i}, dnm(i), recStart(i), ref);
    isFirst = ~isempty(rf.session) && strcmp(fileparts(fa{i}), rf.session);
    if isempty(rf.session)                                  % never-ASV animal
        if ~INCLUDE_NO_ASV, continue; end
        W = {'PreASV', @(t, dur) t >= 0 & t <= dur};
    else
        modes{end+1} = mode; %#ok<AGROW>
        f2 = @(t) ts2FromReference(t, ani{i}, dnm(i), recStart(i), ref);
        later = isfinite(inf_t(i)) && ~isFirst && strcmp(LATER_INFUSIONS, 'preonly');
        cut = Inf; if later, cut = inf_t(i); end
        W = {'PreASV',   @(t, dur) f2(t) < 0;
             'ASVshort', @(t, dur) f2(t) >= SHORT_WIN(1) & f2(t) < SHORT_WIN(2) & t < cut;
             'ASVlong',  @(t, dur) f2(t) >= LONG_AFTER & t < cut};
    end
    % ---- saline windows (first SAL before first ASV only) ----
    [o3, ~, ~, rfS] = ts2FromReference(0, ani{i}, dnm(i), recStart(i), refS);
    if ~isempty(rfS.session) && isfinite(o3)
        o2 = NaN; if ~isempty(rf.session), o2 = ts2FromReference(0, ani{i}, dnm(i), recStart(i), ref); end
        if ~isfinite(o2) || o2 < o3                      % SAL infusion at t = -o3 has ts2 = o2 - o3 < 0
            if isfinite(o2), preA = @(t) t + o2 < 0; else, preA = @(t) true(size(t)); end
            isFirstS = strcmp(fileparts(fa{i}), rfS.session);
            cut3 = Inf;
            if isfinite(inf_t(i)) && ~isFirstS && strcmp(LATER_INFUSIONS, 'preonly'), cut3 = inf_t(i); end
            W = [W; {'PreSAL',   @(t, dur) t + o3 < 0 & preA(t);
                     'SALshort', @(t, dur) t + o3 >= SHORT_WIN(1) & t + o3 < SHORT_WIN(2) & t < cut3 & preA(t);
                     'SALlong',  @(t, dur) t + o3 >= LONG_AFTER & t < cut3 & preA(t)}]; %#ok<AGROW>
        end
    end
    % ---- SAL sessions after the first ASV: own pre-infusion part, and 2-8 h after its infusion ----
    if strcmp(cnd{i}, 'SAL') && isfinite(inf_t(i)) && ~isempty(rf.session)
        o2s = ts2FromReference(0, ani{i}, dnm(i), recStart(i), ref);
        if isfinite(o2s) && inf_t(i) + o2s > 0
            ti = inf_t(i);
            W = [W; {'SALaPre',   @(t, dur) t < ti;
                     'SALaShort', @(t, dur) t - ti >= SHORT_WIN(1) & t - ti < SHORT_WIN(2)}]; %#ok<AGROW>
        end
    end
    for w = 1:size(W, 1)
        [s, u] = local_window(fa{i}, W{w,2}, MIN_EPOCH, MIN_STRAT);
        if isempty(s), continue; end
        S(end+1,:) = [{W{w,1}, ani{i}, fa{i}}, s]; U{end+1,1} = u; %#ok<AGROW>
    end
end
T = cell2table(S, 'VariableNames', [{'Group','Animal','Session','Hours'}, strcat('rate_', ST), strcat('n_', ST), {'fracINT','fracDeep'}]);
fprintf('ts2 mode over ASV-animal sessions: %d clock, %d date\n', sum(strcmp(modes, 'clock')), sum(strcmp(modes, 'date')));
for g = G
    k = strcmp(T.Group, g{1});
    fprintf('%-9s %3d session-windows, %2d animals, %5.0f h | median units: PYR %g, INT %g, PYRdeep %g, PYRsup %g\n', g{1}, sum(k), ...
        numel(unique(T.Animal(k))), sum(T.Hours(k)), median(T.n_PYR(k)), median(T.n_INT(k)), median(T.n_PYRdeep(k)), median(T.n_PYRsup(k)));
end

%% ---- tests ------------------------------------------------------------------------------------
prs = {'Control','PreASV'; 'Control','ASVshort'; 'Control','ASVlong'; 'PreASV','ASVshort'; 'PreASV','ASVlong';
       'PreSAL','SALshort'; 'PreSAL','SALlong'; 'SALshort','ASVshort'; 'SALlong','ASVlong'; 'SALshort','SALaShort'};
pn  = {'C-Pre','C-Ashort','C-Along','Pre-Ashort','Pre-Along','PreS-Sshort','PreS-Slong','Sshort-Ashort','Slong-Along','Sshort-SAshort'};
ppr = {'PreASV','ASVshort'; 'PreASV','ASVlong'; 'PreSAL','SALshort'; 'PreSAL','SALlong'; 'SALaPre','SALaShort'};
ppn = {'Pre-Ashort','Pre-Along','PreS-Sshort','PreS-Slong','SAPre-SAshort'};
dd  = {'ASVshort','PreASV','SALshort','PreSAL'; 'ASVlong','PreASV','SALlong','PreSAL'; 'SALaShort','SALaPre','SALshort','PreSAL'};
ddn = {'short (A-Pre)-(S-PreS)','long (A-Pre)-(S-PreS)','SAL response after - before ASV'};
fprintf(['perm / LME: session-windows, unpaired.  paired: sign-flip on animal means.\n' ...
         'DiD: sign-flip on per-animal (ASV change from PreASV) - (SAL change from PreSAL).\n']);
fprintf('\nrate: geometric mean across session-windows of the session median (tests on log10); composition: mean\n');
fprintf(['%-12s' repmat(' | %-14s', 1, numel(G)) '\n'], 'metric', G{:});
mets = [strcat('rate_', ST), {'fracINT','fracDeep'}];
OUT = struct();
for mm = mets
    nm = mm{1}; isRate = startsWith(nm, 'rate_');
    t = T; if isRate, t.(nm) = log10(t.(nm)); end
    v = cellfun(@(g) t.(nm)(strcmp(t.Group, g) & isfinite(t.(nm))), G, 'UniformOutput', false);
    if isRate, s = cellfun(@(x) sprintf('%.2f Hz (%d)', 10^mean(x), numel(x)), v, 'UniformOutput', false);
    else,      s = cellfun(@(x) sprintf('%.3f (%d)', mean(x), numel(x)), v, 'UniformOutput', false); end
    fprintf(['%-12s' repmat(' | %-14s', 1, numel(G)) '\n'], nm, s{:});
    pr = nan(1, size(prs, 1)); pl = pr;
    for j = 1:size(prs, 1)
        pr(j) = local_perm(v{strcmp(G, prs{j,1})}, v{strcmp(G, prs{j,2})}, NPERM);
        pl(j) = local_lme(t, nm, prs{j,1}, prs{j,2});
    end
    ps = nan(1, size(ppr, 1));
    for j = 1:size(ppr, 1), ps(j) = local_paired(t, nm, ppr{j,1}, ppr{j,2}); end
    pd = nan(size(dd, 1), 2);
    for j = 1:size(dd, 1)
        q = local_did(t, nm, dd{j,:}); if numel(q) == 2, pd(j,:) = q; end
    end
    fprintf('%12s   perm   %s\n', '', strjoin(arrayfun(@(j) sprintf('%s %.3f', pn{j}, pr(j)), 1:numel(pr), 'UniformOutput', false), ' | '));
    fprintf('%12s   LME    %s\n', '', strjoin(arrayfun(@(j) sprintf('%s %.3f', pn{j}, pl(j)), 1:numel(pl), 'UniformOutput', false), ' | '));
    fprintf('%12s   paired %s\n', '', strjoin(arrayfun(@(j) sprintf('%s %.3f', ppn{j}, ps(j)), 1:numel(ps), 'UniformOutput', false), ' | '));
    fprintf('%12s   DiD    %s\n', '', strjoin(arrayfun(@(j) sprintf('%s %.3f (%d animals)', ddn{j}, pd(j,1), pd(j,2)), 1:size(dd,1), 'UniformOutput', false), ' | '));
    OUT.(nm) = struct('v', {v}, 'pPerm', pr, 'pLME', pl, 'pPaired', ps, 'pDiD', pd);
end
save(fullfile(CACHE_ROOT, 'compareRateWindowsByCellType.mat'), 'T', 'U', 'OUT', 'SHORT_WIN', 'LONG_AFTER', 'LATER_INFUSIONS');

%% ---- figure ----------------------------------------------------------------------------------
col = {[0.4 0.4 0.4], [0 0 0], [0.85 0.2 0.1], [0.55 0.1 0.6], [0.45 0.6 0.75], [0.2 0.6 1], [0 0.25 0.6], [0.6 0.8 0.5], [0.2 0.55 0.2]};
fig = figure('Color','w','Position',[60 60 1600 340]);
tiledlayout(1, numel(ST) + 1, 'TileSpacing','compact');
for s = 1:numel(ST)
    nexttile; hold on
    for j = 1:numel(G)
        x = []; wg = [];
        for q = find(strcmp(T.Group, G{j}))'
            u = U{q}; r = u.rate(u.(ST{s}) & u.rate > 0);
            if numel(r) < MIN_STRAT, continue; end
            x = [x; r]; wg = [wg; repmat(1/numel(r), numel(r), 1)]; %#ok<AGROW>
        end
        if isempty(x), plot(NaN, NaN, '-', 'Color', col{j}); continue; end
        [x, o] = sort(x); c = cumsum(wg(o)); c = c / c(end);
        plot(x, c, '-', 'Color', col{j}, 'LineWidth', 2);
    end
    set(gca, 'XScale', 'log'); xlabel('firing rate in window (Hz)'); title(ST{s});
    if s == 1, ylabel('cumulative fraction (session-windows equal)'); legend(G, 'Location','southeast', 'Box','off'); end
end
nexttile; hold on
for j = 1:numel(G)
    x = T.fracINT(strcmp(T.Group, G{j}) & isfinite(T.fracINT));
    plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 3);
    if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', col{j}, 'LineWidth', 2); end
end
xlim([0.5 numel(G)+0.5]); xticks(1:numel(G)); xticklabels(G); xtickangle(30); title('fraction INT');
exportgraphics(fig, fullfile(CACHE_ROOT, 'compareRateWindowsByCellType.pdf'), 'ContentType','vector');


%% =========================================================================
function [s, u] = local_window(spkFile, inWin, minEpoch, minStrat)
% rates of every unit inside the window defined by inWin(t, dur) (t = s from recording start)
s = {}; u = table();
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
% window length on a 1 s grid (inWin may need ts2 evaluation, so sample it)
g = (0.5:1:dur)';
winLen = sum(inWin(g, dur));
if winLen < minEpoch, return; end
r = cellfun(@(x) sum(inWin(x, dur)), ts) / winLen;

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
u = table(r, true(nU, 1), isP, isI, isP & isD, isP & isS, 'VariableNames', {'rate','ALL','PYR','INT','PYRdeep','PYRsup'});
st = u.Properties.VariableNames(2:end);
rs = nan(1, 5); ns = zeros(1, 5);
for k = 1:5
    m = u.(st{k}) & r > 0; ns(k) = sum(m);
    if ns(k) >= minStrat, rs(k) = median(r(m)); end
end
fI = NaN; if ns(2) + ns(3) > 0, fI = ns(3) / (ns(2) + ns(3)); end
fD = NaN; if ns(4) + ns(5) > 0, fD = ns(4) / (ns(4) + ns(5)); end
s = [{winLen / 3600}, num2cell(rs), num2cell(ns), {fI, fD}];
end


function p = local_paired(t, nm, a, b)
% sign-flip on per-animal (b - a), animals with both
p = NaN;
A = t(strcmp(t.Group, a) & isfinite(t.(nm)), :); B = t(strcmp(t.Group, b) & isfinite(t.(nm)), :);
ids = intersect(unique(A.Animal), unique(B.Animal));
if numel(ids) < 2, return; end
d = arrayfun(@(k) mean(B.(nm)(strcmp(B.Animal, ids{k}))) - mean(A.(nm)(strcmp(A.Animal, ids{k}))), 1:numel(ids))';
n = numel(d); obs = abs(mean(d));
if n <= 14
    Sg = 2*(dec2bin(0:2^n-1, n) == '1') - 1; p = mean(abs(Sg*d/n) >= obs - 1e-12);
else
    nul = arrayfun(@(k) mean((2*(rand(n,1) > 0.5) - 1) .* d), 1:1e4); p = (1 + sum(abs(nul) >= obs)) / (1 + 1e4);
end
end


function p = local_perm(x, y, nperm)
p = NaN;
if numel(x) >= 2 && numel(y) >= 2, p = permtest(x, y, nperm); end
end


function p = local_lme(T, nm, a, b)
p = NaN;
t = T(ismember(T.Group, {a, b}) & isfinite(T.(nm)), :);
if height(t) < 6 || numel(unique(t.Group)) < 2, return; end
t.y = t.(nm); t.Group = categorical(t.Group, {a, b}); t.Animal = categorical(t.Animal);
try
    lme = fitlme(t, 'y ~ Group + (1|Animal)');
    [~, ~, C] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
    if isa(C, 'dataset'), C = dataset2table(C); end
    k = find(startsWith(string(C.Name), 'Group'), 1);
    if ~isempty(k), p = C.pValue(k); end
catch ME
    warning('cmp:lme', '%s: %s', nm, ME.message);
end
end


function s = local_ab(g)
s = strrep(strrep(strrep(strrep(g, 'Control', 'C'), 'PreASV', 'Pre'), 'ASVshort', 'Short'), 'ASVlong', 'Long');
end


function a = local_animal(f)
tok = regexp(f, '(CP\d+)', 'tokens', 'once');
if ~isempty(tok), a = tok{1}; return; end
[d, ~] = fileparts(f); [~, a] = fileparts(d);
a = regexprep(a, '_\d{6}(_\d{6})?$', '');
end


function p = local_did(t, nm, a1, a0, b1, b0)
% sign-flip on per-animal (a1 - a0) - (b1 - b0), animals with all four windows
p = NaN;
mu = @(g, id) mean(t.(nm)(strcmp(t.Group, g) & strcmp(t.Animal, id) & isfinite(t.(nm))));
ids = unique(t.Animal);
d = nan(numel(ids), 1);
for k = 1:numel(ids)
    d(k) = (mu(a1, ids{k}) - mu(a0, ids{k})) - (mu(b1, ids{k}) - mu(b0, ids{k}));
end
d = d(isfinite(d)); n = numel(d);
if n < 2, return; end
obs = abs(mean(d));
if n <= 14
    Sg = 2*(dec2bin(0:2^n-1, n) == '1') - 1; p = mean(abs(Sg*d/n) >= obs - 1e-12);
else
    nul = arrayfun(@(k) mean((2*(rand(n,1) > 0.5) - 1) .* d), 1:1e4); p = (1 + sum(abs(nul) >= obs)) / (1 + 1e4);
end
p = [p, n];
end
