%% compareRateByCellType.m
% Firing rate by putative cell type and layer, from cell_metrics alone
% (cell_metrics.firingRate, .putativeCellType, .deepSuperficial).
%   Control   never-infused recordings under CONTROL_ROOT
%   PreASV    TransPlasticity days before the animal's first ASV (+ never-ASV animals)
%   PostASV   TransPlasticity days after the first ASV
% First-ASV day left out. firingRate is a whole-session rate, so it cannot
% be restricted to the pre-infusion part: INFUSION_DAYS = 'exclude' drops
% every other session with an Infusion_start marker; 'whole' keeps them
% (whole-session rate, includes the post-infusion period).
%
% STRATA  ALL, PYR ('Pyramidal'), INT ('Interneuron'), PYRdeep, PYRsup
% PER SESSION  rate_<s> = median firingRate over the stratum (>= MIN_STRAT units),
%              n_<s>, fracINT = INT/(PYR+INT), fracDeep = PYRdeep/(PYRdeep+PYRsup)
% TESTS  session level: permutation + LME y ~ Group + (1|Animal) (Satterthwaite)
%        unit level:    log10 rate ~ Group + (1|Animal) + (1|Session)

ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
INFUSION_DAYS = 'exclude';      % 'exclude' | 'whole'
INCLUDE_NO_ASV = true;
MIN_STRAT = 3; NPERM = 1e4;
CM = 'cell_metrics.cellinfo.mat';

%% ---- sessions ------------------------------------------------------------------------
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, CM));
cnd = cellfun(@(f) conditionFromPath(f, 'TransPlasticity'), fa, 'UniformOutput', false);
ani = cellfun(@(f) animalFromPath(f, 'TransPlasticity'),    fa, 'UniformOutput', false);
dnm = cellfun(@sessionDateFromPath, fa);
k = ismember(cnd, CONDS) & ~cellfun(@isempty, ani) & isfinite(dnm);
fa = fa(k); cnd = cnd(k); ani = ani(k); dnm = dnm(k);
inf_t = cellfun(@infusionTime, fa);          % looks for the .evt.sti in the file's folder
recStart = cellfun(@recordingStartFromPath, fa);
ref = firstASVReference(fa, ani, cnd, dnm, inf_t, recStart, 'ASV');

fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, CM));
fc = fc(~isfinite(cellfun(@infusionTime, fc)));

S = {}; U = {};
for i = 1:numel(fc)
    [s, u] = local_session(fc{i}, MIN_STRAT);
    if isempty(s), continue; end
    S(end+1,:) = [{'Control', controlAnimal(fc{i}), fc{i}}, s]; U{end+1,1} = u; %#ok<AGROW>
end
nFirst = 0; nInf = 0; nNo = 0;
for i = 1:numel(fa)
    [~, ~, dayOff, rf] = ts2FromReference(0, ani{i}, dnm(i), recStart(i), ref);
    if isempty(rf.session)
        if ~INCLUDE_NO_ASV, continue; end
        g = 'PreASV'; nNo = nNo + 1;
    elseif dayOff == 0
        nFirst = nFirst + 1; continue
    elseif dayOff < 0, g = 'PreASV';
    else,              g = 'PostASV';
    end
    if isfinite(inf_t(i)) && strcmp(INFUSION_DAYS, 'exclude'), nInf = nInf + 1; continue; end
    [s, u] = local_session(fa{i}, MIN_STRAT);
    if isempty(s), continue; end
    S(end+1,:) = [{g, ani{i}, fa{i}}, s]; U{end+1,1} = u; %#ok<AGROW>
end
ST = {'ALL','PYR','INT','PYRdeep','PYRsup'};
T = cell2table(S, 'VariableNames', [{'Group','Animal','Session'}, strcat('rate_', ST), strcat('n_', ST), {'fracINT','fracDeep'}]);
G = {'Control','PreASV','PostASV'};
fprintf('left out: %d first-ASV-day session(s), %d other infusion-day session(s) (INFUSION_DAYS = %s); %d never-ASV session(s) in PreASV\n', ...
    nFirst, nInf, INFUSION_DAYS, nNo);
for g = G
    k = strcmp(T.Group, g{1});
    fprintf('%-8s %3d sessions, %2d animals | median units/session: ALL %g, PYR %g, INT %g, PYRdeep %g, PYRsup %g\n', g{1}, sum(k), ...
        numel(unique(T.Animal(k))), median(T.n_ALL(k)), median(T.n_PYR(k)), median(T.n_INT(k)), median(T.n_PYRdeep(k)), median(T.n_PYRsup(k)));
end

%% ---- unit table for the unit-level LME ------------------------------------------------------
UT = table();
for q = 1:height(T)
    u = U{q}; u.Group = repmat(T.Group(q), height(u), 1); u.Animal = repmat(T.Animal(q), height(u), 1);
    u.Session = repmat(T.Session(q), height(u), 1); UT = [UT; u]; %#ok<AGROW>
end

%% ---- tests ------------------------------------------------------------------------------------
prs = {'Control','PreASV'; 'Control','PostASV'; 'PreASV','PostASV'};
fprintf('\nrate: geometric mean across sessions of the session median (tests on log10); composition: mean fraction\n');
fprintf('%-12s | %-15s | %-15s | %-15s || %-24s || %-24s || %-24s\n', 'metric', G{:}, ...
    'perm C-Pre/C-Post/Pre-Post', 'LME session level', 'LME unit level (b C-Pre)');
mets = [strcat('rate_', ST), {'fracINT','fracDeep'}];
OUT = struct();
for mm = mets
    nm = mm{1}; isRate = startsWith(nm, 'rate_');
    t = T; if isRate, t.(nm) = log10(t.(nm)); end
    v = cellfun(@(g) t.(nm)(strcmp(t.Group, g) & isfinite(t.(nm))), G, 'UniformOutput', false);
    pr = [local_perm(v{1}, v{2}, NPERM), local_perm(v{1}, v{3}, NPERM), local_perm(v{2}, v{3}, NPERM)];
    pl = nan(1, 3); pu = nan(1, 3); bu = nan(1, 3);
    for j = 1:3
        [~, pl(j)] = local_lme(t, nm, prs{j,1}, prs{j,2}, 'y ~ Group + (1|Animal)');
        if isRate
            sn = extractAfter(nm, 'rate_');
            ut = UT(UT.(sn) & UT.rate > 0, :); ut.lr = log10(ut.rate);
            [bu(j), pu(j)] = local_lme(ut, 'lr', prs{j,1}, prs{j,2}, 'y ~ Group + (1|Animal) + (1|Session)');
        end
    end
    if isRate, s = cellfun(@(x) sprintf('%.2f Hz (n=%d)', 10^mean(x), numel(x)), v, 'UniformOutput', false);
    else,      s = cellfun(@(x) sprintf('%.3f (n=%d)', mean(x), numel(x)), v, 'UniformOutput', false); end
    fu = '';
    if isRate, fu = sprintf('%.3f / %.3f / %.3f (x%.2f)', pu, 10^bu(1)); end
    fprintf('%-12s | %-15s | %-15s | %-15s || %.3f / %.3f / %.3f  || %.3f / %.3f / %.3f  || %s\n', nm, s{:}, pr, pl, fu);
    OUT.(nm) = struct('v', {v}, 'pPerm', pr, 'pLME', pl, 'pUnit', pu, 'bUnit', bu);
end
fprintf('unit-level x = PreASV / Control rate ratio from the unit LME.\n');
save(fullfile(CACHE_ROOT, 'compareRateByCellType.mat'), 'T', 'UT', 'OUT', 'INFUSION_DAYS');

%% ---- figure: session-weighted rate CDFs per stratum + INT fraction ------------------------------
col = {[0.4 0.4 0.4], [0 0 0], [0.55 0.1 0.6]};
fig = figure('Color','w','Position',[60 60 1600 340]);
tiledlayout(1, numel(ST) + 1, 'TileSpacing','compact');
for s = 1:numel(ST)
    nexttile; hold on
    for j = 1:3
        x = []; wg = [];
        for q = find(strcmp(T.Group, G{j}))'
            u = U{q}; r = u.rate(u.(ST{s}) & u.rate > 0);
            if numel(r) < MIN_STRAT, continue; end
            x = [x; r]; wg = [wg; repmat(1/numel(r), numel(r), 1)]; %#ok<AGROW>
        end
        if isempty(x), continue; end
        [x, o] = sort(x); c = cumsum(wg(o)); c = c / c(end);
        plot(x, c, '-', 'Color', col{j}, 'LineWidth', 2);
    end
    set(gca, 'XScale', 'log'); xlabel('firing rate (Hz)'); title(ST{s});
    if s == 1, ylabel('cumulative fraction (sessions equal)'); legend(G, 'Location','southeast', 'Box','off'); end
end
nexttile; hold on
for j = 1:3
    x = T.fracINT(strcmp(T.Group, G{j}) & isfinite(T.fracINT));
    plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 3);
    if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', col{j}, 'LineWidth', 2); end
end
xlim([0.5 3.5]); xticks(1:3); xticklabels(G); title('fraction INT');
exportgraphics(fig, fullfile(CACHE_ROOT, 'compareRateByCellType.pdf'), 'ContentType','vector');


%% =========================================================================
function [s, u] = local_session(f, minStrat)
% per-session stratum medians, counts, composition; u = per-unit table
s = {}; u = table();
L = load(f, 'cell_metrics'); cm = L.cell_metrics;
r = double(cm.firingRate(:));
n = numel(r);
ct = cellstr(string(cm.putativeCellType(:)));
ly = repmat({''}, n, 1);
if isfield(cm, 'deepSuperficial') && numel(cm.deepSuperficial) == n, ly = cellstr(string(cm.deepSuperficial(:))); end
if numel(ct) ~= n, warning('ct:len', '%s: firingRate %d vs putativeCellType %d', f, n, numel(ct)); return; end
isP = contains(ct, 'Pyramidal'); isI = contains(ct, 'Interneuron');
isD = strcmp(ly, 'Deep'); isS = strcmp(ly, 'Superficial');
u = table(r, true(n, 1), isP, isI, isP & isD, isP & isS, 'VariableNames', {'rate','ALL','PYR','INT','PYRdeep','PYRsup'});
st = {'ALL','PYR','INT','PYRdeep','PYRsup'};
rs = nan(1, 5); ns = zeros(1, 5);
for k = 1:5
    m = u.(st{k}) & isfinite(r) & r > 0; ns(k) = sum(m);
    if ns(k) >= minStrat, rs(k) = median(r(m)); end
end
fI = NaN; if ns(2) + ns(3) > 0, fI = ns(3) / (ns(2) + ns(3)); end
fD = NaN; if ns(4) + ns(5) > 0, fD = ns(4) / (ns(4) + ns(5)); end
s = [num2cell(rs), num2cell(ns), {fI, fD}];
end


function a = local_animal(f)
tok = regexp(f, '(CP\d+)', 'tokens', 'once');
if ~isempty(tok), a = tok{1}; return; end
[d, ~] = fileparts(f); [~, a] = fileparts(d);
a = regexprep(a, '_\d{6}(_\d{6})?$', '');
end


function p = local_perm(x, y, nperm)
p = NaN;
if numel(x) >= 2 && numel(y) >= 2, p = permtest(x, y, nperm); end
end


function [b, p] = local_lme(T, nm, a, bName, form)
% Group coefficient (bName minus a) and its Satterthwaite p
b = NaN; p = NaN;
t = T(ismember(T.Group, {a, bName}) & isfinite(T.(nm)), :);
if height(t) < 6 || numel(unique(t.Group)) < 2, return; end
t.y = t.(nm);
t.Group = categorical(t.Group, {a, bName}); t.Animal = categorical(t.Animal);
if ismember('Session', t.Properties.VariableNames), t.Session = categorical(t.Session); end
try
    lme = fitlme(t, form);
    [~, ~, C] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
    if isa(C, 'dataset'), C = dataset2table(C); end
    k = find(startsWith(string(C.Name), 'Group'), 1);
    if isempty(k), return; end
    b = C.Estimate(k); p = C.pValue(k);
catch ME
    warning('cmp:lme', '%s: %s', nm, ME.message);
end
end
