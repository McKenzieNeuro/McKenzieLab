%% compareRatePrePostControl.m
% Firing-rate difference between Control, PreASV and PostASV, and whether
% the Control vs PreASV dS difference survives adjustment for rate.
% Groups, windows and sessions as in compareSyncPrePostControl.m (first-ASV
% day left out; other infusion days pre-infusion only).
%
% PER-SESSION RATE METRICS (spikes in the session's window / window length),
% over EVERY unit in spikes.times (no rate band, not the dS pair subset)
%   rateAll    median unit rate
%   rateLog    mean log10 unit rate (units with >= 1 spike in the window); equals the
%              mean log10 geometric-mean rate over all unit pairs
% Rates are tested on log10.
%
% dS ADJUSTMENT (LME, sessions nested in animals)
%   m0  dS ~ Group + (1|Animal)
%   m1  dS ~ Group + rateLog + (1|Animal)
%   m2  dS ~ Group + rateLog + logN + (1|Animal)      logN = log10(all units)
%   mI  dS ~ Group * rateLog + (1|Animal)     (does the rate slope differ by group?)
%   cs  m0 restricted to sessions inside the shared rateLog range (common support)
% dS itself is still the SDIST value (its own unit / pair selection).

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
MIN_EPOCH = 1800; NPERM = 1e4;
INCLUDE_NO_ASV = true;
MIN_UNITS = 0;               % e.g. 10 to drop small sessions

fa = getAllExtFiles(ROOT, 'mat', 1);         fa = fa(endsWith(fa, 'spikeDistanceProfile.mat'));
fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, 'spikeDistanceProfile_control.mat'));

rows = {}; UR = {};
for i = 1:numel(fc)
    L = load(fc{i}, 'SDIST'); X = L.SDIST;
    [y, ur] = local_session(X, [0 X.bin.edges(end)], fc{i});
    rows(end+1,:) = [{'Control', controlAnimal(fc{i}), fc{i}, X.N}, y]; UR{end+1,1} = ur; %#ok<AGROW>
end
for i = 1:numel(fa)
    L = load(fa{i}, 'SDIST'); X = L.SDIST; m = X.meta;
    noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
    if noASV
        if ~INCLUDE_NO_ASV, continue; end
        g = 'PreASV';
    elseif ~isfinite(m.dayFromFirstASV) || m.dayFromFirstASV == 0
        continue
    elseif m.dayFromFirstASV < 0
        g = 'PreASV';
    else
        g = 'PostASV';
    end
    w = [0, X.bin.edges(end)];
    if isfinite(m.infusionSec), w = [0, m.infusionSec]; end
    if diff(w) < MIN_EPOCH, continue; end
    [y, ur] = local_session(X, w, fa{i});
    rows(end+1,:) = [{g, m.animal, fa{i}, X.N}, y]; UR{end+1,1} = ur; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Group','Animal','Session','N','Hours','dS', ...
    'rateAll','rateLog','nAll'});
keep = T.nAll >= MIN_UNITS; T = T(keep, :); UR = UR(keep);
T.logN = log10(T.nAll);
G = {'Control','PreASV','PostASV'};
for g = G
    k = strcmp(T.Group, g{1});
    fprintf('%-8s %3d sessions, %2d animals, units %s\n', g{1}, sum(k), ...
        numel(unique(T.Animal(k))), mat2str(prctile(T.nAll(k), [0 50 100])));
end

%% ---- rate tests --------------------------------------------------------------------
fprintf('\nrates (session values; Hz shown as geometric mean across sessions, tests on log10)\n');
fprintf('%-9s | %-16s | %-16s | %-16s || %-26s || %-26s\n', 'metric', G{:}, 'perm: C-Pre / C-Post / Pre-Post', 'LME: same three');
rm = {'rateAll','rateLog'};
OUT = struct();
for mm = rm
    nm = mm{1}; isLog = true;
    t = T; if ~strcmp(nm, 'rateLog'), t.(nm) = log10(t.(nm)); end     % rateLog is already log10
    v = cellfun(@(g) t.(nm)(strcmp(t.Group, g) & isfinite(t.(nm))), G, 'UniformOutput', false);
    pr = [local_perm(v{1}, v{2}, NPERM), local_perm(v{1}, v{3}, NPERM), local_perm(v{2}, v{3}, NPERM)];
    pl = [local_lmeP(t, nm, 'Control', 'PreASV'), local_lmeP(t, nm, 'Control', 'PostASV'), local_lmeP(t, nm, 'PreASV', 'PostASV')];
    if isLog, s = cellfun(@(x) sprintf('%.2f Hz (n=%d)', 10^mean(x), numel(x)), v, 'UniformOutput', false);
    else,     s = cellfun(@(x) sprintf('%.2f (n=%d)', mean(x), numel(x)), v, 'UniformOutput', false); end
    fprintf('%-9s | %-16s | %-16s | %-16s || %.3f / %.3f / %.3f    || %.3f / %.3f / %.3f\n', nm, s{:}, pr, pl);
    OUT.(nm) = struct('v', {v}, 'pPerm', pr, 'pLME', pl);
end

%% ---- dS adjusted for rate ------------------------------------------------------------
fprintf('\ndS adjusted for rate (Group coefficient = second group minus first; dS units)\n');
fprintf('%-18s | %-28s | %-28s | %-28s\n', 'model', 'Control vs PreASV', 'Control vs PostASV', 'PreASV vs PostASV');
mods = {'m0  Group',               'y ~ Group + (1|Animal)';
        'm1  + rateLog',          'y ~ Group + rateLog + (1|Animal)';
        'm2  + rateLog + logN',   'y ~ Group + rateLog + logN + (1|Animal)'};
prs = {'Control','PreASV'; 'Control','PostASV'; 'PreASV','PostASV'};
ADJ = struct();
for q = 1:size(mods, 1)
    c = cell(1, 3);
    for j = 1:3
        [b, p] = local_lmeCoef(T, 'dS', prs{j,1}, prs{j,2}, mods{q,2}, 'Group');
        c{j} = sprintf('%+.4f  p=%.3f', b, p); ADJ.(sprintf('m%d', q-1))(j,:) = [b p];
    end
    fprintf('%-18s | %-28s | %-28s | %-28s\n', mods{q,1}, c{:});
end
% rate slope and group x rate interaction
c = cell(1, 3); c2 = cell(1, 3);
for j = 1:3
    [b, p] = local_lmeCoef(T, 'dS', prs{j,1}, prs{j,2}, 'y ~ Group + rateLog + (1|Animal)', 'rateLog');
    c{j} = sprintf('%+.4f  p=%.3f', b, p);
    [b, p] = local_lmeCoef(T, 'dS', prs{j,1}, prs{j,2}, 'y ~ Group * rateLog + (1|Animal)', ':');
    c2{j} = sprintf('%+.4f  p=%.3f', b, p);
end
fprintf('%-18s | %-28s | %-28s | %-28s\n', 'rateLog slope', c{:});
fprintf('%-18s | %-28s | %-28s | %-28s\n', 'Group x rate', c2{:});
% dS-rate slope within each group alone: y ~ rateLog + (1|Animal)
c = cell(1, 3);
for j = 1:3
    t = T(strcmp(T.Group, G{j}) & isfinite(T.rateLog) & isfinite(T.dS), :);
    c{j} = 'n/a';
    if height(t) < 6, continue; end
    t.y = t.dS; t.Animal = categorical(t.Animal);
    try
        if numel(unique(t.Animal)) == height(t), lme = fitlm(t, 'y ~ rateLog'); C = lme.Coefficients; b = C{'rateLog','Estimate'}; p = C{'rateLog','pValue'};
        else
            lme = fitlme(t, 'y ~ rateLog + (1|Animal)');
            [~, ~, C] = fixedEffects(lme, 'DFMethod', 'satterthwaite'); if isa(C, 'dataset'), C = dataset2table(C); end
            k = find(string(C.Name) == "rateLog"); b = C.Estimate(k); p = C.pValue(k);
        end
        c{j} = sprintf('%+.4f  p=%.3f', b, p);
    catch ME
        warning('cmp:lme', 'slope %s: %s', G{j}, ME.message);
    end
end
fprintf('%-18s | %-28s | %-28s | %-28s\n', 'slope in Control / PreASV / PostASV alone', c{:});
% common support: sessions within the overlap of the two groups' rateLog ranges
c = cell(1, 3);
for j = 1:3
    a = T.rateLog(strcmp(T.Group, prs{j,1})); bb = T.rateLog(strcmp(T.Group, prs{j,2}));
    lo = max(min(a), min(bb)); hi = min(max(a), max(bb));
    t = T(T.rateLog >= lo & T.rateLog <= hi, :);
    [b, p] = local_lmeCoef(t, 'dS', prs{j,1}, prs{j,2}, 'y ~ Group + (1|Animal)', 'Group');
    c{j} = sprintf('%+.4f  p=%.3f (n=%d/%d)', b, p, sum(strcmp(t.Group, prs{j,1})), sum(strcmp(t.Group, prs{j,2})));
end
fprintf('%-18s | %-28s | %-28s | %-28s\n', 'cs  common support', c{:});
fprintf('rateLog = session mean of log10 unit rate over all units.\n');
save(fullfile(CACHE_ROOT, 'compareRatePrePostControl.mat'), 'T', 'UR', 'OUT', 'ADJ');

%% ---- figure ---------------------------------------------------------------------------
col = {[0.4 0.4 0.4], [0 0 0], [0.55 0.1 0.6]};
fig = figure('Color','w','Position',[60 60 1400 380]);
tiledlayout(1, 3, 'TileSpacing','compact');
% unit rate CDF, every session weighted equally (band units)
nexttile; hold on
for j = 1:3
    idx = find(strcmp(T.Group, G{j})); x = []; wg = [];
    for q = idx(:)'
        u = UR{q}; u = u(isfinite(u) & u > 0);
        if isempty(u), continue; end
        x = [x; u(:)]; wg = [wg; repmat(1/numel(u), numel(u), 1)]; %#ok<AGROW>
    end
    [x, o] = sort(x); cc = cumsum(wg(o)); cc = cc / cc(end);
    plot(x, cc, '-', 'Color', col{j}, 'LineWidth', 2);
end
set(gca, 'XScale', 'log'); xlabel('unit rate in window (Hz), all units'); ylabel('cumulative fraction (sessions equal)');
legend(G, 'Location','southeast', 'Box','off');
% per-session pair rate
nexttile; hold on
for j = 1:3
    x = 10.^T.rateLog(strcmp(T.Group, G{j}) & isfinite(T.rateLog));
    plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 3);
    plot(j + [-0.25 0.25], 10^mean(log10(x))*[1 1], '-', 'Color', col{j}, 'LineWidth', 2);
end
set(gca, 'YScale', 'log'); xlim([0.5 3.5]); xticks(1:3); xticklabels(G); ylabel('unit rate (Hz, geometric mean)');
title('all units', 'FontWeight','normal');
% dS vs pair rate
nexttile; hold on
for j = 1:3
    k = strcmp(T.Group, G{j}) & isfinite(T.rateLog); x = T.rateLog(k); y = T.dS(k);
    plot(10.^x, y, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 4);
    if sum(k) >= 3
        b = polyfit(x, y, 1); xx = linspace(min(x), max(x), 2);
        plot(10.^xx, polyval(b, xx), '-', 'Color', col{j}, 'LineWidth', 1.5);
    end
end
set(gca, 'XScale', 'log'); yline(0, 'k:');
xlabel('unit rate (Hz, geometric mean, all units)'); ylabel('dS = Snull - S'); title('dS vs rate (lines: OLS per group)', 'FontWeight','normal');
exportgraphics(fig, fullfile(CACHE_ROOT, 'compareRatePrePostControl.pdf'), 'ContentType','vector');


%% =========================================================================
function [y, uRate] = local_session(X, w, sdFile)
% {hours, dS, rateAll, rateLog, nAll}; uRate = rates of all units in the window
ed = X.bin.edges(:); ctr = 0.5 * (ed(1:end-1) + ed(2:end));
in = ctr >= w(1) & ctr < w(2);
S = X.bin.S(in); Sn = double(X.bin.Snull(in, :));
hours = sum(in) * diff(ed(1:2)) / 3600;
dS = mean(mean(Sn, 1, 'omitnan')) - mean(S, 'omitnan');

% spikes file: next to the SDIST file (meta.file can be stale if folders moved)
d = dir(fullfile(fileparts(sdFile), '*spikes.cellinfo.mat'));
if ~isempty(d), sf = fullfile(d(1).folder, d(1).name);
elseif exist(X.meta.file, 'file'), sf = X.meta.file;
else
    warning('rate:nospikes', 'no spikes.cellinfo.mat next to %s', sdFile);
    uRate = []; y = {hours, dS, NaN, NaN, NaN}; return
end
if numel(d) > 1, warning('rate:multi', '%d spikes files in %s; using %s', numel(d), d(1).folder, d(1).name); end
L = load(sf, 'spikes');
ts = L.spikes.times(:);
ts = cellfun(@(x) unique(reshape(x(~isnan(x)), [], 1)), ts, 'UniformOutput', false);
T = diff(w);
r = cellfun(@(x) sum(x >= w(1) & x < w(2)), ts) / T;     % every unit, rate in window
act = r > 0;
uRate = r;
y = {hours, dS, median(r(act)), mean(log10(r(act))), numel(ts)};
end


function p = local_perm(x, y, nperm)
p = NaN;
if numel(x) >= 2 && numel(y) >= 2, p = permtest(x, y, nperm); end
end


function p = local_lmeP(T, nm, a, b)
[~, p] = local_lmeCoef(T, nm, a, b, 'y ~ Group + (1|Animal)', 'Group');
end


function [b, p] = local_lmeCoef(T, nm, a, bName, form, term)
% coefficient and p of the first fixed-effect row whose name contains TERM
% ('Group' matches the group dummy; ':' matches the interaction)
b = NaN; p = NaN;
t = T(ismember(T.Group, {a, bName}) & isfinite(T.(nm)), :);
if height(t) < 6 || numel(unique(t.Group)) < 2, return; end
t.y = t.(nm);
t.Group = categorical(t.Group, {a, bName}); t.Animal = categorical(t.Animal);
try
    lme = fitlme(t, form);
    [~, ~, C] = fixedEffects(lme, 'DFMethod', 'satterthwaite');   % same df as the anova p in compareSync*
    if isa(C, 'dataset'), C = dataset2table(C); end
    nmC = string(C.Name);
    if strcmp(term, ':'), k = find(contains(nmC, ':'), 1);
    elseif strcmp(term, 'Group'), k = find(startsWith(nmC, 'Group') & ~contains(nmC, ':'), 1);
    else, k = find(nmC == term, 1); end
    if isempty(k), return; end
    b = C.Estimate(k); p = C.pValue(k);
catch ME
    warning('cmp:lme', '%s: %s', nm, ME.message);
end
end