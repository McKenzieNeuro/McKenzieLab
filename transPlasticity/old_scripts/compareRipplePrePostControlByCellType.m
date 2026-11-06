%% compareRipplePrePostControlByCellType.m
% compareRipplePrePostControl split by putative cell type and pyramidal layer.
%   Control   never-infused recordings, whole recording
%   PreASV    days before the first ASV (+ never-ASV animals, INCLUDE_NO_ASV)
%   PostASV   days after the first ASV
% First-ASV day left out; other infusion days pre-infusion only. One value
% per session; tests: permutation and LME y ~ Group + (1|Animal), Satterthwaite.
%
% Labels: cell_metrics.putativeCellType / deepSuperficial in the session
% folder, matched to RP.units by spikes index (checked against cell_metrics.cluID; see matchCellMetrics).
% STRATA  PYR, INT, PYRdeep, PYRsup
%
% PER SESSION AND STRATUM s
%   frac_s     recruitment: mean over ripples of the fraction of s units with >= 1 spike
%   rRip_s     median over s units of in-ripple rate (spikes / total ripple time, Hz)
%   gain_s     median over s units of in-ripple rate / unit rate (RP.units.Rate,
%              whole recording: on pre-infusion-only windows the denominator
%              includes the rest of the session)
%   ppc_s      mean unit PPC over s units with >= MINSPK in-ripple spikes
%   pop_s      PPC of all in-ripple spikes of s pooled
%   ph_s       circular mean phase of pooled s spikes (printed, not tested; phase 0
%              is defined differently in the two datasets, so compare within TransPlasticity)
% Ripple-level (not by type): rate, logDur, freq.

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
MIN_EPOCH = 1800; MINSPK = 20; MIN_STRAT = 3; NPERM = 1e4;
INCLUDE_NO_ASV = true;
ST = {'PYR','INT','PYRdeep','PYRsup'};
G  = {'Control','PreASV','PostASV'};

fa = getAllExtFiles(ROOT, 'mat', 1);         fa = fa(endsWith(fa, 'rippleEvents.mat'));
fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, 'rippleEvents_control.mat'));

rows = {}; UP = {};
for i = 1:numel(fc)
    L = load(fc{i}, 'RIPPLES'); RP = L.RIPPLES;
    [y, up] = local_epoch(RP, [0 RP.meta.recordingSec], fc{i}, ST, MINSPK, MIN_STRAT);
    if isempty(y), continue; end
    rows(end+1,:) = [{'Control', controlAnimal(fc{i}), fc{i}}, y]; UP{end+1,1} = up; %#ok<AGROW>
end
nSkip = 0;
for i = 1:numel(fa)
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
    if noASV
        if ~INCLUDE_NO_ASV, continue; end
        g = 'PreASV';
    elseif ~isfinite(m.dayFromFirstASV) || m.dayFromFirstASV == 0
        nSkip = nSkip + 1; continue
    elseif m.dayFromFirstASV < 0, g = 'PreASV';
    else,                         g = 'PostASV';
    end
    w = [0, m.recordingSec];
    if isfinite(m.infusionSec), w = [0, m.infusionSec]; end
    if diff(w) < MIN_EPOCH, continue; end
    [y, up] = local_epoch(RP, w, fa{i}, ST, MINSPK, MIN_STRAT);
    if isempty(y), continue; end
    rows(end+1,:) = [{g, m.animal, fa{i}}, y]; UP{end+1,1} = up; %#ok<AGROW>
end
pm = {'frac','rRip','gain','ppc','pop','ph','n'};
vn = {'Group','Animal','Session','Event','Hours','logDur','freq','labelled'};
for s = ST, vn = [vn, strcat(pm, ['_' s{1}])]; end %#ok<AGROW>
T = cell2table(rows, 'VariableNames', vn);
T.rate = T.Event ./ T.Hours;
fprintf('first-ASV day left out: %d session(s)\n', nSkip);
bad = find(T.labelled < 0.95);
if ~isempty(bad), fprintf('sessions with < 95%% of units labelled:\n'); end
for q = bad(:)', fprintf('   %.0f%%  %s\n', 100*T.labelled(q), T.Session{q}); end
for g = G
    k = strcmp(T.Group, g{1});
    fprintf('%-8s %3d sessions, %2d animals, %.0f h | median units: PYR %g, INT %g, PYRdeep %g, PYRsup %g\n', g{1}, sum(k), ...
        numel(unique(T.Animal(k))), sum(T.Hours(k)), median(T.n_PYR(k)), median(T.n_INT(k)), median(T.n_PYRdeep(k)), median(T.n_PYRsup(k)));
end

%% ---- tests ------------------------------------------------------------------------------
mets = {'rate','logDur','freq'};
for s = ST, mets = [mets, strcat({'frac','rRip','gain','ppc','pop'}, ['_' s{1}])]; end %#ok<AGROW>
fprintf('\nsession means (mean +- SEM); rate, rRip, gain shown as geometric means and tested on log10\n');
fprintf('%-12s | %-20s | %-20s | %-20s || %-26s || %-26s\n', 'metric', G{:}, 'perm: C-Pre / C-Post / Pre-Post', 'LME: same three');
OUT = struct();
for mm = mets
    nm = mm{1}; isLog = any(startsWith(nm, {'rate','rRip','gain'}));
    t = T; if isLog, t.(nm) = log10(t.(nm)); t.(nm)(~isfinite(t.(nm))) = NaN; end
    v = cellfun(@(g) t.(nm)(strcmp(t.Group, g) & isfinite(t.(nm))), G, 'UniformOutput', false);
    pr = [local_perm(v{1}, v{2}, NPERM), local_perm(v{1}, v{3}, NPERM), local_perm(v{2}, v{3}, NPERM)];
    pl = [local_lme(t, nm, 'Control', 'PreASV'), local_lme(t, nm, 'Control', 'PostASV'), local_lme(t, nm, 'PreASV', 'PostASV')];
    if isLog, s = cellfun(@(x) sprintf('%.3g (n=%d)', 10^mean(x), numel(x)), v, 'UniformOutput', false);
    else,     s = cellfun(@local_ms, v, 'UniformOutput', false); end
    fprintf('%-12s | %-20s | %-20s | %-20s || %.3f / %.3f / %.3f    || %.3f / %.3f / %.3f\n', nm, s{:}, pr, pl);
    OUT.(nm) = struct('v', {v}, 'pPerm', pr, 'pLME', pl);
end
fprintf('\npooled preferred phase (deg, circular mean over sessions; within-dataset comparison only)\n');
for s = ST
    nm = ['ph_' s{1}];
    c = cellfun(@(g) rad2deg(angle(mean(exp(1i*T.(nm)(strcmp(T.Group, g) & isfinite(T.(nm))))))), G);
    fprintf('%-10s  %s\n', s{1}, strjoin(arrayfun(@(j) sprintf('%s %+6.1f', G{j}, c(j)), 1:3, 'UniformOutput', false), ' | '));
end
save(fullfile(CACHE_ROOT, 'compareRipplePrePostControlByCellType.mat'), 'T', 'UP', 'OUT');

%% ---- figure: recruitment per session, and session-weighted unit-PPC CDFs, per stratum ----------
col = {[0.4 0.4 0.4], [0 0 0], [0.55 0.1 0.6]};
fig = figure('Color','w','Position',[60 60 1400 640]);
tl = tiledlayout(2, numel(ST), 'TileSpacing','compact');
for s = 1:numel(ST)
    nexttile(s); hold on
    for j = 1:3
        x = OUT.(['frac_' ST{s}]).v{j};
        plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 3);
        if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', col{j}, 'LineWidth', 2); end
    end
    xlim([0.5 3.5]); xticks(1:3); xticklabels(G); title(['recruitment ' ST{s}]);
    if s == 1, ylabel('fraction of units active per ripple'); end
    ax = nexttile(numel(ST) + s); hold(ax, 'on');
    for j = 1:3
        x = []; wg = [];
        for q = find(strcmp(T.Group, G{j}))'
            u = UP{q}.ppc(UP{q}.(ST{s}) & isfinite(UP{q}.ppc));
            if numel(u) < MIN_STRAT, continue; end
            x = [x; u]; wg = [wg; repmat(1/numel(u), numel(u), 1)]; %#ok<AGROW>
        end
        if isempty(x), plot(ax, NaN, NaN, '-', 'Color', col{j}); continue; end
        [x, o] = sort(x); c = cumsum(wg(o)); c = c / c(end);
        plot(ax, x, c, '-', 'Color', col{j}, 'LineWidth', 2);
    end
    xline(ax, 0, 'k:'); xlabel(ax, ['unit PPC, ' ST{s}]);
    if s == 1, ylabel(ax, 'cumulative fraction (sessions equal)'); legend(ax, G, 'Location','southeast', 'Box','off'); end
end
title(tl, 'Ripples by cell type: controls vs TransPlasticity before / after first ASV (infusion day excluded)');
exportgraphics(fig, fullfile(CACHE_ROOT, 'compareRipplePrePostControlByCellType.pdf'), 'ContentType','vector');


%% =========================================================================
function [y, up] = local_epoch(RP, w, ripFile, ST, minSpk, minStrat)
y = {}; up = table();
% ---- labels for RP.units ----
fld = fileparts(ripFile);
dc = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(dc), warning('rc:cm', 'no cell_metrics in %s', fld); return; end
C = load(fullfile(dc(1).folder, dc(1).name), 'cell_metrics'); cm = C.cell_metrics;
nC = numel(cm.putativeCellType);
cmType = cellstr(string(cm.putativeCellType(:)));
cmLay = repmat({''}, nC, 1);
if isfield(cm, 'deepSuperficial') && numel(cm.deepSuperficial) == nC, cmLay = cellstr(string(cm.deepSuperficial(:))); end
nU = height(RP.units);
uid = []; if ismember('UID', RP.units.Properties.VariableNames), uid = RP.units.UID; end
nSp = [];
sf = dir(fullfile(fld, '*spikes.cellinfo.mat'));
if ~isempty(sf)
    Ls = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); nSp = numel(Ls.spikes.times);
end
loc = matchCellMetrics(cm, RP.units.SpikesIdx, uid, nSp);
ok = loc > 0;
ct = repmat({''}, nU, 1); ly = ct;
ct(ok) = cmType(loc(ok)); ly(ok) = cmLay(loc(ok));
isP = contains(ct, 'Pyramidal'); isI = contains(ct, 'Interneuron');
isD = strcmp(ly, 'Deep'); isS = strcmp(ly, 'Superficial');
M = [isP, isI, isP & isD, isP & isS];

% ---- ripples and in-ripple spikes in the window ----
in = find(RP.ts1 >= w(1) & RP.ts1Off <= w(2));
ev = RP.ev(in, :);
nR = numel(in);
if nR == 0, return; end
sp = RP.spikes(ismember(RP.spikes.Ripple, in) & isfinite(RP.spikes.Phase), :);
spAll = RP.spikes(ismember(RP.spikes.Ripple, in), :);       % counts include spikes without phase
ripT = sum(RP.ts1Off(in) - RP.ts1(in));
[~, rIdx] = ismember(spAll.Ripple, in);
A = sparse(rIdx, spAll.Unit, 1, nR, nU) > 0;                  % ripple x unit active
nSpk = accumarray(spAll.Unit, 1, [nU 1]);
rRip = nSpk / ripT;
gain = rRip ./ RP.units.Rate(:);
th = sp.Phase; u = sp.Unit;
n = accumarray(u, 1, [nU 1]);
Cs = accumarray(u, cos(th), [nU 1]); Sn = accumarray(u, sin(th), [nU 1]);
ppc = (Cs.^2 + Sn.^2 - n) ./ (n .* (n - 1)); ppc(n < minSpk) = NaN;
up = table(ppc, rRip, gain, M(:,1), M(:,2), M(:,3), M(:,4), 'VariableNames', [{'ppc','rRip','gain'}, ST]);

y = {nR, diff(w) / 3600, mean(log(ev.Dur), 'omitnan'), mean(ev.Freq, 'omitnan'), mean(ok)};
for k = 1:numel(ST)
    mk = M(:, k); nk = sum(mk);
    v = nan(1, 6);
    if nk >= minStrat
        v(1) = mean(full(sum(A(:, mk), 2)) / nk);
        v(2) = median(rRip(mk));
        g = gain(mk); v(3) = median(g(isfinite(g)));
        v(4) = mean(ppc(mk), 'omitnan');
        nP = sum(n(mk)); CP = sum(Cs(mk)); SP = sum(Sn(mk));
        if nP >= 2, v(5) = (CP^2 + SP^2 - nP) / (nP * (nP - 1)); v(6) = atan2(SP, CP); end
    end
    y = [y, num2cell(v), {nk}]; %#ok<AGROW>
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


function s = local_ms(v)
if isempty(v), s = 'n/a'; return; end
s = sprintf('%.3g +- %.2g (n=%d)', mean(v), std(v) / sqrt(numel(v)), numel(v));
end
