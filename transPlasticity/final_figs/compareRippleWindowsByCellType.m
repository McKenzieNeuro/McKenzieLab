%% compareRippleWindowsByCellType.m
% Ripple properties by putative cell type / pyramidal layer in four groups:
%   Control    never-infused recordings, whole recording
%   PreASV     any time before the animal's first ASV infusion (ts2 < 0),
%              including SAL days; never-ASV animals added (INCLUDE_NO_ASV)
%   ASVshort   SHORT_WIN after the first ASV infusion (default 2-8 h)
%   ASVlong    >= LONG_AFTER after the first ASV infusion (default 24 h)
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
% LATER_INFUSIONS = 'preonly': on infusion days after the first ASV only the
% time before that day's infusion counts toward ASVlong.
% Windows are intervals in ts1 from the session's ts2 - ts1 offset
% (RIPPLES.ts2, exact within a session; across days exact in 'clock' mode,
% +-24 h in 'date' mode). Same windows as compareRateWindowsByCellType.m.
%
% Labels: cell_metrics in the session folder, matched with matchCellMetrics.
% Metrics per session-window and stratum (PYR, INT, PYRdeep, PYRsup):
%   frac_ (recruitment), xfrac_ (recruitment minus the Poisson chance of >= 1 spike
%   in a ripple of that duration at the unit's out-of-ripple rate in the same window),
%   gainW_ (in-ripple rate / out-of-ripple rate in the same window),
%   rRip_ (in-ripple rate), gain_ (in-ripple / unit
%   rate; denominator RIPPLES.units.Rate = whole recording), ppc_ (mean unit
%   PPC, units >= MINSPK in-ripple spikes), pop_ (pooled PPC), ph_ (pooled
%   phase, printed only). Ripple level: rate, logDur, freq.
% One replicate = one session's portion of a window. Tests: permutation and
% LME y ~ Group + (1|Animal) (Satterthwaite) for C-Pre, C-ASVshort, C-ASVlong,
% Pre-ASVshort, Pre-ASVlong, PreSAL-SALshort, PreSAL-SALlong, SALshort-ASVshort,
% SALlong-ASVlong; paired sign-flip on animal means for the four within-animal
% changes; DiD = sign-flip on per-animal (ASV change) - (SAL change).

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
if exist('SHARED', 'var')     % set by runAllAnalyses.m
    CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root; CONTROL_ROOT = SHARED.controlRoot;
end
SHORT_WIN  = [2 8] * 3600;
LONG_AFTER = 24 * 3600;
LATER_INFUSIONS = 'preonly';
INCLUDE_NO_ASV = true;
MIN_EPOCH = 1800; MINSPK = 20; MIN_STRAT = 3; NPERM = 1e4;
ST = {'PYR','INT','PYRdeep','PYRsup'};
G  = {'Control','PreASV','ASVshort','ASVlong','PreSAL','SALshort','SALlong','SALaPre','SALaShort','BaseAll'};
% BaseAll = cohort baseline: EVERY session, before that session's own infusion (whole recording
% when it has none). Overlaps the other groups; used only for Control vs iPKR (C-Base).

fa = getAllExtFiles(ROOT, 'mat', 1);         fa = fa(endsWith(fa, 'rippleEvents.mat'));
fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, 'rippleEvents_control.mat'));

rows = {}; UP = {};
for i = 1:numel(fc)
    L = load(fc{i}, 'RIPPLES'); RP = L.RIPPLES;
    [y, up] = local_epoch(RP, [0 RP.meta.recordingSec], fc{i}, ST, MINSPK, MIN_STRAT);
    if isempty(y), continue; end
    rows(end+1,:) = [{'Control', controlAnimal(fc{i}), fc{i}}, y]; UP{end+1,1} = up; %#ok<AGROW>
end
for i = 1:numel(fa)
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    rec = m.recordingSec;
    noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
    off = RP.ts2(:) - RP.ts1(:); off = median(off(isfinite(off)));     % ts2 = ts1 + off
    off3 = NaN;
    if isfield(RP, 'ts3'), off3 = RP.ts3(:) - RP.ts1(:); off3 = median(off3(isfinite(off3))); end
    if isempty(off3), off3 = NaN; end
    if noASV
        if ~INCLUDE_NO_ASV, continue; end
        W = {'PreASV', [0 rec]}; pre = [0 rec];
    elseif isempty(off) || ~isfinite(off)
        W = cell(0, 2); pre = [0 rec];                   % no ts2: drug windows unknown, BaseAll still counts
    else
        pre = [0, min(rec, -off)];
        isFirst = strcmp(m.folder, m.firstASV.session);
        cut = rec;
        if isfinite(m.infusionSec) && ~isFirst && strcmp(LATER_INFUSIONS, 'preonly'), cut = m.infusionSec; end
        W = {'PreASV',   [0, min(rec, -off)];
             'ASVshort', [max(0, SHORT_WIN(1) - off), min(cut, SHORT_WIN(2) - off)];
             'ASVlong',  [max(0, LONG_AFTER - off), cut]};
    end
    % ---- saline windows (first SAL before first ASV only) ----
    hasSAL = isfield(m, 'firstSAL') && ~isempty(m.firstSAL.session) && isfinite(off3);
    if hasSAL && (noASV || off - off3 < 0)               % SAL infusion (t = -off3) has ts2 = off - off3
        isFirstS = strcmp(m.folder, m.firstSAL.session);
        cut3 = rec;
        if isfinite(m.infusionSec) && ~isFirstS && strcmp(LATER_INFUSIONS, 'preonly'), cut3 = m.infusionSec; end
        andPre = @(a) [max([a(1), pre(1), 0]), min([a(2), pre(2), rec])];
        W = [W; {'PreSAL',   andPre([0, -off3]);
                 'SALshort', andPre([SHORT_WIN(1) - off3, min(cut3, SHORT_WIN(2) - off3)]);
                 'SALlong',  andPre([LONG_AFTER - off3, cut3])}]; %#ok<AGROW>
    end
    % ---- SAL sessions after the first ASV: own pre-infusion part, and 2-8 h after its infusion ----
    if strcmp(m.condition, 'SAL') && isfinite(m.infusionSec) && ~noASV && isfinite(off) && m.infusionSec + off > 0
        ti = m.infusionSec;
        W = [W; {'SALaPre',   [0, ti];
                 'SALaShort', [ti + SHORT_WIN(1), min(rec, ti + SHORT_WIN(2))]}]; %#ok<AGROW>
    end
    tiB = rec; if isfield(m, 'infusionSec') && isfinite(m.infusionSec), tiB = min(rec, m.infusionSec); end
    W = [W; {'BaseAll', [0, tiB]}]; %#ok<AGROW>
    for k = 1:size(W, 1)
        w = W{k,2};
        if diff(w) < MIN_EPOCH, continue; end
        [y, up] = local_epoch(RP, w, fa{i}, ST, MINSPK, MIN_STRAT);
        if isempty(y), continue; end
        rows(end+1,:) = [{W{k,1}, m.animal, fa{i}}, y]; UP{end+1,1} = up; %#ok<AGROW>
    end
end
pm = {'frac','rRip','gain','ppc','pop','ph','xfrac','gainW','n'};
vn = {'Group','Animal','Session','Event','Hours','logDur','freq','labelled'};
for s = ST, vn = [vn, strcat(pm, ['_' s{1}])]; end %#ok<AGROW>
T = cell2table(rows, 'VariableNames', vn);
T.rate = T.Event ./ T.Hours;
bad = find(T.labelled < 0.95);
if ~isempty(bad), fprintf('session-windows with < 95%% of units labelled:\n'); end
for q = bad(:)', fprintf('   %.0f%%  %s\n', 100*T.labelled(q), T.Session{q}); end
for g = G
    k = strcmp(T.Group, g{1});
    fprintf('%-9s %3d session-windows, %2d animals, %5.0f h | median units: PYR %g, INT %g, PYRdeep %g, PYRsup %g\n', g{1}, sum(k), ...
        numel(unique(T.Animal(k))), sum(T.Hours(k)), median(T.n_PYR(k)), median(T.n_INT(k)), median(T.n_PYRdeep(k)), median(T.n_PYRsup(k)));
end

%% ---- tests -------------------------------------------------------------------------------
prs = {'Control','PreASV'; 'Control','ASVshort'; 'Control','ASVlong'; 'PreASV','ASVshort'; 'PreASV','ASVlong';
       'PreSAL','SALshort'; 'PreSAL','SALlong'; 'SALshort','ASVshort'; 'SALlong','ASVlong'; 'SALshort','SALaShort'; 'Control','BaseAll'};
pn  = {'C-Pre','C-Ashort','C-Along','Pre-Ashort','Pre-Along','PreS-Sshort','PreS-Slong','Sshort-Ashort','Slong-Along','Sshort-SAshort','C-Base'};
ppr = {'PreASV','ASVshort'; 'PreASV','ASVlong'; 'PreSAL','SALshort'; 'PreSAL','SALlong'; 'SALaPre','SALaShort'};
ppn = {'Pre-Ashort','Pre-Along','PreS-Sshort','PreS-Slong','SAPre-SAshort'};
dd  = {'ASVshort','PreASV','SALshort','PreSAL'; 'ASVlong','PreASV','SALlong','PreSAL'; 'SALaShort','SALaPre','SALshort','PreSAL'};
ddn = {'short (A-Pre)-(S-PreS)','long (A-Pre)-(S-PreS)','SAL response after - before ASV'};
fprintf(['perm / LME: session-windows, unpaired.  paired: sign-flip on animal means.\n' ...
         'DiD: sign-flip on per-animal (ASV change from PreASV) - (SAL change from PreSAL).\n']);
mets = {'rate','logDur','freq'};
for s = ST, mets = [mets, strcat({'frac','xfrac','rRip','gain','gainW','ppc','pop'}, ['_' s{1}])]; end %#ok<AGROW>
fprintf('\nmean +- SEM across session-windows; rate, rRip, gain as geometric means, tested on log10\n');
fprintf(['%-12s' repmat(' | %-20s', 1, numel(G)) '\n'], 'metric', G{:});
OUT = struct();
for mm = mets
    nm = mm{1}; isLog = any(startsWith(nm, {'rate','rRip','gain'}));
    t = T; if isLog, t.(nm) = log10(t.(nm)); t.(nm)(~isfinite(t.(nm))) = NaN; end
    v = cellfun(@(g) t.(nm)(strcmp(t.Group, g) & isfinite(t.(nm))), G, 'UniformOutput', false);
    if isLog, s = cellfun(@(x) sprintf('%.3g (n=%d)', 10^mean(x), numel(x)), v, 'UniformOutput', false);
    else,     s = cellfun(@local_ms, v, 'UniformOutput', false); end
    fprintf(['%-12s' repmat(' | %-20s', 1, numel(G)) '\n'], nm, s{:});
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
fprintf('\npooled preferred phase (deg, circular mean over session-windows; within-dataset only)\n');
for s = ST
    nm = ['ph_' s{1}];
    c = cellfun(@(g) rad2deg(angle(mean(exp(1i*T.(nm)(strcmp(T.Group, g) & isfinite(T.(nm))))))), G);
    fprintf('%-10s  %s\n', s{1}, strjoin(arrayfun(@(j) sprintf('%s %+6.1f', G{j}, c(j)), 1:numel(G), 'UniformOutput', false), ' | '));
end
save(fullfile(CACHE_ROOT, 'compareRippleWindowsByCellType.mat'), 'T', 'UP', 'OUT', 'SHORT_WIN', 'LONG_AFTER', 'LATER_INFUSIONS');

%% ---- figure --------------------------------------------------------------------------------
col = {[0.4 0.4 0.4], [0 0 0], [0.85 0.2 0.1], [0.55 0.1 0.6], [0.45 0.6 0.75], [0.2 0.6 1], [0 0.25 0.6], [0.6 0.8 0.5], [0.2 0.55 0.2], [0.75 0.55 0.1]};
fig = figure('Color','w','Position',[60 60 1400 640]);
tl = tiledlayout(2, numel(ST), 'TileSpacing','compact');
for s = 1:numel(ST)
    nexttile(s); hold on
    for j = 1:numel(G)
        x = OUT.(['frac_' ST{s}]).v{j};
        plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 3);
        if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', col{j}, 'LineWidth', 2); end
    end
    xlim([0.5 numel(G)+0.5]); xticks(1:numel(G)); xticklabels(G); xtickangle(30); title(['recruitment ' ST{s}]);
    if s == 1, ylabel('fraction of units active per ripple'); end
    ax = nexttile(numel(ST) + s); hold(ax, 'on');
    for j = 1:numel(G)
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
    if s == 1, ylabel(ax, 'cumulative fraction (session-windows equal)'); legend(ax, G, 'Location','southeast', 'Box','off'); end
end
title(tl, 'Ripples by cell type: controls; before / 2-8 h / > 24 h after first ASV; same for first SAL (before ASV)');
exportgraphics(fig, fullfile(CACHE_ROOT, 'compareRippleWindowsByCellType.pdf'), 'ContentType','vector');


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
% out-of-ripple rate in this window (spikes file), for rate-corrected recruitment and gain
rOut = nan(nU, 1);
if ~isempty(sf)
    tt = Ls.spikes.times(:); si = RP.units.SpikesIdx(:);
    if max(si) <= numel(tt)
        nWin = cellfun(@(x) sum(x >= w(1) & x < w(2)), tt(si));
        rOut = max(nWin(:) - nSpk, 0) / max(diff(w) - ripT, eps);
    end
end
gainW = rRip ./ rOut;
dur = RP.ts1Off(in) - RP.ts1(in);
th = sp.Phase; u = sp.Unit;
n = accumarray(u, 1, [nU 1]);
Cs = accumarray(u, cos(th), [nU 1]); Sn = accumarray(u, sin(th), [nU 1]);
ppc = (Cs.^2 + Sn.^2 - n) ./ (n .* (n - 1)); ppc(n < minSpk) = NaN;
up = table(ppc, rRip, gain, M(:,1), M(:,2), M(:,3), M(:,4), 'VariableNames', [{'ppc','rRip','gain'}, ST]);

y = {nR, diff(w) / 3600, mean(log(ev.Dur), 'omitnan'), mean(ev.Freq, 'omitnan'), mean(ok)};
for k = 1:numel(ST)
    mk = M(:, k); nk = sum(mk);
    v = nan(1, 8);
    if nk >= minStrat
        v(1) = mean(full(sum(A(:, mk), 2)) / nk);
        if all(isfinite(rOut(mk)))
            pExp = 1 - exp(-dur(:) * rOut(mk)');                    % ripples x units, Poisson chance of >= 1 spike
            v(7) = v(1) - mean(mean(pExp, 2));                     % recruitment beyond chance at the window's own rate
            gw = gainW(mk); v(8) = median(gw(isfinite(gw)));
        end
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