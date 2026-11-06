%% plotUnitPPCPolar.m
% Intuition for unit PPC (pairwise phase consistency of a unit's spikes to the
% ripple oscillation), in three figures:
%
%   1. What PPC values look like: von Mises phase distributions whose PPC is
%      0.02 / 0.09 / 0.16 / 0.35 (the range seen in the data), each with 200
%      sampled spikes, plus the same for a uniform (PPC 0) unit sampled with 30
%      spikes to show why the resultant length R is biased and PPC is not.
%   2. Real units: in-ripple spike phases of the best, median and worst-locked
%      well-sampled units in one control session, one TransPlasticity drug-naive
%      window and one post-exposure window.
%   3. Distribution of unit PPC across all well-sampled units by group.
%
% PPC = mean over all spike pairs of cos(theta_i - theta_j)
%     = (|sum exp(i*theta)|^2 - n) / (n(n-1)),
% an unbiased estimate of R^2 (squared resultant length of the unit's true
% phase distribution), so it does not grow just because a unit fired few spikes.
% Phase 0 is the detector's zero: Hilbert peak of the band-passed LFP for
% TransPlasticity, the stored phase_ts convention for the controls. Compare
% WIDTH across datasets, not where the peak sits.

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
MINSPK = 20;
LAT = 4*3600; POSTLEN = 2*3600;
NB = 24;                                   % polar histogram bins

%% ---- 1. synthetic: what a given PPC looks like --------------------------------
ppcT = [0.02 0.09 0.16 0.35];
f1 = figure('Color','w','Position',[60 60 1400 330]);
t1 = tiledlayout(f1, 1, numel(ppcT) + 1, 'TileSpacing','compact');
th = linspace(-pi, pi, 361);
rng(1);
for k = 1:numel(ppcT)
    kap = fzero(@(x) (besseli(1, x) / besseli(0, x))^2 - ppcT(k), [1e-6 50]);
    dens = exp(kap * cos(th)) / (2*pi*besseli(0, kap));
    x = local_vmrand(kap, 200);
    ax = polaraxes(t1); ax.Layout.Tile = k; hold(ax, 'on');
    polarhistogram(ax, x, NB, 'Normalization','pdf', 'FaceColor',[0.6 0.6 0.6], 'EdgeColor','none');
    polarplot(ax, th, dens, 'k-', 'LineWidth', 1.5);
    [p, R] = local_ppc(x);
    title(ax, sprintf('PPC %.2f  (R %.2f, circ SD %.0f%s)\nsample n=200: PPC %.2f, R %.2f', ppcT(k), sqrt(ppcT(k)), ...
        rad2deg(sqrt(-2*log(sqrt(ppcT(k))))), char(176), p, R), 'FontWeight','normal');
    ax.RTickLabel = {};
end
x = 2*pi*rand(30, 1) - pi;                  % uniform, few spikes
ax = polaraxes(t1); ax.Layout.Tile = numel(ppcT) + 1;
polarhistogram(ax, x, NB, 'Normalization','pdf', 'FaceColor',[0.85 0.4 0.3], 'EdgeColor','none');
[p, R] = local_ppc(x);
title(ax, sprintf('no locking, n=30 spikes\nR %.2f (biased up)  PPC %.2f', R, p), 'FontWeight','normal');
ax.RTickLabel = {};
title(t1, 'What unit PPC values mean (von Mises phase distributions)');

%% ---- gather units by group -------------------------------------------------------
fa = getAllExtFiles(ROOT, 'mat', 1);         fa = fa(endsWith(fa, 'rippleEvents.mat'));
fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, 'rippleEvents_control.mat'));
U = struct('group',{},'session',{},'unit',{},'phases',{},'ppc',{},'R',{},'pref',{});
for i = 1:numel(fc)
    L = load(fc{i}, 'RIPPLES'); RP = L.RIPPLES;
    U = [U, local_units(RP, [0 RP.meta.recordingSec], 'Control', fc{i}, MINSPK)]; %#ok<AGROW>
end
for i = 1:numel(fa)
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    dur = m.recordingSec; infS = m.infusionSec;
    drugFree = [0 dur]; if isfinite(infS), drugFree = [0 infS]; end
    noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
    if noASV || (isfinite(m.dayFromFirstASV) && m.dayFromFirstASV < 0)
        U = [U, local_units(RP, drugFree, 'Naive', fa{i}, MINSPK)]; %#ok<AGROW>
    elseif strcmp(m.folder, m.firstASV.session)
        U = [U, local_units(RP, [0 infS], 'Naive', fa{i}, MINSPK), ...
               local_units(RP, [infS + LAT, min(dur, infS + LAT + POSTLEN)], 'ASVacute', fa{i}, MINSPK)]; %#ok<AGROW>
    elseif m.dayFromFirstASV >= 1
        U = [U, local_units(RP, drugFree, 'Post', fa{i}, MINSPK)]; %#ok<AGROW>
    end
end
G = {'Control','Naive','ASVacute','Post'};
col = {[0.4 0.4 0.4], [0 0 0], [0.85 0.2 0.1], [0.55 0.1 0.6]};
for g = G, fprintf('%-9s %4d well-sampled units (>= %d in-ripple spikes)\n', g{1}, sum(strcmp({U.group}, g{1})), MINSPK); end

%% ---- 2. real units: best / median / worst in one example window per group ------------
show = {'Control','Naive','Post'};
f2 = figure('Color','w','Position',[60 420 1200 900]);
t2 = tiledlayout(f2, numel(show), 3, 'TileSpacing','compact');
lab = {'most locked','median','least locked'};
for r = 1:numel(show)
    ug = U(strcmp({U.group}, show{r}));
    if isempty(ug), continue; end
    % example window: the session with the most well-sampled units
    [ss, ~, j] = unique({ug.session}); [~, b] = max(accumarray(j(:), 1)); ue = ug(strcmp({ug.session}, ss{b}));
    [~, o] = sort([ue.ppc], 'descend'); pick = o([1, ceil(numel(o)/2), numel(o)]);
    [~, sname] = fileparts(fileparts(ss{b}));
    for c = 1:3
        u = ue(pick(c));
        ax = polaraxes(t2); ax.Layout.Tile = (r-1)*3 + c;
        polarhistogram(ax, u.phases, NB, 'Normalization','pdf', 'FaceColor',col{strcmp(G, show{r})}, 'EdgeColor','none');
        hold(ax, 'on'); rl = ax.RLim(2);
        polarplot(ax, [u.pref u.pref], [0 rl], 'k-', 'LineWidth', 1.5);    % preferred phase
        title(ax, sprintf('%s: %s unit %d\nPPC %.3f  R %.2f  n=%d', show{r}, lab{c}, u.unit, u.ppc, u.R, numel(u.phases)), ...
            'FontWeight','normal', 'Interpreter','none');
        ax.RTickLabel = {};
    end
    fprintf('example %-8s %s (%d units)\n', show{r}, sname, numel(ue));
end
title(t2, 'In-ripple spike phase of single units (line = preferred phase)');

%% ---- 3. distribution of unit PPC by group ------------------------------------------------
f3 = figure('Color','w','Position',[1280 420 520 420]); hold on
for k = 1:numel(G)
    v = sort([U(strcmp({U.group}, G{k})).ppc]);
    if isempty(v), continue; end
    plot(v, (1:numel(v)) / numel(v), '-', 'Color', col{k}, 'LineWidth', 2);
end
xline(0, 'k:');
xlabel('unit PPC'); ylabel('cumulative fraction of units');
legend(G, 'Location','southeast', 'Box','off');
title('Unit PPC across all well-sampled units (units pooled; not a test)', 'FontWeight','normal');

exportgraphics(f1, fullfile(CACHE_ROOT, 'unitPPC_intuition.pdf'), 'ContentType','vector');
exportgraphics(f2, fullfile(CACHE_ROOT, 'unitPPC_examples.pdf'), 'ContentType','vector');
exportgraphics(f3, fullfile(CACHE_ROOT, 'unitPPC_cdf.pdf'), 'ContentType','vector');


%% =========================================================================
function U = local_units(RP, w, grp, file, minSpk)
% well-sampled units in window w: their in-ripple phases, PPC, R, preferred phase
U = struct('group',{},'session',{},'unit',{},'phases',{},'ppc',{},'R',{},'pref',{});
if diff(w) <= 0, return; end
in = find(RP.ts1 >= w(1) & RP.ts1Off <= w(2));
sp = RP.spikes(ismember(RP.spikes.Ripple, in) & isfinite(RP.spikes.Phase), :);
for u = unique(sp.Unit)'
    x = sp.Phase(sp.Unit == u);
    if numel(x) < minSpk, continue; end
    [p, R] = local_ppc(x);
    U(end+1) = struct('group',grp, 'session',file, 'unit',u, 'phases',x, 'ppc',p, 'R',R, ...
        'pref',angle(sum(exp(1i*x)))); %#ok<AGROW>
end
end


function [p, R] = local_ppc(x)
n = numel(x); z = sum(exp(1i * x(:)));
R = abs(z) / n;
p = (abs(z)^2 - n) / (n * (n - 1));
end


function x = local_vmrand(kap, n)
% von Mises samples by inverse CDF on a fine grid
g = linspace(-pi, pi, 4001); pdf = exp(kap * cos(g)); cdf = cumtrapz(g, pdf); cdf = cdf / cdf(end);
[cdf, iu] = unique(cdf); x = interp1(cdf, g(iu), rand(n, 1));
end
