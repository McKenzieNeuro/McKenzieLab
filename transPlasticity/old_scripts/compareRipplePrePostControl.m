%% compareRipplePrePostControl.m
% Ripple properties, session by session:
%   Control   never-infused recordings (R:\HuszarR), whole recording
%   PreASV    TransPlasticity sessions on days BEFORE the animal's first ASV
%   PostASV   TransPlasticity sessions on days AFTER the animal's first ASV
% The first-ASV day itself is left out. On any other infusion day (SAL before
% ASV, later ASV/SAL) only the pre-infusion part of the recording is used, so
% every window is drug- and procedure-free. INCLUDE_NO_ASV adds the sessions
% of animals that never received ASV to PreASV (drug-naive throughout).
%
% REPLICATE = SESSION in every group (each session one value). The LME rows
% add a random intercept per animal, i.e. the same comparison with sessions
% nested in animals, because sessions of one mouse share probe and units.
%
% METRICS: rate (ripples/h), logDur, freq, fracActive, unitPPC (mean over
% units with >= MINSPK in-ripple spikes), popPPC. Figure: per-session values
% by group, and the unit-PPC distribution with every SESSION weighted equally
% (pooling units lets sessions with many units dominate).
%
% Between datasets, lab / sorting / probes / detector / phase method differ
% (see compareRippleLevels.m); PreASV vs PostASV is within dataset.

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
MIN_EPOCH = 1800; MINSPK = 20; NPERM = 1e4;
INCLUDE_NO_ASV = true;

fa = getAllExtFiles(ROOT, 'mat', 1);         fa = fa(endsWith(fa, 'rippleEvents.mat'));
fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, 'rippleEvents_control.mat'));

rows = {}; UP = {};
for i = 1:numel(fc)
    L = load(fc{i}, 'RIPPLES'); RP = L.RIPPLES;
    [y, up] = local_epoch(RP, [0 RP.meta.recordingSec], MINSPK);
    rows(end+1,:) = [{'Control', controlAnimal(fc{i}), fc{i}}, y]; UP{end+1,1} = up; %#ok<AGROW>
end
nSkip = 0; nNoASV = 0;
for i = 1:numel(fa)
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
    if noASV
        if ~INCLUDE_NO_ASV, continue; end
        g = 'PreASV'; nNoASV = nNoASV + 1;
    elseif ~isfinite(m.dayFromFirstASV) || m.dayFromFirstASV == 0
        nSkip = nSkip + 1; continue                % first-ASV day: left out
    elseif m.dayFromFirstASV < 0
        g = 'PreASV';
    else
        g = 'PostASV';
    end
    w = [0, m.recordingSec];
    if isfinite(m.infusionSec), w = [0, m.infusionSec]; end   % other infusion days: pre-infusion only
    if diff(w) < MIN_EPOCH, continue; end
    [y, up] = local_epoch(RP, w, MINSPK);
    rows(end+1,:) = [{g, m.animal, fa{i}}, y]; UP{end+1,1} = up; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Group','Animal','Session','Event','Hours','logDur','freq', ...
    'fracActive','unitPPC','popPPC'});
T.rate = T.Event ./ T.Hours;
G = {'Control','PreASV','PostASV'};
fprintf('first-ASV day left out: %d session(s); PreASV includes %d session(s) of never-ASV animals\n', nSkip, nNoASV);
for g = G
    k = strcmp(T.Group, g{1});
    fprintf('%-8s %3d sessions, %2d animals, %.0f h\n', g{1}, sum(k), numel(unique(T.Animal(k))), sum(T.Hours(k)));
end

%% ---- tests -------------------------------------------------------------------------------
metrics = {'rate','logDur','freq','fracActive','unitPPC','popPPC'};
fprintf('\nsession means (mean +- SEM across sessions)\n');
fprintf('%-10s | %-20s | %-20s | %-20s || %-30s || %-30s\n', 'metric', G{:}, ...
    'perm p: C-Pre / C-Post / Pre-Post', 'LME (animal RE): same three');
OUT = struct();
for mm = metrics
    nm = mm{1};
    v = cellfun(@(g) T.(nm)(strcmp(T.Group, g) & isfinite(T.(nm))), G, 'UniformOutput', false);
    pr = [local_perm(v{1}, v{2}, NPERM), local_perm(v{1}, v{3}, NPERM), local_perm(v{2}, v{3}, NPERM)];
    pl = [local_lme(T, nm, 'Control', 'PreASV'), local_lme(T, nm, 'Control', 'PostASV'), local_lme(T, nm, 'PreASV', 'PostASV')];
    OUT.(nm) = struct('v', {v}, 'pPerm', pr, 'pLME', pl);
    s = cellfun(@local_ms, v, 'UniformOutput', false);
    fprintf('%-10s | %-20s | %-20s | %-20s || %.3f / %.3f / %.3f    || %.3f / %.3f / %.3f\n', nm, s{:}, pr, pl);
end
fprintf(['perm: sessions as independent replicates. LME: y ~ Group + (1|Animal), sessions nested in\n' ...
         'animals; with only a few control mice its p for Control contrasts is the more honest one.\n']);
save(fullfile(CACHE_ROOT, 'compareRipplePrePostControl.mat'), 'T', 'UP', 'OUT');

%% ---- figure ------------------------------------------------------------------------------------
col = {[0.4 0.4 0.4], [0 0 0], [0.55 0.1 0.6]};
fig = figure('Color','w','Position',[60 60 1500 560]);
tl = tiledlayout(2, numel(metrics), 'TileSpacing','compact');
for k = 1:numel(metrics)
    nexttile(k); hold on
    for j = 1:3
        x = OUT.(metrics{k}).v{j};
        plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 3);
        if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', col{j}, 'LineWidth', 2); end
    end
    xlim([0.5 3.5]); xticks(1:3); xticklabels(G); title(metrics{k}, 'Interpreter','none');
end
% unit PPC distribution, every session weighted equally
ax = nexttile(numel(metrics) + 1, [1 numel(metrics)]); hold(ax, 'on');
for j = 1:3
    idx = find(strcmp(T.Group, G{j}));
    x = []; wgt = [];
    for q = idx(:)'
        u = UP{q}; u = u(isfinite(u));
        if isempty(u), continue; end
        x = [x; u(:)]; wgt = [wgt; repmat(1/numel(u), numel(u), 1)]; %#ok<AGROW>
    end
    [x, o] = sort(x); c = cumsum(wgt(o)); c = c / c(end);
    plot(ax, x, c, '-', 'Color', col{j}, 'LineWidth', 2);
end
xline(ax, 0, 'k:'); xlabel(ax, 'unit PPC'); ylabel(ax, 'cumulative fraction (sessions weighted equally)');
legend(ax, G, 'Location','southeast', 'Box','off');
title(tl, 'Ripples: never-infused controls vs TransPlasticity before / after first ASV (infusion day excluded)');
exportgraphics(fig, fullfile(CACHE_ROOT, 'compareRipplePrePostControl.pdf'), 'ContentType','vector');


%% =========================================================================
function [y, ppcU] = local_epoch(RP, w, minSpk)
% {nRip, hours, logDur, freq, fracActive, unitPPC, popPPC} and per-unit PPCs, for one window
in = RP.ts1 >= w(1) & RP.ts1Off <= w(2);
ev = RP.ev(in, :);
sp = RP.spikes(ismember(RP.spikes.Ripple, find(in)) & isfinite(RP.spikes.Phase), :);
nU = height(RP.units); th = sp.Phase; u = sp.Unit;
n = accumarray(u, 1, [nU 1]);
Cs = accumarray(u, cos(th), [nU 1]); Sn = accumarray(u, sin(th), [nU 1]);
ppc = (Cs.^2 + Sn.^2 - n) ./ (n .* (n - 1));
ppcU = ppc(n >= minSpk);
nP = sum(n); popPPC = NaN;
if nP >= 2, popPPC = (sum(Cs)^2 + sum(Sn)^2 - nP) / (nP * (nP - 1)); end
y = {height(ev), diff(w) / 3600, mean(log(ev.Dur), 'omitnan'), mean(ev.Freq, 'omitnan'), ...
     mean(ev.FracActive, 'omitnan'), mean(ppcU), popPPC};
end


function p = local_perm(x, y, nperm)
p = NaN;
if numel(x) >= 2 && numel(y) >= 2, p = permtest(x, y, nperm); end
end


function p = local_lme(T, nm, a, b)
% y ~ Group + (1|Animal) on sessions of groups a and b
p = NaN;
t = T(ismember(T.Group, {a, b}) & isfinite(T.(nm)), :);
if height(t) < 6 || numel(unique(t.Group)) < 2, return; end
t.y = t.(nm); if strcmp(nm, 'rate'), t.y = log(t.y); end
t.Group = categorical(t.Group, {a, b}); t.Animal = categorical(t.Animal);
try
    lme = fitlme(t, 'y ~ Group + (1|Animal)');
    av = anova(lme, 'DFMethod', 'Satterthwaite');
    p = av.pValue(strcmp(av.Term, 'Group'));
catch ME
    warning('cmp:lme', '%s: %s', nm, ME.message);
end
end


function s = local_ms(v)
if isempty(v), s = 'n/a'; return; end
s = sprintf('%.3g +- %.2g (n=%d)', mean(v), std(v) / sqrt(numel(v)), numel(v));
end
