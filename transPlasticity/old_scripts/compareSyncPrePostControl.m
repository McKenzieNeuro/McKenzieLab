%% compareSyncPrePostControl.m
% SPIKE-distance synchrony, one value per session, in the same groups as
% compareRipplePrePostControl:
%   Control   never-infused recordings      spikeDistanceProfile_control.mat
%   PreASV    days BEFORE the first ASV     spikeDistanceProfile.mat
%   PostASV   days AFTER the first ASV      spikeDistanceProfile.mat
% First-ASV day left out; other infusion days use only the pre-infusion part.
% INCLUDE_NO_ASV adds never-ASV animals to PreASV.
%
% PER-SESSION METRICS (all from bin.S / bin.Snull over the session's window)
%   S        mean SPIKE-distance (dissimilarity; lower = more synchronous)
%   Snull    mean of the interval-jitter nulls (mean over surrogates)
%   dS       Snull - S: excess synchrony in S units (> 0 = more synchronous than null)
%   zS       dS / SD over surrogates of their window means
%   eBlk     mean of per-block e over WIN-s blocks, e = (med(null) - S)/(med - p2.5)
%   fracSig  fraction of WIN-s blocks with block p <= 0.05 (observed below >=95% of nulls)
% S itself depends on the units' rates and number (and the pair subset), which
% differ between datasets and sessions; dS / zS / eBlk / fracSig are each
% relative to the session's own null and are the comparable quantities.
%
% TESTS: permutation with sessions as replicates, and LME y ~ Group + (1|Animal)
% (sessions nested in animals: use this one for the Control contrasts).

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
MIN_EPOCH = 1800; WIN = 60; NPERM = 1e4;
INCLUDE_NO_ASV = true;

fa = getAllExtFiles(ROOT, 'mat', 1);         fa = fa(endsWith(fa, 'spikeDistanceProfile.mat'));
fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, 'spikeDistanceProfile_control.mat'));
fprintf('%d infusion-study files, %d control files\n', numel(fa), numel(fc));

rows = {};
for i = 1:numel(fc)
    L = load(fc{i}, 'SDIST'); X = L.SDIST;
    y = local_session(X, [0 X.bin.edges(end)], WIN);
    rows(end+1,:) = [{'Control', controlAnimal(fc{i}), fc{i}, X.N}, y]; %#ok<AGROW>
end
nSkip = 0; nNoASV = 0;
for i = 1:numel(fa)
    L = load(fa{i}, 'SDIST'); X = L.SDIST; m = X.meta;
    noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
    if noASV
        if ~INCLUDE_NO_ASV, continue; end
        g = 'PreASV'; nNoASV = nNoASV + 1;
    elseif ~isfinite(m.dayFromFirstASV) || m.dayFromFirstASV == 0
        nSkip = nSkip + 1; continue
    elseif m.dayFromFirstASV < 0
        g = 'PreASV';
    else
        g = 'PostASV';
    end
    w = [0, X.bin.edges(end)];
    if isfinite(m.infusionSec), w = [0, m.infusionSec]; end
    if diff(w) < MIN_EPOCH, continue; end
    y = local_session(X, w, WIN);
    rows(end+1,:) = [{g, m.animal, fa{i}, X.N}, y]; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'Group','Animal','Session','N','Hours','S','Snull','dS','zS','eBlk','fracSig'});
G = {'Control','PreASV','PostASV'};
fprintf('first-ASV day left out: %d session(s); PreASV includes %d session(s) of never-ASV animals\n', nSkip, nNoASV);
for g = G
    k = strcmp(T.Group, g{1});
    fprintf('%-8s %3d sessions, %2d animals, %.0f h, units per session %s\n', g{1}, sum(k), ...
        numel(unique(T.Animal(k))), sum(T.Hours(k)), mat2str(prctile(T.N(k), [0 50 100])));
end

metrics = {'S','Snull','dS','zS','eBlk','fracSig'};
fprintf('\nsession means (mean +- SEM across sessions)\n');
fprintf('%-8s | %-22s | %-22s | %-22s || %-30s || %-30s\n', 'metric', G{:}, ...
    'perm p: C-Pre / C-Post / Pre-Post', 'LME (animal RE): same three');
OUT = struct();
for mm = metrics
    nm = mm{1};
    v = cellfun(@(g) T.(nm)(strcmp(T.Group, g) & isfinite(T.(nm))), G, 'UniformOutput', false);
    pr = [local_perm(v{1}, v{2}, NPERM), local_perm(v{1}, v{3}, NPERM), local_perm(v{2}, v{3}, NPERM)];
    pl = [local_lme(T, nm, 'Control', 'PreASV'), local_lme(T, nm, 'Control', 'PostASV'), local_lme(T, nm, 'PreASV', 'PostASV')];
    OUT.(nm) = struct('v', {v}, 'pPerm', pr, 'pLME', pl);
    s = cellfun(@local_ms, v, 'UniformOutput', false);
    fprintf('%-8s | %-22s | %-22s | %-22s || %.3f / %.3f / %.3f    || %.3f / %.3f / %.3f\n', nm, s{:}, pr, pl);
end
fprintf('S lower = more synchronous; dS, zS, eBlk, fracSig higher = more excess synchrony vs the jitter null.\n');
save(fullfile(CACHE_ROOT, 'compareSyncPrePostControl.mat'), 'T', 'OUT', 'WIN');

col = {[0.4 0.4 0.4], [0 0 0], [0.55 0.1 0.6]};
fig = figure('Color','w','Position',[60 60 1500 320]);
tiledlayout(1, numel(metrics), 'TileSpacing','compact');
for k = 1:numel(metrics)
    nexttile; hold on
    for j = 1:3
        x = OUT.(metrics{k}).v{j};
        plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 3);
        if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', col{j}, 'LineWidth', 2); end
    end
    xlim([0.5 3.5]); xticks(1:3); xticklabels(G); title(metrics{k});
end
exportgraphics(fig, fullfile(CACHE_ROOT, 'compareSyncPrePostControl.pdf'), 'ContentType','vector');


%% =========================================================================
function y = local_session(X, w, W)
% {hours, S, Snull, dS, zS, eBlk, fracSig} over window w
ed = X.bin.edges(:); ctr = 0.5 * (ed(1:end-1) + ed(2:end));
in = ctr >= w(1) & ctr < w(2);
S = X.bin.S(in); Sn = double(X.bin.Snull(in, :));
bw = diff(ed(1:2));
hours = sum(in) * bw / 3600;
mS = mean(S, 'omitnan'); mN = mean(Sn, 1, 'omitnan');           % window mean, per surrogate
dS = mean(mN) - mS; zS = dS / std(mN);
nb = max(1, round(W / bw)); nW = floor(numel(S) / nb);
eB = NaN; fs = NaN;
if nW >= 1
    Sw = mean(reshape(S(1:nW*nb), nb, nW), 1)';
    Snw = reshape(mean(reshape(Sn(1:nW*nb, :), nb, nW, []), 1), nW, []);
    med = median(Snw, 2); lo = prctile(Snw, 2.5, 2);
    e = (med - Sw) ./ (med - lo);
    pb = (1 + sum(Snw <= Sw, 2)) ./ (size(Snw, 2) + 1);
    eB = mean(e, 'omitnan'); fs = mean(pb <= 0.05);
end
y = {hours, mS, mean(mN), dS, zS, eB, fs};
end


function p = local_perm(x, y, nperm)
p = NaN;
if numel(x) >= 2 && numel(y) >= 2, p = permtest(x, y, nperm); end
end


function p = local_lme(T, nm, a, b)
p = NaN;
t = T(ismember(T.Group, {a, b}) & isfinite(T.(nm)), :);
if height(t) < 6 || numel(unique(t.Group)) < 2, return; end
t.y = t.(nm);
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
