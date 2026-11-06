%% plotTaggedSynchrony.m
% Figure for the OptoRAM tagged-synchrony result. Run after runTaggedSynchrony.m
% (reads taggedSynchronyResults.mat; recomputes only what that file does not
% store, for ONE representative session).
%
%   A  raster, representative session: the modulated group (k units) above a
%      random same-size draw of untagged units, on the stitched light-off
%      timeline, with the SPIKE-Synchronization profile C(t) of each group
%      underneath (shared x axis)
%   B  that session's test: observed C of the modulated group against its
%      null of size-matched random draws (the stored Cnull)
%   C  pairwise SPIKE-Sync in that session: modulated-modulated,
%      modulated-untagged and untagged-untagged pairs (ECDFs)
%   D  light-off firing rates, modulated vs untagged, same session -- the
%      confound check for C: a rate difference between the groups would need
%      reporting even though the coincidence window is rate-adaptive
%   E  every session's z, grouped by animal (dot size = k), animal means as
%      bars, animal-level sign-flip p in the title
%   F  every session's observed C against the mean of its own null (+-1 SD);
%      points above the unity line are sessions where the modulated group is
%      more synchronous than same-size random groups
%
% CHOICES THAT ARE NOT COSMETIC
%   Representative session (REP_IDX = []): the session whose z is closest to
%   the MEDIAN z across sessions, not the best one. Set REP_IDX to override.
%   Raster window: the WIN-s window whose local C of the modulated group
%   (spike-weighted mean of the profile, i.e. the same estimator as C itself)
%   is closest to that session's overall C, among windows with at least
%   median spiking. So the raster shows what the session's C looks like on
%   average, not its most synchronous stretch.
%   Untagged draw in A: drawn from the UNTAGGED pool only, for a clean visual
%   contrast. The stored null in B/E/F drew from whatever pool
%   taggedSynchrony used ('NullPool', default 'all' = tagged + untagged).
%   Panel C gives no p-value: pairs share units, so pairs are not independent
%   observations. The inferential test is B (per session) and E (animals).
%
% The stitched timeline is rebuilt from what the result itself stores
% (r.rawPulses, r.params.Pad), so it cannot drift from the run. As a check,
% the recomputed light-off duration and C are compared with the stored ones.

CACHE_ROOT  = 'C:\Users\samckenzie\syncCache';
SAVE_TAG    = '';       % must match runTaggedSynchrony's SAVE_TAG
REP_IDX     = [];       % index into U; [] = session with median z
WIN         = 10;       % s of stitched light-off time in the raster
PROF_BIN    = 0.5;      % s, bin for the C(t) trace
PAIR_MAX_UNTAGGED = 80; % cap on untagged units in panel C (pairs grow as n^2)
SEED        = 1;

colT = [0.85 0.37 0.01];     % modulated group
colU = [0.45 0.45 0.45];     % untagged
colX = [0.46 0.44 0.70];     % modulated-untagged pairs

%% ---- load session-level results -----------------------------------------
resFile = fullfile(CACHE_ROOT, ['taggedSynchronyResults' SAVE_TAG '.mat']);
if ~isfile(resFile)
    error('plotTaggedSynchrony:noResults', '%s not found -- run runTaggedSynchrony first.', resFile);
end
L = load(resFile, 'U', 'anim');
U = L.U(:); anim = L.anim(:);
if isempty(U), error('plotTaggedSynchrony:empty', 'no usable sessions in %s', resFile); end
if ~isfield(U, 'params') || isempty(fieldnames(U(1).params))
    error('plotTaggedSynchrony:oldResults', ...
        '%s predates the current taggedSynchrony (no stored params/pulses) -- rerun runTaggedSynchrony.', resFile);
end
switch lower(U(1).params.ModulationSign)
    case 'decrease', GROUP_LABEL = 'suppressed';
    case 'increase', GROUP_LABEL = 'activated';
    otherwise,       GROUP_LABEL = 'modulated';
end
PAD = U(1).params.Pad;
z    = [U.z]';
Cobs = [U.Cobs]';
nT   = [U.nTagged]';
Cnm  = arrayfun(@(r) mean(r.Cnull), U);
Cnsd = arrayfun(@(r) std(r.Cnull),  U);

if isempty(REP_IDX), [~, REP_IDX] = min(abs(z - median(z))); end
r = U(REP_IDX);
[~, sessName] = fileparts(fileparts(r.file));
fprintf('representative session %d/%d: %s (%s), k = %d, z = %+.2f\n', ...
    REP_IDX, numel(U), sessName, anim{REP_IDX}, r.nTagged, r.z);

%% ---- rebuild the stitched light-off timeline exactly as taggedSynchrony --
S = load(r.file, 'spikes');
tsAll = cellfun(@(x) sort(x(:)), S.spikes.times(:), 'UniformOutput', false);
lastSpk = cellfun(@(x) max([x; -Inf]), tsAll);
dur = max(lastSpk(isfinite(lastSpk)));
[tsOff, durOff] = stitchOutIntervals(tsAll, dur, r.rawPulses, r.params.Pad);
if abs(durOff - r.lightOffDur) > 1e-6
    warning('plotTaggedSynchrony:dur', ['stitched light-off duration %.2f s ~= stored %.2f s ', ...
        '-- the spikes file has changed since the run.'], durOff, r.lightOffDur);
end

iT = r.taggedIdx(:); iU = r.untaggedIdx(:); k = numel(iT);

[~, ~, profT] = spikeSync(tsOff(iT), 'Interval', [0 durOff]);
if abs(mean(profT.C) - r.Cobs) > 1e-9
    warning('plotTaggedSynchrony:C', 'recomputed C %.6f ~= stored %.6f', mean(profT.C), r.Cobs);
end

rng(SEED);
iUd = iU(randperm(numel(iU), k));
[~, ~, profU] = spikeSync(tsOff(iUd), 'Interval', [0 durOff]);

%% ---- choose the raster window (representative, not best) ----------------
WIN = min(WIN, durOff);
wStart = 0:WIN/2:(durOff - WIN);
locC = nan(size(wStart)); nSp = zeros(size(wStart));
for w = 1:numel(wStart)
    m = profT.t >= wStart(w) & profT.t < wStart(w) + WIN;
    nSp(w) = sum(m);
    if nSp(w) > 0, locC(w) = mean(profT.C(m)); end
end
cost = abs(locC - r.Cobs);
cost(nSp < max(2*k, median(nSp)) | ~isfinite(cost)) = Inf;
[~, w] = min(cost);
win = wStart(w) + [0 WIN];
mU = profU.t >= win(1) & profU.t < win(2);
locCU = mean(profU.C(mU));

%% ---- pairwise C, representative session ---------------------------------
iUp = iU;
if numel(iUp) > PAIR_MAX_UNTAGGED
    iUp = iUp(sort(randperm(numel(iUp), PAIR_MAX_UNTAGGED)));
end
allIdx = [iT; iUp]; nAll = numel(allIdx);
isT = [true(k,1); false(numel(iUp),1)];
Cp = nan(nAll);
fprintf('pairwise SPIKE-Sync on %d units (%d pairs)...\n', nAll, nAll*(nAll-1)/2);
for a = 1:nAll-1
    for b = a+1:nAll
        Cp(a,b) = spikeSync(tsOff(allIdx([a b])), 'Interval', [0 durOff]);
    end
end
up = triu(true(nAll), 1);
cTT = Cp(up & (isT & isT'));
cUU = Cp(up & (~isT & ~isT'));
cTU = Cp(up & xor(isT, isT'));

FR = cellfun(@numel, tsOff) / durOff;

%% ================================ FIGURE ==================================
fig = figure('Color','w', 'Position',[40 40 1450 1050]);
tl = tiledlayout(fig, 4, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf(['Light-off synchrony of the %s group vs size-matched random units (OptoRAM; ', ...
    'light-on excluded with %g ms pad)'], GROUP_LABEL, 1000*PAD), 'FontWeight','bold');

% ---- A: raster ------------------------------------------------------------
axR = nexttile(tl, 1, [2 2]); hold(axR,'on');
[~, oU] = sort(FR(iUd), 'descend');
[~, oT] = sort(FR(iT),  'descend');
local_ticks(axR, tsOff(iUd(oU)), 1:k,         win, colU);
local_ticks(axR, tsOff(iT(oT)),  k+1+(1:k),   win, colT);
yline(axR, k+1, '-', 'Color',[0.85 0.85 0.85]);
set(axR, 'YTick',[ (1+k)/2, k+1+(1+k)/2 ], ...
    'YTickLabel',{sprintf('untagged (random %d)',k), sprintf('%s (k = %d)',GROUP_LABEL,k)}, ...
    'YDir','normal', 'XTickLabel',[], 'TickLength',[0 0], 'Box','off');
ylim(axR, [0.3 2*k+1.7]); xlim(axR, win);
title(axR, sprintf(['A   %s (%s), %g s of stitched light-off time.  Local C: %s %.3f, ', ...
    'untagged %.3f  (session C = %.3f)'], strrep(sessName,'_','\_'), anim{REP_IDX}, WIN, ...
    GROUP_LABEL, locC(w), locCU, r.Cobs), 'FontWeight','normal', 'HorizontalAlignment','left', ...
    'Units','normalized', 'Position',[0 1.01 0]);

% ---- B: session null ------------------------------------------------------
axB = nexttile(tl, 3); hold(axB,'on');
histogram(axB, r.Cnull, 30, 'FaceColor',colU, 'EdgeColor','none', 'FaceAlpha',0.6);
xline(axB, r.Cobs, '-', 'Color',colT, 'LineWidth',2);
xlabel(axB, 'C (light off)'); ylabel(axB, 'null draws');
title(axB, sprintf('B   observed vs %d size-matched draws\nz = %+.2f, p = %.3f', ...
    numel(r.Cnull), r.z, r.p), 'FontWeight','normal');
legend(axB, {'random k units', GROUP_LABEL}, 'Location','best', 'Box','off');

% ---- C: pairwise ECDFs ----------------------------------------------------
axC = nexttile(tl, 6); hold(axC,'on');
hC = [local_ecdf(axC, cUU, colU), local_ecdf(axC, cTU, colX), local_ecdf(axC, cTT, colT)];
xlabel(axC, 'pairwise SPIKE-Sync'); ylabel(axC, 'cumulative fraction');
legend(axC, hC, { sprintf('untagged-untagged  n=%d, med %.3f', numel(cUU), median(cUU)), ...
                  sprintf('%s-untagged  n=%d, med %.3f', GROUP_LABEL, numel(cTU), median(cTU)), ...
                  sprintf('%s-%s  n=%d, med %.3f', GROUP_LABEL, GROUP_LABEL, numel(cTT), median(cTT)) }, ...
    'Location','southeast', 'Box','off', 'FontSize',8);
title(axC, 'C   pairs, same session (pairs share units: descriptive)', 'FontWeight','normal');

% ---- A (lower): C(t) profile ----------------------------------------------
axP = nexttile(tl, 7, [1 2]); hold(axP,'on');
[tcU, ccU] = local_profileTrace(profU, win, PROF_BIN);
[tcT, ccT] = local_profileTrace(profT, win, PROF_BIN);
stairs(axP, tcU - PROF_BIN/2, ccU, 'Color',colU, 'LineWidth',1.2);
stairs(axP, tcT - PROF_BIN/2, ccT, 'Color',colT, 'LineWidth',1.5);
yline(axP, r.Cobs, ':', 'Color',colT);
ylim(axP, [0 1]); xlim(axP, win);
xlabel(axP, 'stitched light-off time (s)');
ylabel(axP, sprintf('C(t), %.1f-s bins', PROF_BIN));
legend(axP, {'untagged draw', GROUP_LABEL, 'session C'}, 'Location','northeast', 'Box','off', ...
    'Orientation','horizontal');
linkaxes([axR axP], 'x');

% ---- D: firing rates ------------------------------------------------------
axD = nexttile(tl, 9); hold(axD,'on');
rng(SEED);
scatter(axD, 1 + 0.3*(rand(numel(iU),1)-0.5), FR(iU), 14, colU, 'filled', 'MarkerFaceAlpha',0.5);
scatter(axD, 2 + 0.3*(rand(k,1)-0.5),         FR(iT), 22, colT, 'filled', 'MarkerFaceAlpha',0.8);
plot(axD, [0.75 1.25], median(FR(iU))*[1 1], 'k-', 'LineWidth',2);
plot(axD, [1.75 2.25], median(FR(iT))*[1 1], 'k-', 'LineWidth',2);
set(axD, 'YScale','log', 'XTick',[1 2], 'XTickLabel',{'untagged', GROUP_LABEL}); xlim(axD, [0.4 2.6]);
ylabel(axD, 'light-off rate (Hz)');
title(axD, sprintf('D   rates, same session (median %.2f vs %.2f Hz)', ...
    median(FR(iU)), median(FR(iT))), 'FontWeight','normal');

% ---- E: per-session z by animal -------------------------------------------
axE = nexttile(tl, 10, [1 2]); hold(axE,'on');
ua = unique(anim); nA = numel(ua); colA = lines(nA);
zA = nan(nA,1); xs = nan(numel(U),1);
rng(SEED);
for a = 1:nA
    m = find(strcmp(anim, ua{a}));
    xs(m) = a + 0.4*(rand(numel(m),1) - 0.5);
    scatter(axE, xs(m), z(m), 10 + 5*nT(m), colA(a,:), 'filled', 'MarkerFaceAlpha',0.6);
    zA(a) = mean(z(m));
    plot(axE, a + [-0.3 0.3], zA(a)*[1 1], 'k-', 'LineWidth',2);
end
plot(axE, xs(REP_IDX), z(REP_IDX), 'ko', 'MarkerSize',13, 'LineWidth',1.2);
yline(axE, 0, ':');
set(axE, 'XTick',1:nA, 'XTickLabel',ua); xlim(axE, [0.4 nA+0.6]);
ylabel(axE, 'z  (C vs size-matched null)');
pA = local_signflip(zA, 1e4);
ttl = sprintf(['E   %d sessions, %d animals (dot size = k; bar = animal mean; ring = panel A).  ', ...
    'Animal-level sign-flip: %d/%d > 0, p = %.3f'], numel(U), nA, sum(zA>0), nA, pA);
if nA < 6, ttl = sprintf('%s (exact floor %.3f)', ttl, 2/2^nA); end
title(axE, ttl, 'FontWeight','normal');

% ---- F: observed vs null mean ---------------------------------------------
axF = nexttile(tl, 12); hold(axF,'on');
for a = 1:nA
    m = strcmp(anim, ua{a});
    errorbar(axF, Cnm(m), Cobs(m), [], [], Cnsd(m), Cnsd(m), 'o', 'Color',colA(a,:), ...
        'MarkerFaceColor',colA(a,:), 'MarkerSize',4, 'CapSize',0);
end
lim = [min([Cnm - Cnsd; Cobs]) max([Cnm + Cnsd; Cobs])];
lim = lim + [-1 1]*0.05*diff(lim);
plot(axF, lim, lim, 'k:');
axis(axF, 'square'); xlim(axF, lim); ylim(axF, lim);
xlabel(axF, 'null mean C (\pm1 SD)'); ylabel(axF, sprintf('observed C, %s', GROUP_LABEL));
title(axF, sprintf('F   %d/%d sessions above unity', sum(Cobs > Cnm), numel(U)), 'FontWeight','normal');

exportgraphics(fig, fullfile(CACHE_ROOT,'taggedSynchrony_fig.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile(CACHE_ROOT,'taggedSynchrony_fig.png'), 'Resolution',300);


%% =========================================================================
function local_ticks(ax, ts, yRow, win, col)
% one line object per block: ticks separated by NaNs
X = []; Y = [];
for n = 1:numel(ts)
    t = ts{n}(ts{n} >= win(1) & ts{n} < win(2)); t = t(:)';
    if isempty(t), continue; end
    X = [X, [t; t; nan(1,numel(t))]]; %#ok<AGROW>
    Y = [Y, repmat([yRow(n)-0.4; yRow(n)+0.4; NaN], 1, numel(t))]; %#ok<AGROW>
end
plot(ax, X(:), Y(:), '-', 'Color', col, 'LineWidth', 0.8);
end

function [tc, cc] = local_profileTrace(prof, win, binW)
% mean C(t_k) of the spikes falling in each bin (NaN where no spikes)
edges = win(1):binW:win(2);
if edges(end) < win(2), edges(end+1) = edges(end) + binW; end
m = prof.t >= edges(1) & prof.t < edges(end);
b = discretize(prof.t(m), edges);
cc = accumarray(b(:), prof.C(m), [numel(edges)-1 1], @mean, NaN);
tc = edges(1:end-1)' + binW/2;
end

function h = local_ecdf(ax, v, col)
v = sort(v(isfinite(v)));
h = stairs(ax, [v(1); v], [0; (1:numel(v))'/numel(v)], 'Color', col, 'LineWidth', 1.6);
end

function p = local_signflip(d, nperm)
d = d(:); d = d(isfinite(d));
n = numel(d); obs = mean(d); null = zeros(nperm,1);
for k = 1:nperm, null(k) = mean((2*(rand(n,1)>0.5)-1).*d); end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end