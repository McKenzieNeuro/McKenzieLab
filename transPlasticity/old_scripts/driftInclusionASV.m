%% driftInclusionASV.m
% Does ASV change what the exploration contributes to the ripples that follow it?
%
% THE ENDPOINT, AND WHY THIS ONE
% Nine never-infused control sessions established that an exploration's binned co-firing
% structure predicts the covariance of later ripples above chance:
%       fastCtx = 0.043 against an analytic chance of k/(n-1) = 0.020, p = 0.0039
% and that this is NOT experience-specific: once temporal proximity is controlled, the template
% from the other exploration does just as well (xFastAdj null). So the quantity that is real and
% measurable here is how much exploration-related covariance structure appears in later ripples,
% not whether it is tied to one particular experience.
%
% That is enough for the drug question, and this is the point worth being explicit about. Every
% confound that defeated the within-session specificity tests -- temporal proximity, expert type,
% redundancy with the preceding home cage, the better-template main effect -- is COMMON TO BOTH
% ARMS. ASV and saline sessions share the layout, the lags and the construction, so a between-
% group contrast differences them out by design. Specificity was needed to interpret the control
% result on its own; it is not needed to ask whether a drug moves the measure.
%
% fastCtx is also the only endpoint with an exact null (a random basis captures k/(n-1) of the
% residual variance by symmetry, no surrogate required) and no dependence on a twin whose
% position confounds it.
%
% THE DESIGN THE INFUSION CREATES
%   HC  ConA1  HC  Infusion  HC  ConB1  HC  ConA2  HC  ConB2  HC
% ConA1 is encoded BEFORE the infusion and ConB1 after it, so within one session ConA1 is the
% undrugged reference and ConB1 the drugged test. ConA1's targets are truncated at the infusion
% so its replay window is undrugged too. The endpoint is then a difference in differences:
%       (ConB1 - ConA1)_ASV  -  (ConB1 - ConA1)_SAL
% which removes the animal, the probe, the cell yield, the session layout and the infusion
% procedure itself, since both arms are infused.
%
% Controls are computed identically and printed as a reference row only. They are never tested
% against the infused arms: they have no infusion, so the within-session contrast does not exist
% for them and too much else differs.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT       = 'R:\TransPlasticity';
CONTROLS = { ...
    'R:\HuszarR\e13_26m1_211019', [15149.6 16749.2], [7605.7 9322]; ...
    'R:\HuszarR\e13_26m1_211116', [13608 13939],     [7469.8 8707.2]; ...
    'R:\HuszarR\e13_26m1_211119', [14911.9 16446],   [6181.3 7572]; ...
    'R:\HuszarR\e15_13f1_220118', [13137 15046],     [5636.7 7700]; ...
    'R:\HuszarR\e15_13f1_220119', [21219 22700],     [6836.6 8860]; ...
    'R:\HuszarR\e15_13f1_220120', [18017.9 19801],   [9610.5 11950]; ...
    'R:\HuszarR\e16_3m2_211211',  [9584.6 11203],    [5294.1 6511]; ...
    'R:\HuszarR\e16_3m2_211212',  [16416.3 17731],   [6109.3 7845]; ...
    'R:\HuszarR\e16_3m2_211213',  [13291 15323],     [7567.9 9359]};   % folder, NOVEL, FAMILIAR

SUB_WIN     = 600;        % s, target window width
GUARD       = 600;        % s after the context ends before targets begin
TGT_MAX     = 4*3600;     % s of home cage after the context used as targets
BIN_H       = 1;          % h, bins for the time course
BEH_BIN     = 0.25;       % s, bin width for the behaviour-time template
GEO_RANK    = 2;          % covariance directions
MIN_BIN_CTX = 200;
MIN_SPK_CTX = 200;
MIN_RIP_SUB = 25;
MIN_UNITS   = 9;
% ---- CELL-COUNT MATCHING, which this comparison cannot do without ------------------------
% Subtracting chance removes the chance LEVEL but not the yield dependence of the excess above
% it: with more cells the basis is better estimated, so the same underlying structure produces a
% larger excess. The arms differ here (SAL about 19 pyramidal cells, ASV about 21) and the
% Huszar controls have about 212, so none of these are comparable as they stand.
% Every session is therefore reduced to the SAME number of randomly chosen cells, repeated over
% several draws and averaged. Cell identity is the only thing the draws destroy, and they
% destroy it equally in both arms.
N_CELL_SUB  = 9;         % cells used per session; [] disables matching
N_SUB_DRAW  = 50;         % random draws averaged over
% N_CELL_SUB costs sessions: everything below it is dropped, and at 15 that removed CP15 SALpre
% (13 units), CP17 ASV1 (10), CP17 SALpre (12) and CP19 ASV1 (10), leaving only 3 animals with
% both arms. Dropping to 10 recovers CP15 and CP17 as paired animals, at the price of chance
% rising to GEO_RANK/9. Lower GEO_RANK with it: 2 of 9 is 0.222, about what 4 of 14 gives now.
N_RIP_MATCH = 50;        % ripples per target window, subsampled to this; [] disables
% A target window with more ripples has a better-estimated covariance, so an external basis
% captures closer to its true value and further above chance. SAL averages 186 ripples per
% window and ASV 216, which biases the comparison toward ASV. Matching removes it.
% With N_CELL_SUB cells the residuals span N_CELL_SUB-1 dimensions, so GEO_RANK of them is a
% large fraction and chance sits high: 4 of 14 is 0.286. That costs dynamic range. Lower
% GEO_RANK if the excess is small relative to its own session-to-session spread.
FREE_GAIN   = true;
EXCL_INF    = false;      % as in compareReactivation.m
INF_EXCL_WIN = [0 1800];  % s [before after] Infusion_start, removed from home cage
TRUNC_NEXT  = true;       % stop a context's targets at the next context's start
SAVE_OUT    = true;

rows = {};

%% ---- TransPlasticity -------------------------------------------------------------------
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
led = cell(0, 3);
for i = 1:numel(fa)
    fld = fileparts(fa{i});
    [~, snL] = fileparts(fld);
    ev = local_events(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConA1_start'))
        led(end+1,:) = {snL, '', 'no event file or no ConA1'}; continue %#ok<AGROW>
    end
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    g = local_group(RP, ev);
    if isempty(g), led(end+1,:) = {snL, '', 'not SALpre or ASV1'}; continue; end %#ok<AGROW>
    HC = local_hcBlocks(ev, RP, EXCL_INF, INF_EXCL_WIN);
    tOf = @(lab) ev.time(find(strcmp(ev.description, lab), 1));
    infOn = tOf('Infusion_start'); if isempty(infOn), infOn = NaN; end
    IV = struct('name', {}, 'iv', {});
    for c = {'ConA1', 'ConB1'}
        on = tOf([c{1} '_start']); off = tOf([c{1} '_end']);
        if isempty(on) || isempty(off), led(end+1,:) = {snL, g, ['no ' c{1} ' markers']}; continue; end %#ok<AGROW>
        IV(end+1) = struct('name', c{1}, 'iv', [on off]); %#ok<AGROW>
    end
    allCtx = zeros(0, 2);
    for c = {'ConA1','ConB1','ConA2','ConB2'}
        on = tOf([c{1} '_start']); off = tOf([c{1} '_end']);
        if ~isempty(on) && ~isempty(off), allCtx(end+1,:) = [on off]; end %#ok<AGROW>
    end
    r = local_session(fld, RP, IV, allCtx, HC, infOn, g, m.animal, snL, ...
        SUB_WIN, GUARD, TGT_MAX, BEH_BIN, GEO_RANK, MIN_BIN_CTX, MIN_SPK_CTX, ...
        MIN_RIP_SUB, MIN_UNITS, FREE_GAIN, TRUNC_NEXT, N_CELL_SUB, N_SUB_DRAW, N_RIP_MATCH);
    if isempty(r), led(end+1,:) = {snL, g, 'no usable context expert and targets'}; end %#ok<AGROW>
    rows = [rows; r]; %#ok<AGROW>
end

%% ---- controls, as a reference row only ---------------------------------------------------
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1}; nov = CONTROLS{i,2}; fam = CONTROLS{i,3};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), continue; end
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES;
    IV = struct('name', {'ConA1','ConB1'}, 'iv', {fam, nov});
    allCtx = sortrows([fam; nov]);
    rec = max([RP.ts1Off(:); 1]);
    r = local_session(fld, RP, IV, allCtx, [0 rec], NaN, 'Control', local_name(fld), local_name(fld), ...
        SUB_WIN, GUARD, TGT_MAX, BEH_BIN, GEO_RANK, MIN_BIN_CTX, MIN_SPK_CTX, ...
        MIN_RIP_SUB, MIN_UNITS, FREE_GAIN, TRUNC_NEXT, N_CELL_SUB, N_SUB_DRAW, N_RIP_MATCH);
    rows = [rows; r]; %#ok<AGROW>
end

if isempty(rows), error('asv:none', 'no session produced a usable context expert and targets'); end
T = cell2table(rows, 'VariableNames', {'Group','Animal','Session','Context','HoursAfter', ...
    'nRip','nUnits','fastCtx','fastChance','fastAbove','slowDev','floorDev','drugged', ...
    'nSpkCtx','ctxMin','nCellUsed'});

% ---- SESSION ORDER, which is confounded with Group by construction ------------------------
% SALpre is by definition the saline sessions BEFORE the first ASV, so ASV1 is always later in
% each animal's sequence. Anything that drifts across days -- electrode position, unit yield,
% familiarity with the contexts, handling -- therefore differs between the arms before any drug
% is given, and that is exactly what the ConA1 baseline shows. DayIdx is days since that
% animal's first included session, so the comparison can be made with order held fixed.
T.DayIdx = nan(height(T), 1);
dn = arrayfun(@(q) local_sessionDate(T.Session{q}), (1:height(T))');
anT = unique(T.Animal);
for q = 1:numel(anT)
    k = strcmp(T.Animal, anT{q});
    if any(isfinite(dn(k))), T.DayIdx(k) = dn(k) - min(dn(k & isfinite(dn))); end
end
fprintf('\nsession order (days since each animal''s first included session):\n');
for q = 1:numel(anT)
    k = strcmp(T.Animal, anT{q});
    [ss, iu] = unique(T.Session(k)); dd = T.DayIdx(k); gg = T.Group(k);
    fprintf('   %-8s', anT{q});
    for z = 1:numel(ss), fprintf('  %s %s d+%.0f;', ss{z}, gg{iu(z)}, dd(iu(z))); end
    fprintf('\n');
end

%% ---------------------------------------------------------------------------------------------
%  REPORT
%% ---------------------------------------------------------------------------------------------
fprintf('\n================ ASV AND THE EXPLORATION CONTRIBUTION TO LATER RIPPLES ================\n');
fprintf('endpoint: fastAbove = captured variance of the exploration basis on post-context ripple\n');
fprintf('residuals, MINUS the exact chance level k/(n-1) for that session''s cell count.\n');
fprintf('%d sessions, %d target windows\n', numel(unique(T.Session)), height(T));
if ~isempty(led)
    fprintf('\nnot contributing:\n');
    for q = 1:size(led, 1), fprintf('   %-34s %-7s %s\n', led{q,:}); end
end

G = {'SALpre', 'ASV1', 'Control'};
G = G(ismember(G, unique(T.Group)));
fprintf('\n--- fastAbove by group and context (session means, animal means in brackets) ---\n');
fprintf('%-9s %-7s %10s %10s %10s %8s %8s\n', 'group', 'context', 'fastCtx', 'chance', 'above', 'n ses', 'n anim');
for gg = G
    for c = {'ConA1', 'ConB1'}
        k = strcmp(T.Group, gg{1}) & strcmp(T.Context, c{1});
        if ~any(k), continue; end
        [xs, ses] = local_sesMean(T.fastAbove(k), T.Session(k)); %#ok<ASGLU>
        xa = local_sesMean(T.fastAbove(k), T.Animal(k));
        fprintf('%-9s %-7s %10.4f %10.4f %10.4f %8d %8d\n', gg{1}, c{1}, ...
            mean(T.fastCtx(k), 'omitnan'), mean(T.fastChance(k), 'omitnan'), ...
            mean(xa, 'omitnan'), numel(xs), numel(xa));
    end
end

fprintf('\n--- is the exploration contribution above chance within each arm? ---\n');
for gg = G
    for c = {'ConA1', 'ConB1'}
        k = strcmp(T.Group, gg{1}) & strcmp(T.Context, c{1});
        if ~any(k), continue; end
        x = local_sesMean(T.fastAbove(k), T.Animal(k)); x = x(isfinite(x));
        if numel(x) < 2, continue; end
        fprintf('  %-9s %-7s above chance = %+.4f, p = %.4f (n = %d animals)\n', ...
            gg{1}, c{1}, mean(x), local_signflip(x), numel(x));
    end
end

fprintf('\n--- PRIMARY: ConB1, encoded and replayed under drug. ASV vs SAL ---\n');
kB = strcmp(T.Context, 'ConB1') & ismember(T.Group, {'SALpre','ASV1'});
local_groupTest(T(kB,:), 'fastAbove');

fprintf('\n--- PRE-DRUG BASELINE CHECK: ConA1, encoded BEFORE the infusion. ASV vs SAL ---\n');
fprintf('ConA1 happens before anything is infused and its targets are truncated at the\n');
fprintf('infusion, so the two arms CANNOT differ here for any pharmacological reason. Whatever\n');
fprintf('appears on this line is a difference between the sessions that happen to be in each\n');
fprintf('arm, and it is present in the ConB1 line above as well. Only the excess of ConB1 over\n');
fprintf('this can be attributed to the drug, which is exactly what the DiD below computes.\n');
kA0 = strcmp(T.Context, 'ConA1') & ismember(T.Group, {'SALpre','ASV1'});
local_groupTest(T(kA0,:), 'fastAbove');

fprintf('\n--- PRIMARY DiD: (ConB1 - ConA1) within session, ASV vs SAL ---\n');
fprintf('ConA1 is encoded before the infusion and its targets are truncated at it, so this is\n');
fprintf('the drugged contrast minus the undrugged one in the same animal, then between arms.\n');
Dd = local_within(T, 'fastAbove', 'ConB1', 'ConA1');
local_groupTest(Dd, 'd');
fprintf('\n  the same difference in differences for the SLOW component:\n');
Ds = local_within(T, 'slowDev', 'ConB1', 'ConA1');
local_groupTest(Ds, 'd');
fprintf('  slowDev is the context expert''s deviance above the window''s own oracle, in floor\n');
fprintf('  units, so SMALLER IS BETTER PREDICTION and a negative ASV-SAL means the exploration\n');
fprintf('  rate vector predicts later ripple rates MORE closely under ASV, not less.\n');

fprintf('\n--- secondary: slow (mean rate) deviance of the context expert, in floor units ---\n');
kB2 = strcmp(T.Context, 'ConB1') & ismember(T.Group, {'SALpre','ASV1'});
local_groupTest(T(kB2,:), 'slowDev');
fprintf('  and its own pre-drug baseline, ConA1:\n');
local_groupTest(T(kA0,:), 'slowDev');
fprintf('  Reported as a negative control. In controls the slow component showed no\n');
fprintf('  experience specificity once lag was removed, so a drug effect here would be a\n');
fprintf('  statement about rate, not about replay content.\n');

fprintf('\n--- time course of fastAbove after ConB1 ---\n');
eb = 0:BIN_H:TGT_MAX/3600;
fprintf('%-9s %10s %10s %8s\n', 'group', 'bin (h)', 'above', 'n anim');
for gg = G
    for b = 1:numel(eb) - 1
        k = strcmp(T.Group, gg{1}) & strcmp(T.Context, 'ConB1') & ...
            T.HoursAfter >= eb(b) & T.HoursAfter < eb(b+1);
        if ~any(k), continue; end
        x = local_sesMean(T.fastAbove(k), T.Animal(k));
        fprintf('%-9s %4.1f-%4.1f %10.4f %8d\n', gg{1}, eb(b), eb(b+1), mean(x, 'omitnan'), numel(x));
    end
end

fprintf('\n--- sampling check (the contrast is only as good as this) ---\n');
fprintf('%-9s %-7s %9s %9s %9s %9s\n', 'group', 'context', 'ctx min', 'spikes', 'nUnits', 'nRip/win');
for gg = G
    for c = {'ConA1', 'ConB1'}
        k = strcmp(T.Group, gg{1}) & strcmp(T.Context, c{1});
        if ~any(k), continue; end
        fprintf('%-9s %-7s %9.1f %9.0f %9.0f %9.0f\n', gg{1}, c{1}, ...
            mean(T.ctxMin(k), 'omitnan'), mean(T.nSpkCtx(k), 'omitnan'), ...
            mean(T.nUnits(k), 'omitnan'), mean(T.nRip(k), 'omitnan'));
    end
end
fprintf('  cells actually used per session: %.0f (recorded: %.0f). Subtracting chance removes\n', ...
    mean(T.nCellUsed, 'omitnan'), mean(T.nUnits, 'omitnan'));
fprintf('  the chance LEVEL but not the yield dependence of the excess above it, so the counts\n');
fprintf('  must be equal across sessions, not merely accounted for. Check that nCellUsed is the\n');
fprintf('  same everywhere; if it is not, N_CELL_SUB was disabled or a session fell below it.\n');
fprintf('  Exploration duration and ripple count per window are NOT matched: if they differ\n');
fprintf('  between arms, add them as covariates before believing a group difference.\n');

%% ---- figure --------------------------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [70 70 1150 430]);
tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
col = containers.Map({'SALpre','ASV1','Control'}, {[0.2 0.6 1], [0.85 0.2 0.1], [0.45 0.45 0.45]});

nexttile(tl); hold on
xt = 0; lbl = {};
for gg = G
    for c = {'ConA1', 'ConB1'}
        k = strcmp(T.Group, gg{1}) & strcmp(T.Context, c{1});
        if ~any(k), continue; end
        xt = xt + 1;
        x = local_sesMean(T.fastAbove(k), T.Animal(k));
        plot(xt + 0.12*randn(size(x)), x, 'o', 'Color', col(gg{1}), 'MarkerSize', 4);
        plot(xt + [-0.25 0.25], [1 1]*mean(x, 'omitnan'), '-', 'LineWidth', 2, 'Color', col(gg{1}));
        lbl{end+1} = sprintf('%s\n%s', gg{1}, c{1}); %#ok<AGROW>
    end
end
yline(0, 'k:'); set(gca, 'XTick', 1:xt, 'XTickLabel', lbl);
ylabel('fastAbove (captured variance over chance)'); title('exploration contribution'); box off

nexttile(tl); hold on
if ~isempty(Dd)
    xt = 0; lbl = {};
    for gg = {'SALpre','ASV1'}
        k = strcmp(Dd.Group, gg{1});
        if ~any(k), continue; end
        xt = xt + 1;
        x = local_sesMean(Dd.d(k), Dd.Animal(k));
        plot(xt + 0.12*randn(size(x)), x, 'o', 'Color', col(gg{1}), 'MarkerSize', 5);
        plot(xt + [-0.25 0.25], [1 1]*mean(x, 'omitnan'), '-', 'LineWidth', 2, 'Color', col(gg{1}));
        lbl{end+1} = gg{1}; %#ok<AGROW>
    end
    yline(0, 'k:'); set(gca, 'XTick', 1:xt, 'XTickLabel', lbl);
end
ylabel('ConB1 - ConA1'); title('within-session drug contrast'); box off

nexttile(tl); hold on
for gg = G
    mu = nan(1, numel(eb)-1); se = nan(1, numel(eb)-1); ct = (eb(1:end-1) + eb(2:end))/2;
    for b = 1:numel(eb) - 1
        k = strcmp(T.Group, gg{1}) & strcmp(T.Context, 'ConB1') & ...
            T.HoursAfter >= eb(b) & T.HoursAfter < eb(b+1);
        if ~any(k), continue; end
        x = local_sesMean(T.fastAbove(k), T.Animal(k));
        mu(b) = mean(x, 'omitnan'); se(b) = std(x, 'omitnan') / max(sqrt(sum(isfinite(x))), 1);
    end
    errorbar(ct, mu, se, '-o', 'LineWidth', 1.4, 'Color', col(gg{1}), ...
        'MarkerFaceColor', col(gg{1}), 'DisplayName', gg{1});
end
yline(0, 'k:'); xlabel('h after ConB1'); ylabel('fastAbove');
title('time course'); legend('Location', 'best', 'Box', 'off'); box off
title(tl, 'Does ASV change what the exploration contributes to later ripples?');

if SAVE_OUT
    save(fullfile(CACHE_ROOT, 'driftInclusionASV.mat'), 'T', 'Dd', 'GEO_RANK', 'SUB_WIN', ...
        'TGT_MAX', 'BEH_BIN', 'INF_EXCL_WIN');
    fprintf('\nsaved driftInclusionASV.mat\n');
end


%% ---------------------------------------------------------------------------------------------
%  per-session worker
%% ---------------------------------------------------------------------------------------------
function r = local_session(fld, RP, IV, allCtx, HC, infOn, grp, animal, sname, ...
    subWin, guard, tgtMax, behBin, geoRank, minBin, minSpk, minRipSub, minUnits, freeGain, ...
    truncNext, nCellSub, nSubDraw, nRipMatch)
r = {};
[X, dur, tOn, nU, ts] = local_ripCounts(fld, RP);
if isempty(X) || nU < minUnits
    fprintf('%-34s skipped: %d PYR units\n', sname, nU); return
end
nDrop = 0; ripSeen = [];
nSub = nU; nDraw = 1;
if ~isempty(nCellSub)
    if nU < nCellSub
        fprintf('%-34s skipped: %d PYR units, below the matched count of %d\n', sname, nU, nCellSub); return
    end
    nSub = nCellSub; nDraw = nSubDraw;
end
inHC = false(size(tOn));
for q = 1:size(HC, 1), inHC = inHC | (tOn >= HC(q,1) & tOn < HC(q,2)); end

for w = 1:numel(IV)
    iv = IV(w).iv; nm = IV(w).name;
    drugged = isfinite(infOn) && iv(1) > infOn;        % encoded after the infusion
    [eC, ZC] = local_behVec(ts, iv, behBin);
    if isempty(eC) || size(ZC, 1) < minBin || sum(ZC(:)) < minSpk
        fprintf('%-34s %s: %d bins, %d spikes, no expert\n', sname, nm, size(ZC,1), round(sum(ZC(:))));
        continue
    end
    % one cell subset per draw, used for the basis AND for the target residuals, so the
    % expert and the thing it predicts always describe the same population
    cellSets = cell(1, nDraw); BCd = cell(1, nDraw);
    for d = 1:nDraw
        if nSub >= nU, cs = (1:nU)'; else, cs = randperm(nU, nSub)'; end
        cellSets{d} = cs;
        BCd{d} = local_binBasis(ZC(:, cs), geoRank);
    end
    keepD = ~cellfun(@isempty, BCd);
    cellSets = cellSets(keepD); BCd = BCd(keepD);
    if isempty(BCd), continue; end

    % targets: home cage after this context, stopping at the next context, and for a
    % pre-infusion context also stopping at the infusion so its replay window is undrugged
    tStop = iv(2) + tgtMax;
    if truncNext
        nx = allCtx(allCtx(:,1) > iv(2), 1);
        if ~isempty(nx), tStop = min(tStop, min(nx)); end
    end
    if isfinite(infOn) && iv(1) < infOn, tStop = min(tStop, infOn); end
    kPost = find(tOn >= iv(2) + guard & tOn < tStop & inHC);
    [~, tT, ripT] = local_subWinsIdx(X, dur, tOn, kPost, subWin, minRipSub);
    if isempty(ripT), continue; end

    for q = 1:numel(ripT)
        k = ripT{q};
        if ~isempty(nRipMatch)
            if numel(k) < nRipMatch
                nDrop = nDrop + 1; ripSeen(end+1) = numel(k); %#ok<AGROW>
                continue                                     % too few to match: drop the window
            end
            k = k(unique(round(linspace(1, numel(k), nRipMatch))));   % spread across the window
        end
        ripSeen(end+1) = numel(k); %#ok<AGROW>
        Xt = X(k,:); dt = dur(k);
        if size(Xt, 1) < 2 * minRipSub, continue; end
        vv = nan(1, numel(BCd)); rr = nan(1, numel(BCd));
        for d = 1:numel(BCd)
            [vv(d), rr(d)] = local_captOne(Xt(:, cellSets{d}), dt, BCd{d}, geoRank);
        end
        vC = mean(vv, 'omitnan'); vR = mean(rr, 'omitnan');
        if ~isfinite(vC), continue; end
        dCa = local_devG(Xt, dt, eC, freeGain);
        mt = sum(Xt, 1)' / max(sum(dt), eps);
        ora = local_poissDev(Xt, dt, mt);
        ha = 1:2:size(Xt,1); hb = 2:2:size(Xt,1);
        mA = sum(Xt(ha,:), 1)' / max(sum(dt(ha)), eps);
        mB = sum(Xt(hb,:), 1)' / max(sum(dt(hb)), eps);
        spl = (local_devG(Xt(hb,:), dt(hb), mA, freeGain) * numel(hb) + ...
               local_devG(Xt(ha,:), dt(ha), mB, freeGain) * numel(ha)) / size(Xt, 1);
        flr = max(spl - ora, eps);
        r(end+1,:) = {grp, char(string(animal)), sname, nm, (tT(q) - iv(2))/3600, ...
            numel(k), nU, vC, vR, vC - vR, (dCa - ora)/flr, flr, drugged, ...
            sum(ZC(:)), diff(iv)/60, nSub}; %#ok<AGROW>
    end
    fprintf('%-34s %-6s %s: %.0f min expert, %d bins, %d target windows%s\n', sname, grp, nm, ...
        diff(iv)/60, size(ZC,1), numel(ripT), local_tern(drugged, '  [drugged]', ''));
end
if nDrop > 0
    % A window below nRipMatch is discarded, and a session whose windows are ALL discarded
    % vanishes from the analysis entirely while still printing its window counts above. That
    % is how an animal silently loses its pair. The median available count says how far
    % nRipMatch would have to come down to keep it.
    fprintf('%-34s %-6s  nRipMatch dropped %d of %d windows (available ripples: median %d, max %d)\n', ...
        sname, grp, nDrop, numel(ripSeen), round(median(ripSeen)), max(ripSeen));
end
end


%% ---------------------------------------------------------------------------------------------
%  statistics
%% ---------------------------------------------------------------------------------------------
function local_groupTest(TT, fld)
% THIS IS A WITHIN-ANIMAL DESIGN. Most animals contribute both a SALpre and an ASV1 session,
% so animals cannot be split into two groups: the contrast is PAIRED inside each animal that has
% both, and the LME below uses every session including the animals that have only one arm.
if isempty(TT) || ~any(strcmp(TT.Group, 'ASV1')) || ~any(strcmp(TT.Group, 'SALpre'))
    fprintf('  not enough of both arms to test\n'); return
end
an = unique(TT.Animal);
pa = nan(numel(an), 1); sA = nan(numel(an), 1); sS = nan(numel(an), 1);
for q = 1:numel(an)
    kA = strcmp(TT.Animal, an{q}) & strcmp(TT.Group, 'ASV1');
    kS = strcmp(TT.Animal, an{q}) & strcmp(TT.Group, 'SALpre');
    if any(kA), sA(q) = mean(TT.(fld)(kA), 'omitnan'); end
    if any(kS), sS(q) = mean(TT.(fld)(kS), 'omitnan'); end
    if any(kA) && any(kS), pa(q) = sA(q) - sS(q); end
end
nPair = sum(isfinite(pa));
fprintf('  %s: %d animals with both arms, %d ASV-only, %d SAL-only\n', fld, nPair, ...
    sum(isfinite(sA) & ~isfinite(sS)), sum(isfinite(sS) & ~isfinite(sA)));
if nPair >= 2
    fprintf('    PAIRED within animal: ASV-SAL = %+.4f, p = %.4f (sign-flip, n = %d)\n', ...
        mean(pa, 'omitnan'), local_signflip(pa), nPair);
    fprintf('    per animal:');
    for q = find(isfinite(pa))'
        fprintf(' %s %+.3f;', an{q}, pa(q));
    end
    fprintf('\n');
else
    fprintf('    too few animals with both arms for the paired test\n');
end
if exist('fitlme', 'file') == 2
    try
        k2 = ismember(TT.Group, {'SALpre','ASV1'});
        hasDay = ismember('DayIdx', TT.Properties.VariableNames) && any(isfinite(TT.DayIdx(k2)));
        dv = zeros(sum(k2), 1);
        if hasDay, dv = TT.DayIdx(k2); dv(~isfinite(dv)) = 0; end
        D = table(TT.(fld)(k2), categorical(TT.Group(k2)), categorical(TT.Animal(k2)), ...
            categorical(TT.Session(k2)), dv, ...
            'VariableNames', {'y','Group','Animal','Session','DayIdx'});
        % SESSION MUST BE A RANDOM EFFECT. Target windows are ~20 per session and are not
        % independent of each other: a model with only (1|Animal) treats each window as a
        % replicate and returns a df in the hundreds for an experiment with ten sessions. That
        % is pseudoreplication and it inflates the test by an order of magnitude. Nesting
        % session inside animal puts the df back where the design puts it.
        lme = fitlme(D, 'y ~ Group + (1|Animal) + (1|Animal:Session)');
        % SATTERTHWAITE, not residual df. lme.Coefficients defaults to residual degrees of
        % freedom, which counts target windows: adding the nested random effect changes the
        % model but not that default, so the df stays in the hundreds and the test stays
        % pseudoreplicated. Satterthwaite takes the df from the variance components.
        try
            [~, ~, cf] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
        catch
            cf = lme.Coefficients;
        end
        j = find(~strcmp(cf.Name, '(Intercept)'), 1);
        % MATLAB orders the categorical alphabetically, so ASV1 is the reference and the
        % coefficient is SALpre minus ASV1. Printed as ASV-SAL to match the paired test above.
        fprintf('    LME y ~ Group + (1|Animal) + (1|Animal:Session): ASV-SAL = %+.4f, t = %.2f, df = %.1f, p = %.4f\n', ...
            -cf.Estimate(j), -cf.tStat(j), cf.DF(j), cf.pValue(j));
        fprintf('    [%s is the reference level; %d windows, %d sessions, %d animals]\n', ...
            char(cf.Name{j}), height(D), numel(unique(D.Session)), numel(unique(D.Animal)));
        fprintf('    (the LME keeps the single-arm animals; the paired test above does not)\n');
        % the same model with session order held fixed. If Group only matters because ASV
        % sessions come later, the coefficient collapses here and the arm difference is order.
        if hasDay && numel(unique(dv)) > 2
            try
                lme2 = fitlme(D, 'y ~ Group + DayIdx + (1|Animal) + (1|Animal:Session)');
                [~, ~, c2] = fixedEffects(lme2, 'DFMethod', 'satterthwaite');
                jg = find(startsWith(c2.Name, 'Group'), 1);
                jd = find(strcmp(c2.Name, 'DayIdx'), 1);
                fprintf('    + DayIdx covariate:  ASV-SAL = %+.4f, t = %.2f, df = %.1f, p = %.4f   [day: %+.4f/day, p = %.4f]\n', ...
                    -c2.Estimate(jg), -c2.tStat(jg), c2.DF(jg), c2.pValue(jg), ...
                    c2.Estimate(jd), c2.pValue(jd));
            catch err
                fprintf('    + DayIdx covariate not fitted (%s)\n', err.message);
            end
        end
    catch err
        fprintf('     LME not fitted (%s)\n', err.message);
    end
end
end


function D = local_within(T, fld, c1, c2)
% per session: value(c1) - value(c2), carrying the group and animal labels
s = unique(T.Session); rows = {};
for q = 1:numel(s)
    k1 = strcmp(T.Session, s{q}) & strcmp(T.Context, c1);
    k2 = strcmp(T.Session, s{q}) & strcmp(T.Context, c2);
    if ~any(k1) || ~any(k2), continue; end
    i1 = find(k1, 1);
    rows(end+1,:) = {T.Group{i1}, T.Animal{i1}, s{q}, ...
        mean(T.(fld)(k1), 'omitnan') - mean(T.(fld)(k2), 'omitnan')}; %#ok<AGROW>
end
if isempty(rows)
    D = table(cell(0,1), cell(0,1), cell(0,1), zeros(0,1), ...
        'VariableNames', {'Group','Animal','Session','d'});
else
    D = cell2table(rows, 'VariableNames', {'Group','Animal','Session','d'});
end
end


function [x, g] = local_sesMean(y, grp)
g = unique(grp); x = nan(numel(g), 1);
for q = 1:numel(g), x(q) = mean(y(strcmp(grp, g{q})), 'omitnan'); end
end


function p = local_signflip(x)
x = x(isfinite(x)); n = numel(x);
if n < 2, p = NaN; return; end
if n > 20
    s = sign(randn(n, 2000)); s(:,1) = 1;
else
    s = double(dec2bin(0:2^n - 1) - '0')' * 2 - 1;
end
m = abs(x(:)' * s / n);
p = max(mean(m >= abs(mean(x)) - 1e-12), 2 / 2^n);
end


function s = local_tern(c, a, b)
if c, s = a; else, s = b; end
end


%% ---------------------------------------------------------------------------------------------
%  measurement helpers (same definitions as driftInclusionTest.m)
%% ---------------------------------------------------------------------------------------------
function [v, vR] = local_captOne(Xt, dt, B, k)
% captured variance of a basis on the window's per-ripple residuals, and the EXACT chance level
% for the same dimension: a random k-dim basis captures k/(n-1) by symmetry.
v = NaN; vR = NaN;
if isempty(B) || size(Xt, 1) < 4, return; end
mt = sum(Xt, 1)' / max(sum(dt), eps);
rho = mt / max(norm(mt), eps);
R = Xt ./ max(dt(:), eps) - mt';
R = R - (R * rho) * rho';
vTot = sum(R(:).^2);
if vTot <= eps, return; end
Bo = local_orthTo(B, rho);
if isempty(Bo), return; end
kE = min(k, size(Bo, 2));
v = sum((R * Bo(:,1:kE)).^2, 'all') / vTot;
vR = kE / max(size(R, 2) - 1, 1);
end


function [lam, Z] = local_behVec(ts, iv, binw)
lam = []; Z = zeros(0, numel(ts));
if iv(2) <= iv(1), return; end
ed = iv(1):binw:iv(2);
if numel(ed) < 3, return; end
Z = zeros(numel(ed) - 1, numel(ts));
for u = 1:numel(ts), Z(:,u) = histcounts(ts{u}, ed)'; end
lam = sum(Z, 1)' / (iv(2) - iv(1));
end


function B = local_binBasis(Z, k)
B = [];
if size(Z, 1) < 5, return; end
m = mean(Z, 1)';
rho = m / max(norm(m), eps);
R = Z - m';
R = R - (R * rho) * rho';
B = local_topV(R, k);
end


function Bo = local_orthTo(B, rho)
Bo = [];
if isempty(B), return; end
B = B - rho * (rho' * B);
[Q, Rq] = qr(B, 0);
keep = abs(diag(Rq)) > 1e-8 * max(abs(diag(Rq)));
Bo = Q(:, keep);
end


function d = local_devG(X, dur, lam, useGain)
lam = max(lam(:), 0);
if useGain
    den = sum(dur(:)) * sum(lam);
    if den > eps, lam = max(sum(X(:)) / den, eps) * lam; end
end
d = local_poissDev(X, dur, lam);
end


function d = local_poissDev(X, dur, lam)
Lam = max(dur(:) * lam(:)', 1e-9);
T = X .* log(max(X, 1e-12) ./ Lam);
T(X == 0) = 0;
d = 2 * sum(T - (X - Lam), 'all') / size(X, 1);
end


function V = local_topV(R, k)
k = max(min([k, size(R, 2) - 1, size(R, 1) - 1]), 1);
[~, ~, V] = svd(R, 'econ');
V = V(:, 1:min(k, size(V, 2)));
end


function [M, tC, ripC] = local_subWinsIdx(X, dur, tOn, keep, subWin, minRip)
M = []; tC = []; ripC = {};
if isempty(keep), return; end
tk = tOn(keep);
ed = min(tk):subWin:max(tk) + subWin;
M = nan(size(X, 2), numel(ed) - 1); tC = nan(1, numel(ed) - 1); ripC = cell(1, numel(ed) - 1);
for b = 1:numel(ed) - 1
    k = keep(tk >= ed(b) & tk < ed(b+1));
    if numel(k) < minRip, continue; end
    M(:,b) = sum(X(k,:), 1)' / sum(dur(k));
    tC(b) = (ed(b) + ed(b+1)) / 2;
    ripC{b} = k;
end
ok = all(isfinite(M), 1) & ~cellfun(@isempty, ripC);
M = M(:, ok); tC = tC(ok); ripC = ripC(ok);
end


function [X, dur, tOn, nU, ts] = local_ripCounts(fld, RP)
X = []; dur = []; tOn = []; nU = 0; ts = {};
sf = dir(fullfile(fld, '*spikes.cellinfo.mat')); cf = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(sf) || isempty(cf), return; end
S = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); sp = S.spikes;
C = load(fullfile(cf(1).folder, cf(1).name), 'cell_metrics'); cm = C.cell_metrics;
ts = cellfun(@(x) sort(reshape(x(~isnan(x)), [], 1)), sp.times(:), 'UniformOutput', false);
n = numel(ts);
clu = []; if isfield(sp, 'UID') && numel(sp.UID) == n, clu = sp.UID(:); end
if isfield(sp, 'cluID') && numel(sp.cluID) == n, clu = sp.cluID(:); end
loc = matchCellMetrics(cm, (1:n)', clu, n); okL = loc > 0;
typ = repmat({''}, n, 1); ct = cellstr(string(cm.putativeCellType(:))); typ(okL) = ct(loc(okL));
idx = find(contains(typ, 'Pyramidal'));
if isempty(idx), ts = {}; return; end
on = RP.ts1(:); off = RP.ts1Off(:);
[on, o] = sort(on); off = off(o);
keep = off > on & [true; on(2:end) > off(1:end-1)];
on = on(keep); off = off(keep);
if numel(on) < 2, ts = {}; return; end
edges = reshape([on'; off'], [], 1);
X = zeros(numel(on), numel(idx));
for u = 1:numel(idx)
    c = histcounts(ts{idx(u)}, edges);
    X(:,u) = c(1:2:end)';
end
dur = off - on; tOn = on; nU = numel(idx); ts = ts(idx);
end


function d = local_sessionDate(sname)
% Recording date from the session folder name, as a datenum. Names look like CP12_250604_075757
% or CP5_230922: the first six-digit field after an underscore is YYMMDD.
d = NaN;
t = regexp(char(sname), '_(\d{6})(?:_|$)', 'tokens', 'once');
if isempty(t), return; end
x = t{1};
try
    d = datenum(2000 + str2double(x(1:2)), str2double(x(3:4)), str2double(x(5:6)));
catch
    d = NaN;
end
end


function s = local_name(fld)
[~, s] = fileparts(fld);
end


%% ---------------------------------------------------------------------------------------------
%  session markers: verbatim copies from compareReactivation.m so this script stands alone
%% ---------------------------------------------------------------------------------------------
function ev = local_events(fld)
ev = []; T = {}; X = []; src = {};
c = dir(fullfile(fld, '*.evt.sti'));
if ~isempty(c)
    try
        ws_ = warning('off', 'all');
        E = LoadEvents(fullfile(c(1).folder, c(1).name));
        warning(ws_);
        if isfield(E, 'time') && isfield(E, 'description')
            T = cellstr(E.description(:)); X = E.time(:); src{end+1} = c(1).name;
        end
    catch
        if exist('ws_', 'var'), warning(ws_); end
    end
end
[b, ~, ~] = cellfun(@local_parseLabel, T, 'UniformOutput', false);
if ~any(startsWith(b, 'Con'))
    [~, sn] = fileparts(fld);
    dirs = {fld, fileparts(fld), fileparts(fileparts(fld))}; found = false;
    for k = 1:numel(dirs)
        a = dir(fullfile(dirs{k}, '*transitions*.mat'));
        for q = 1:numel(a)
            try
                L = load(fullfile(a(q).folder, a(q).name), 'all_events'); A = L.all_events;
            catch
                continue
            end
            m = true(height(A), 1);
            if ismember('video', A.Properties.VariableNames)
                hit = contains(string(A.video), sn);
                if any(hit), m = hit; elseif k > 1, continue; end
            end
            T = [T(~startsWith(b, {'HC','Con'})); cellstr(string(A.event(m)))]; %#ok<AGROW>
            X = [X(~startsWith(b, {'HC','Con'})); double(A.Intan_ts(m))]; %#ok<AGROW>
            src{end+1} = fullfile(a(q).folder, a(q).name); %#ok<AGROW>
            found = true; break
        end
        if found, break; end
    end
end
if isempty(T), return; end
[X, o] = sort(X); T = T(o);
[b, num, edge] = cellfun(@local_parseLabel, T, 'UniformOutput', false);
keep = ~cellfun(@isempty, b);
unparsed = unique(T(~keep));
T = T(keep); X = X(keep); b = b(keep); num = num(keep); edge = edge(keep);
cnt = struct('ConA', 0, 'ConB', 0); lab = cell(size(T));
for k = 1:numel(T)
    if startsWith(b{k}, 'Con')
        if isempty(num{k})
            if strcmp(edge{k}, 'start'), cnt.(b{k}) = cnt.(b{k}) + 1; end
            num{k} = sprintf('%d', max(1, cnt.(b{k})));
        elseif strcmp(edge{k}, 'start')
            cnt.(b{k}) = str2double(num{k});
        end
        lab{k} = sprintf('%s%s_%s', b{k}, num{k}, edge{k});
    else
        lab{k} = sprintf('%s_%s', b{k}, edge{k});
    end
end
starts = find(endsWith(lab, '_start') & startsWith(lab, 'Con'));
for k = starts(:)'
    nm = erase(lab{k}, '_start');
    if ~any(strcmp(lab, [nm '_end']) & X > X(k))
        nx = find(X > X(k), 1);
        tEnd = X(k) + 1800; if ~isempty(nx), tEnd = X(nx); end
        lab{end+1, 1} = [nm '_end']; X(end+1, 1) = tEnd; %#ok<AGROW>
        unparsed{end+1, 1} = sprintf('(%s_end inferred)', nm); %#ok<AGROW>
    end
end
[X, o] = sort(X); lab = lab(o);
ev = struct('time', X, 'description', {lab}, 'source', strjoin(src, ' + '), 'unparsed', {unparsed});
end


function [b, num, edge] = local_parseLabel(s)
b = ''; num = ''; edge = '';
x = lower(strtrim(char(s)));
x = regexprep(x, '[\s_\-\.:]+', ' ');
t = regexp(x, ['^(con(?:text)? ?[ab]|familiar|novel|hc|home ?cage|infusion|inf) ?(\d*) ?' ...
               '(start|begin|onset|on|end|stop|offset|off)$'], 'tokens', 'once');
if isempty(t), return; end
w = t{1};
if startsWith(w, 'con'), b = ['Con' upper(w(end))];
elseif strcmp(w, 'familiar'), b = 'ConA';
elseif strcmp(w, 'novel'), b = 'ConB';
elseif startsWith(w, 'h'), b = 'HC';
else, b = 'Infusion';
end
num = t{2};
if any(strcmp(t{3}, {'start','begin','onset','on'})), edge = 'start'; else, edge = 'end'; end
end


function HC = local_hcBlocks(ev, RP, exclInf, infWin)
rec = NaN;
if isfield(RP, 'meta') && isfield(RP.meta, 'recordingSec'), rec = RP.meta.recordingSec; end
if ~isfinite(rec), rec = max([RP.ts1Off(:); ev.time(:)]) + 1; end
if nargin < 3, exclInf = false; end
if nargin < 4, infWin = [0 0]; end
lab = ev.description; t = ev.time;
ks = find(startsWith(lab, 'Con') & endsWith(lab, '_start'));
if exclInf, ks = [ks; find(strcmp(lab, 'Infusion_start'))]; end
C = zeros(0, 2);
for k = ks(:)'
    nm = erase(lab{k}, '_start');
    j = find(strcmp(lab, [nm '_end']) & t > t(k), 1);
    if ~isempty(j), C(end+1, :) = [t(k) t(j)]; end %#ok<AGROW>
end
if any(infWin > 0)
    ti = t(strcmp(lab, 'Infusion_start'));
    C = [C; ti(:) - infWin(1), ti(:) + infWin(2)];
end
C = sortrows(C);
HC = zeros(0, 2); t0 = 0;
for k = 1:size(C, 1)
    if C(k, 1) > t0, HC(end+1, :) = [t0 C(k, 1)]; end %#ok<AGROW>
    t0 = max(t0, C(k, 2));
end
if rec > t0, HC(end+1, :) = [t0 rec]; end
end


function g = local_group(RP, ev)
m = RP.meta; g = '';
inf_ = m.infusionSec;
if ~isfinite(inf_) && nargin > 1 && ~isempty(ev)
    k = find(strcmp(ev.description, 'Infusion_start'), 1);
    if isempty(k), k = find(strcmp(ev.description, 'ConA1_end'), 1); end
    if ~isempty(k), inf_ = ev.time(k); end
end
if ~isfinite(inf_), return; end
noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
off = RP.ts2(:) - RP.ts1(:); off = median(off(isfinite(off)));
if noASV || isempty(off) || ~isfinite(off), after = false; else, after = inf_ + off > 0; end
isFirst = ~noASV && strcmp(m.folder, m.firstASV.session);
switch m.condition
    case 'SAL'
        if ~after, g = 'SALpre'; end
    case 'ASV'
        if isFirst, g = 'ASV1'; end
end
end