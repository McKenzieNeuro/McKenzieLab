%% driftReverseControl.m
% Stage 2 of the drift-prediction model, on the never-infused control sessions only.
% Three controls, all built on the reverse-explained-variance idea: a measure that is driven by
% time-invariant structure, or by a monotone trend any smoother would capture, gives the same
% answer with time running backwards. Differencing the two directions removes it.
%
%   (1) FORECAST DIRECTION.  dDE = DE(A -> B) - DE(B -> A), two blocks of EQUAL width separated
%       by a gap. A directed flow predicts forward better than backward; a static smooth
%       gradient predicts both equally.
%   (2) GEOMETRY DIRECTION.  The fast-timescale content, i.e. the ripple-to-ripple covariance
%       subspace, is FORECAST rather than merely checked for persistence: a Grassmann geodesic
%       fitted to the subspace trajectory inside the fit block and extrapolated across the gap.
%       Scored against the frozen basis, and run in both directions.
%   (3) ACQUISITION DIRECTION.  Delta(tau, t) = D(L_pre(tau); t) - D(L_pre(tau) + e_tau; t):
%       what one expert fitted at tau buys for windows after tau, against a library that holds
%       only windows before it. Mirrored backwards. Adding any expert to any library helps a
%       little, and it helps both directions equally, so the forward-minus-backward difference
%       is what isolates content the session did not already contain.
%
% WHY THE BLOCKS ARE SYMMETRIC
% driftHorizonControl.m fits a 2 h anchor and tests a 30 min window. That design has no reverse:
% 30 min holds three sub-windows, too few for a flow. Here both blocks are BLK_WID wide, so
% either can fit and either can be the target, and the two directions share cells, ripple count
% (MATCH_RIP) and sub-window count.
%
% THE FREE GAIN, AND WHY IT IS NOT OPTIONAL HERE
% EV/REV is a correlation, so it is invariant to uniform rate scaling. Poisson deviance is not.
% These sessions have rate trends, and predicting a low-rate block from a high-rate one costs
% differently than the reverse, which yields dDE > 0 with no flow present at all. Every predictor
% therefore carries one free multiplicative gain, fitted on the target block in closed form
% (g = sum(X) / (sum(dur) * sum(lambda))). One parameter, and it restores exactly the invariance
% that makes dEV fair. Both versions are reported: de*_g is gain-free, de*_raw is not, and the
% difference between them is the size of the artefact the control removes.
% The gain is applied to every rung of the ladder, D_stat and D_split included, so DE measures
% prediction of PATTERN at both ends of its ratio. The oracle is unchanged by it (its ML gain
% is 1 by construction).
%
% WHY THE BASIS IS NOT PUSHED THROUGH THE MEAN-RATE FLOW
% The obvious forecast of the covariance subspace is B_hat = orth(L * B), with L = U A U' the
% linear part of the fitted DMD. It is wrong here: L has rank DMD_RANK in the subspace of the
% sub-window MEAN trajectory, which is not the subspace spanned by the per-ripple residual PCs,
% so L * B collapses a 9-dimensional basis onto 3 directions it knows nothing about. The flow
% fitted on mean rates carries no information about where the covariance is going. So the
% subspace gets its own model: bases from the first and second half of the fit block define a
% geodesic on the Grassmann manifold, extrapolated across the gap.
%
% UNIT OF DATA: per-ripple spike-count vectors, counts scaled by ripple duration, PYR only.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';  %#ok<NASGU>
CONTROLS = { ...
    'R:\HuszarR\e13_26m1_211019'; 'R:\HuszarR\e13_26m1_211116'; 'R:\HuszarR\e13_26m1_211119'; ...
    'R:\HuszarR\e15_13f1_220118'; 'R:\HuszarR\e15_13f1_220119'; 'R:\HuszarR\e15_13f1_220120'; ...
    'R:\HuszarR\e16_3m2_211211';  'R:\HuszarR\e16_3m2_211212';  'R:\HuszarR\e16_3m2_211213'};

BLK_WID     = 3600;            % s, width of BOTH blocks (fit and target swap roles)
SUB_WIN     = 600;             % s, sub-window for the slow trajectory inside a block
GAPS        = [0.5 1 2 3]*3600;   % s between the end of one block and the start of the other
K_RANK      = 10;              % 1 radial + K_RANK-1 covariance directions, for the rate ladder
GEO_RANK    = 4;               % covariance directions FORECAST by the geodesic. Deliberately
                               % smaller than K_RANK-1: a half-block holds a few hundred ripples,
                               % and at rank 9 the estimated subspace sits near the random-
                               % subspace angle, where the step and its sampling floor are both
                               % saturated and their difference is noise. The random-subspace
                               % angle is printed so this is visible rather than assumed.
TT_MAX      = 6;               % cap on the geodesic lever arm (see the tt computation below)
DMD_RANK    = 3;               % truncation; BLK_WID/SUB_WIN sub-windows cannot support more
CLAMP_DMD   = true;
MIN_RIP_SUB = 25;              % ripples a sub-window needs
MIN_RIP_BLK = 60;              % ripples a block needs
MIN_RIP_HLF = 30;              % ripples each half of a block needs for its own basis
MIN_SUB_BLK = 5;               % usable sub-windows a block needs
MIN_UNITS   = 20;
MATCH_RIP   = true;            % subsample both blocks to a common ripple count
FREE_GAIN   = true;            % report the gain-free ladder as the primary one
SAVE_OUT    = true;

% ---- acquisition scan: which epoch supplies training data the rest of the session needs ----
ACQ_ON       = true;
ACQ_GAP      = 1;              % sub-windows skipped either side of tau, against adjacency leakage
ACQ_MIN_SIDE = 3;              % windows needed on each side for tau to be scorable
ACQ_RIP_SUB  = 25;

rows = {};
ACQ = struct('name', {}, 'tau_h', {}, 'A', {}, 'fwd', {}, 'bwd', {}, 'ripRate', {});
for i = 1:numel(CONTROLS)
    fld = CONTROLS{i};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), warning('rev:file', 'missing %s', rf); continue; end
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES;
    [X, dur, tOn, nU] = local_ripCounts(fld, RP);
    if isempty(X) || nU < MIN_UNITS
        fprintf('%-34s skipped: %d PYR units\n', local_name(fld), nU); continue
    end
    t0 = min(tOn); rec = max(tOn) - t0;
    need = 2 * BLK_WID + max(GAPS);
    if rec < need
        fprintf('%-34s skipped: %.1f h recording, needs %.1f h\n', ...
            local_name(fld), rec/3600, need/3600); continue
    end

    % ---- symmetric block pairs, one per gap -------------------------------------------------
    for g = GAPS
        a0 = t0; a1 = a0 + BLK_WID;
        b0 = a1 + g; b1 = b0 + BLK_WID;
        if b1 > t0 + rec, continue; end
        kA = find(tOn >= a0 & tOn < a1);
        kB = find(tOn >= b0 & tOn < b1);
        if numel(kA) < MIN_RIP_BLK || numel(kB) < MIN_RIP_BLK, continue; end
        if MATCH_RIP
            nM = min(numel(kA), numel(kB));
            kA = kA(unique(round(linspace(1, numel(kA), nM))));
            kB = kB(unique(round(linspace(1, numel(kB), nM))));
        end

        BK = struct();
        [BK(1).X, BK(1).d, BK(1).t] = local_pick(X, dur, tOn, kA);
        [BK(2).X, BK(2).d, BK(2).t] = local_pick(X, dur, tOn, kB);
        for q = 1:2
            [BK(q).M, BK(q).tM] = local_subMeans(BK(q).X, BK(q).d, BK(q).t, SUB_WIN, MIN_RIP_SUB);
            BK(q).ctr = mean([min(BK(q).t), max(BK(q).t)]);
        end
        if size(BK(1).M, 2) < MIN_SUB_BLK || size(BK(2).M, 2) < MIN_SUB_BLK, continue; end

        % ---- (1) forecast, both directions ---------------------------------------------------
        % dir 1: fit block 1, predict block 2.  dir 2: fit block 2 with its sub-windows in
        % reverse order and time negated, predict block 1. The mirror is exact: same number of
        % extrapolation steps, same model class, same ripple count.
        DE = struct('g', nan(2, 3), 'raw', nan(2, 3), 'drift', nan(1, 2), 'floor', nan(1, 2));
        for d12 = 1:2
            fitB = BK(d12); tgtB = BK(3 - d12);
            sgn = 1; if d12 == 2, sgn = -1; end
            Mf = fitB.M; tf = fitB.tM;
            if d12 == 2, Mf = fliplr(Mf); tf = -fliplr(tf); end
            S = local_slowFit(Mf, tf, DMD_RANK, CLAMP_DMD);
            s = max(round((sgn * tgtB.ctr - tf(end)) / SUB_WIN), 1);
            pred = local_slowPred(S, s, sgn * tgtB.ctr);
            for useG = [true false]
                Lad = local_ladder(tgtB.X, tgtB.d, pred, useG);
                drift = Lad.stat - Lad.oracle;
                de = [(Lad.stat - Lad.dmd), (Lad.stat - Lad.ar1), (Lad.stat - Lad.ramp)] / ...
                     max(drift, eps);
                if useG
                    DE.g(d12,:) = de; DE.drift(d12) = drift; DE.floor(d12) = Lad.split - Lad.oracle;
                else
                    DE.raw(d12,:) = de;
                end
            end
        end

        % ---- (2) geometry, forecast and frozen, both directions -------------------------------
        GE = nan(2, 4);     % [forecast frozen oracle split] captured variance
        angR = nan(1, 2); angF = nan(1, 2); angN = nan(1, 2);   % step, sampling floor, saturation
        for d12 = 1:2
            fitB = BK(d12); tgtB = BK(3 - d12);
            nh = floor(size(fitB.X, 1) / 2);
            if nh < MIN_RIP_HLF, continue; end
            if d12 == 1
                i1 = 1:nh; i2 = nh+1:size(fitB.X, 1);          % early half, late half
            else
                i2 = 1:nh; i1 = nh+1:size(fitB.X, 1);          % reversed: late is "earlier"
            end
            B1 = local_blockBasis(fitB.X(i1,:), fitB.d(i1), GEO_RANK);
            B2 = local_blockBasis(fitB.X(i2,:), fitB.d(i2), GEO_RANK);
            if isempty(B1) || isempty(B2), continue; end
            % VELOCITY FLOOR. The geodesic is extrapolated to tt = 1 + (gap + BLK/2)/(BLK/2),
            % which is 4 at the shortest gap and 8 at the longest, so any angle between B1 and B2
            % is multiplied by that lever arm. Two independent samples of the SAME subspace are
            % already separated by a sampling angle, and extrapolating that angle is pure noise
            % amplification. The floor is estimated by splitting one half into odd and even
            % ripples, which holds the subspace fixed and the exposure matched, and is subtracted
            % from every principal angle before stepping - the same floor logic used for the rate
            % and variance endpoints elsewhere.
            oa = i1(1:2:end); ob = i1(2:2:end);
            Ba = local_blockBasis(fitB.X(oa,:), fitB.d(oa), GEO_RANK);
            Bb = local_blockBasis(fitB.X(ob,:), fitB.d(ob), GEO_RANK);
            thFlr = local_angles(Ba, Bb);
            angR(d12) = rad2deg(median(local_angles(B1, B2), 'omitnan'));
            angF(d12) = rad2deg(median(thFlr, 'omitnan'));
            angN(d12) = rad2deg(local_randAngle(nU, size(B1, 2)));
            tt = min(1 + (g + BLK_WID/2) / (BLK_WID/2), TT_MAX);
            Bf = local_grassStep(B1, B2, tt, median(thFlr, 'omitnan'));
            GE(d12,:) = local_geomScore(tgtB.X, tgtB.d, Bf, B2, GEO_RANK, MIN_RIP_HLF);
        end

        rrA = numel(kA) / BLK_WID * 3600; rrB = numel(kB) / BLK_WID * 3600;
        popA = sum(sum(BK(1).X, 1)' / sum(BK(1).d)); popB = sum(sum(BK(2).X, 1)' / sum(BK(2).d));
        rows(end+1,:) = {local_name(fld), nU, numel(kA), numel(kB), g/3600, ...
            DE.g(1,1), DE.g(2,1), DE.g(1,1) - DE.g(2,1), ...
            DE.raw(1,1), DE.raw(2,1), DE.raw(1,1) - DE.raw(2,1), ...
            DE.g(1,2), DE.g(2,2), DE.g(1,3), DE.g(2,3), ...
            DE.drift(1), DE.drift(2), DE.floor(1), DE.floor(2), ...
            GE(1,1), GE(1,2), GE(1,3), GE(1,4), ...
            GE(2,1), GE(2,2), GE(2,3), GE(2,4), angR(1), angF(1), angN(1), ...
            log10(max(popA,eps)), log10(max(popB,eps)/max(popA,eps)), ...
            log10(max(rrA,eps)), log10(max(rrB,eps)/max(rrA,eps))}; %#ok<AGROW>
    end

    % ---- (3) acquisition scan over the whole session ----------------------------------------
    if ACQ_ON
        [Ms, tCs, ripS] = local_subWins(X, dur, tOn, SUB_WIN, ACQ_RIP_SUB);
        ok = all(isfinite(Ms), 1) & ~cellfun(@isempty, ripS);
        Ms = Ms(:, ok); tCs = tCs(ok); ripS = ripS(ok);
        nS = size(Ms, 2);
        if nS >= 2*(ACQ_MIN_SIDE + ACQ_GAP) + 1
            [Dev, flr, rrS] = local_devMatrix(X, dur, Ms, ripS, FREE_GAIN, SUB_WIN);
            for tau = 1:nS
                pre  = 1:(tau - 1 - ACQ_GAP);
                post = (tau + 1 + ACQ_GAP):nS;
                if numel(pre) < ACQ_MIN_SIDE || numel(post) < ACQ_MIN_SIDE, continue; end
                % forward: library = pre, targets = post.  backward: library = post, targets = pre.
                fwd = local_delta(Dev, pre,  post, tau, flr);
                bwd = local_delta(Dev, post, pre,  tau, flr);
                ACQ(end+1) = struct('name', local_name(fld), 'tau_h', (tCs(tau) - t0)/3600, ...
                    'A', fwd - bwd, 'fwd', fwd, 'bwd', bwd, 'ripRate', rrS(tau)); %#ok<AGROW>
            end
        else
            fprintf('%-34s acquisition skipped: %d usable sub-windows\n', local_name(fld), nS);
        end
    end
end
if isempty(rows), error('rev:none', 'no control session produced a usable symmetric block pair'); end

T = cell2table(rows, 'VariableNames', {'Session','nUnits','nA','nB','Gap', ...
    'deFwd','deRev','dDE','deFwdRaw','deRevRaw','dDERaw', ...
    'ar1Fwd','ar1Rev','rampFwd','rampRev', ...
    'driftFwd','driftRev','floorFwd','floorRev', ...
    'gFcFwd','gFzFwd','gOrFwd','gSpFwd','gFcRev','gFzRev','gOrRev','gSpRev', ...
    'subAngle','subAngleFloor','subAngleRand', ...
    'logRateA','dLogRate','logRipRateA','dLogRipRate'});
T.gGainFwd = T.gFcFwd - T.gFzFwd;      % what forecasting the subspace buys over freezing it
T.gGainRev = T.gFcRev - T.gFzRev;
T.dGeom    = T.gGainFwd - T.gGainRev;  % and the direction control on that
A = struct2table(ACQ);

%% ---------------------------------------------------------------------------------------------
%  REPORT
%% ---------------------------------------------------------------------------------------------
fprintf('\n================ REVERSE-DIRECTION CONTROLS, CONTROL SESSIONS ================\n');
fprintf('%d sessions, %d block pairs, %g h blocks, rank %d, DMD rank %d, free gain %d\n\n', ...
    numel(unique(T.Session)), height(T), BLK_WID/3600, K_RANK, DMD_RANK, FREE_GAIN);

fprintf('--- (1) forecast direction: DE forward vs reverse (gain-free) ---\n');
fprintf('%6s %8s %8s %8s %10s %8s\n', 'gap(h)', 'DEfwd', 'DErev', 'dDE', 'p(signfl)', 'n');
gv = unique(T.Gap);
for q = 1:numel(gv)
    k = T.Gap == gv(q);
    x = local_animalMean(T.dDE(k), T.Session(k));
    fprintf('%6.1f %8.3f %8.3f %8.3f %10.4f %8d\n', gv(q), ...
        mean(T.deFwd(k), 'omitnan'), mean(T.deRev(k), 'omitnan'), ...
        mean(x, 'omitnan'), local_signflip(x), numel(x));
end
fprintf('\nsame, WITHOUT the free gain (the artefact the control removes):\n');
for q = 1:numel(gv)
    k = T.Gap == gv(q);
    fprintf('%6.1f %8.3f %8.3f %8.3f  (gain-free dDE %6.3f)\n', gv(q), ...
        mean(T.deFwdRaw(k), 'omitnan'), mean(T.deRevRaw(k), 'omitnan'), ...
        mean(T.dDERaw(k), 'omitnan'), mean(T.dDE(k), 'omitnan'));
end
fprintf('\nDE is reported only where the forward drift clears its split-half floor:\n');
okF = T.driftFwd > T.floorFwd; okR = T.driftRev > T.floorRev;
fprintf('  forward above floor %d/%d block pairs, reverse %d/%d\n', ...
    sum(okF), height(T), sum(okR), height(T));

fprintf('\n--- (2) geometry: does forecasting the subspace beat freezing it ---\n');
fprintf('%6s %9s %9s %9s %9s %9s %10s\n', 'gap(h)', 'fcast', 'frozen', 'gain', 'oracle', 'split', 'p(gain)');
for q = 1:numel(gv)
    k = T.Gap == gv(q);
    x = local_animalMean(T.gGainFwd(k), T.Session(k));
    fprintf('%6.1f %9.3f %9.3f %9.3f %9.3f %9.3f %10.4f\n', gv(q), ...
        mean(T.gFcFwd(k), 'omitnan'), mean(T.gFzFwd(k), 'omitnan'), mean(x, 'omitnan'), ...
        mean(T.gOrFwd(k), 'omitnan'), mean(T.gSpFwd(k), 'omitnan'), local_signflip(x));
end
fprintf('subspace step %.1f deg, sampling floor %.1f deg, random-subspace angle %.1f deg (medians)\n', ...
    median(T.subAngle, 'omitnan'), median(T.subAngleFloor, 'omitnan'), median(T.subAngleRand, 'omitnan'));
fprintf('  step at the floor: no velocity survives shrinkage, the forecast IS the frozen basis,\n');
fprintf('  which is the honest answer rather than a failure. Floor near the random-subspace\n');
fprintf('  angle: the subspace is not estimable at rank %d from a half block, so neither the\n', GEO_RANK);
fprintf('  step nor its floor means anything and GEO_RANK must come down before this is read.\n');
fprintf('direction control on the geometry forecast (forward gain minus reverse gain):\n');
xg = local_animalMean(T.dGeom, T.Session);
fprintf('  dGeom = %.4f, p = %.4f, n = %d\n', mean(xg, 'omitnan'), local_signflip(xg), numel(xg));
fprintf('  read: a subspace forecast that works equally well backwards is a static gradient,\n');
fprintf('        not a trajectory, however much it beats the frozen basis.\n');

if ~isempty(ACQ)
    fprintf('\n--- (3) acquisition: which epoch supplies needed training data ---\n');
    fprintf('%-34s %8s %8s %8s %8s %10s\n', 'session', 'tau(h)', 'A', 'fwd', 'bwd', 'rip/h');
    sn = unique(A.name);
    for q = 1:numel(sn)
        k = strcmp(A.name, sn{q});
        [~, j] = max(A.A(k));
        sub = A(k,:);
        fprintf('%-34s %8.2f %8.3f %8.3f %8.3f %10.1f\n', sn{q}, sub.tau_h(j), ...
            sub.A(j), sub.fwd(j), sub.bwd(j), sub.ripRate(j));
    end
    fprintf('\nA = (what an expert at tau buys for LATER windows) - (what it buys for EARLIER ones),\n');
    fprintf('both in units of the target window''s split-half floor. Adding any expert to any\n');
    fprintf('library helps, and it helps both directions alike, so A near 0 means the epoch holds\n');
    fprintf('nothing the session did not already contain. A > 1 means one epoch carries more than\n');
    fprintf('resampling noise for the rest of the session. The ripple rate is printed because a\n');
    fprintf('state change is also an acquisition in these terms.\n');
end

%% ---------------------------------------------------------------------------------------------
%  FIGURES
%% ---------------------------------------------------------------------------------------------
f1 = figure('Color', 'w', 'Position', [80 80 1180 420]);
tl = tiledlayout(f1, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile(tl); hold on
mf = arrayfun(@(h) mean(T.deFwd(T.Gap == h), 'omitnan'), gv);
mr = arrayfun(@(h) mean(T.deRev(T.Gap == h), 'omitnan'), gv);
plot(gv, mf, '-o', 'LineWidth', 1.4, 'Color', [0.15 0.35 0.70], 'DisplayName', 'forward');
plot(gv, mr, '-s', 'LineWidth', 1.4, 'Color', [0.75 0.35 0.15], 'DisplayName', 'reverse');
yline(0, 'k:'); xlabel('gap (h)'); ylabel('drift explained (DMD)');
title('forecast, both directions'); legend('Location', 'best'); box off

nexttile(tl); hold on
md = arrayfun(@(h) mean(T.dDE(T.Gap == h), 'omitnan'), gv);
mdr = arrayfun(@(h) mean(T.dDERaw(T.Gap == h), 'omitnan'), gv);
plot(gv, md, '-o', 'LineWidth', 1.6, 'Color', [0.15 0.35 0.70], 'DisplayName', 'gain-free');
plot(gv, mdr, '--^', 'LineWidth', 1.2, 'Color', [0.55 0.55 0.55], 'DisplayName', 'raw');
yline(0, 'k:'); xlabel('gap (h)'); ylabel('dDE = forward - reverse');
title('direction control'); legend('Location', 'best'); box off

nexttile(tl); hold on
gf = arrayfun(@(h) mean(T.gFcFwd(T.Gap == h), 'omitnan'), gv);
gz = arrayfun(@(h) mean(T.gFzFwd(T.Gap == h), 'omitnan'), gv);
go = arrayfun(@(h) mean(T.gOrFwd(T.Gap == h), 'omitnan'), gv);
gs = arrayfun(@(h) mean(T.gSpFwd(T.Gap == h), 'omitnan'), gv);
plot(gv, go, '-k', 'LineWidth', 1.0, 'DisplayName', 'oracle');
plot(gv, gs, ':k', 'LineWidth', 1.0, 'DisplayName', 'split floor');
plot(gv, gf, '-o', 'LineWidth', 1.5, 'Color', [0.15 0.35 0.70], 'DisplayName', 'forecast');
plot(gv, gz, '-s', 'LineWidth', 1.3, 'Color', [0.75 0.35 0.15], 'DisplayName', 'frozen');
xlabel('gap (h)'); ylabel('captured variance'); title('fast-content geometry');
legend('Location', 'best'); box off
title(tl, sprintf('Reverse-direction controls, %d control sessions', numel(unique(T.Session))));

if ~isempty(ACQ)
    f2 = figure('Color', 'w', 'Position', [80 540 1180 420]);
    t2 = tiledlayout(f2, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(t2); hold on
    sn = unique(A.name);
    for q = 1:numel(sn)
        k = strcmp(A.name, sn{q});
        plot(A.tau_h(k), A.A(k), '-', 'LineWidth', 1.0);
    end
    yline(0, 'k:'); yline(1, 'k--');
    xlabel('tau (h into session)'); ylabel('A(tau) = forward - backward gain');
    title('acquisition scan, one line per session'); box off
    nexttile(t2); hold on
    plot(A.ripRate, A.A, '.', 'MarkerSize', 8, 'Color', [0.35 0.35 0.35]);
    xlabel('ripple rate at tau (per h)'); ylabel('A(tau)');
    title('state confound: A against ripple rate'); box off
    if numel(A.A) > 3
        cc = corr(A.ripRate, A.A, 'rows', 'complete', 'type', 'Spearman');
        subtitle(sprintf('Spearman r = %.2f', cc));
    end
    title(t2, 'Does one epoch supply training data the rest of the session needs?');
end

if SAVE_OUT
    save(fullfile(pwd, 'driftReverseControl_out.mat'), 'T', 'A', ...
        'BLK_WID', 'SUB_WIN', 'GAPS', 'K_RANK', 'DMD_RANK', 'FREE_GAIN', ...
        'ACQ_GAP', 'ACQ_MIN_SIDE', 'MATCH_RIP');
    fprintf('\nsaved driftReverseControl_out.mat\n');
end


%% ---------------------------------------------------------------------------------------------
%  helpers
%% ---------------------------------------------------------------------------------------------
function [Xs, ds, ts] = local_pick(X, dur, tOn, k)
Xs = X(k,:); ds = dur(k); ts = tOn(k);
end


function [M, tM] = local_subMeans(X, dur, tOn, subWin, minRip)
% sub-window mean rate vectors inside one block, dropping thin sub-windows
ed = min(tOn):subWin:max(tOn) + subWin;
M = nan(size(X, 2), numel(ed) - 1); tM = nan(1, numel(ed) - 1);
for b = 1:numel(ed) - 1
    k = tOn >= ed(b) & tOn < ed(b+1);
    if sum(k) < minRip, continue; end
    M(:,b) = sum(X(k,:), 1)' / sum(dur(k));
    tM(b) = (ed(b) + ed(b+1)) / 2;
end
ok = all(isfinite(M), 1) & isfinite(tM);
M = M(:, ok); tM = tM(ok);
end


function S = local_slowFit(M, tM, r, clampIt)
% the three slow models, fitted on a sub-window trajectory given in its own time direction
S.M = M; S.tM = tM;
[S.A, S.U, S.mu] = local_fitDMD(M, r, clampIt);
[S.ar, S.arMu] = local_fitAR1(M);
S.ramp = local_fitRamp(M, tM);
end


function pred = local_slowPred(S, s, tTarget)
% s extrapolation steps ahead IN THE FITTED DIRECTION, tTarget in the fitted time coordinate
pred.stat  = S.M(:, end);
pred.block = mean(S.M, 2);
pred.dmd   = local_stepDMD(S.A, S.U, S.mu, S.M(:, end), s);
pred.ar1   = S.arMu + (S.ar .^ max(s, 0)) .* (S.M(:, end) - S.arMu);
pred.ramp  = S.ramp(:,1) + S.ramp(:,2) * tTarget;
fn = fieldnames(pred);
for q = 1:numel(fn), pred.(fn{q}) = max(pred.(fn{q}), 0); end
end


function Lad = local_ladder(XW, dW, pred, useGain)
% Poisson deviance ladder on the target block. Every rung that is a PREDICTION carries the free
% gain when useGain is true, so DE is a ratio of pattern errors at both ends. The oracle is
% invariant to it.
fn = fieldnames(pred);
for q = 1:numel(fn)
    Lad.(fn{q}) = local_devG(XW, dW, pred.(fn{q}), useGain);
end
mW = sum(XW, 1)' / sum(dW);
Lad.oracle = local_poissDev(XW, dW, mW);
ha = 1:2:size(XW,1); hb = 2:2:size(XW,1);
mA = sum(XW(ha,:), 1)' / max(sum(dW(ha)), eps);
mB = sum(XW(hb,:), 1)' / max(sum(dW(hb)), eps);
Lad.split = (local_devG(XW(hb,:), dW(hb), mA, useGain) * numel(hb) + ...
             local_devG(XW(ha,:), dW(ha), mB, useGain) * numel(ha)) / size(XW, 1);
end


function d = local_devG(X, dur, lam, useGain)
% Poisson deviance with one free multiplicative gain fitted on the target in closed form.
% lambda_ij = g * dur_i * lam_j, so the ML gain is total counts over total predicted exposure.
lam = max(lam(:), 0);
if useGain
    den = sum(dur(:)) * sum(lam);
    if den > eps
        g = sum(X(:)) / den;
        lam = max(g, eps) * lam;
    end
end
d = local_poissDev(X, dur, lam);
end


function B = local_blockBasis(X, dur, k)
% covariance directions of per-ripple residuals in one block, radial direction projected out
B = [];
if size(X, 1) < 3, return; end
m = sum(X, 1)' / max(sum(dur), eps);
rho = m / max(norm(m), eps);
R = X ./ max(dur(:), eps) - m';
R = R - (R * rho) * rho';
B = local_topV(R, k);
end


function th = local_angles(Ba, Bb)
% principal angles, in radians, between two subspaces of possibly different dimension
th = NaN;
if isempty(Ba) || isempty(Bb), return; end
k = min(size(Ba, 2), size(Bb, 2));
s = svd(Ba(:, 1:k)' * Bb(:, 1:k));
th = acos(min(max(s, -1), 1));
end


function th = local_randAngle(n, k)
% median principal angle between two independent random k-dim subspaces of R^n: the saturation
% reference. A sampling floor that reaches this means the subspace carries no information at
% rank k, and any step measured against it is an artefact of the rank, not of the data.
th = NaN;
if k < 1 || n <= k, return; end
[Qa, ~] = qr(randn(n, k), 0); [Qb, ~] = qr(randn(n, k), 0);
th = median(acos(min(max(svd(Qa' * Qb), -1), 1)));
end


function Bt = local_grassStep(B1, B2, tt, thFloor)
% Grassmann geodesic through B1 (tt = 0) and B2 (tt = 1), evaluated at tt, with each principal
% angle shrunk by thFloor first so that only motion exceeding the sampling angle is extrapolated.
% Velocity from the CS decomposition: Z = (I - B1 B1') B2 (B1' B2)^-1 = U tan(Theta) V'.
k = min(size(B1, 2), size(B2, 2));
B1 = B1(:, 1:k); B2 = B2(:, 1:k);
Mi = B1' * B2;
if k < 1 || rcond(Mi) < 1e-10, Bt = B2; return; end
Z = (B2 - B1 * Mi) / Mi;
[U, S, V] = svd(Z, 'econ');
th = atan(diag(S));
if nargin > 3 && isfinite(thFloor)
    th = max(th - thFloor, 0);          % a step at the sampling floor is not a step
end
Bt = B1 * V * diag(cos(th * tt)) * V' + U(:, 1:k) * diag(sin(th * tt)) * V';
[Bt, ~] = qr(Bt, 0);
Bt = Bt(:, 1:k);
end


function v = local_geomScore(XW, dW, Bfc, Bfz, k, minHalf)
% captured variance of the target block's per-ripple residuals by [forecast frozen oracle split].
% The target's OWN radial direction is removed from both the residuals and every candidate basis,
% so none of the four numbers can be inflated by a level difference between blocks.
v = nan(1, 4);
mW = sum(XW, 1)' / max(sum(dW), eps);
rho = mW / max(norm(mW), eps);
R = XW ./ max(dW(:), eps) - mW';
R = R - (R * rho) * rho';
vTot = sum(R(:).^2);
if vTot <= eps, return; end
Bor = local_topV(R, k);
ha = 1:2:size(R,1); hb = 2:2:size(R,1);
if numel(ha) < minHalf || numel(hb) < minHalf, return; end
Bsp = local_topV(R(ha,:), k);
cand = {local_orthTo(Bfc, rho), local_orthTo(Bfz, rho), local_orthTo(Bor, rho), local_orthTo(Bsp, rho)};
kEff = min(cellfun(@(b) size(b, 2), cand));
if kEff < 1, return; end
for q = 1:4
    Bq = cand{q}(:, 1:kEff);
    if q == 4
        v(q) = sum((R(hb,:) * Bq).^2, 'all') / max(sum(R(hb,:).^2, 'all'), eps);
    else
        v(q) = sum((R * Bq).^2, 'all') / vTot;
    end
end
end


function Bo = local_orthTo(B, rho)
% B with the rho component removed and re-orthonormalised, dropping any collapsed column
Bo = [];
if isempty(B), return; end
B = B - rho * (rho' * B);
[Q, Rq] = qr(B, 0);
keep = abs(diag(Rq)) > 1e-8 * max(abs(diag(Rq)));
Bo = Q(:, keep);
end


function [Dev, flr, rr] = local_devMatrix(X, dur, M, ripS, useGain, subWin)
% Dev(e,t) = deviance of expert e (a sub-window rate vector) on sub-window t.
% flr(t) = split-half floor of window t, in the same per-ripple units.
nS = size(M, 2);
Dev = nan(nS, nS); flr = nan(1, nS); rr = nan(1, nS);
for t = 1:nS
    k = ripS{t};
    Xt = X(k,:); dt = dur(k);
    for e = 1:nS
        Dev(e,t) = local_devG(Xt, dt, M(:,e), useGain);
    end
    mt = sum(Xt, 1)' / max(sum(dt), eps);
    ora = local_poissDev(Xt, dt, mt);
    ha = 1:2:size(Xt,1); hb = 2:2:size(Xt,1);
    if numel(ha) > 1 && numel(hb) > 1
        mA = sum(Xt(ha,:), 1)' / max(sum(dt(ha)), eps);
        mB = sum(Xt(hb,:), 1)' / max(sum(dt(hb)), eps);
        sp = (local_devG(Xt(hb,:), dt(hb), mA, useGain) * numel(hb) + ...
              local_devG(Xt(ha,:), dt(ha), mB, useGain) * numel(ha)) / size(Xt, 1);
        flr(t) = max(sp - ora, eps);
    end
    rr(t) = numel(k) / subWin * 3600;             % ripples per hour in this sub-window
end
end


function v = local_delta(Dev, lib, tgt, tau, flr)
% mean over target windows of what adding expert tau to library lib buys, in floor units.
% Nested by construction, so the gain is non-negative; the floor is what makes it readable.
d = nan(1, numel(tgt));
for q = 1:numel(tgt)
    t = tgt(q);
    base = min(Dev(lib, t));
    aug  = min(base, Dev(tau, t));
    if isfinite(flr(t)), d(q) = (base - aug) / flr(t); end
end
v = mean(d, 'omitnan');
end


function x = local_animalMean(y, grp)
% one value per session, so the sign-flip test has independent replicates
g = unique(grp);
x = nan(numel(g), 1);
for q = 1:numel(g)
    x(q) = mean(y(strcmp(grp, g{q})), 'omitnan');
end
end


function p = local_signflip(x)
% exact two-sided sign-flip test on the mean, floored at 2/2^n
x = x(isfinite(x)); n = numel(x);
if n < 2, p = NaN; return; end
if n > 20
    s = sign(randn(n, 2000)); s(:,1) = 1;
else
    s = double(dec2bin(0:2^n - 1) - '0')' * 2 - 1;
end
m = abs(x(:)' * s / n);                 % one flipped mean per sign pattern
p = max(mean(m >= abs(mean(x)) - 1e-12), 2 / 2^n);
end


function [M, tC, ripC] = local_subWins(X, dur, tOn, subWin, minRip)
% rate vector, centre time and ripple indices for every sub-window of the whole recording
ed = min(tOn):subWin:max(tOn) + subWin;
M = nan(size(X, 2), numel(ed) - 1); tC = nan(1, numel(ed) - 1); ripC = cell(1, numel(ed) - 1);
for b = 1:numel(ed) - 1
    k = find(tOn >= ed(b) & tOn < ed(b+1));
    if numel(k) < minRip, continue; end
    M(:,b) = sum(X(k,:), 1)' / sum(dur(k));
    tC(b) = (ed(b) + ed(b+1)) / 2;
    ripC{b} = k;
end
end


function [X, dur, tOn, nU] = local_ripCounts(fld, RP)
% per-ripple spike counts for pyramidal cells, from spike times and ripple edges
X = []; dur = []; tOn = []; nU = 0;
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
if isempty(idx), return; end
on = RP.ts1(:); off = RP.ts1Off(:);
[on, o] = sort(on); off = off(o);
keep = off > on & [true; on(2:end) > off(1:end-1)];
on = on(keep); off = off(keep);
if numel(on) < 2, return; end
edges = reshape([on'; off'], [], 1);
X = zeros(numel(on), numel(idx));
for u = 1:numel(idx)
    c = histcounts(ts{idx(u)}, edges);
    X(:,u) = c(1:2:end)';
end
dur = off - on; tOn = on; nU = numel(idx);
end


function V = local_topV(R, k)
k = max(min([k, size(R, 2) - 1, size(R, 1) - 1]), 1);
[~, ~, V] = svd(R, 'econ');
V = V(:, 1:min(k, size(V, 2)));
end


function [A, U, mu] = local_fitDMD(M, r, clampIt)
mu = mean(M, 2);
Mc = M - mu;
r = min([r, size(M, 2) - 1, size(M, 1)]);
[U, ~, ~] = svd(Mc, 'econ'); U = U(:, 1:max(r, 1));
Y = U' * Mc;
A = Y(:, 2:end) / Y(:, 1:end-1);
if clampIt
    [W, D] = eig(A); d = diag(D);
    bad = abs(d) > 1;
    if any(bad) && rank(W) == size(W, 1)
        d(bad) = d(bad) ./ abs(d(bad));
        A = real(W * diag(d) / W);
    end
end
end


function m = local_stepDMD(A, U, mu, mLast, s)
y = U' * (mLast - mu);
m = mu + U * (A^max(s, 0) * y);
end


function [a, mu] = local_fitAR1(M)
mu = mean(M, 2);
Mc = M - mu;
num = sum(Mc(:, 1:end-1) .* Mc(:, 2:end), 2);
den = sum(Mc(:, 1:end-1).^2, 2);
a = num ./ max(den, eps);
a = max(min(a, 0.999), -0.999);
end


function B = local_fitRamp(M, t)
G = [ones(numel(t), 1), t(:)];
B = (G \ M')';
end


function d = local_poissDev(X, dur, lam)
Lam = max(dur(:) * lam(:)', 1e-9);
T = X .* log(max(X, 1e-12) ./ Lam);
T(X == 0) = 0;
d = 2 * sum(T - (X - Lam), 'all') / size(X, 1);
end


function s = local_name(fld)
[~, s] = fileparts(fld);
end
