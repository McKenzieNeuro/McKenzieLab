%% driftNoveltyScan.m
% When does something enter the population that no earlier model expresses?
%
% Two traces are computed for every sub-window t of a session, both from a library fitted
% STRICTLY on sub-windows before t (anything else leaks the answer into the question):
%
%   N(t)  model novelty.  Best past expert against the window's own oracle, scaled by the
%         window's split-half floor:
%             N(t) = [ min_e D_t(expert_e) - D_t(oracle) ] / [ D_t(split) - D_t(oracle) ]
%         N ~ 1: the best the past offers is as good as resampling the window itself.
%         N >> 1: no past model expresses this window.
%
%   B(t)  basis novelty.  How much of the window's residual variance the best past basis
%         captures, against what an independent contemporaneous basis captures:
%             B(t) = v_selfCV(t) - max_p v_pastBasis(p -> t)
%         B > 0: the population is in directions no earlier basis reaches.
%
% N and B dissociate. N up with B flat is the same directions under a new flow; B up with N
% flat is a new direction that the old dynamics could still carry.
%
% DECLARING AN EPOCH is not enough on its own, because a one-off anomaly also scores high.
% The second criterion is prequential: fit an expert on the candidate epoch, add it to the
% library, and ask whether LATER windows are predicted better with it than without. Transient
% surprise fails that; a new model that has joined the set passes it.
%
% CALIBRATION, TWO NULLS, ANSWERING TWO DIFFERENT QUESTIONS.
%   Order permutation keeps each window's own statistics and destroys all temporal structure, so
%   it asks: is there ANY structure here? The null is the distribution of the maximum over time,
%   since that is the statistic being read off the trace. This sets the within-session threshold.
%   Time reversal runs the whole scan on the session played backwards. It preserves the
%   autocorrelation, the drift magnitude and every window's content, and flips only the direction,
%   so it asks the sharper question: is this an ACQUISITION rather than merely a change? A window
%   that is novel because something arrived and stayed is novel going forwards and unremarkable
%   going backwards. One reversal gives one number per session, so the test is the paired
%   forward-minus-reversed maximum across sessions, by sign-flip. This is the same logic as
%   reverse explained variance: difference the two directions and what is symmetric cancels.
%
% WHAT THIS CANNOT SEPARATE. Entering a new brain state is, in these terms, a genuinely new
% model. Novelty that coincides with a change in ripple rate is a state event until independent
% state labels say otherwise; the ripple rate is reported alongside every epoch for that reason.
%
% Helpers are local copies of those in driftHorizonControl.m so this script runs standalone.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
CONTROLS = { ...
    'R:\HuszarR\e13_26m1_211019'; 'R:\HuszarR\e13_26m1_211116'; 'R:\HuszarR\e13_26m1_211119'; ...
    'R:\HuszarR\e15_13f1_220118'; 'R:\HuszarR\e15_13f1_220119'; 'R:\HuszarR\e15_13f1_220120'; ...
    'R:\HuszarR\e16_3m2_211211';  'R:\HuszarR\e16_3m2_211212';  'R:\HuszarR\e16_3m2_211213'};

% Context exploration intervals, for marking on the figure ONLY: nothing in the scan uses them,
% which is the point of an unsupervised detector. They are drawn so that the traces can be read
% against the times the design actually prescribes, instead of against nothing.
% Columns: folder, ConB1 (novel), ConA1 (familiar). Same values as compareReactivation.m.
CTXMARK = { ...
    'e13_26m1_211019', [15149.6 16749.2], [7605.7 9322]; ...
    'e13_26m1_211116', [13608 13939],     [7469.8 8707.2]; ...
    'e13_26m1_211119', [14911.9 16446],   [6181.3 7572]; ...
    'e15_13f1_220118', [13137 15046],     [5636.7 7700]; ...
    'e15_13f1_220119', [21219 22700],     [6836.6 8860]; ...
    'e15_13f1_220120', [18017.9 19801],   [9610.5 11950]; ...
    'e16_3m2_211211',  [9584.6 11203],    [5294.1 6511]; ...
    'e16_3m2_211212',  [16416.3 17731],   [6109.3 7845]; ...
    'e16_3m2_211213',  [13291 15323],     [7567.9 9359]};

SCAN_WIN    = 600;        % s, sub-window width: also the resolution of an onset
MIN_RIP_SUB = 25;         % ripples a sub-window needs
MIN_UNITS   = 20;
K_RANK      = 10;         % rank of each basis (1 radial + K_RANK-1 residual PCs)
WARMUP      = 6;          % sub-windows of library before scanning starts
MAX_BASIS   = 25;         % most recent past bases compared for B(t) (cost control)
MIN_RUN     = 2;          % consecutive windows above threshold to declare an epoch
SURR_N      = 10;         % order-permutation surrogates; 0 skips calibration
SURR_PCT    = 95;         % percentile of the surrogate maximum used as the threshold
SURR_REV    = true;       % also scan the session time-reversed, as the directedness control
DMD_RANK    = 5;
SEED        = 11;
SAVE_OUT    = true;

rng(SEED, 'twister');
scanRows = {}; epoRows = {}; revRows = {};
for i = 1:numel(CONTROLS)
    fld = CONTROLS{i};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), warning('nov:file', 'missing %s', rf); continue; end
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES;
    [X, dur, tOn, nU] = local_ripCounts(fld, RP);
    if isempty(X) || nU < MIN_UNITS
        fprintf('%-34s skipped: %d PYR units\n', local_name(fld), nU); continue
    end
    [M, tC, ripC] = local_subWins(X, dur, tOn, SCAN_WIN, MIN_RIP_SUB);
    nT = size(M, 2);
    if nT < WARMUP + 4
        fprintf('%-34s skipped: %d usable sub-windows\n', local_name(fld), nT); continue
    end
    opt = struct('kRank', K_RANK, 'warmup', WARMUP, 'maxBasis', MAX_BASIS, 'dmdRank', DMD_RANK, ...
        'subWin', SCAN_WIN);

    % ---- the observed scan, in recording order -------------------------------------------
    [Nt, Bt, dMin] = local_scan(X, dur, M, tC, ripC, 1:nT, opt);

    % ---- surrogates: the same windows, their order permuted -------------------------------
    sN = nan(SURR_N, 1); sB = nan(SURR_N, 1);
    for q = 1:SURR_N
        o = randperm(nT);
        [n2, b2] = local_scan(X, dur, M(:, o), tC(o), ripC(o), 1:nT, opt);
        sN(q) = max(n2, [], 'omitnan'); sB(q) = max(b2, [], 'omitnan');
    end
    thrN = inf; thrB = inf;
    if SURR_N > 0
        thrN = prctile(sN, SURR_PCT); thrB = prctile(sB, SURR_PCT);
    end
    % POSITION-MATCHED THRESHOLD. The global surrogate maximum above is biased upward, badly.
    % A permutation reshuffles which window lands at each position but every position still
    % scores against "everything before it", so whichever window falls just past WARMUP is
    % judged by a library of WARMUP experts. If any window in the session is idiosyncratic,
    % some permutation drops it there and the surrogate maximum records a novelty the observed
    % trace can never produce, because in recording order that window has its real library.
    % The library-size effect is a function of POSITION, so the null is taken per position and
    % the observed trace is standardised against it before any maximum is read.
    if SURR_N >= 5
        PN = nan(SURR_N, nT); PB = nan(SURR_N, nT);
        for q = 1:SURR_N
            o = randperm(nT);
            [PN(q,:), PB(q,:)] = local_scan(X, dur, M(:, o), tC(o), ripC(o), 1:nT, opt);
        end
        zN = local_posZ(Nt, PN); zB = local_posZ(Bt, PB);
        thrZ = prctile(max(local_posZloo(PN), [], 2, 'omitnan'), SURR_PCT);
        thrZB = prctile(max(local_posZloo(PB), [], 2, 'omitnan'), SURR_PCT);
    else
        zN = nan(size(Nt)); zB = nan(size(Bt)); thrZ = inf; thrZB = inf;
    end

    % ---- the directedness control: the same session played backwards ----------------------
    % Only the direction changes. The window contents, the autocorrelation and the total drift
    % are all preserved, so a difference between the two maxima cannot be explained by any of
    % them. Times are negated as well as reversed so the ramp expert still sees increasing time.
    if SURR_REV
        oR = nT:-1:1;
        [nR, bR] = local_scan(X, dur, M(:, oR), -tC(oR), ripC(oR), 1:nT, opt);
        revRows(end+1,:) = {local_name(fld), max(Nt, [], 'omitnan'), max(nR, [], 'omitnan'), ...
            max(Bt, [], 'omitnan'), max(bR, [], 'omitnan')}; %#ok<AGROW>
    end

    % ---- declare epochs ------------------------------------------------------------------
    rr = cellfun(@numel, ripC) / SCAN_WIN * 3600;          % ripples per hour: the state proxy
    if isfinite(thrZ)
        epN = local_runs(zN > thrZ, MIN_RUN);              % position-matched, the primary rule
        epB = local_runs(zB > thrZB, MIN_RUN);
    else
        epN = local_runs(Nt > thrN, MIN_RUN);
        epB = local_runs(Bt > thrB, MIN_RUN);
    end
    for e = [{epN, 'model'}; {epB, 'basis'}]'
        ep = e{1}; kind = e{2};
        for q = 1:size(ep, 1)
            w = ep(q,1):ep(q,2);
            if ep(q,2) >= nT - 1, continue; end             % no future left to test the payoff
            [gain, frac] = local_payoff(X, dur, M, tC, ripC, w, opt);
            dRR = mean(rr(w)) / max(mean(rr(max(1, ep(q,1)-3):ep(q,1)-1)), eps);
            epoRows(end+1,:) = {local_name(fld), kind, tC(ep(q,1))/3600, numel(w), ...
                max(Nt(w)), max(Bt(w)), gain, frac, mean(rr(w)), dRR}; %#ok<AGROW>
        end
    end
    for t = 1:nT
        scanRows(end+1,:) = {local_name(fld), tC(t)/3600, numel(ripC{t}), Nt(t), Bt(t), ...
            dMin(t), rr(t), thrN, thrB, zN(t), zB(t), thrZ, thrZB}; %#ok<AGROW>
    end
    fprintf('%-34s %3d sub-win, global thr N %.2f (max %.2f), position-matched z thr %.2f (max %.2f), epochs: %d model, %d basis\n', ...
        local_name(fld), nT, thrN, max(Nt, [], 'omitnan'), thrZ, max(zN, [], 'omitnan'), ...
        size(epN, 1), size(epB, 1));
end
if isempty(scanRows), error('nov:none', 'no session produced a usable scan'); end
SC = cell2table(scanRows, 'VariableNames', {'Session','Hour','nRip','N','B','dMin','ripPerHour', ...
    'thrN','thrB','zN','zB','thrZ','thrZB'});
if isempty(epoRows), epoRows = cell(0, 10); end
EP = cell2table(epoRows, 'VariableNames', {'Session','Kind','Hour','nWin','maxN','maxB', ...
    'payoffDev','payoffFrac','ripPerHour','ripRateRatio'});

%% ---- report -----------------------------------------------------------------------------
fprintf('\n%d sessions scanned at %g min resolution\n', numel(unique(SC.Session)), SCAN_WIN/60);
fprintf('\nEPOCHS  (payoffDev > 0 means later windows are predicted better once the epoch''s\n');
fprintf('expert joins the library; payoffFrac is the share of later windows improved)\n');
fprintf('%-30s %-6s %7s %5s %7s %8s %10s %9s %9s\n', 'session', 'kind', 'hour', 'nWin', ...
    'maxN', 'maxB', 'payoffDev', 'payoffFrac', 'ripRatio');
for q = 1:height(EP)
    fprintf('%-30s %-6s %7.2f %5d %7.2f %+8.3f %+10.4f %9.2f %9.2f\n', EP.Session{q}, EP.Kind{q}, ...
        EP.Hour(q), EP.nWin(q), EP.maxN(q), EP.maxB(q), EP.payoffDev(q), EP.payoffFrac(q), EP.ripRateRatio(q));
end
if isempty(epoRows)
    fprintf('   none: no window exceeded its surrogate threshold for %d consecutive sub-windows\n', MIN_RUN);
else
    kept = EP.payoffDev > 0 & EP.payoffFrac > 0.5;
    fprintf('\n%d of %d epochs pay forward; %d of those sit within 20%% of the local ripple rate\n', ...
        sum(kept), height(EP), sum(kept & abs(EP.ripRateRatio - 1) < 0.2));
    fprintf('(an epoch whose ripple rate also jumps is a state event until a scorer says otherwise)\n');
end
k = isfinite(SC.N) & isfinite(SC.ripPerHour);
if sum(k) > 10
    fprintf('\nstate check: corr(N, ripples/h) = %+.2f, corr(B, ripples/h) = %+.2f over %d windows\n', ...
        corr(SC.N(k), SC.ripPerHour(k)), corr(SC.B(k), SC.ripPerHour(k)), sum(k));
end

if ~isempty(revRows)
    RV = cell2table(revRows, 'VariableNames', {'Session','maxNfwd','maxNrev','maxBfwd','maxBrev'});
    dN = RV.maxNfwd - RV.maxNrev; dB = RV.maxBfwd - RV.maxBrev;
    fprintf('\nDIRECTEDNESS (forward vs the same session reversed, paired across sessions)\n');
    fprintf('%-30s %8s %8s %8s %8s\n', 'session', 'Nfwd', 'Nrev', 'Bfwd', 'Brev');
    for q = 1:height(RV)
        fprintf('%-30s %8.2f %8.2f %+8.3f %+8.3f\n', RV.Session{q}, ...
            RV.maxNfwd(q), RV.maxNrev(q), RV.maxBfwd(q), RV.maxBrev(q));
    end
    fprintf('  max N: forward - reversed = %+.3f, p = %.4f (n = %d)\n', ...
        mean(dN, 'omitnan'), local_signflip(dN), sum(isfinite(dN)));
    fprintf('  max B: forward - reversed = %+.3f, p = %.4f (n = %d)\n', ...
        mean(dB, 'omitnan'), local_signflip(dB), sum(isfinite(dB)));
    fprintf('  a difference near zero means the trace marks a change, not an acquisition: the\n');
    fprintf('  window is equally inexpressible from either side, which is what a transition\n');
    fprintf('  between two long-lived states looks like.\n');
else
    RV = cell2table(cell(0, 5), 'VariableNames', {'Session','maxNfwd','maxNrev','maxBfwd','maxBrev'});
end

%% ---- figure ------------------------------------------------------------------------------
ses = unique(SC.Session); ns = numel(ses);
fig = figure('Color', 'w', 'Position', [60 60 1250 210*ns]);
tl = tiledlayout(ns, 1, 'TileSpacing', 'compact');
for q = 1:ns
    ax = nexttile; hold on
    k = strcmp(SC.Session, ses{q});

    % contexts first, so the traces draw over them. ConA1 familiar, ConB1 novel.
    jm = find(strcmp(CTXMARK(:,1), ses{q}), 1);
    if ~isempty(jm)
        yl = [0 1];
        for cc = [3 2]                               % 3 = familiar, 2 = novel
            iv = CTXMARK{jm, cc} / 3600;
            if cc == 2, fc = [0.95 0.80 0.80]; lb = 'ConB1'; else, fc = [0.85 0.88 0.95]; lb = 'ConA1'; end
            patch(ax, [iv(1) iv(2) iv(2) iv(1)], [yl(1) yl(1) yl(2) yl(2)], fc, ...
                'EdgeColor', 'none', 'FaceAlpha', 0.75, 'HandleVisibility', 'off');
            text(mean(iv), yl(2), lb, 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
        end
        ax.YLimMode = 'auto';
    end

    yyaxis left
    plot(SC.Hour(k), SC.N(k), '-o', 'LineWidth', 1.4, 'MarkerSize', 3);
    plot(SC.Hour(k), SC.B(k) * 10, '-s', 'LineWidth', 1.1, 'MarkerSize', 3);
    yline(SC.thrN(find(k, 1)), 'k--');
    if any(isfinite(SC.zN(k)))
        yyaxis left
        plot(SC.Hour(k), SC.zN(k), '-', 'Color', [0.1 0.5 0.2], 'LineWidth', 1.2);
        yline(SC.thrZ(find(k, 1)), '--', 'Color', [0.1 0.5 0.2]);
    end
    ylabel('N  (B x10, z)');
    yyaxis right
    plot(SC.Hour(k), SC.ripPerHour(k), '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
    ylabel('ripples / h');
    ke = strcmp(EP.Session, ses{q});
    for e = find(ke)'
        xline(EP.Hour(e), ':', EP.Kind{e}, 'LabelVerticalAlignment', 'bottom');
    end
    xlabel('h into the recording'); title(ses{q}, 'FontWeight', 'normal', 'Interpreter', 'none');
    if q == 1, legend({'model novelty N', 'basis novelty B (x10)', 'global threshold', ...
            'position-matched z', 'z threshold', 'ripple rate'}, 'Location', 'northwest', 'Box', 'off'); end
end
title(tl, 'Novelty scan, with the prescribed context times marked (shading is not used by the scan)');
if SAVE_OUT
    exportgraphics(fig, fullfile(CACHE_ROOT, 'driftNoveltyScan.pdf'), 'ContentType', 'vector');
    save(fullfile(CACHE_ROOT, 'driftNoveltyScan.mat'), 'SC', 'EP', 'RV', 'SCAN_WIN', 'K_RANK', ...
        'WARMUP', 'MIN_RUN', 'SURR_N', 'SURR_PCT', 'SURR_REV');
end


%% =========================================================================
function [Nt, Bt, dMin] = local_scan(X, dur, M, tC, ripC, ord, o)
% causal scan: for every window t, a library fitted only on windows before it
nT = numel(ord);
Nt = nan(1, nT); Bt = nan(1, nT); dMin = nan(1, nT);
Bas = cell(1, nT);                                          % each window's basis, computed once
for p = 1:nT
    rp = ripC{p};
    if numel(rp) < 5, continue; end
    Xp = X(rp,:); dp = dur(rp); mp = sum(Xp, 1)' / sum(dp);
    Rp = Xp ./ max(dp, eps) - mp';
    rhp = mp / max(norm(mp), eps);
    Bas{p} = local_basis(Rp - (Rp * rhp) * rhp', o.kRank);
end
for t = (o.warmup + 1):nT
    past = 1:(t - 1);
    r = ripC{t};
    if numel(r) < 5, continue; end
    XW = X(r,:); dW = dur(r);
    mW = sum(XW, 1)' / sum(dW);

    % --- model novelty -----------------------------------------------------------------
    dBest = inf;
    for p = past                                            % every past window as a memory expert
        dBest = min(dBest, local_poissDev(XW, dW, M(:, p)));
    end
    E = local_slowExperts(M(:, past), tC(past), o);         % and the parametric slow models
    for q = 1:numel(E)
        dBest = min(dBest, local_poissDev(XW, dW, local_predictAt(E(q), tC(t))));
    end
    dOra = local_poissDev(XW, dW, mW);
    ha = 1:2:size(XW,1); hb = 2:2:size(XW,1);
    if numel(ha) < 2 || numel(hb) < 2, continue; end
    mA = sum(XW(ha,:), 1)' / sum(dW(ha)); mB = sum(XW(hb,:), 1)' / sum(dW(hb));
    dSpl = (local_poissDev(XW(hb,:), dW(hb), mA) * numel(hb) + ...
            local_poissDev(XW(ha,:), dW(ha), mB) * numel(ha)) / size(XW, 1);
    Nt(t) = (dBest - dOra) / max(dSpl - dOra, eps);
    dMin(t) = dBest;

    % --- basis novelty -----------------------------------------------------------------
    R = XW ./ max(dW, eps) - mW';
    rho = mW / max(norm(mW), eps);
    Rr = R - (R * rho) * rho';
    vTot = max(sum(Rr(:).^2), eps);
    Va = local_basis(Rr(ha,:), o.kRank);                    % contemporaneous, cross-validated
    vSelf = sum((Rr(hb,:) * Va).^2, 'all') / max(sum(Rr(hb,:).^2, 'all'), eps);
    vPast = 0;
    pb = past(max(1, numel(past) - o.maxBasis + 1):end);
    for p = pb
        if isempty(Bas{p}), continue; end
        kE = min(size(Bas{p}, 2), size(Va, 2));
        vPast = max(vPast, sum((Rr * Bas{p}(:, 1:kE)).^2, 'all') / vTot);
    end
    Bt(t) = vSelf - vPast;
end
end


function [gain, frac] = local_payoff(X, dur, M, tC, ripC, w, o)
% does adding an expert fitted on the epoch improve prediction of LATER windows?
gain = NaN; frac = NaN;
nT = size(M, 2);
fut = (w(end) + 1):nT;
if isempty(fut), return; end
mNew = mean(M(:, w), 2);
d0 = nan(1, numel(fut)); d1 = d0;
for q = 1:numel(fut)
    t = fut(q); r = ripC{t};
    if numel(r) < 5, continue; end
    XW = X(r,:); dW = dur(r);
    base = inf;
    for p = 1:(w(1) - 1), base = min(base, local_poissDev(XW, dW, M(:, p))); end
    E = local_slowExperts(M(:, 1:(w(1) - 1)), tC(1:(w(1) - 1)), o);
    for e = 1:numel(E), base = min(base, local_poissDev(XW, dW, local_predictAt(E(e), tC(t)))); end
    d0(q) = base;
    d1(q) = min(base, local_poissDev(XW, dW, mNew));
end
k = isfinite(d0) & isfinite(d1);
if ~any(k), return; end
gain = mean(d0(k) - d1(k));                                  % positive = better with the new expert
frac = mean(d1(k) < d0(k) - 1e-12);
end


function E = local_slowExperts(M, t, o)
% the parametric slow models of driftHorizonControl, fitted on the past only
E = struct('kind', {}, 'p', {});
if size(M, 2) < 4, return; end
[A, U, mu] = local_fitDMD(M, o.dmdRank);
E(end+1) = struct('kind', 'dmd', 'p', struct('A', A, 'U', U, 'mu', mu, 'mLast', M(:,end), ...
    'tLast', t(end), 'sub', o.subWin));
[aa, am] = local_fitAR1(M);
E(end+1) = struct('kind', 'ar1', 'p', struct('a', aa, 'mu', am, 'mLast', M(:,end), ...
    'tLast', t(end), 'sub', o.subWin));
E(end+1) = struct('kind', 'ramp', 'p', struct('B', local_fitRamp(M, t)));
end


function lam = local_predictAt(e, t)
switch e.kind
    case 'ramp', lam = e.p.B(:,1) + e.p.B(:,2) * t;
    case 'ar1'
        s = max(round((t - e.p.tLast) / e.p.sub), 0);
        lam = e.p.mu + (e.p.a .^ s) .* (e.p.mLast - e.p.mu);
    otherwise
        s = max(round((t - e.p.tLast) / e.p.sub), 0);
        y = e.p.U' * (e.p.mLast - e.p.mu);
        lam = e.p.mu + e.p.U * (e.p.A^s * y);
end
lam = max(lam, 0);
end


function B = local_basis(R, k)
% orthonormal rank-(k-1) basis of the residual rows of R
k = max(min([k - 1, size(R, 2) - 1, size(R, 1) - 1]), 1);
[~, ~, V] = svd(R, 'econ');
B = V(:, 1:min(k, size(V, 2)));
end


function ep = local_runs(flag, minRun)
% [start end] of every run of true at least minRun long
ep = zeros(0, 2);
flag = flag(:)'; flag(isnan(flag)) = false;
d = diff([false, flag, false]);
s = find(d == 1); e = find(d == -1) - 1;
keep = (e - s + 1) >= minRun;
if any(keep), ep = [s(keep)', e(keep)']; end
end


function [M, tC, ripC] = local_subWins(X, dur, tOn, subWin, minRip)
ed = min(tOn):subWin:max(tOn);
M = nan(size(X, 2), numel(ed) - 1); tC = nan(1, numel(ed) - 1); ripC = cell(1, numel(ed) - 1);
for b = 1:numel(ed) - 1
    k = find(tOn >= ed(b) & tOn < ed(b+1));
    if numel(k) < minRip, continue; end
    M(:,b) = sum(X(k,:), 1)' / sum(dur(k)); tC(b) = (ed(b) + ed(b+1))/2; ripC{b} = k;
end
ok = all(isfinite(M), 1) & isfinite(tC);
M = M(:, ok); tC = tC(ok); ripC = ripC(ok);
end


function [X, dur, tOn, nU] = local_ripCounts(fld, RP)
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


function [A, U, mu] = local_fitDMD(M, r)
mu = mean(M, 2);
Mc = M - mu;
r = min([r, size(M, 2) - 1, size(M, 1)]);
[U, ~, ~] = svd(Mc, 'econ'); U = U(:, 1:max(r, 1));
Y = U' * Mc;
A = Y(:, 2:end) / Y(:, 1:end-1);
[W, D] = eig(A); d = diag(D); bad = abs(d) > 1;
if any(bad) && rank(W) == size(W, 1)
    d(bad) = d(bad) ./ abs(d(bad));
    A = real(W * diag(d) / W);
end
end


function [a, mu] = local_fitAR1(M)
mu = mean(M, 2);
Mc = M - mu;
a = sum(Mc(:, 1:end-1) .* Mc(:, 2:end), 2) ./ max(sum(Mc(:, 1:end-1).^2, 2), eps);
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


function z = local_posZ(y, P)
% standardise the observed trace against the surrogate distribution AT THE SAME POSITION, using
% the median and a robust scale so one wild surrogate cannot set the yardstick
z = nan(size(y));
for t = 1:numel(y)
    c = P(:,t); c = c(isfinite(c));
    if numel(c) < 3 || ~isfinite(y(t)), continue; end
    s = 1.4826 * median(abs(c - median(c)));
    z(t) = (y(t) - median(c)) / max(s, eps);
end
end


function Z = local_posZloo(P)
% the same standardisation applied to each surrogate against the others, so the threshold is
% read off the null distribution of the statistic actually being used
Z = nan(size(P));
for q = 1:size(P, 1)
    k = true(size(P, 1), 1); k(q) = false;
    Z(q,:) = local_posZ(P(q,:), P(k,:));
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
m = abs(x(:)' * s / n);
p = max(mean(m >= abs(mean(x)) - 1e-12), 2 / 2^n);
end


function s = local_name(fld)
[~, s] = fileparts(fld);
end