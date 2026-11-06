%% noveltyEventAligned.m
% Is the novelty trace elevated AFTER a prescribed context exploration, across sessions?
%
% WHY THIS INSTEAD OF DECLARING EPOCHS
% driftNoveltyScan.m asks each session, on its own, whether some window is more novel than a
% max-statistic null allows. That is localisation from about thirty windows, and it is the wrong
% place to spend the power in this design: with SURR_N surrogates the threshold is set by the
% extreme-value behaviour of a small null, which is why the per-session z thresholds come back at
% 5 to 21 and nothing is ever declared. The one test on that page that does reach significance is
% the directedness contrast, and it reaches it because it is paired across sessions with a single
% number each. This script applies the same structure to the prescribed times: one number per
% session per context, tested across sessions.
%
% THE TEST
% Align each session's novelty trace to a context, take the mean over a window before it starts
% and over a window after it ends, and difference them. Then sign-flip across sessions.
%
%       delta = mean N(post) - mean N(pre),     post = [0, POST_WIN] h after the context ends
%                                               pre  = [-PRE_WIN, 0] h before it starts
%
% ConA1 (familiar) and ConB1 (novel) are run separately. A rise after both is a statement about
% exploration; a rise after ConB1 only is a statement about novelty, and the within-session
% difference ConB1 - ConA1 is what separates them.
%
% THE STATE CONFOUND, HANDLED RATHER THAN NOTED
% N correlates +0.62 with ripple rate across windows, and ripple rate rises during the rest that
% follows exploration, so a raw post-minus-pre rise is expected even with nothing learned. Every
% endpoint is therefore computed three ways:
%       raw        N as it comes
%       resid      N after regressing out log ripple rate WITHIN session, so the common
%                  state-driven component is removed before the alignment is applied
%       z          the position-matched standardised trace, which already divides out the
%                  library-size effect but not the state effect
% The honest result is the residualised one. If the raw effect survives and the residualised one
% does not, the finding is post-exploration sleep, not acquisition.
%
% INPUT: driftNoveltyScan.mat, written by driftNoveltyScan.m. Nothing is recomputed here.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
PRE_WIN  = 2.0;      % h before the context STARTS
POST_WIN = 2.0;      % h after the context ENDS
BIN      = 0.5;      % h, bins for the time course
BINS     = -3:BIN:4; % h relative to the context end
MIN_WIN  = 2;        % windows a session needs on each side to contribute

% ---- RIPPLE-RATE MATCHED PRE-WINDOWS -----------------------------------------------------
% Regressing N on log ripple rate removes the average linear effect of state. It does not make
% the two periods comparable: in these sessions ConA1 is always the FIRST exploration, so its
% pre-window sits in a low-ripple awake baseline and its post-window is the first descent into
% sleep, while ConB1 is already preceded by sleep-like activity. That is a difference in state
% occupancy, not in slope, and a regression cannot repair it.
% So each post window is instead paired with a pre window of the same ripple rate, drawn from
% anywhere before the context, and the contrast is computed over matched pairs only.
% The number of pairs is itself the diagnostic: if the pre period never reaches the post
% period's ripple rates, few pairs form, and that means the original comparison was between
% two states that do not overlap rather than between two times.
MATCH_PRE = true;
CALIPER   = 0.05;    % log10 ripples/h: the widest rate difference a pair may have
MIN_PAIR  = 3;       % pairs a session needs for the matched endpoint
% The caliper has to be tight to be worth anything. Simulating a post period whose ripple rate
% is 1.5x the pre period (the ConA1 case here) with N driven ONLY by rate, so the true effect
% is zero, the residual artefact left by matching is:
%       caliper   artefact   median pairs   sessions usable   achieved rate ratio
%         0.20     +0.100         8               99%              1.34
%         0.15     +0.071         7               99%              1.24
%         0.10     +0.044         6               96%              1.13
%         0.05     +0.021         4               84%              1.04
%         0.03     +0.002         4               71%              1.02
% A loose caliper matches in name only: at 0.15 the pairs still differ by 24% in rate and more
% than half the artefact survives. 0.05 is the compromise; 0.03 is the honest setting if enough
% sessions still form pairs. Read ripRatioMatch to see which of these actually happened.

CTXMARK = { ...
    'e13_26m1_211019', [15149.6 16749.2], [7605.7 9322]; ...
    'e13_26m1_211116', [13608 13939],     [7469.8 8707.2]; ...
    'e13_26m1_211119', [14911.9 16446],   [6181.3 7572]; ...
    'e15_13f1_220118', [13137 15046],     [5636.7 7700]; ...
    'e15_13f1_220119', [21219 22700],     [6836.6 8860]; ...
    'e15_13f1_220120', [18017.9 19801],   [9610.5 11950]; ...
    'e16_3m2_211211',  [9584.6 11203],    [5294.1 6511]; ...
    'e16_3m2_211212',  [16416.3 17731],   [6109.3 7845]; ...
    'e16_3m2_211213',  [13291 15323],     [7567.9 9359]};   % folder, NOVEL, FAMILIAR

L = load(fullfile(CACHE_ROOT, 'driftNoveltyScan.mat'), 'SC');
SC = L.SC;
if ~ismember('zN', SC.Properties.VariableNames)
    warning('align:noz', 'no zN column: rerun driftNoveltyScan.m, the z endpoint will be skipped');
    SC.zN = nan(height(SC), 1); SC.zB = nan(height(SC), 1);
end

ses = unique(SC.Session);
rows = {}; tcRows = {};
for i = 1:numel(ses)
    k = strcmp(SC.Session, ses{i});
    S = SC(k,:);
    jm = find(strcmp(CTXMARK(:,1), ses{i}), 1);
    if isempty(jm), fprintf('%-30s no context times\n', ses{i}); continue; end
    novH = CTXMARK{jm,2} / 3600; famH = CTXMARK{jm,3} / 3600;

    % ---- residualise the trace on log ripple rate, within session -----------------------
    % This is the state control. b is reported so the strength of the confound is visible per
    % session rather than only as the pooled +0.62.
    ok = isfinite(S.N) & isfinite(S.ripPerHour) & S.ripPerHour > 0;
    bRate = NaN; resN = nan(height(S), 1); resB = nan(height(S), 1);
    if sum(ok) > 4
        G = [ones(sum(ok), 1), log10(S.ripPerHour(ok))];
        cN = G \ S.N(ok); resN(ok) = S.N(ok) - G * cN; bRate = cN(2);
        cB = G \ S.B(ok); resB(ok) = S.B(ok) - G * cB;
    end

    % windows inside any context are excluded everywhere: the scan scores them, but they are
    % the event, not a measurement of what the event left behind
    inCtx = (S.Hour >= famH(1) & S.Hour < famH(2)) | (S.Hour >= novH(1) & S.Hour < novH(2));

    for c = 1:2
        if c == 1, iv = famH; nm = 'ConA1'; else, iv = novH; nm = 'ConB1'; end
        % the ConA1 post window is truncated at ConB1 onset, otherwise in several sessions it
        % runs straight into the novel exposure and measures that instead
        pEnd = iv(2) + POST_WIN;
        if c == 1, pEnd = min(pEnd, novH(1)); end
        kPre  = S.Hour >= iv(1) - PRE_WIN & S.Hour < iv(1) & ~inCtx;
        kPost = S.Hour >= iv(2) & S.Hour < pEnd & ~inCtx;
        if sum(kPre) < MIN_WIN || sum(kPost) < MIN_WIN
            fprintf('%-30s %s: %d pre, %d post windows, skipped\n', ses{i}, nm, sum(kPre), sum(kPost));
            continue
        end
        % matched pre-windows: candidates are every pre-context window, not just PRE_WIN
        dNm = NaN; dBm = NaN; rrM = NaN; nPair = 0;
        if MATCH_PRE
            iPre  = find(S.Hour < iv(1) & ~inCtx);
            iPost = find(kPost);
            [ip, iq] = local_matchPairs(S.ripPerHour(iPre), S.ripPerHour(iPost), CALIPER);
            nPair = numel(ip);
            if nPair >= MIN_PAIR
                a = iPre(ip); b = iPost(iq);
                dNm = mean(S.N(b), 'omitnan') - mean(S.N(a), 'omitnan');
                dBm = mean(S.B(b), 'omitnan') - mean(S.B(a), 'omitnan');
                rrM = mean(S.ripPerHour(b), 'omitnan') / max(mean(S.ripPerHour(a), 'omitnan'), eps);
            end
        end
        rows(end+1,:) = {ses{i}, nm, ...
            mean(S.N(kPost), 'omitnan')    - mean(S.N(kPre), 'omitnan'), ...
            mean(resN(kPost), 'omitnan')   - mean(resN(kPre), 'omitnan'), ...
            mean(S.zN(kPost), 'omitnan')   - mean(S.zN(kPre), 'omitnan'), ...
            mean(S.B(kPost), 'omitnan')    - mean(S.B(kPre), 'omitnan'), ...
            mean(resB(kPost), 'omitnan')   - mean(resB(kPre), 'omitnan'), ...
            mean(S.ripPerHour(kPost), 'omitnan') / max(mean(S.ripPerHour(kPre), 'omitnan'), eps), ...
            dNm, dBm, rrM, nPair, ...
            sum(kPre), sum(kPost), bRate, (iv(2)-iv(1))*60}; %#ok<AGROW>

        for b = 1:numel(BINS) - 1
            kb = S.Hour - iv(2) >= BINS(b) & S.Hour - iv(2) < BINS(b+1) & ~inCtx;
            if ~any(kb), continue; end
            tcRows(end+1,:) = {ses{i}, nm, (BINS(b) + BINS(b+1))/2, ...
                mean(S.N(kb), 'omitnan'), mean(resN(kb), 'omitnan'), ...
                mean(S.B(kb), 'omitnan'), mean(S.ripPerHour(kb), 'omitnan')}; %#ok<AGROW>
        end
    end
end
if isempty(rows), error('align:none', 'no session had usable pre and post windows'); end
T = cell2table(rows, 'VariableNames', {'Session','Context','dN','dNresid','dZ','dB','dBresid', ...
    'ripRatio','dNmatch','dBmatch','ripRatioMatch','nPair','nPre','nPost','bRate','ctxMin'});
TC = cell2table(tcRows, 'VariableNames', {'Session','Context','Bin','N','Nresid','B','ripPerHour'});

%% ---------------------------------------------------------------------------------------------
%  REPORT
%% ---------------------------------------------------------------------------------------------
fprintf('\n======== NOVELTY ALIGNED TO THE PRESCRIBED CONTEXTS ========\n');
fprintf('pre = %g h before the context starts, post = %g h after it ends, %d sessions\n', ...
    PRE_WIN, POST_WIN, numel(unique(T.Session)));
fprintf('one number per session per context; the test is a sign-flip over sessions\n');

CX = unique(T.Context);
fprintf('\n%-8s %-10s %9s %9s %9s %7s\n', 'context', 'endpoint', 'mean', 'p', 'ripRatio', 'n');
for c = 1:numel(CX)
    for f = {'dN','dNresid','dZ','dB','dBresid','dNmatch','dBmatch'}
        k = strcmp(T.Context, CX{c});
        x = T.(f{1})(k); x = x(isfinite(x));
        if numel(x) < 2, continue; end
        rr = mean(T.ripRatio(k), 'omitnan');
        if endsWith(f{1}, 'match'), rr = mean(T.ripRatioMatch(k), 'omitnan'); end
        fprintf('%-8s %-10s %9.3f %9.4f %9.2f %7d\n', CX{c}, f{1}, mean(x), local_signflip(x), ...
            rr, numel(x));
    end
end
fprintf('\ndN is raw, dNresid is after log ripple rate is regressed out within session, dZ is the\n');
fprintf('position-matched trace. ripRatio is post/pre ripple rate: when it exceeds 1 the raw\n');
fprintf('endpoint is partly state, and dNresid is the one to read.\n');
fprintf('per-session slope of N on log ripple rate: median %+.2f, range %+.2f to %+.2f\n', ...
    median(T.bRate, 'omitnan'), min(T.bRate), max(T.bRate));

fprintf('\n--- NOVEL vs FAMILIAR, paired within session ---\n');
for f = {'dN','dNresid','dZ','dB','dBresid','dNmatch','dBmatch'}
    d = local_pairCtx(T, f{1}, 'ConB1', 'ConA1');
    if sum(isfinite(d)) < 2, continue; end
    fprintf('  %-8s ConB1 - ConA1 = %+.3f, p = %.4f (n = %d)\n', f{1}, ...
        mean(d, 'omitnan'), local_signflip(d), sum(isfinite(d)));
end
fprintf('  A rise after both contexts is about exploration. A rise after ConB1 only, surviving\n');
fprintf('  the rate regression, is the claim that the novel experience added something.\n');

fprintf('\n--- per session ---\n');
fprintf('%-30s %-8s %8s %8s %8s %8s %8s %6s\n', 'session', 'context', 'dN', 'dNresid', ...
    'dNmatch', 'ripRat', 'rrMatch', 'nPair');
for q = 1:height(T)
    fprintf('%-30s %-8s %8.3f %8.3f %8.3f %8.2f %8.2f %6d\n', T.Session{q}, T.Context{q}, ...
        T.dN(q), T.dNresid(q), T.dNmatch(q), T.ripRatio(q), T.ripRatioMatch(q), T.nPair(q));
end
fprintf('\nnPair is how many post windows found a pre window at the same ripple rate within a\n');
fprintf('%.2f log10 caliper. A context whose pair count collapses had no overlapping pre state\n', CALIPER);
fprintf('at all, and for that context the unmatched contrast was comparing wake with sleep.\n');
fprintf('rrMatch should sit near 1.00 wherever pairs formed; if it does not, widen MIN_PAIR or\n');
fprintf('narrow CALIPER and accept fewer sessions.\n');

%% ---------------------------------------------------------------------------------------------
%  FIGURE
%% ---------------------------------------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [70 70 1150 430]);
tl = tiledlayout(fig, 1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
cols = [0.25 0.45 0.75; 0.80 0.30 0.20];
for p = 1:3
    nexttile(tl); hold on
    if p == 1, fv = 'N'; ttl = 'novelty, raw'; elseif p == 2, fv = 'Nresid'; ttl = 'novelty, rate-residualised';
    else, fv = 'ripPerHour'; ttl = 'ripple rate'; end
    for c = 1:numel(CX)
        bv = unique(TC.Bin); mu = nan(size(bv)); se = nan(size(bv));
        for b = 1:numel(bv)
            k = strcmp(TC.Context, CX{c}) & TC.Bin == bv(b);
            x = TC.(fv)(k); x = x(isfinite(x));
            if numel(x) < 2, continue; end
            mu(b) = mean(x); se(b) = std(x) / sqrt(numel(x));
        end
        errorbar(bv, mu, se, '-o', 'LineWidth', 1.4, 'Color', cols(min(c,2),:), ...
            'MarkerFaceColor', cols(min(c,2),:), 'DisplayName', CX{c});
    end
    xline(0, 'k--', 'context ends');
    xlabel('h relative to the context ending'); ylabel(ttl); title(ttl);
    if p == 1, legend('Location', 'best', 'Box', 'off'); end
    box off
end
title(tl, 'Novelty aligned to the prescribed exposures, across sessions');

save(fullfile(CACHE_ROOT, 'noveltyEventAligned.mat'), 'T', 'TC', 'PRE_WIN', 'POST_WIN', 'BINS');
fprintf('\nsaved noveltyEventAligned.mat\n');


%% ---------------------------------------------------------------------------------------------
function [ip, iq] = local_matchPairs(ratePre, ratePost, caliper)
% greedy 1:1 nearest-neighbour matching of post windows to pre windows on log ripple rate,
% without replacement and within a caliper. Post windows are taken in order of how extreme
% their rate is, so the hardest ones are matched while candidates remain.
ip = []; iq = [];
a = log10(max(ratePre(:), eps)); b = log10(max(ratePost(:), eps));
if isempty(a) || isempty(b), return; end
[~, ord] = sort(abs(b - median(a)), 'descend');
used = false(numel(a), 1);
for q = ord(:)'
    d = abs(a - b(q)); d(used) = inf;
    [dm, j] = min(d);
    if dm <= caliper
        used(j) = true; ip(end+1) = j; iq(end+1) = q; %#ok<AGROW>
    end
end
end


function d = local_pairCtx(T, fld, c1, c2)
s = unique(T.Session); d = nan(numel(s), 1);
for q = 1:numel(s)
    k1 = strcmp(T.Session, s{q}) & strcmp(T.Context, c1);
    k2 = strcmp(T.Session, s{q}) & strcmp(T.Context, c2);
    if any(k1) && any(k2), d(q) = mean(T.(fld)(k1), 'omitnan') - mean(T.(fld)(k2), 'omitnan'); end
end
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