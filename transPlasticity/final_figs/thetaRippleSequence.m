%% thetaRippleSequence.m
% Does the theta-timescale order of cell pairs during exploration recur as within-ripple
% spike timing afterwards (and not before)?  Dragoi & Buzsaki 2006; Diba & Buzsaki 2007.
%
% FOR EACH CONTEXT EXPOSURE X (ConA1 familiar, drug-free; ConB1 novel, after infusion)
%   running / theta part of X: identical to compareReactivation (tracking speed > SPEED_MIN
%   in 1 s windows, else LFP theta/delta > THETA_RATIO, else the whole exposure). Spikes are
%   kept inside those windows and the windows are separated in time, so no CCG lag spans a gap.
%   units: PYR, >= MIN_SPK_RUN spikes in the running part, present in rippleEvents.
%   CCG    of j relative to i (positive lag = j after i), CCG_BIN bins, +-CCG_HALF, smoothed
%          with a Gaussian of SD CCG_SMOOTH.
%   theta-modulated pair: >= MIN_CCG_N counts, spectral peak of the (mean-removed, Hann
%          windowed) CCG between 4 and 14 Hz lies in THETA_BAND, and mean power in THETA_BAND
%          is >= THETA_CCG_RATIO x the mean power in FLANKS. THETA_ONLY = false keeps every pair.
%   theta phase lag: least-squares fit CCG(t) = b0 + b1 t + b2 t^2 + A cos(w t - phi) over
%          +-CCG_HALF with w = 2 pi f_peak; phi (deg, -180..180) is the phase lag and
%          tauTheta = phi / w (s) the lag of the theta peak nearest zero. phi > 0: j after i.
%   ripple lag: in the PRE and POST home-cage blocks of compareReactivation (PRE = the hour
%          before X; for ConB1 all home cage since ConA1 minus the 30 min after infusion onset;
%          POST = up to 30 min after X), the mean of t_j - t_i over all spike pairs in the same
%          ripple with |t_j - t_i| <= RIP_LAG (centre of mass of the ripple CCG); a pair needs
%          >= MIN_RIP_CO ripples in which both fire, in BOTH blocks (same pairs before and after).
%   per exposure, over pairs: rPre = corr(phi, ripple lag PRE), rPost = corr(phi, ripple lag
%          POST), dR = rPost - rPre, rPostP = partial corr(phi, POST | PRE); Pearson (also
%          Spearman, *_s). rPost > 0 = forward replay of the theta order, < 0 = reverse.
% TESTS: within group, sign-flip over animal means (controls: sessions too); between groups,
%   permutation over sessions, LME y ~ Group + (1|Animal), paired sign-flip; ConB1 - ConA1 within
%   session. Pooled-pair rows are descriptive (pairs are not independent).
% GROUPS (GROUPING = 'all', default): Control (Huszar), SAL = every SAL session, ASV = every ASV
%   session (the session list prints SALpre / ASV1 / later for each). 'first': SALpre / ASV1 only.
% RANK ORDER (primary): per ripple, Spearman correlation of the active template cells' theta order
%   with their order in the ripple; forward (rho > 0) and backward (rho < 0) both count. See the block.
% POOLING (DO_PAIRS): pair-level phase lag vs ripple lag, pairs merged across sessions of a group.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT       = 'R:\TransPlasticity';
CELLS = 'PYR';
MIN_SPK_RUN = 100;          % spikes in the running / theta part of the exposure
MIN_UNITS = 5; MIN_PAIRS = 10;
SPEED_MIN = 3; THETA_RATIO = 2;           % running / theta mask, as compareReactivation
PAIR_EXCLUDE = 'none';                    % 'none' | 'sameChannel' | 'sameShank'
CCG_BIN = 0.005; CCG_HALF = 0.5; CCG_SMOOTH = 0.010;   % s
MIN_CCG_N = 200;                          % counts in the CCG
THETA_BAND = [6 10]; FLANKS = [3 5; 11 14]; THETA_CCG_RATIO = 2;   % Hz; theta-modulation criterion
THETA_ONLY = true;                        % false: every pair (sensitivity)
RIP_LAG = 0.05; MIN_RIP_CO = 5;           % s; co-firing ripples per block
EXCL_INF = false; INF_EXCL_WIN = [0 1800]; PRE_B_ALLGAP = true;
PRE_WIN = 3600; POST_WIN = 1800;
CONTEXTS = {'ConA1','ConB1'};
NPERM = 1e4; NPERM_POOL = 2000; NBOOT = 2000;
DO_PAIRS = false;           % pair-level pooled phase-lag vs ripple-lag analysis (slow; superseded by RANK ORDER)
MIN_TPL = 5;                % template: connected theta-modulated cells
MIN_ACTIVE = 5;             % template cells active in a ripple (4 cannot reach p < 0.05)
TIME_STAT = 'com';          % 'com' mean spike time in the ripple | 'first' first spike
ALPHA_RIP = 0.05; NSHUF = 500; MIN_RIP_EL = 20;
GROUPING = 'all';           % 'all': every SAL / every ASV session with the contexts (pairs pooled within
                            % group, primary); 'first': SALpre / ASV1 only, as compareReactivation
CONTROLS = { ...
    'R:\HuszarR\e13_26m1_211019', [15149.6 16749.2], [7605.7 9322];
    'R:\HuszarR\e13_26m1_211116', [13608 13939],     [7469.8 8707.2];
    'R:\HuszarR\e13_26m1_211119', [14911.9 16446],   [6181.3 7572];
    'R:\HuszarR\e15_13f1_220118', [13137 15046],     [5636.7 7700];
    'R:\HuszarR\e15_13f1_220119', [21219 22700],     [6836.6 8860];
    'R:\HuszarR\e15_13f1_220120', [18017.9 19801],   [9610.5 11950];
    'R:\HuszarR\e16_3m2_211211',  [9584.6 11203],    [5294.1 6511];
    'R:\HuszarR\e16_3m2_211212',  [16416.3 17731],   [6109.3 7845];
    'R:\HuszarR\e16_3m2_211213',  [13291 15323],     [7567.9 9359]};   % folder, novel, familiar
if exist('SHARED', 'var')
    CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root; PAIR_EXCLUDE = SHARED.pairExclude;
end
o = struct('cells', CELLS, 'minSpk', MIN_SPK_RUN, 'minUnits', MIN_UNITS, 'minPairs', MIN_PAIRS, ...
    'speedMin', SPEED_MIN, 'thetaRatio', THETA_RATIO, 'pairExclude', PAIR_EXCLUDE, ...
    'bin', CCG_BIN, 'half', CCG_HALF, 'smooth', CCG_SMOOTH, 'minN', MIN_CCG_N, 'band', THETA_BAND, ...
    'flanks', FLANKS, 'ratio', THETA_CCG_RATIO, 'thetaOnly', THETA_ONLY, 'ripLag', RIP_LAG, ...
    'minCo', MIN_RIP_CO, 'preWin', PRE_WIN, 'postWin', POST_WIN, 'minTpl', MIN_TPL, 'minActive', MIN_ACTIVE, ...
    'timeStat', TIME_STAT, 'alpha', ALPHA_RIP, 'nShuf', NSHUF, 'minRipEl', MIN_RIP_EL);
[o.EN, o.AN, o.exact] = local_nullTables(400, ALPHA_RIP, 2000);

rows = {};
%% ---- controls ----------------------------------------------------------------------------
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), warning('trs:file', 'missing %s', rf); continue; end
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES;
    nov = CONTROLS{i,2}; fam = CONTROLS{i,3}; rec = RP.meta.recordingSec;
    E = struct('name', {'ConA1','ConB1'}, 'on', {fam(1), nov(1)}, 'off', {fam(2), nov(2)}, 'preIv', {[], []});
    if fam(1) < nov(1)
        E(1).pre = [0 fam(1)];      E(1).post = [fam(2) nov(1)];
        E(2).pre = [fam(2) nov(1)]; E(2).post = [nov(2) rec];
        if PRE_B_ALLGAP, E(2).preIv = [fam(2) nov(1)]; end
    else
        E(2).pre = [0 nov(1)];      E(2).post = [nov(2) fam(1)];
        E(1).pre = [nov(2) fam(1)]; E(1).post = [fam(2) rec];
    end
    for e = E
        [y, P, R] = local_seq(RP, fld, e, o);
        rows(end+1,:) = [{'Control', 'Control', controlAnimal(rf), rf, e.name}, y, {{P}}, {{R}}]; %#ok<AGROW>
    end
end

%% ---- TransPlasticity -----------------------------------------------------------------------
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
for i = 1:numel(fa)
    fld = fileparts(fa{i});
    ev = local_events(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConA1_start')), continue; end
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    sub = local_group(RP, ev);                       % SALpre | ASV1 | '' (later SAL / repeat ASV)
    if isempty(sub), sub = [m.condition 'later']; end
    switch GROUPING
        case 'first', g = local_group(RP, ev);
        otherwise,   g = ''; if any(strcmp(m.condition, {'SAL','ASV'})), g = m.condition; end
    end
    if isempty(g), continue; end
    HC = local_hcBlocks(ev, RP, EXCL_INF, INF_EXCL_WIN);
    for c = 1:numel(CONTEXTS)
        on = ev.time(strcmp(ev.description, [CONTEXTS{c} '_start']));
        off = ev.time(strcmp(ev.description, [CONTEXTS{c} '_end']));
        if isempty(on) || isempty(off), continue; end
        e = struct('name', CONTEXTS{c}, 'on', on(1), 'off', off(1), 'pre', [], 'post', [], 'preIv', []);
        kp = find(abs(HC(:,2) - e.on) <= 1, 1); kq = find(abs(HC(:,1) - e.off) <= 1, 1);
        if PRE_B_ALLGAP && strcmp(CONTEXTS{c}, 'ConB1')
            e.preIv = local_gapIv(HC, ev, e.on); kp = [];
            if ~isempty(e.preIv), kp = 0; end
        end
        if isempty(kp) || isempty(kq), continue; end
        if kp > 0, e.pre = HC(kp,:); else, e.pre = [e.preIv(1,1) e.preIv(end,2)]; end
        e.post = HC(kq,:);
        [y, P, R] = local_seq(RP, fld, e, o);
        rows(end+1,:) = [{g, sub, m.animal, fa{i}, e.name}, y, {{P}}, {{R}}]; %#ok<AGROW>
    end
end

VN = {'Group','Sub','Animal','File','Context','nUnits','nPairs','nTheta','nUsed','thetaFrac','src', ...
      'rPre','rPost','dR','rPostP','rPre_s','rPost_s','dR_s','slopePost','why','P','R'};
T = cell2table(rows, 'VariableNames', VN);
fprintf('\nexposures scored (>= %d pairs): %d of %d\n', MIN_PAIRS, sum(isfinite(T.rPost)), height(T));
w = T(~isfinite(T.rPost), :);
for k = 1:height(w), [~, sn] = fileparts(fileparts(w.File{k})); fprintf('   %-8s %-34s %s %s\n', w.Group{k}, sn, w.Context{k}, w.why{k}); end
fprintf('criteria: PYR >= %d spikes running; CCG %g ms bins +-%g s, >= %d counts; theta pair: peak in %g-%g Hz and band/flank power >= %g (THETA_ONLY %d); ripple lag |dt| <= %g ms, >= %d co-firing ripples in PRE and POST\n', ...
    MIN_SPK_RUN, 1e3*CCG_BIN, CCG_HALF, MIN_CCG_N, THETA_BAND, THETA_CCG_RATIO, THETA_ONLY, 1e3*RIP_LAG, MIN_RIP_CO);

G = {'Control', 'SAL', 'ASV'};
if strcmp(GROUPING, 'first'), G = {'Control', 'SALpre', 'ASV1'}; end

%% ---- RANK ORDER (primary): every ripple tested against the theta order, forward OR backward ----
% template: each theta-modulated cell's position in the theta sequence of the exposure, from the
%   pairwise theta phase lags by least squares (x_j - x_i = phi_ij over all theta-modulated pairs of
%   the largest connected set; fitR2 = how well one order explains the pairwise lags).
% ripple: template cells active in the ripple (>= MIN_ACTIVE; with 4 cells no ripple can reach
%   p < 0.05), ordered by TIME_STAT spike time; rho = Spearman(template position, ripple order).
%   rho > 0 forward, < 0 backward. Per-ripple p: exact for <= 8 cells, t-approximation above.
% per block (PRE / POST): exAbs = mean(|rho| - E|rho| under the null for that ripple's cell count),
%   exSig = fraction significant (|rho|, two-sided) minus the null rate; fwd / bwd = fraction
%   significant forward / backward. Block p: template positions shuffled across cells (one shuffle
%   for all ripples of the block), NSHUF times. POST - PRE: ripple labels permuted within session.
%   Both excess measures are direction-free, so forward and backward replay add, never cancel.
RO = table();
for k = 1:height(T)
    R = T.R{k};
    r = struct2table(rmfield(R, {'rip', 'why'}), 'AsArray', true);
    r.Group = T.Group(k); r.Sub = T.Sub(k); r.Animal = T.Animal(k); r.File = T.File(k); r.Context = T.Context(k);
    r.why = {R.why}; r.rip = {R.rip};
    RO = [RO; r]; %#ok<AGROW>
end
okR = isfinite(RO.dExAbs);
fprintf('\n================ RANK ORDER: ripple order vs theta order (forward and backward) ================\n');
fprintf('template: >= %d connected theta-modulated cells; ripple: >= %d template cells active (%s spike time), alpha %.2f two-sided; %d template shuffles; >= %d eligible ripples per block\n', ...
    MIN_TPL, MIN_ACTIVE, TIME_STAT, ALPHA_RIP, NSHUF, MIN_RIP_EL);
fprintf('exposures scored: %d of %d; not scored:\n', sum(okR), height(RO));
w = RO(~okR, :);
for k = 1:height(w), [~, sn] = fileparts(fileparts(w.File{k})); fprintf('   %-8s %-34s %s %s\n', w.Group{k}, sn, w.Context{k}, w.why{k}); end
ROk = RO(okR, :);
for c = CONTEXTS
    fprintf('\n---- %s ----\n', c{1});
    for g = G
        t = ROk(strcmp(ROk.Group, g{1}) & strcmp(ROk.Context, c{1}), :);
        if isempty(t), continue; end
        fprintf('%s: %d sessions (%d animals) | template cells median %g (fit R2 %.2f) | eligible ripples PRE %g POST %g (medians)\n', ...
            g{1}, height(t), numel(unique(t.Animal)), median(t.nTpl), median(t.fitR2), median(t.nPre), median(t.nPost));
        for k = 1:height(t)
            [~, sn] = fileparts(fileparts(t.File{k}));
            fprintf('      %-34s %-8s tpl %3d  PRE %4d rip exAbs %+.3f (p=%.3f) sig %4.1f%% [F %4.1f B %4.1f] | POST %4d rip exAbs %+.3f (p=%.3f) sig %4.1f%% [F %4.1f B %4.1f] | dExAbs %+.3f (p=%.3f)\n', ...
                sn, t.Sub{k}, t.nTpl(k), t.nPre(k), t.exAbsPre(k), t.pAbsPre(k), 100*(t.fwdPre(k) + t.bwdPre(k)), 100*t.fwdPre(k), 100*t.bwdPre(k), ...
                t.nPost(k), t.exAbsPost(k), t.pAbsPost(k), 100*(t.fwdPost(k) + t.bwdPost(k)), 100*t.fwdPost(k), 100*t.bwdPost(k), t.dExAbs(k), t.pDAbs(k));
        end
        s = cellfun(@(nm) local_cell(t, nm), {'exAbsPre','exAbsPost','dExAbs','exSigPre','exSigPost','dExSig'}, 'UniformOutput', false);
        fprintf('   mean over sessions (p animals/sessions): exAbs PRE %s | POST %s | POST-PRE %s\n', s{1:3});
        fprintf('                                            exSig PRE %s | POST %s | POST-PRE %s\n', s{4:6});
        fprintf('   blocks p < 0.05 (template shuffle): PRE %d of %d, POST %d of %d\n', sum(t.pAbsPre < 0.05), height(t), sum(t.pAbsPost < 0.05), height(t));
        Rp = vertcat(t.rip{:});
        for b = 0:1
            q = Rp(Rp.post == b, :); nm = {'PRE', 'POST'};
            fprintf('   pooled ripples %-4s %6d: forward %5.2f%%, backward %5.2f%% (null %4.2f%% each)\n', nm{b+1}, height(q), ...
                100*mean(q.sig & q.rho > 0), 100*mean(q.sig & q.rho < 0), 50*mean(o.AN(min(q.n, numel(o.AN)))));
        end
    end
    fprintf('  between groups (perm over sessions | LME | paired animals)\n');
    for nm = {'exAbsPost', 'dExAbs', 'dExSig'}
        fprintf('    %-9s', nm{1});
        for pr = {G{2}, G{3}; G{1}, G{3}; G{1}, G{2}}'
            a = ROk(strcmp(ROk.Group, pr{1}) & strcmp(ROk.Context, c{1}), :); b = ROk(strcmp(ROk.Group, pr{2}) & strcmp(ROk.Context, c{1}), :);
            fprintf(' | %s vs %s %+.3f vs %+.3f: %.3f, %.3f, %s', pr{1}, pr{2}, mean(a.(nm{1})), mean(b.(nm{1})), ...
                local_perm(a.(nm{1}), b.(nm{1}), NPERM), local_lme([a; b], nm{1}, pr{1}, pr{2}), local_pairedStr([a; b], nm{1}, pr{1}, pr{2}));
        end
        fprintf('\n');
    end
end
fprintf('\n---- ConB1 - ConA1 within session (dExAbs, exAbsPost) ----\n');
for nm = {'dExAbs', 'exAbsPost'}
    for g = G
        t = ROk(strcmp(ROk.Group, g{1}), :);
        S = unique(t.File); d = nan(numel(S), 1); an = cell(numel(S), 1);
        for k = 1:numel(S)
            a = t(strcmp(t.File, S{k}) & strcmp(t.Context, 'ConA1'), :); b = t(strcmp(t.File, S{k}) & strcmp(t.Context, 'ConB1'), :);
            if height(a) == 1 && height(b) == 1, d(k) = b.(nm{1}) - a.(nm{1}); an{k} = a.Animal{1}; end
        end
        ok = isfinite(d); if ~any(ok), continue; end
        ua = unique(an(ok)); dm = cellfun(@(x) mean(d(ok & strcmp(an, x))), ua);
        fprintf('  %-9s %-8s %+.3f +- %.3f (%d sessions, %d animals)  animals p = %.3f, sessions p = %.3f\n', nm{1}, g{1}, ...
            mean(d(ok)), std(d(ok)) / sqrt(sum(ok)), sum(ok), numel(ua), local_signflip(dm), local_signflip(d(ok)));
    end
end
fr = figure('Color', 'w', 'Position', [100 100 1100 600]);
tr = tiledlayout(fr, numel(CONTEXTS), numel(G), 'TileSpacing', 'compact');
ed = -1:0.1:1;
for ci = 1:numel(CONTEXTS)
    for gi = 1:numel(G)
        ax = nexttile(tr); hold(ax, 'on');
        t = ROk(strcmp(ROk.Group, G{gi}) & strcmp(ROk.Context, CONTEXTS{ci}), :);
        if isempty(t), axis(ax, 'off'); continue; end
        Rp = vertcat(t.rip{:});
        histogram(ax, Rp.rho(Rp.post == 0), ed, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.2);
        histogram(ax, Rp.rho(Rp.post == 1), ed, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', [0.85 0.33 0.1], 'LineWidth', 1.2);
        title(ax, sprintf('%s %s (%d sessions)', G{gi}, CONTEXTS{ci}, height(t)), 'FontWeight', 'normal');
        if ci == numel(CONTEXTS), xlabel(ax, 'ripple rank correlation with theta order'); end
    end
end
title(tr, 'per-ripple rank correlation: PRE (grey), POST (orange); both tails are replay');
try, exportgraphics(fr, fullfile(CACHE_ROOT, 'thetaRippleRankOrder.pdf')); catch, end

if DO_PAIRS
%% ---- POOLED: theta-modulated pairs merged across the sessions of a group (primary) ----------
% Pair orientation (which cell is i) is arbitrary, so every correlation is taken through the
% origin: r = sum(w x y) / sqrt(sum(w x^2) sum(w y^2)), which is unchanged when a pair's phase
% lag and ripple lags all change sign. Computed from per-session sums, so a session can be
% weighted: 'pairs' = every pair counts once (large sessions dominate), 'sessions' = every
% session carries equal total weight (primary). Tests:
%   r > 0          shuffle ripple lags among pairs WITHIN each session (keeps session structure)
%   dR             swap PRE / POST labels per session (exact when <= 14 sessions); per pair (liberal)
%   between groups permute sessions between the two groups (exact when <= 2e4 assignments)
%   CI             bootstrap over sessions
Pall = table();
for k = 1:height(T)
    P = T.P{k};
    if isempty(P) || ~any(P.used), continue; end
    P = P(P.used, :); n = height(P);
    P.Group = repmat(T.Group(k), n, 1); P.Sub = repmat(T.Sub(k), n, 1); P.Context = repmat(T.Context(k), n, 1);
    P.Animal = repmat(T.Animal(k), n, 1); P.File = repmat(T.File(k), n, 1);
    Pall = [Pall; P]; %#ok<AGROW>
end
fprintf('\n%d usable pairs in total (THETA_ONLY %d; >= %d co-firing ripples in PRE and POST)\n', height(Pall), THETA_ONLY, MIN_RIP_CO);
WM = {'sessions', 'pairs'};
STAT = {'rPre', 'rPost', 'dR', 'rPost|PRE', 'slope ms/100deg'};
POOL = struct();
for c = CONTEXTS
    fprintf('\n================ POOLED %s: theta phase lag vs within-ripple lag ================\n', c{1});
    for g = G
        k = strcmp(Pall.Group, g{1}) & strcmp(Pall.Context, c{1});
        if sum(k) < 3, continue; end
        Q = Pall(k, :);
        [us, ~, sx] = unique(Q.File); nS = numel(us);
        fprintf('%s: %d pairs from %d sessions, %d animals\n', g{1}, height(Q), nS, numel(unique(Q.Animal)));
        for s = 1:nS
            [~, sn] = fileparts(fileparts(us{s}));
            fprintf('      %-34s %-8s %4d pairs (%4.1f%%)\n', sn, Q.Sub{find(sx == s, 1)}, sum(sx == s), 100*mean(sx == s));
        end
        M = local_sums(Q.phi, Q.ripPre, Q.ripPost, sx, nS);
        for wm = 1:2
            st = local_pool(M, WM{wm});
            [pPre, pPost, pPart] = local_shuffle(Q, sx, nS, WM{wm}, NPERM_POOL);
            pS = local_swapSess(M, WM{wm}, NPERM_POOL); pP = local_swapPair(Q, sx, nS, WM{wm}, NPERM_POOL);
            ci = local_boot(M, WM{wm}, NBOOT);
            fprintf('   weight %-8s rPre %+.3f (p=%.3f)  rPost %+.3f [%+.3f %+.3f] (p=%.3f)  dR %+.3f [%+.3f %+.3f] (sessions p=%.3f, pairs p=%.3f)  rPost|PRE %+.3f (p=%.3f)  slope %+.2f\n', ...
                WM{wm}, st(1), pPre, st(2), ci(1,2), ci(2,2), pPost, st(3), ci(1,3), ci(2,3), pS, pP, st(4), pPart, st(5));
            POOL.(c{1}).(g{1}).(WM{wm}) = struct('stat', st, 'pPre', pPre, 'pPost', pPost, 'pPartial', pPart, 'pDRsess', pS, 'pDRpair', pP, 'ci', ci, 'nPairs', height(Q), 'nSess', nS);
        end
    end
    fprintf('  between groups (permute sessions between the two groups), weight = sessions | pairs\n');
    for pr = {G{2}, G{3}; G{1}, G{3}; G{1}, G{2}}'
        ka = strcmp(Pall.Group, pr{1}) & strcmp(Pall.Context, c{1}); kb = strcmp(Pall.Group, pr{2}) & strcmp(Pall.Context, c{1});
        if sum(ka) < 3 || sum(kb) < 3, continue; end
        fprintf('     %-8s vs %-8s', pr{1}, pr{2});
        for wm = 1:2
            [d, p] = local_between(Pall(ka, :), Pall(kb, :), WM{wm}, NPERM_POOL);
            fprintf(' | %s: dRpost %+.3f p=%.3f, dDR %+.3f p=%.3f, dPartial %+.3f p=%.3f', WM{wm}, d(2), p(2), d(3), p(3), d(4), p(4));
        end
        fprintf('\n');
    end
end
fprintf('\n================ POOLED ConB1 - ConA1 (swap context labels per session) ================\n');
for g = G
    ka = strcmp(Pall.Group, g{1}) & strcmp(Pall.Context, 'ConA1'); kb = strcmp(Pall.Group, g{1}) & strcmp(Pall.Context, 'ConB1');
    if sum(ka) < 3 || sum(kb) < 3, continue; end
    for wm = 1:2
        [d, p] = local_ctxSwap(Pall(ka, :), Pall(kb, :), WM{wm}, NPERM_POOL);
        fprintf('  %-8s %-8s rPost B-A %+.3f p=%.3f | dR B-A %+.3f p=%.3f | partial B-A %+.3f p=%.3f\n', g{1}, WM{wm}, d(2), p(2), d(3), p(3), d(4), p(4));
    end
end

%% ---- per session (secondary: only exposures with >= MIN_PAIRS pairs) ----------------------
Ts = T(isfinite(T.rPost), :);
fprintf('\n================ PER SESSION (>= %d pairs; %d exposures) ================\n', MIN_PAIRS, height(Ts));
for c = CONTEXTS
    for g = G
        t = Ts(strcmp(Ts.Group, g{1}) & strcmp(Ts.Context, c{1}), :);
        if isempty(t), continue; end
        s = cellfun(@(nm) local_cell(t, nm), {'rPre','rPost','dR','rPostP'}, 'UniformOutput', false);
        fprintf('  %s %-8s %2d sess (%2d an), median %4.0f pairs | rPre %s | rPost %s | dR %s | rPost|PRE %s\n', c{1}, g{1}, height(t), numel(unique(t.Animal)), median(t.nUsed), s{:});
    end
end

%% ---- figure: pooled pairs, PRE grey, POST coloured ----------------------------------------
fh = figure('Color', 'w', 'Position', [80 80 1100 640]);
tl = tiledlayout(fh, numel(CONTEXTS), numel(G), 'TileSpacing', 'compact');
for ci = 1:numel(CONTEXTS)
    for gi = 1:numel(G)
        ax = nexttile(tl); hold(ax, 'on');
        Pp = Pall(strcmp(Pall.Group, G{gi}) & strcmp(Pall.Context, CONTEXTS{ci}), :);
        if isempty(Pp), axis(ax, 'off'); continue; end
        x = [Pp.phi; -Pp.phi];                                  % both orientations of every pair
        scatter(ax, x, 1e3*[Pp.ripPre; -Pp.ripPre], 6, [0.7 0.7 0.7], 'filled');
        scatter(ax, x, 1e3*[Pp.ripPost; -Pp.ripPost], 6, [0.85 0.33 0.1], 'filled');
        xline(ax, 0, ':'); yline(ax, 0, ':');
        xlim(ax, [-180 180]); ylim(ax, 1e3*RIP_LAG*[-1 1]);
        title(ax, sprintf('%s %s (%d pairs, %d sessions)', G{gi}, CONTEXTS{ci}, height(Pp), numel(unique(Pp.File))), 'FontWeight', 'normal');
        if gi == 1, ylabel(ax, 'ripple lag (ms)'); end
        if ci == numel(CONTEXTS), xlabel(ax, 'theta phase lag (deg)'); end
    end
end
title(tl, 'theta phase lag during exploration vs within-ripple lag, pooled pairs: PRE (grey), POST (orange)');
try, exportgraphics(fh, fullfile(CACHE_ROOT, 'thetaRippleSequence.pdf')); catch, end
end
T.P = []; T.R = []; Tsave = T; RO.rip = []; ROsave = RO;
if ~exist('Pall', 'var'), Pall = table(); POOL = struct(); end
save(fullfile(CACHE_ROOT, 'thetaRippleSequence.mat'), 'Tsave', 'ROsave', 'Pall', 'POOL', 'o');
fprintf('saved thetaRippleSequence.mat / .pdf\n');


%% ==== helpers ===============================================================================
function [y, P, R] = local_seq(RP, fld, e, o)
% y = {nUnits, nPairs, nTheta, nUsed, thetaFrac, src, rPre, rPost, dR, rPostP, rPre_s, rPost_s, dR_s, slopePost, why}
y = {0, 0, 0, 0, NaN, '', NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ''};
P = table(); R = local_rankDefault();
sf = dir(fullfile(fld, '*spikes.cellinfo.mat')); cf = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(sf) || isempty(cf), y{end} = 'no spikes or cell_metrics'; return; end
S = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); sp = S.spikes;
C = load(fullfile(cf(1).folder, cf(1).name), 'cell_metrics'); cm = C.cell_metrics;
ts = cellfun(@(x) unique(x(~isnan(x))), sp.times(:), 'UniformOutput', false);
nU = numel(ts);
clu = []; if isfield(sp, 'UID') && numel(sp.UID) == nU, clu = sp.UID(:); end
if isfield(sp, 'cluID') && numel(sp.cluID) == nU, clu = sp.cluID(:); end
loc = matchCellMetrics(cm, (1:nU)', clu, nU);
typ = repmat({''}, nU, 1); okL = loc > 0;
ct = cellstr(string(cm.putativeCellType(:))); typ(okL) = ct(loc(okL));
grp = nan(nU, 1); gf = '';
switch o.pairExclude
    case 'sameChannel', if isfield(cm, 'maxWaveformCh1'), gf = 'maxWaveformCh1'; elseif isfield(cm, 'maxWaveformCh'), gf = 'maxWaveformCh'; end
    case 'sameShank',   if isfield(cm, 'shankID'), gf = 'shankID'; end
end
if ~isempty(gf) && numel(cm.(gf)) == numel(ct), gv = double(cm.(gf)(:)); grp(okL) = gv(loc(okL)); end
% ---- running / theta part, intervals separated so no lag spans a gap ----
[mask, src] = local_runMask(fld, e.on, e.off, o.speedMin);
if isempty(src), [mask, src] = local_thetaMask(RP, e.on, e.off, o.thetaRatio); end
y{5} = sum(diff(mask, 1, 2)) / (e.off - e.on); y{6} = src;
if isempty(mask), y{end} = 'no running / theta'; return; end
iv = local_merge(mask);
tRun = cellfun(@(t) local_blockTimes(t, iv, 10 + 2*o.half), ts, 'UniformOutput', false);
nRun = cellfun(@numel, tRun);
keep = nRun >= o.minSpk;
if strcmp(o.cells, 'PYR'), keep = keep & contains(typ, 'Pyramidal'); end
si = RP.units.SpikesIdx(:);
if max(si) > nU, y{end} = 'RIPPLES / spikes sort mismatch'; return; end
uMap = zeros(nU, 1); uMap(si) = 1:numel(si);
keep = keep & uMap > 0;
idx = find(keep);
if numel(idx) < o.minUnits, y{end} = sprintf('< %d units', o.minUnits); return; end
% ---- ripple blocks (as compareReactivation) ----
preIv = local_preIv(e, o.preWin);
post = [max(e.post(1), e.off), min(e.post(2), e.off + o.postWin)];
inPre = find(local_ripInIv(RP, preIv)); inPost = find(RP.ts1 >= post(1) & RP.ts1Off <= post(2));
tA = local_ripTimes(RP, inPre, uMap(idx)); tB = local_ripTimes(RP, inPost, uMap(idx));
QA = local_ripCounts(RP, inPre, uMap(idx)) > 0; QB = local_ripCounts(RP, inPost, uMap(idx)) > 0;
% ---- pairs ----
pm = triu(true(numel(idx)), 1);
if any(isfinite(grp(idx))), pm = pm & (grp(idx) ~= grp(idx)'); end
[pI, pJ] = find(pm);
nP = numel(pI);
edges = (-o.half:o.bin:o.half)'; edges = edges - o.bin/2; edges(end+1) = edges(end) + o.bin;
lag = edges(1:end-1) + o.bin/2;
kern = exp(-0.5 * ((-4*o.smooth:o.bin:4*o.smooth)' / o.smooth).^2); kern = kern / sum(kern);
fs = 1 / o.bin; nfft = 2^nextpow2(8 * numel(lag)); f = (0:nfft-1)' * fs / nfft;
hw = 0.5 - 0.5*cos(2*pi*(0:numel(lag)-1)' / (numel(lag)-1));
nC = zeros(nP, 1); ratio = nan(nP, 1); fPk = nan(nP, 1); phi = nan(nP, 1); tauT = nan(nP, 1);
rPreL = nan(nP, 1); rPostL = nan(nP, 1); coA = zeros(nP, 1); coB = zeros(nP, 1);
for q = 1:nP
    a = tRun{idx(pI(q))}; b = tRun{idx(pJ(q))};
    d = local_diffs(a, b, o.half + o.bin);
    h = histcounts(d, edges)'; nC(q) = sum(h);
    coA(q) = sum(QA(:, pI(q)) & QA(:, pJ(q))); coB(q) = sum(QB(:, pI(q)) & QB(:, pJ(q)));
    if coA(q) >= o.minCo
        dd = local_diffs(tA{pI(q)}, tA{pJ(q)}, o.ripLag); if ~isempty(dd), rPreL(q) = mean(dd); end
    end
    if coB(q) >= o.minCo
        dd = local_diffs(tB{pI(q)}, tB{pJ(q)}, o.ripLag); if ~isempty(dd), rPostL(q) = mean(dd); end
    end
    if nC(q) < o.minN, continue; end
    hs = conv(h, kern, 'same');
    x = (hs - mean(hs)) .* hw; F = abs(fft(x, nfft)).^2;
    sb = f >= 4 & f <= 14; [~, k] = max(F .* sb); fPk(q) = f(k);
    inB = f >= o.band(1) & f <= o.band(2);
    fl = false(size(f)); for r = 1:size(o.flanks, 1), fl = fl | (f >= o.flanks(r,1) & f <= o.flanks(r,2)); end
    ratio(q) = mean(F(inB)) / mean(F(fl));
    wv = 2*pi*fPk(q);
    X = [ones(size(lag)), lag, lag.^2, cos(wv*lag), sin(wv*lag)];
    bc = X \ hs;
    phi(q) = atan2(bc(5), bc(4)) * 180 / pi;           % A cos(w t - phi): peak at t = phi / w
    tauT(q) = (phi(q) * pi / 180) / wv;
end
isTh = nC >= o.minN & fPk >= o.band(1) & fPk <= o.band(2) & ratio >= o.ratio;
used = isfinite(phi) & isfinite(rPreL) & isfinite(rPostL);
if o.thetaOnly, used = used & isTh; end
P = table(pI, pJ, nC, ratio, fPk, phi, tauT, rPreL, rPostL, coA, coB, isTh, used, ...
    'VariableNames', {'i','j','nCCG','thetaRatio','fPeak','phi','tauTheta','ripPre','ripPost','coPre','coPost','theta','used'});
y{2} = nP; y{3} = sum(isTh); y{4} = sum(used); y{1} = numel(idx);
R = local_rankOrder(RP, P, idx, uMap, inPre, inPost, o);   % per-ripple rank order (all theta-modulated pairs)
if sum(used) < o.minPairs, y{end} = sprintf('%d usable pairs (%d theta-modulated)', sum(used), sum(isTh)); return; end
u = P(used, :);
rPre = corr(u.phi, u.ripPre); rPost = corr(u.phi, u.ripPost);
rs1 = corr(u.phi, u.ripPre, 'Type', 'Spearman'); rs2 = corr(u.phi, u.ripPost, 'Type', 'Spearman');
rPP = corr(u.ripPre, u.ripPost);
rPostP = (rPost - rPre*rPP) / sqrt((1 - rPre^2) * (1 - rPP^2));
bq = polyfit(u.phi, 1e3*u.ripPost, 1);
y(7:14) = {rPre, rPost, rPost - rPre, rPostP, rs1, rs2, rs2 - rs1, 100*bq(1)};
end


function tb = local_blockTimes(t, iv, gap)
% spikes inside the intervals iv; each interval shifted by gap x its index so lags never span gaps
tb = zeros(0, 1);
if isempty(t) || isempty(iv), return; end
k = discretize(t, reshape(iv', [], 1));
in = mod(k, 2) == 1;
tb = t(in) + gap * (k(in) + 1) / 2;
end


function d = local_diffs(a, b, w)
% all b - a with |b - a| <= w; a, b sorted column vectors
d = zeros(0, 1);
if isempty(a) || isempty(b), return; end
a = a(:); b = b(:);
lo = local_countLE(a - w - eps, b) + 1; hi = local_countLE(a + w, b);
n = max(hi - lo + 1, 0);
if ~any(n), return; end
ra = repelem((1:numel(a))', n);
off = (1:sum(n))' - repelem(cumsum([0; n(1:end-1)]), n);
d = b(repelem(lo, n) + off - 1) - a(ra);
end


function s = local_cell(t, nm)
v = t.(nm); [am, ~] = local_animalMean(t, nm);
s = sprintf('%+.3f +- %.3f p=%.3f/%.3f', mean(v, 'omitnan'), std(v, 'omitnan') / sqrt(sum(isfinite(v))), local_signflip(am), local_signflip(v));
end


function M = local_sums(x, a, b, sx, nS)
% per-session sums: [n, xx, aa, bb, xa, xb, ab]
f = @(v) accumarray(sx, v, [nS 1]);
M = [f(ones(size(x))), f(x.^2), f(a.^2), f(b.^2), f(x.*a), f(x.*b), f(a.*b)];
end


function st = local_pool(M, wm)
% [rPre rPost dR rPost|PRE slope(ms per 100 deg)], correlations through the origin
w = ones(size(M, 1), 1);
if strcmp(wm, 'sessions'), w = 1 ./ max(M(:,1), 1); end
S = sum(M .* w, 1);
rxa = S(5) / sqrt(S(2) * S(3)); rxb = S(6) / sqrt(S(2) * S(4)); rab = S(7) / sqrt(S(3) * S(4));
pp = (rxb - rxa*rab) / sqrt((1 - rxa^2) * (1 - rab^2));
st = [rxa, rxb, rxb - rxa, pp, 1e5 * S(6) / S(2)];     % lags in s: 1e3 ms x 100 deg
end


function [pPre, pPost, pPart] = local_shuffle(Q, sx, nS, wm, nperm)
% shuffle each ripple-lag column among the pairs of the same session
[~, o0] = sort(sx); x = Q.phi(o0); a = Q.ripPre(o0); b = Q.ripPost(o0); s = sx(o0);
st0 = local_pool(local_sums(x, a, b, s, nS), wm);
cnt = zeros(1, 3);
for k = 1:nperm
    [~, pa] = sortrows([s, rand(numel(s), 1)]); [~, pb] = sortrows([s, rand(numel(s), 1)]);
    sA = local_pool(local_sums(x, a(pa), b, s, nS), wm);
    sB = local_pool(local_sums(x, a, b(pb), s, nS), wm);
    cnt = cnt + [abs(sA(1)) >= abs(st0(1)), abs(sB(2)) >= abs(st0(2)), abs(sB(4)) >= abs(st0(4))];
end
p = (1 + cnt) / (1 + nperm); pPre = p(1); pPost = p(2); pPart = p(3);
end


function p = local_swapSess(M, wm, nperm)
% dR: swap PRE / POST sums per session; exact over 2^nS when nS <= 14
st0 = local_pool(M, wm); nS = size(M, 1);
exact = nS <= 14;
if exact, Sg = dec2bin(0:2^nS-1, nS) == '1'; else, Sg = rand(nperm, nS) > 0.5; end
c = 0;
for k = 1:size(Sg, 1)
    Mk = M; sw = Sg(k, :)';
    Mk(sw, [3 4 5 6]) = M(sw, [4 3 6 5]);
    st = local_pool(Mk, wm); c = c + (abs(st(3)) >= abs(st0(3)) - 1e-12);
end
p = c / size(Sg, 1); if ~exact, p = (1 + c) / (1 + size(Sg, 1)); end
end


function p = local_swapPair(Q, sx, nS, wm, nperm)
% dR: swap PRE / POST lags per pair (pairs are not independent: liberal)
st0 = local_pool(local_sums(Q.phi, Q.ripPre, Q.ripPost, sx, nS), wm); c = 0;
for k = 1:nperm
    sw = rand(height(Q), 1) > 0.5; a = Q.ripPre; b = Q.ripPost; a(sw) = Q.ripPost(sw); b(sw) = Q.ripPre(sw);
    st = local_pool(local_sums(Q.phi, a, b, sx, nS), wm); c = c + (abs(st(3)) >= abs(st0(3)));
end
p = (1 + c) / (1 + nperm);
end


function ci = local_boot(M, wm, nboot)
% 95% bootstrap interval over sessions, rows = [lo; hi]
nS = size(M, 1); B = nan(nboot, 5);
for k = 1:nboot, B(k, :) = local_pool(M(randi(nS, nS, 1), :), wm); end
ci = prctile(B, [2.5 97.5], 1);
end


function [d, p] = local_between(QA, QB, wm, nperm)
% difference (B - A) of pooled stats; sessions permuted between groups
[ua, ~, sa] = unique(QA.File); [ub, ~, sb] = unique(QB.File);
MA = local_sums(QA.phi, QA.ripPre, QA.ripPost, sa, numel(ua)); MB = local_sums(QB.phi, QB.ripPre, QB.ripPost, sb, numel(ub));
d = local_pool(MB, wm) - local_pool(MA, wm);
M = [MA; MB]; nA = size(MA, 1); n = size(M, 1);
exact = nchoosek(n, nA) <= 2e4;
if exact
    C = nchoosek(1:n, nA); L = false(size(C, 1), n);
    for k = 1:size(C, 1), L(k, C(k, :)) = true; end
else
    L = false(nperm, n);
    for k = 1:nperm, L(k, randperm(n, nA)) = true; end
end
c = zeros(1, 5);
for k = 1:size(L, 1)
    dk = local_pool(M(~L(k, :), :), wm) - local_pool(M(L(k, :), :), wm);
    c = c + (abs(dk) >= abs(d) - 1e-12);
end
p = c / size(L, 1); if ~exact, p = (1 + c) / (1 + size(L, 1)); end
end


function [d, p] = local_ctxSwap(QA, QB, wm, nperm)
% ConB1 - ConA1 within a group: swap the two contexts' labels per session
us = unique([QA.File; QB.File]); nS = numel(us);
[~, sa] = ismember(QA.File, us); [~, sb] = ismember(QB.File, us);
MA = local_sums(QA.phi, QA.ripPre, QA.ripPost, sa, nS); MB = local_sums(QB.phi, QB.ripPre, QB.ripPost, sb, nS);
d = local_pool(MB, wm) - local_pool(MA, wm);
exact = nS <= 14;
if exact, Sg = dec2bin(0:2^nS-1, nS) == '1'; else, Sg = rand(nperm, nS) > 0.5; end
c = zeros(1, 5);
for k = 1:size(Sg, 1)
    sw = Sg(k, :)'; A = MA; Bm = MB; A(sw, :) = MB(sw, :); Bm(sw, :) = MA(sw, :);
    A = A(A(:,1) > 0, :); Bm = Bm(Bm(:,1) > 0, :);
    if isempty(A) || isempty(Bm), continue; end
    dk = local_pool(Bm, wm) - local_pool(A, wm); c = c + (abs(dk) >= abs(d) - 1e-12);
end
p = c / size(Sg, 1); if ~exact, p = (1 + c) / (1 + size(Sg, 1)); end
end


function R = local_rankDefault()
R = struct('nTpl', 0, 'fitR2', NaN, 'nPre', 0, 'nPost', 0, ...
    'exAbsPre', NaN, 'exAbsPost', NaN, 'dExAbs', NaN, 'pAbsPre', NaN, 'pAbsPost', NaN, 'pDAbs', NaN, ...
    'exSigPre', NaN, 'exSigPost', NaN, 'dExSig', NaN, 'pSigPre', NaN, 'pSigPost', NaN, 'pDSig', NaN, ...
    'fwdPre', NaN, 'bwdPre', NaN, 'fwdPost', NaN, 'bwdPost', NaN, 'rip', {table()}, 'why', '');
end


function R = local_rankOrder(RP, P, idx, uMap, inPre, inPost, o)
% theta-order template from pairwise phase lags, then per-ripple rank correlation in PRE / POST
R = local_rankDefault();
if isempty(P), R.why = 'no pairs'; return; end
th = P(P.theta & isfinite(P.phi), :);
if height(th) < 3, R.why = sprintf('%d theta-modulated pairs', height(th)); return; end
gr = graph(th.i, th.j, [], numel(idx)); bins = conncomp(gr);
cnt = accumarray(bins(:), 1); [~, big] = max(cnt); cc = find(bins == big);
if numel(cc) < o.minTpl, R.why = sprintf('template: %d connected theta-modulated cells', numel(cc)); return; end
k = ismember(th.i, cc) & ismember(th.j, cc); th = th(k, :); h = height(th);
[~, ii] = ismember(th.i, cc); [~, jj] = ismember(th.j, cc);
D = full(sparse([(1:h)'; (1:h)'], [ii(:); jj(:)], [-ones(h, 1); ones(h, 1)], h, numel(cc)));
x = lsqminnorm(D, th.phi);                       % theta position (deg): later = larger
R.nTpl = numel(cc); R.fitR2 = 1 - sum((D*x - th.phi).^2) / sum(th.phi.^2);
cols = uMap(idx(cc));
A = local_rankBlock(RP, inPre, cols, x, o); B = local_rankBlock(RP, inPost, cols, x, o);
R.nPre = numel(A.rho); R.nPost = numel(B.rho);
R.rip = table([A.rho; B.rho], [A.n; B.n], [zeros(R.nPre, 1); ones(R.nPost, 1)], [A.sig; B.sig], 'VariableNames', {'rho', 'n', 'post', 'sig'});
if R.nPre < o.minRipEl || R.nPost < o.minRipEl
    R.why = sprintf('%d / %d eligible ripples (PRE / POST), %d template cells', R.nPre, R.nPost, R.nTpl); return
end
R.exAbsPre = A.exAbs; R.exAbsPost = B.exAbs; R.pAbsPre = A.pAbs; R.pAbsPost = B.pAbs;
R.exSigPre = A.exSig; R.exSigPost = B.exSig; R.pSigPre = A.pSig; R.pSigPost = B.pSig;
R.fwdPre = A.fwd; R.bwdPre = A.bwd; R.fwdPost = B.fwd; R.bwdPost = B.bwd;
e = [A.e; B.e]; s = [A.es; B.es]; lab = [false(R.nPre, 1); true(R.nPost, 1)];
R.dExAbs = mean(e(lab)) - mean(e(~lab)); R.dExSig = mean(s(lab)) - mean(s(~lab));
c = [0 0];
for q = 1:o.nShuf
    L = lab(randperm(numel(lab)));
    c = c + [abs(mean(e(L)) - mean(e(~L))) >= abs(R.dExAbs), abs(mean(s(L)) - mean(s(~L))) >= abs(R.dExSig)];
end
R.pDAbs = (1 + c(1)) / (1 + o.nShuf); R.pDSig = (1 + c(2)) / (1 + o.nShuf);
end


function B = local_rankBlock(RP, rk, cols, x, o)
B = struct('rho', zeros(0, 1), 'n', zeros(0, 1), 'sig', false(0, 1), 'e', zeros(0, 1), 'es', zeros(0, 1), ...
    'exAbs', NaN, 'exSig', NaN, 'pAbs', NaN, 'pSig', NaN, 'fwd', NaN, 'bwd', NaN);
if isempty(rk), return; end
sp = RP.spikes(ismember(RP.spikes.Ripple, rk) & ismember(RP.spikes.Unit, cols), :);
if isempty(sp), return; end
[~, ri] = ismember(sp.Ripple, rk); [~, ci] = ismember(sp.Unit, cols);
nc = numel(cols);
if strcmp(o.timeStat, 'first'), fn = @min; else, fn = @mean; end
Tm = accumarray([ri(:) ci(:)], sp.ts1(:), [numel(rk) nc], fn, NaN);
Act = ~isnan(Tm); n = sum(Act, 2); el = n >= o.minActive;
if ~any(el), return; end
Act = Act(el, :); Tm = Tm(el, :); n = n(el);
Tm(~Act) = Inf; [~, ord] = sort(Tm, 2);
Rt = zeros(size(Tm));
for r = 1:size(Tm, 1), Rt(r, ord(r, 1:n(r))) = 1:n(r); end
rho = local_rhoFast(Act, Rt, x, n);
sig = local_pRho(rho, n, o) < o.alpha;
nm = numel(o.EN); nn = min(n, nm);
B.rho = rho; B.n = n; B.sig = sig;
B.e = abs(rho) - o.EN(nn); B.es = double(sig) - o.AN(nn);
B.exAbs = mean(B.e); B.exSig = mean(B.es); B.fwd = mean(sig & rho > 0); B.bwd = mean(sig & rho < 0);
obsA = mean(abs(rho)); obsS = mean(sig); cA = 0; cS = 0;
for q = 1:o.nShuf                                 % one shuffle of template positions for every ripple
    rs = local_rhoFast(Act, Rt, x(randperm(nc)), n);
    cA = cA + (mean(abs(rs)) >= obsA); cS = cS + (mean(local_pRho(rs, n, o) < o.alpha) >= obsS);
end
B.pAbs = (1 + cA) / (1 + o.nShuf); B.pSig = (1 + cS) / (1 + o.nShuf);
end


function rho = local_rhoFast(Act, Rt, x, n)
% Spearman rho per ripple (no ties): rank of each active cell's template position among the
% ripple's active cells, against its rank in time
C = double(x(:) <= x(:)');
Rx = (double(Act) * C) .* Act;
rho = 1 - 6 * sum((Rx - Rt).^2, 2) ./ (n .* (n.^2 - 1));
end


function p = local_pRho(rho, n, o)
% two-sided p of Spearman rho with n cells: exact for n <= numel(o.exact), t-approximation above
p = ones(size(rho));
for m = unique(n(:))'
    k = n == m;
    if m <= numel(o.exact) && ~isempty(o.exact{m})
        E = o.exact{m}; aq = abs(rho(k)) - 1e-9;
        pk = interp1(E.u, E.ge, aq, 'next');
        pk(aq < E.u(1)) = 1; pk(isnan(pk)) = 0;
        p(k) = pk;
    else
        t = rho(k) .* sqrt((m - 2) ./ max(1 - rho(k).^2, eps));
        p(k) = 2 * tcdf(-abs(t), m - 2);
    end
end
end


function [EN, AN, EX] = local_nullTables(nMax, alpha, nSim)
% null of Spearman rho for n cells: E|rho|, P(p < alpha); exact tables (n <= 8) for the p values
EN = nan(nMax, 1); AN = nan(nMax, 1); EX = cell(8, 1);
for m = 3:8
    Pm = perms(1:m); a = abs(1 - 6 * sum((Pm - (1:m)).^2, 2) / (m * (m^2 - 1)));
    u = unique(a); ge = arrayfun(@(v) mean(a >= v - 1e-12), u);
    EX{m} = struct('u', u, 'ge', ge);
    pa = interp1(u, ge, a - 1e-9, 'next'); pa(isnan(pa)) = 0; pa(a - 1e-9 < u(1)) = 1;
    EN(m) = mean(a); AN(m) = mean(pa < alpha);
end
for m = 9:nMax
    [~, Pm] = sort(rand(nSim, m), 2);
    r = 1 - 6 * sum((Pm - (1:m)).^2, 2) / (m * (m^2 - 1));
    t = r .* sqrt((m - 2) ./ max(1 - r.^2, eps));
    EN(m) = mean(abs(r)); AN(m) = mean(2 * tcdf(-abs(t), m - 2) < alpha);
end
end


%% ---- copied unchanged from compareReactivation.m ----
function iv = local_merge(mask)
% union of [start end] rows
iv = zeros(0, 2);
if isempty(mask), return; end
mask = sortrows(mask); iv = mask(1,:);
for k = 2:size(mask, 1)
    if mask(k,1) <= iv(end,2), iv(end,2) = max(iv(end,2), mask(k,2)); else, iv(end+1,:) = mask(k,:); end %#ok<AGROW>
end
end


function m = local_inIntervals(t, iv)
m = false(size(t));
for k = 1:size(iv, 1), m = m | (t >= iv(k,1) & t < iv(k,2)); end
end


function tc = local_ripTimes(RP, rk, cols)
% in-ripple spike times per unit, each ripple moved to its own block so lags never cross ripples
gap = 10;
sp = RP.spikes(ismember(RP.spikes.Ripple, rk) & ismember(RP.spikes.Unit, cols), :);
[~, ri] = ismember(sp.Ripple, rk); [~, ci] = ismember(sp.Unit, cols);
t = sp.ts1 - RP.ts1(sp.Ripple) + ri * gap;
tc = accumarray(ci, t, [numel(cols) 1], @(x) {sort(x)}, {zeros(0, 1)});
end


function k = local_countLE(x, c)
% number of elements of sorted c that are <= x, for each x
[~, ~, k] = histcounts(x, [-Inf; c; Inf]); k = k(:) - 1;
end


function Cn = local_ripCounts(RP, rk, cols)
% ripples x units spike counts, units = rows of RP.units given by cols
sp = RP.spikes(ismember(RP.spikes.Ripple, rk) & ismember(RP.spikes.Unit, cols), :);
[~, ri] = ismember(sp.Ripple, rk); [~, ci] = ismember(sp.Unit, cols);
Cn = full(sparse(ri, ci, 1, numel(rk), numel(cols)));
end


function [mask, src] = local_runMask(fld, a, b, vmin)
% 1 s windows inside [a b] whose mean speed > vmin; src = '' when there is no tracking file
mask = zeros(0, 2); src = '';
d = dir(fullfile(fld, '*.Behavior.mat'));
if isempty(d), return; end
[~, o_] = sort(endsWith({d.name}, 'Tracking.Behavior.mat')); d = d(o_);   % <session>.Behavior.mat first
try
    L = load(fullfile(d(1).folder, d(1).name), 'behavior'); bh = L.behavior;
    t = bh.timestamps(:); x = bh.position.x(:); y = bh.position.y(:);
    lin = []; if isfield(bh.position, 'lin'), lin = bh.position.lin(:); end
catch ME
    warning('react:track', '%s: %s', d(1).name, ME.message); return
end
inExp = t >= a & t <= b;
useLin = numel(lin) == numel(t) && sum(isfinite(lin(inExp))) > 10;   % linearised position (cm) when available
if useLin, ok = isfinite(t) & isfinite(lin); else, ok = isfinite(t) & isfinite(x) & isfinite(y); end
t = t(ok); x = x(ok); y = y(ok); if useLin, lin = lin(ok); end
if numel(t) < 3, return; end
[t, iu] = unique(t); x = x(iu); y = y(iu); if useLin, lin = lin(iu); end
if useLin
    v = [0; abs(diff(lin)) ./ diff(t)];
    v(v > 200) = NaN;                          % jumps where the linearisation wraps (track ends / arms)
else
    v = [0; hypot(diff(x), diff(y)) ./ diff(t)];
end
dt = median(diff(t)); v = movmean(v, max(1, round(0.25 / dt)), 'omitnan');   % 250 ms smoothing
w0 = (floor(a):1:ceil(b) - 1)'; w0 = w0(w0 >= a & w0 + 1 <= b);
if isempty(w0), src = sprintf('tracking speed > %g', vmin); return; end
[~, ~, bi] = histcounts(t, [w0; w0(end) + 1]);
okb = bi > 0;
ms = accumarray(bi(okb), v(okb), [numel(w0) 1], @mean, NaN);
% units: vmin is cm/s. Position stored in metres (or normalised track units) gives
% speeds < ~1; detect that from this exposure's own speeds and rescale.
scale = 1;
if ~useLin                                     % x/y: detect metres from this exposure's speeds
    v95 = prctile(ms(isfinite(ms)), 95);
    if isfinite(v95) && v95 < 1, scale = 100; end
end
k = find(ms * scale > vmin);
mask = [w0(k), w0(k) + 1];
if useLin, src = sprintf('tracking position.lin speed > %g cm/s', vmin);
else, src = sprintf('tracking x/y speed > %g cm/s', vmin); end
if scale ~= 1, src = [src ' (x/y x100: metres)']; end
end


function [mask, src] = local_thetaMask(RP, a, b, ratio)
% 1 s windows inside [a b] with theta/delta power ratio > ratio; whole [a b] when no LFP
mask = [a b]; src = 'whole exposure (no LFP)';
m = RP.meta;
if ~isfield(m, 'lfpFile') || ~exist(m.lfpFile, 'file') || ~isfield(m, 'channel1based'), return; end
[fd, bn] = fileparts(m.lfpFile); xf = fullfile(fd, [bn '.xml']);
if ~exist(xf, 'file'), return; end
try
    xml = LoadXml(xf); fs = m.fs;
    x = double(LoadBinary(m.lfpFile, 'nChannels', xml.nChannels, 'channels', m.channel1based, ...
        'frequency', fs, 'start', a, 'duration', b - a));
catch ME
    warning('react:lfp', '%s: %s', m.lfpFile, ME.message); return
end
fs = round(fs); nw = floor(numel(x) / fs);
if nw < 1, mask = zeros(0, 2); return; end
X = reshape(x(1:nw*fs), fs, nw);
F = abs(fft(X - mean(X, 1))).^2; f = (0:fs-1)';      % 1 s windows: 1 Hz resolution
th = sum(F(f >= 6 & f <= 10, :), 1); de = sum(F(f >= 1 & f <= 4, :), 1);
k = find(th ./ de > ratio);
mask = [a + (k(:) - 1), a + k(:)];
src = sprintf('LFP theta/delta > %g', ratio);
end


function ev = local_events(fld)
% Session markers, normalised to HC_start/_end, Infusion_start/_end, ConA1_start, ConB2_end, ...
% Sources: the .evt.sti (LoadEvents); if it has no context markers, an annotations file
% (*transitions*.mat with table all_events: event, Intan_ts, ..., video) in the session folder
% or up to two folders above, rows matched to this session by folder name in the video path.
% Labels like 'ConA start' / 'ConA-end' / 'conb1_Start' are accepted; contexts without a number
% are numbered by order of their starts. ev.source says which file(s) were used.
ev = []; T = {}; X = []; src = {};
c = dir(fullfile(fld, '*.evt.sti'));
if ~isempty(c)
    try
        ws_ = warning('off', 'all');                    % LoadEvents prints an FMAToolbox notice on every call
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
if ~any(startsWith(b, 'Con'))                         % no contexts in the .evt.sti: annotations file
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
                if any(hit), m = hit; elseif k > 1, continue; end   % file above the session: must name it
            end
            T = [T(~startsWith(b, {'HC','Con'})); cellstr(string(A.event(m)))]; %#ok<AGROW>  keep evt infusion markers
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
% a context start without its end: end = the next marker of any kind, else start + 30 min
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


function IV = local_preIv(e, preWin)
% PRE intervals (n x 2): e.preIv when given (used as is), else the last preWin of e.pre before e.on
if isfield(e, 'preIv') && ~isempty(e.preIv)
    IV = e.preIv;
else
    IV = [max(e.pre(1), e.on - preWin), min(e.pre(2), e.on)];
end
IV = IV(diff(IV, 1, 2) > 0, :);
end


function k = local_ripInIv(RP, IV)
% ripples lying wholly inside any of the intervals IV (n x 2)
k = false(numel(RP.ts1), 1);
for q = 1:size(IV, 1), k = k | (RP.ts1(:) >= IV(q,1) & RP.ts1Off(:) <= IV(q,2)); end
end


function IV = local_gapIv(HC, ev, on)
% home cage (HC, already minus contexts and infusion exclusion) between the end of the last context
% before 'on' and 'on'
lab = ev.description;
pe = ev.time(startsWith(lab, 'Con') & endsWith(lab, '_end') & ev.time <= on + 1);
prevEnd = 0; if ~isempty(pe), prevEnd = max(pe); end
IV = [max(HC(:,1), prevEnd), min(HC(:,2), on)];
IV = IV(diff(IV, 1, 2) > 0, :);
end


function HC = local_hcBlocks(ev, RP, exclInf, infWin)
% Home-cage intervals = [0, end of recording] minus every ConA*/ConB* interval, minus
% [Infusion_start - infWin(1), Infusion_start + infWin(2)] for every infusion start, and minus the
% marked Infusion_start..Infusion_end interval when exclInf.
% HC markers are not used (often missing or unpaired, e.g. two HC_start in a row).
% Each context start is paired with the first matching _end after it.
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


function [b, num, edge] = local_parseLabel(s)
% Context / home-cage / infusion markers in the formats seen so far, e.g.
%   'ConA1_start'  'ConA start'  'Con A start'  'ContextB-2 end'  'conb_Start'  'Familiar start'
%   'Novel_end'  'HC_start'  'Home cage end'  'homecage_start'  'Infusion_start'  'Inf end'
% -> base (ConA | ConB | HC | Infusion), number ('' if none), edge (start | end)
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


function g = local_group(RP, ev)
% SALpre | ASV1 | '' (SAL sessions after the first ASV and repeat ASV sessions are left out).
% The infusion time comes from RIPPLES.meta, else the events' Infusion_start, else ConA1_end
% (only its side of the first ASV matters here).
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


function [v, ids] = local_animalMean(T, nm)
ids = unique(T.Animal); v = cellfun(@(a) mean(T.(nm)(strcmp(T.Animal, a))), ids);
end


function p = local_signflip(d)
d = d(isfinite(d)); n = numel(d); p = NaN;
if n < 2, return; end
obs = abs(mean(d));
if n <= 14
    Sg = 2*(dec2bin(0:2^n-1, n) == '1') - 1; p = mean(abs(Sg*d(:)/n) >= obs - 1e-12);
else
    nul = arrayfun(@(k) mean((2*(rand(n,1) > 0.5) - 1) .* d(:)), 1:1e4); p = (1 + sum(abs(nul) >= obs)) / (1 + 1e4);
end
end


function s = local_pairedStr(T, nm, a, b)
A = T(strcmp(T.Group, a), :); Bt = T(strcmp(T.Group, b), :);
ids = intersect(unique(A.Animal), unique(Bt.Animal));
if numel(ids) < 2, s = sprintf('n/a (%d animals)', numel(ids)); return; end
d = cellfun(@(x) mean(Bt.(nm)(strcmp(Bt.Animal, x))) - mean(A.(nm)(strcmp(A.Animal, x))), ids);
s = sprintf('%.3f (%d animals)', local_signflip(d), numel(ids));
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
    warning('react:lme', '%s: %s', nm, ME.message);
end
end