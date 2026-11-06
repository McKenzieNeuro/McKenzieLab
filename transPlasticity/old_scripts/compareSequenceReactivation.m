%% compareSequenceReactivation.m
% Do ordered L-cell sequences expressed during a context exposure recur in the ripples that follow
% it more than in the ripples that preceded it?  Session layout, groups and PRE / POST blocks are
% those of compareReactivation.m; helpers marked =CR= are copied from it unchanged.
%
% TEMPLATES (from the exposure itself)
%   theta part of exposure X: running windows from <session>.Tracking.Behavior.mat when it exists,
%   else theta/delta > THETA_RATIO in 1 s windows, else the whole exposure (thetaSrc says which);
%   then every ripple detected inside the exposure, +- RIP_PAD, is removed (thetaFrac is what is left).
%   Units: CELLS, >= MIN_SPK spikes in that part, present in RIPPLES.units.
%   BURST_ISI > 0: a spike of a cell < BURST_ISI after that cell's previous spike is dropped (theta
%             and ripples alike), so a burst enters chains once.
%   chain     L = SEQ_LEN spikes of L distinct units, successive lags in (THETA_WIN(1), THETA_WIN(2)]
%   bin       THETA_BIN windows within the theta part; a chain belongs to the bin of its first spike.
%             O_b = orderings of cell set C with >= 1 chain in bin b (L! possible orderings).
%   specificity of C: n(pi) = #bins with pi in O_b. Relabelling C's cells within a bin gives
%             P(pi in O_b) = |O_b|/L! for every pi:  E = sum_b |O_b|/L!, V = sum_b p_b(1 - p_b).
%             template order fwd = argmax n(pi);  z = (n(fwd) - E) / sqrt(V)
%             kept: >= MIN_THETA_BINS informative bins (1 <= |O_b| <= L!-1), L! * p(z) < SEL_P,
%             no two units on the same channel (PAIR_EXCLUDE); at most MAX_TEMPLATES, by z.
%             Selection uses only the exposure, so it is independent of PRE / POST.
% SCORING IN RIPPLES
%   same chain definition, lags in SEQ_RIP_WIN, chains never cross ripples
%   stratum = (ripple r, template C);  O_rC = orderings of C present in r
%   2 x 2 per stratum: tested ordering present / absent vs reference orderings present / absent;
%   reference = orderings of C other than fwd, its reverse, the random and the rate order
%   (L = 2: the other ordering). Under relabelling of C's cells E[ad - bc] = 0 in every stratum,
%   so the log Mantel-Haenszel odds ratio psi = 0 is "no order preference".
%   Per exposure the MH numerator / denominator sums are kept (mh), PRE and POST, for 5 orderings:
%     1 fwd   template order                      2 rev   its reverse
%     3 rnd   a fixed random ordering of C (none at L = 2)
%     4 rate  C ordered by exposure rate, highest first (gain change without sequence content)
%     5 fwdK  fwd with PRE and POST strata subsampled to the same distribution of |O_rC|
%     6 fwd~R fwd on the templates whose order is not the rate order
%   dPsi = psi(POST) - psi(PRE)  is the estimand; psi(PRE) > 0 = order present before the exposure
% INFERENCE: pooled MH (sparse strata: MH stays consistent as strata grow in number, not size)
%   group x context: psi pooled over every exposure in the cell (sum of numerators / denominators);
%     each exposure's sums are scaled by w = min(1, STRATA_CAP / nInv), nInv = mean over PRE and
%     POST of its strata with 1 <= |O_rC| <= L!-1. nInv is invariant under relabelling, so the
%     weighted psi keeps its null of 0; the cap stops 150-unit sessions from outweighing 10-unit
%     sessions by 1000x. STRATA_CAP = Inf is plain MH (printed alongside as dPsi(MH)).
%     nEff = Kish effective number of animals, (sum s_a)^2 / sum s_a^2, s_a = animal's weighted strata.
%     CI = bootstrap over animals; p = exact permutation swapping PRE and POST within animal
%   B1 - A1 within group: context labels swapped within animal (all its sessions together)
%   between groups (per context, and on B1 - A1): animals in both groups -> group labels swapped
%     within animal; otherwise animals permuted between groups
%   Per-exposure psi / dPsi are printed only when both PRE and POST have >= MIN_STRATA strata.
% External: matchCellMetrics, controlAnimal, getAllExtFiles, LoadEvents, LoadXml, LoadBinary.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT       = 'R:\TransPlasticity';
CELLS = 'PYR';
THETA_RATIO = 2;            % theta (6-10 Hz) / delta (1-4 Hz) power, 1 s windows
PAIR_EXCLUDE = 'sameChannel';   % 'sameChannel' | 'sameShank' | 'none': sets with such a pair are dropped
SPEED_MIN = 3;              % cm/s, 1 s windows (tracking file, when present)
MIN_SPK = 100;               % spikes in the exposure (theta part)
MIN_UNITS = 10;
MIN_RIP = 50; MAX_RIP = Inf;
MATCH_N = false;            % false: every ripple in the PRE / POST blocks; true: the same number in both
EXCL_INF = false;
INF_EXCL_WIN = [0 1800];    % s [before after] Infusion_start, removed from every PRE/POST block
PRE_B_ALLGAP = true;        % ConB1 PRE = all home cage between ConA1 and ConB1 (minus the infusion window)
PRE_WIN = 3600;             % s: ripples from this long before the exposure (PRE)
POST_WIN = 1800;            % s: ripples from this long after it (POST)
CONTEXTS = {'ConA1','ConB1','ConA2','ConB2'};
NPERM = 1e4;                % random permutations when exhaustive enumeration is too large
% ---- sequences ----
SEQ_LEN = 3;                % L: 2 = pairwise order, 3, or 4 (4 needs large ensembles)
THETA_WIN = [0.002 0.030];  % s: (min, max] lag between successive spikes of a chain in theta
SEQ_RIP_WIN = [0.002 0.017];% s: the same in ripples
RIP_PAD = 0.050;            % s: ripples in the exposure, +- this, are removed from the theta part
BURST_ISI = 0;              % s: > 0 collapses each burst to its first spike (theta and ripples)
THETA_BIN = 0.125;          % s: theta-cycle-scale bins for template presence
MIN_THETA_BINS = 20;        % informative bins a cell set needs in the exposure
SEL_P = 0.01;               % template selection: one-sided p(z) x L! orderings
MAX_TEMPLATES = 500;
MIN_STRATA = 30;            % per-exposure psi printed only above this (pooling uses every stratum)
STRATA_CAP = 300;           % pooled weight per exposure saturates at this many informative strata (Inf = plain MH)
NBOOT = 1000;               % ripple bootstrap (per exposure) and animal bootstrap (pooled)
ROW_CAP = 5e6;              % chain rows built per block of start spikes (blocks shrink to stay under it)
OUT_CAP = 2e8;              % distinct (set, bin, ordering) rows held (16 B each, ~3x that during unique)
SEED = 1;

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
opt = struct('preWin', PRE_WIN, 'postWin', POST_WIN, 'pairExclude', PAIR_EXCLUDE, 'cells', CELLS, ...
    'thetaRatio', THETA_RATIO, 'speedMin', SPEED_MIN, 'minSpk', MIN_SPK, 'minUnits', MIN_UNITS, ...
    'minRip', MIN_RIP, 'maxRip', MAX_RIP, 'matchN', MATCH_N, 'seqLen', SEQ_LEN, 'thetaWin', THETA_WIN, ...
    'ripWin', SEQ_RIP_WIN, 'ripPad', RIP_PAD, 'burstIsi', BURST_ISI, 'thetaBin', THETA_BIN, ...
    'minThetaBins', MIN_THETA_BINS, 'selP', SEL_P, 'maxTemplates', MAX_TEMPLATES, 'minStrata', MIN_STRATA, ...
    'nBoot', NBOOT, 'rowCap', ROW_CAP, 'outCap', OUT_CAP, 'gap', 10, 'seed', SEED);
TY = {'fwd','rev','rnd','rate','fwd |O|-matched','fwd ~= rate'};

meta = cell(0, 4); R = repmat(local_row(), 0, 1); led = cell(0, 4);
fprintf('L = %d\n%-8s %-36s %-6s %5s %5s | %6s %6s | %7s %7s %8s | %s\n', SEQ_LEN, 'group', 'session', 'ctx', ...
    'units', 'tpl', 'strPRE', 'strPOS', 'psiPRE', 'psiPOST', 'dPsi', 'note');
%% ---- controls ------------------------------------------------------------------------
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), warning('seq:file', 'missing %s', rf); continue; end
    L_ = load(rf, 'RIPPLES'); RP = L_.RIPPLES;
    nov = CONTROLS{i,2}; fam = CONTROLS{i,3}; rec = RP.meta.recordingSec;
    E = struct('name', {'ConA1','ConB1'}, 'on', {fam(1), nov(1)}, 'off', {fam(2), nov(2)}, 'preIv', {[], []});
    if fam(1) < nov(1)
        E(1).pre = [0 fam(1)];     E(1).post = [fam(2) nov(1)];
        E(2).pre = [fam(2) nov(1)]; E(2).post = [nov(2) rec];
        if PRE_B_ALLGAP, E(2).preIv = [fam(2) nov(1)]; end
    else
        E(2).pre = [0 nov(1)];     E(2).post = [nov(2) fam(1)];
        E(1).pre = [nov(2) fam(1)]; E(1).post = [fam(2) rec];
    end
    [~, snL] = fileparts(fld);
    for e = E
        s = local_seq(RP, fld, e, opt);
        meta(end+1,:) = {'Control', controlAnimal(rf), rf, e.name}; R(end+1,1) = s; %#ok<AGROW,SAGROW>
        local_progress('Control', snL, e.name, s);
    end
end

%% ---- TransPlasticity ----------------------------------------------------------------------
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
for i = 1:numel(fa)
    fld = fileparts(fa{i});
    ev = local_events(fld);
    [~, snL] = fileparts(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConA1_start'))
        why_ = 'no event file or annotations'; if ~isempty(ev), why_ = ['no ConA1 in ' ev.source]; end
        led(end+1,:) = {snL, '', 'all', why_}; continue %#ok<AGROW>
    end
    L_ = load(fa{i}, 'RIPPLES'); RP = L_.RIPPLES; m = RP.meta;
    g = local_group(RP, ev);
    if isempty(g), led(end+1,:) = {snL, '', 'all', 'not the first SAL before ASV or the first ASV (not in SALpre/ASV1)'}; continue; end %#ok<AGROW>
    HC = local_hcBlocks(ev, RP, EXCL_INF, INF_EXCL_WIN);
    for c = 1:numel(CONTEXTS)
        on = ev.time(strcmp(ev.description, [CONTEXTS{c} '_start']));
        off = ev.time(strcmp(ev.description, [CONTEXTS{c} '_end']));
        if isempty(on) || isempty(off)
            if c <= 2, led(end+1,:) = {snL, g, CONTEXTS{c}, 'no start/end marker'}; end %#ok<AGROW>
            continue
        end
        e = struct('name', CONTEXTS{c}, 'on', on(1), 'off', off(1), 'pre', [], 'post', [], 'preIv', []);
        kp = find(abs(HC(:,2) - e.on) <= 1, 1); kq = find(abs(HC(:,1) - e.off) <= 1, 1);
        if PRE_B_ALLGAP && strcmp(CONTEXTS{c}, 'ConB1')
            e.preIv = local_gapIv(HC, ev, e.on); kp = [];
            if ~isempty(e.preIv), kp = 0; end
        end
        if isempty(kp) || isempty(kq)
            led(end+1,:) = {snL, g, CONTEXTS{c}, sprintf('no home cage before (%d) / after (%d) the exposure', ~isempty(kp), ~isempty(kq))}; %#ok<AGROW>
            continue
        end
        if kp > 0, e.pre = HC(kp,:); else, e.pre = [e.preIv(1,1) e.preIv(end,2)]; end
        e.post = HC(kq,:);
        s = local_seq(RP, fld, e, opt);
        meta(end+1,:) = {g, m.animal, fa{i}, e.name}; R(end+1,1) = s; %#ok<AGROW,SAGROW>
        local_progress(g, snL, e.name, s);
    end
end
if isempty(R), error('seq:none', 'no context exposures'); end

%% ---- table, ledger --------------------------------------------------------------------------
T = [cell2table(meta, 'VariableNames', {'Group','Animal','Session','Context'}), struct2table(R)];
T.why = cellstr(T.why); T.thetaSrc = cellstr(T.thetaSrc);
fprintf('\n%d context exposures: %d pooled (%d also scored alone)\n', height(T), sum(T.pool), sum(isfinite(T.dPsi)));
for r_ = find(~T.pool)'
    if strcmp(T.Group{r_}, 'Control'), continue; end
    [~, snL] = fileparts(fileparts(T.Session{r_}));
    led(end+1,:) = {snL, T.Group{r_}, T.Context{r_}, T.why{r_}}; %#ok<AGROW>
end
fprintf('\nSESSION LEDGER: exports not contributing (session | group | context | reason)\n');
if isempty(led), fprintf('   none\n'); else
    [~, io] = sortrows(led(:, [2 1 3])); led = led(io, :);
    for r_ = 1:size(led, 1), fprintf('   %-34s | %-6s | %-5s | %s\n', led{r_,:}); end
end
P = T(T.pool, :);
if isempty(P), error('seq:none', 'no exposure reached the ripple scoring'); end
P.w = min(1, STRATA_CAP ./ max(1, (P.nInvPre + P.nInvPost) / 2));
P.w1 = ones(height(P), 1);                                  % plain MH
G = intersect({'Control','SALpre','ASV1','SALpost','ASVrep'}, unique(P.Group), 'stable'); G = G(:)';
fprintf('theta template source: %s\n', strjoin(cellfun(@(s) sprintf('%s %d', s, sum(strcmp(P.thetaSrc, s))), ...
    unique(P.thetaSrc), 'UniformOutput', false), ', '));

%% ---- pooled MH per group x context ------------------------------------------------------------
fprintf(['\nPOOLED MH (strata cap %g) by context and group: psi PRE -> POST, dPsi [95%% CI, bootstrap over animals], ' ...
    'p = PRE/POST swapped within animal (exact when <= 14 animals; min p shown) | dPsi plain MH\n'], STRATA_CAP);
RES = cell(0, 9);
for c = CONTEXTS
    for g = G
        k = strcmp(P.Group, g{1}) & strcmp(P.Context, c{1});
        if ~any(k), continue; end
        fprintf(['%-6s %-8s %d exposures, %d animals | nEff animals %.1f (plain MH %.1f) | strata %d -> %d | ' ...
            'templates %d (fwd = rate order %.2f) | all-L participation %.3f -> %.3f\n'], ...
            c{1}, g{1}, sum(k), numel(unique(P.Animal(k))), local_neff(P(k,:), 'w'), local_neff(P(k,:), 'w1'), ...
            sum(P.nStrPre(k)), sum(P.nStrPost(k)), sum(P.nTpl(k)), sum(P.fwdIsRate(k) .* P.nTpl(k)) / sum(P.nTpl(k)), ...
            mean(P.copPre(k)), mean(P.copPost(k)));
        for q = 1:6
            r = local_poolWithin(P(k,:), q, NPERM, NBOOT);
            if ~isfinite(r.d), continue; end
            fprintf('      %-16s psi %+.3f -> %+.3f   dPsi %+.3f [%+.3f %+.3f]   p = %.3f (min %.3f) | %+.3f\n', ...
                TY{q}, r.pre, r.post, r.d, r.lo, r.hi, r.p, r.pmin, local_D(local_X(P(k,:), q, 'w1')));
            RES(end+1,:) = {c{1}, g{1}, TY{q}, r.pre, r.post, r.d, r.lo, r.hi, r.p}; %#ok<AGROW>
        end
    end
end
RES = cell2table(RES, 'VariableNames', {'Context','Group','Order','psiPre','psiPost','dPsi','lo','hi','p'});

%% ---- novel minus familiar within group (A1 is before the infusion) ---------------------------
fprintf('\nB1 - A1 within group (dPsi; context swapped within animal; CI over animals)\n');
BA = cell(0, 6);
for g = G
    kA = strcmp(P.Group, g{1}) & strcmp(P.Context, 'ConA1'); kB = strcmp(P.Group, g{1}) & strcmp(P.Context, 'ConB1');
    if ~any(kA) || ~any(kB), continue; end
    for q = [1 3 4 5 6]
        r = local_poolBetween(P, kA, kB, q, 'Animal', [], NPERM, NBOOT);
        if ~isfinite(r.d), continue; end
        fprintf('   %-8s %-16s %+.3f [%+.3f %+.3f]  p = %.3f  (%s)\n', g{1}, TY{q}, r.d, r.lo, r.hi, r.p, r.design);
        BA(end+1,:) = {g{1}, TY{q}, r.d, r.lo, r.hi, r.p}; %#ok<AGROW>
    end
end
BA = cell2table(BA, 'VariableNames', {'Group','Order','d','lo','hi','p'});

%% ---- between groups ---------------------------------------------------------------------------
prs = {'SALpre','ASV1'; 'SALpre','SALpost'; 'Control','SALpre'; 'ASV1','ASVrep'};
fprintf('\nbetween groups: dPsi(second) - dPsi(first) [95%% CI over animals], p\n');
for j = 1:size(prs, 1)
    if ~all(ismember(prs(j,:), P.Group)), continue; end
    for q = [1 3 4 5 6]
        for c = {'ConA1','ConB1'}
            kA = strcmp(P.Group, prs{j,1}) & strcmp(P.Context, c{1}); kB = strcmp(P.Group, prs{j,2}) & strcmp(P.Context, c{1});
            if ~any(kA) || ~any(kB), continue; end
            r = local_poolBetween(P, kA, kB, q, 'Animal', [], NPERM, NBOOT);
            fprintf('   %-8s vs %-8s %-6s %-16s %+.3f [%+.3f %+.3f]  p = %.3f  (%s)\n', prs{j,:}, c{1}, TY{q}, ...
                r.d, r.lo, r.hi, r.p, r.design);
        end
        isAB = ismember(P.Context, {'ConA1','ConB1'});
        kA = strcmp(P.Group, prs{j,1}) & isAB; kB = strcmp(P.Group, prs{j,2}) & isAB;
        sub = 1 + strcmp(P.Context, 'ConB1');
        r = local_poolBetween(P, kA, kB, q, 'Animal', sub, NPERM, NBOOT);
        fprintf('   %-8s vs %-8s %-6s %-16s %+.3f [%+.3f %+.3f]  p = %.3f  (%s)\n', prs{j,:}, 'B1-A1', TY{q}, ...
            r.d, r.lo, r.hi, r.p, r.design);
    end
end

%% ---- figure ---------------------------------------------------------------------------------
col = containers.Map({'Control','SALpre','ASV1','SALpost','ASVrep'}, ...
    {[0.4 0.4 0.4], [0.2 0.6 1], [0.85 0.2 0.1], [0 0.25 0.6], [0.55 0.1 0.6]});
mk = containers.Map({'fwd','rnd','rate'}, {'o','s','^'}); off_ = containers.Map({'fwd','rnd','rate'}, {-0.2, 0, 0.2});
fig = figure('Color','w','Position',[60 60 1300 380]);
tiledlayout(1, 3, 'TileSpacing','compact');
for c = 1:3
    nexttile; hold on
    for j = 1:numel(G)
        for o_ = {'fwd','rnd','rate'}
            if c < 3, k = strcmp(RES.Context, CONTEXTS{c}) & strcmp(RES.Group, G{j}) & strcmp(RES.Order, o_{1}); X_ = RES;
            else, k = strcmp(BA.Group, G{j}) & strcmp(BA.Order, o_{1}); X_ = BA; end
            if ~any(k), continue; end
            if c < 3, y = X_.dPsi(k); else, y = X_.d(k); end
            x = j + off_(o_{1}); fc = col(G{j}); if ~strcmp(o_{1}, 'fwd'), fc = 'w'; end
            plot([x x], [X_.lo(k) X_.hi(k)], '-', 'Color', col(G{j}));
            plot(x, y, mk(o_{1}), 'Color', col(G{j}), 'MarkerFaceColor', fc, 'MarkerSize', 6);
        end
    end
    yline(0, 'k:'); xlim([0.5 numel(G)+0.5]); xticks(1:numel(G)); xticklabels(G);
    if c < 3, title(sprintf('%s: pooled \\Delta\\psi (o template, s random, ^ rate)', CONTEXTS{c})); ylabel('\psi_{POST} - \psi_{PRE}');
    else, title('B1 - A1, same session'); ylabel('\Delta\Delta\psi'); end
end
exportgraphics(fig, fullfile(CACHE_ROOT, sprintf('compareSequenceReactivation_L%d.pdf', SEQ_LEN)), 'ContentType','vector');
save(fullfile(CACHE_ROOT, sprintf('compareSequenceReactivation_L%d.mat', SEQ_LEN)), 'T', 'RES', 'BA', 'opt', 'CONTROLS');


%% =========================================================================
function s = local_row()
s = struct('nUnits', 0, 'nTpl', 0, 'tplZ', NaN, 'thetaFrac', NaN, 'thetaSrc', '', 'thetaHz', NaN, 'nRipExp', 0, ...
    'nRipPre', 0, 'nRipPost', 0, 'nStrPre', 0, 'nStrPost', 0, 'psiPre', NaN, 'psiPost', NaN, 'dPsi', NaN, ...
    'dPsiLo', NaN, 'dPsiHi', NaN, 'dPsiRev', NaN, 'dPsiRnd', NaN, 'dPsiRate', NaN, 'dPsiK', NaN, ...
    'kPre', NaN, 'kPost', NaN, 'satPre', NaN, 'satPost', NaN, 'copPre', NaN, 'copPost', NaN, ...
    'nInvPre', 0, 'nInvPost', 0, 'fwdIsRate', NaN, ...
    'mh', zeros(1, 24), 'pool', false, 'tplCells', {{zeros(0, 4)}}, 'why', '');
end


function local_progress(g, sn, ctx, s)
fprintf('%-8s %-36s %-6s %5d %5d | %6d %6d | %+7.2f %+7.2f %+8.3f | %s\n', g, sn, ctx, s.nUnits, s.nTpl, ...
    s.nStrPre, s.nStrPost, s.psiPre, s.psiPost, s.dPsi, s.why);
end


function s = local_seq(RP, fld, e, o)
% one context exposure: theta templates, PRE / POST ripple scoring, MH sums
s = local_row();
rng(o.seed);
L = o.seqLen;
sf = dir(fullfile(fld, '*spikes.cellinfo.mat')); cf = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(sf) || isempty(cf), s.why = 'no spikes or cell_metrics'; return; end
S_ = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); sp = S_.spikes;
C_ = load(fullfile(cf(1).folder, cf(1).name), 'cell_metrics'); cm = C_.cell_metrics;
ts = cellfun(@(x) sort(reshape(x(~isnan(x)), [], 1)), sp.times(:), 'UniformOutput', false);
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
si = RP.units.SpikesIdx(:);
if max(si) > nU, s.why = sprintf('RIPPLES indexes unit %d, spikes file has %d (different sort)', max(si), nU); return; end
uMap = zeros(nU, 1); uMap(si) = 1:numel(si);

% ---- theta part of the exposure, minus the ripples inside it ----
[mask, src] = local_runMask(fld, e.on, e.off, o.speedMin);
if isempty(src), [mask, src] = local_thetaMask(RP, e.on, e.off, o.thetaRatio); end
s.thetaSrc = src;
if isempty(mask), s.why = 'no theta in exposure'; return; end
iv = local_merge(mask);
kx = RP.ts1(:) < e.off & RP.ts1Off(:) > e.on; s.nRipExp = sum(kx);
if any(kx), iv = local_subtract(iv, [RP.ts1(kx) - o.ripPad, RP.ts1Off(kx) + o.ripPad]); end
s.thetaFrac = sum(diff(iv, 1, 2)) / (e.off - e.on);
if isempty(iv), s.why = 'no theta left after removing ripples'; return; end
tRun = cellfun(@(t) t(local_inIv(t, iv)), ts, 'UniformOutput', false);
nRun = cellfun(@numel, tRun);
isC = true(nU, 1); if strcmp(o.cells, 'PYR'), isC = contains(typ, 'Pyramidal'); end
keep = nRun >= o.minSpk & uMap > 0 & isC;
idx = find(keep); n = numel(idx);
if n < max(o.minUnits, L)
    s.why = sprintf('< %d units (%d sorted, %d in RIPPLES, %d %s, %d with >= %d theta spikes)', max(o.minUnits, L), ...
        nU, sum(uMap > 0), sum(uMap > 0 & isC), o.cells, n, o.minSpk);
    return
end
s.nUnits = n;
lab = zeros(nU, 1); lab(idx) = 1:n;                       % spikes-file unit -> compact label
rate = nRun(idx) / sum(diff(iv, 1, 2));
tt = vertcat(tRun{idx}); uu = repelem((1:n)', nRun(idx));
[tt, so] = sort(tt); uu = uu(so);
kb = local_deburst(tt, uu, o.burstIsi); tt = tt(kb); uu = uu(kb);
s.thetaHz = numel(tt) / sum(diff(iv, 1, 2));
bid = local_binId(tt, iv, o.thetaBin);

% ---- templates ----
try
    [tpl, why] = local_templates(tt, uu, bid, n, grp(idx), rate, o);
catch ME
    s.why = sprintf('theta chains (%.0f Hz population): %s', s.thetaHz, ME.message); return
end
if ~isempty(why), s.why = why; return; end
s.nTpl = size(tpl.S, 1); s.tplZ = median(tpl.z);
s.tplCells = {reshape(idx(tpl.cellsFwd), size(tpl.cellsFwd))};   % spikes-file units, template order

% ---- PRE / POST ripples (as compareReactivation) ----
preIv = local_preIv(e, o.preWin);
if ~isfield(e, 'postFixed') || ~e.postFixed
    e.post = [max(e.post(1), e.off), min(e.post(2), e.off + o.postWin)];
end
inPre = find(local_ripInIv(RP, preIv));
inPost = find(RP.ts1(:) >= e.post(1) & RP.ts1Off(:) <= e.post(2));
if min(numel(inPre), numel(inPost)) < o.minRip, s.why = sprintf('< %d ripples in pre or post', o.minRip); return; end
if o.matchN
    N = min([numel(inPre), numel(inPost), o.maxRip]);
    inPre = inPre(end-N+1:end); inPost = inPost(1:N);
end
nA = numel(inPre); nB = numel(inPost); s.nRipPre = nA; s.nRipPost = nB;
[MA, copA, why] = local_ripMasks(RP, inPre, si, lab, tpl, n, o);
if isempty(why), [MB, copB, why] = local_ripMasks(RP, inPost, si, lab, tpl, n, o); end
if ~isempty(why), s.why = ['ripple chains: ' why]; return; end
nP = tpl.nP;
s.copPre = mean(copA); s.copPost = mean(copB);
kA = local_popcount(MA.mask, nP); kB = local_popcount(MB.mask, nP);
s.kPre = mean(kA); s.kPost = mean(kB); s.satPre = mean(kA == nP); s.satPost = mean(kB == nP);
s.nInvPre = sum(kA <= nP - 1); s.nInvPost = sum(kB <= nP - 1);     % relabelling-invariant weight basis
s.fwdIsRate = mean(tpl.fwd == tpl.rate);

% ---- Mantel-Haenszel sums: fwd, rev, rnd, rate, |O|-matched fwd, fwd ~= rate ----
selA = false(size(kA)); selB = false(size(kB));
for k = 1:nP-1
    ia = find(kA == k); ib = find(kB == k); mk = min(numel(ia), numel(ib));
    if mk == 0, continue; end
    selA(ia(randperm(numel(ia), mk))) = true; selB(ib(randperm(numel(ib), mk))) = true;
end
ord = {tpl.fwd, tpl.rev, tpl.rnd, tpl.rate, tpl.fwd, tpl.fwd .* (tpl.fwd ~= tpl.rate)};
mh = zeros(1, 24); parts = cell(6, 4);
for q = 1:6
    mA = MA; mB = MB;
    if q == 5, mA = local_sub(MA, selA); mB = local_sub(MB, selB); end
    [parts{q,1}, parts{q,2}, nIa] = local_mhParts(mA, ord{q}, tpl, nA);
    [parts{q,3}, parts{q,4}, nIb] = local_mhParts(mB, ord{q}, tpl, nB);
    mh((q-1)*4 + (1:4)) = cellfun(@sum, parts(q, :));
    if q == 1, s.nStrPre = nIa; s.nStrPost = nIb; end
end
s.mh = mh; s.pool = true;
if min(s.nStrPre, s.nStrPost) < o.minStrata
    s.why = sprintf('pooled only (< %d strata alone)', o.minStrata); return
end
lor = @(nu, de) log(sum(nu) / sum(de));
s.psiPre = lor(parts{1,1}, parts{1,2}); s.psiPost = lor(parts{1,3}, parts{1,4}); s.dPsi = s.psiPost - s.psiPre;
dq = @(q) log(mh((q-1)*4 + 3) / mh((q-1)*4 + 4)) - log(mh((q-1)*4 + 1) / mh((q-1)*4 + 2));
s.dPsiRev = dq(2); s.dPsiRnd = dq(3); s.dPsiRate = dq(4); s.dPsiK = dq(5);
bs = nan(o.nBoot, 1);
for b = 1:o.nBoot
    ia = randi(nA, nA, 1); ib = randi(nB, nB, 1);
    bs(b) = lor(parts{1,3}(ib), parts{1,4}(ib)) - lor(parts{1,1}(ia), parts{1,2}(ia));
end
bs = bs(isfinite(bs));
if ~isempty(bs), s.dPsiLo = prctile(bs, 2.5); s.dPsiHi = prctile(bs, 97.5); end
end


function [tpl, why] = local_templates(t, u, bid, n, grpC, rate, o)
% order-specific L-cell sets in the theta part (see header)
tpl = []; why = '';
L = o.seqLen; [PT, LUT, REV, nP] = local_permTables(L);
nBin = max(bid);
[key, oi] = local_chainOrders(t, u, bid, o.thetaWin, [], o, n, nBin, LUT);
if isempty(key), why = 'no theta chains'; return; end
[kk, ~, gk] = unique(key);
msk = accumarray(gk, 2 .^ (oi - 1));                      % (set, bin, ordering) unique: sum = OR
[uc, ~, gs] = unique(floor((kk - 1) / nBin) + 1);
nS = numel(uc); Sset = local_setDecode(uc, n, L);
k = local_popcount(msk, nP); inf_ = k >= 1 & k <= nP - 1;
g2 = gs(inf_); m2 = msk(inf_); pb = k(inf_) / nP;
Nord = zeros(nS, nP);
for p = 1:nP, Nord(:, p) = accumarray(g2, bitget(m2, p), [nS 1]); end
Ex = accumarray(g2, pb, [nS 1]); Vx = accumarray(g2, pb .* (1 - pb), [nS 1]); nb = accumarray(g2, 1, [nS 1]);
[nMax, fwd] = max(Nord, [], 2);
z = (nMax - Ex) ./ sqrt(Vx); pz = 0.5 * erfc(z / sqrt(2));
sel = nb >= o.minThetaBins & Vx > 0 & nP * pz < o.selP;
g = reshape(grpC(Sset), size(Sset)); pr = nchoosek(1:L, 2); bad = false(nS, 1);
for q = 1:size(pr, 1), bad = bad | (g(:, pr(q,1)) == g(:, pr(q,2)) & isfinite(g(:, pr(q,1)))); end
sel = sel & ~bad;
if ~any(sel), why = sprintf('no order-specific theta sets (%d sets with >= %d informative bins)', ...
        sum(nb >= o.minThetaBins), o.minThetaBins); return; end
j = find(sel); [~, so] = sort(z(j), 'descend'); j = j(so(1:min(end, o.maxTemplates)));
nT = numel(j);
tpl.S = Sset(j, :); tpl.z = z(j); tpl.nb = nb(j); tpl.code = uc(j);
tpl.fwd = fwd(j); tpl.rev = REV(tpl.fwd); tpl.nP = nP; tpl.LUT = LUT;
rL = reshape(rate(tpl.S), size(tpl.S)); [~, ro] = sort(rL, 2, 'descend');   % highest exposure rate first
tpl.rate = LUT(local_pcode(ro, L));
tpl.rnd = zeros(nT, 1); tpl.cellsFwd = zeros(nT, L);
for q = 1:nT
    cand = setdiff(1:nP, [tpl.fwd(q) tpl.rev(q) tpl.rate(q)]);
    if ~isempty(cand), tpl.rnd(q) = cand(randi(numel(cand))); end
    tpl.cellsFwd(q, :) = tpl.S(q, PT(tpl.fwd(q), :));
end
tpl.R = []; tpl.NR = [];                                   % L = 2: reference = the other ordering
if L > 2
    tpl.R = zeros(nT, 1); tpl.NR = zeros(nT, 1);
    for q = 1:nT
        ex = unique([tpl.fwd(q) tpl.rev(q) tpl.rnd(q) tpl.rate(q)]); ex(ex == 0) = [];
        tpl.R(q) = (2^nP - 1) - sum(2 .^ (ex - 1)); tpl.NR(q) = nP - numel(ex);
    end
end
end


function [M, cop, why] = local_ripMasks(RP, rk, si, lab, tpl, n, o)
% (template, ripple, bitmask of the L! orderings present) for ripples rk; cop = per template, fraction
% of ripples in which all L cells fire
L = o.seqLen;
M = struct('tpl', zeros(0, 1), 'rip', zeros(0, 1), 'mask', zeros(0, 1)); why = '';
nT = size(tpl.S, 1); nR = numel(rk); cop = zeros(nT, 1);
labRP = lab(si);                                            % RIPPLES unit row -> compact label
qrows = find(ismember(labRP, unique(tpl.S(:))));
ts1 = RP.ts1(:);
sp = RP.spikes(ismember(RP.spikes.Ripple, rk) & ismember(RP.spikes.Unit, qrows), :);
if height(sp) == 0, return; end
[~, ri] = ismember(sp.Ripple, rk);
uu = labRP(sp.Unit);
t = sp.ts1 - ts1(sp.Ripple) + ri * o.gap;                   % each ripple in its own block: no chain crosses
Pm = sparse(ri, uu, 1, nR, n) > 0;
allL = Pm(:, tpl.S(:,1));
for k = 2:L, allL = allL & Pm(:, tpl.S(:,k)); end
cop = full(sum(allL, 1))' / nR;
[t, so] = sort(t); uu = uu(so); ri = ri(so);
kb = local_deburst(t, uu, o.burstIsi); t = t(kb); uu = uu(kb); ri = ri(kb);
allowed = false(n);
for j = 1:nT, c = tpl.S(j, :); allowed(c, c) = true; end
try
    [key, oi] = local_chainOrders(t, uu, ri, o.ripWin, allowed, o, n, nR, tpl.LUT);
catch ME
    why = ME.message; return
end
if isempty(key), return; end
setc = floor((key - 1) / nR) + 1; rip = mod(key - 1, nR) + 1;
[tf, ti] = ismember(setc, tpl.code);
if ~any(tf), return; end
[kk, ~, g] = unique([ti(tf), rip(tf)], 'rows');
M.tpl = kk(:, 1); M.rip = kk(:, 2); M.mask = accumarray(g, 2 .^ (oi(tf) - 1));
end


function [num, den, nInf] = local_mhParts(M, tbit, tpl, nRip)
% Mantel-Haenszel numerator / denominator summed per ripple, strata = (template, ripple):
% a = tested ordering present, c = reference orderings present, N = #reference + 1
num = zeros(nRip, 1); den = zeros(nRip, 1); nInf = 0;
if isempty(M.tpl), return; end
tb = tbit(M.tpl); ok = tb > 0;                              % tbit 0: ordering undefined (rnd at L = 2)
if ~any(ok), return; end
tb = tb(ok); msk = M.mask(ok); rip = M.rip(ok); ti = M.tpl(ok);
if isempty(tpl.R), rm = (2^tpl.nP - 1) - 2 .^ (tb - 1); nr = (tpl.nP - 1) * ones(size(tb));
else, rm = tpl.R(ti); nr = tpl.NR(ti); end
a = bitget(msk, tb); c = local_popcount(bitand(msk, rm), tpl.nP);
N = nr + 1; b = 1 - a; d = nr - c;
num = accumarray(rip, a .* d ./ N, [nRip 1]); den = accumarray(rip, b .* c ./ N, [nRip 1]);
nInf = sum(a + c > 0 & b + d > 0);
end


function M = local_sub(M, sel)
M.tpl = M.tpl(sel); M.rip = M.rip(sel); M.mask = M.mask(sel);
end


function [key, oi] = local_chainOrders(t, u, bid, win, allowed, o, n, nBin, LUT)
% unique (key, ordering) over all L-chains of distinct cells with successive lags in (win(1), win(2)];
% key = (set code - 1) * nBin + bin of the first spike. t sorted; u compact labels 1..n;
% allowed = n x n pairs that may share a chain ([] = all)
L = o.seqLen; key = zeros(0, 1); oi = zeros(0, 1);
nS = numel(t); if nS < L, return; end
[lo, hi] = local_succ(t, win(1), win(2));
blk = 4096; s0 = 1; acc = {}; accN = 0;
while s0 <= nS
    s1 = min(nS, s0 + blk - 1);
    [P, ok] = local_chainBlock((s0:s1)', lo, hi, u, allowed, o.rowCap, L);
    if ~ok
        if blk == 1, error('seq:cap', 'one start spike exceeds %g chain rows; shrink the window', o.rowCap); end
        blk = max(1, floor(blk / 2)); continue
    end
    if ~isempty(P)
        U = zeros(size(P)); for k = 1:L, U(:, k) = u(P(:, k)); end
        [S, q] = local_setOrder(U, LUT);
        acc{end+1, 1} = unique([(local_setCode(S, n) - 1) * nBin + bid(P(:, 1)), q], 'rows'); %#ok<AGROW>
        accN = accN + size(acc{end}, 1);
        if accN > o.outCap
            X = unique(vertcat(acc{:}), 'rows'); acc = {X}; accN = size(X, 1);
            if accN > o.outCap
                error('seq:cap', '> %.3g distinct (set, bin, ordering) rows by spike %d of %d', o.outCap, s1, nS);
            end
        end
    end
    s0 = s1 + 1;
end
if isempty(acc), return; end
X = unique(vertcat(acc{:}), 'rows'); key = X(:, 1); oi = X(:, 2);
end


function [P, ok] = local_chainBlock(st, lo, hi, u, allowed, rowCap, L)
% all L-chains of distinct cells starting at spikes st; rows = spike indices
useA = ~isempty(allowed); nA = size(allowed, 1);
P = st(:); ok = true;
for lev = 2:L
    [r, c, ok] = local_expand(P(:, end), lo, hi, rowCap);
    if ~ok, P = zeros(0, L); return; end
    P = [P(r, :), c]; un = u(c); k = true(size(c));
    for j = 1:lev-1
        uj = u(P(:, j)); k = k & un ~= uj;
        if useA, k = k & allowed((un - 1)*nA + uj); end
    end
    P = P(k, :);
    if isempty(P), P = zeros(0, L); return; end
end
end


function [r, c, ok] = local_expand(p, lo, hi, rowCap)
% every successor c of every parent spike p(r)
p = p(:); r = zeros(0, 1); c = zeros(0, 1);
ns = max(hi(p) - lo(p) + 1, 0); tot = sum(ns); ok = tot <= rowCap;
if ~ok || tot == 0, return; end
r = repelem((1:numel(p))', ns);
first = cumsum(ns) - ns;                                    % 0-based start of each parent's block
c = lo(p(r)) + ((1:tot)' - 1 - first(r));
end


function [lo, hi] = local_succ(t, a, b)
% successors of spike s: indices lo(s)..hi(s), the spikes with t(s) + a < t <= t(s) + b
n = numel(t); lo = zeros(n, 1); hi = zeros(n, 1); j1 = 1; j2 = 0;
for s = 1:n
    j1 = max(j1, s + 1);
    while j1 <= n && t(j1) - t(s) <= a, j1 = j1 + 1; end
    j2 = max(j2, j1 - 1);
    while j2 < n && t(j2 + 1) - t(s) <= b, j2 = j2 + 1; end
    lo(s) = j1; hi(s) = j2;
end
end


function keep = local_deburst(t, u, isi)
% false for a spike < isi after the same cell's previous spike
keep = true(size(t));
if isi <= 0 || numel(t) < 2, return; end
[~, o] = sortrows([u(:) t(:)]); us = u(o); tsr = t(o);
dup = [false; us(2:end) == us(1:end-1) & diff(tsr) < isi];
keep(o(dup)) = false;
end


function b = local_binId(t, iv, w)
% global index of the w-wide bin of t within the merged intervals iv (t inside iv)
q = discretize(t, [iv(:, 1); Inf]);
nb = ceil(diff(iv, 1, 2) / w); off = [0; cumsum(nb(1:end-1))];
b = off(q) + floor((t - iv(q, 1)) / w) + 1;
end


function iv = local_subtract(iv, cut)
% intervals iv (merged, sorted) minus the union of cut
cut = local_merge(cut); out = zeros(0, 2);
for k = 1:size(iv, 1)
    a = iv(k, 1); b = iv(k, 2);
    c = cut(cut(:, 2) > a & cut(:, 1) < b, :);
    for j = 1:size(c, 1)
        if c(j, 1) > a, out(end+1, :) = [a c(j, 1)]; end %#ok<AGROW>
        a = max(a, c(j, 2));
    end
    if b > a, out(end+1, :) = [a b]; end %#ok<AGROW>
end
iv = out;
end


function m = local_inIv(t, iv)
% t in [iv(k,1), iv(k,2)) for some k; iv sorted, disjoint
m = false(size(t));
if isempty(iv) || isempty(t), return; end
e = reshape(iv', [], 1); q = discretize(t, e);
m = mod(q, 2) == 1; m(t >= e(end)) = false;
end


function [S, oi] = local_setOrder(U, LUT)
% S = cells sorted; oi = ordering index: U(k,:) = S(k, PT(oi(k), :))
L = size(U, 2); S = sort(U, 2); P = zeros(size(U));
for k = 1:L, P(:, k) = sum(S < U(:, k), 2) + 1; end
oi = LUT(local_pcode(P, L));
end


function c = local_pcode(P, L)
c = zeros(size(P, 1), 1); for k = 1:L, c = c * L + (P(:, k) - 1); end; c = c + 1;
end


function c = local_setCode(S, n)
c = zeros(size(S, 1), 1); for k = 1:size(S, 2), c = c * n + (S(:, k) - 1); end; c = c + 1;
end


function S = local_setDecode(c, n, L)
c = c(:) - 1; S = zeros(numel(c), L);
for k = L:-1:1, S(:, k) = mod(c, n) + 1; c = floor(c / n); end
end


function [PT, LUT, REV, nP] = local_permTables(L)
PT = perms(1:L); nP = size(PT, 1);
LUT = zeros(L^L, 1); LUT(local_pcode(PT, L)) = 1:nP;
REV = LUT(local_pcode(fliplr(PT), L));
end


function c = local_popcount(m, nb)
c = zeros(size(m)); for b = 1:nb, c = c + bitget(m, b); end
end


%% ---- pooled inference ------------------------------------------------------------------------
function d = local_D(Y)
% pooled dPsi from rows [numPre denPre numPost denPost]
s = sum(Y, 1); d = log(s(3) / s(4)) - log(s(1) / s(2));
end


function X = local_X(T, q, wf)
% exposure MH sums [numPre denPre numPost denPost] for ordering q, scaled by weight column wf
X = T.mh(:, (q-1)*4 + (1:4)) .* T.(wf);
end


function ne = local_neff(T, wf)
% Kish effective number of animals: s_a = weighted informative strata of animal a
[~, ~, ia] = unique(T.Animal);
s = accumarray(ia, T.(wf) .* (T.nInvPre + T.nInvPost));
ne = sum(s)^2 / sum(s.^2);
end


function [Sg, exact] = local_flips(n, nperm)
exact = n <= 14;
if exact, Sg = dec2bin(0:2^n-1, n) == '1'; else, Sg = rand(nperm, n) > 0.5; end
end


function p = local_p(null, obs, exact)
tol = 1e-12 * max(1, abs(obs));
if exact, p = mean(abs(null) >= abs(obs) - tol);
else, p = (1 + sum(abs(null) >= abs(obs) - tol)) / (1 + numel(null)); end
end


function r = local_poolWithin(Tk, q, nperm, nboot)
% pooled dPsi of the exposures Tk; p: PRE / POST swapped within animal; CI: animals resampled
r = struct('pre', NaN, 'post', NaN, 'd', NaN, 'lo', NaN, 'hi', NaN, 'p', NaN, 'pmin', NaN);
X = local_X(Tk, q, 'w');
[~, ~, ia] = unique(Tk.Animal); nA = max(ia);
A = zeros(nA, 4); for j = 1:nA, A(j, :) = sum(X(ia == j, :), 1); end
s = sum(A, 1); r.pre = log(s(1) / s(2)); r.post = log(s(3) / s(4)); r.d = local_D(A);
if ~isfinite(r.d), return; end
[Sg, exact] = local_flips(nA, nperm); null = nan(size(Sg, 1), 1);
for k = 1:size(Sg, 1)
    Aw = A; sw = Sg(k, :); Aw(sw, :) = A(sw, [3 4 1 2]); null(k) = local_D(Aw);
end
r.p = local_p(null, r.d, exact); r.pmin = local_p(null, max(abs(null)), exact);
bs = nan(nboot, 1);
for b = 1:nboot, bs(b) = local_D(A(randi(nA, nA, 1), :)); end
bs = bs(isfinite(bs));
if ~isempty(bs), r.lo = prctile(bs, 2.5); r.hi = prctile(bs, 97.5); end
end


function r = local_poolBetween(P, kA, kB, q, unitField, sub, nperm, nboot)
% stat = D(side B) - D(side A); D = pooled dPsi, or with sub (1 / 2 per row) D(sub 2) - D(sub 1).
% units (unitField) on both sides: sides swapped within unit; otherwise units permuted between sides.
% CI: animals resampled within {both sides, A only, B only}.
r = struct('d', NaN, 'lo', NaN, 'hi', NaN, 'p', NaN, 'design', '');
X = local_X(P, q, 'w');
if isempty(sub), sub = ones(height(P), 1); nsub = 1; else, nsub = 2; end
un = P.(unitField); uu = unique(un(kA | kB)); nU = numel(uu);
S = zeros(nU, 4, 2, nsub); has = false(nU, 2); ani = cell(nU, 1); kS = {kA, kB};
for j = 1:nU
    m = strcmp(un, uu{j}); ani{j} = P.Animal{find(m & (kA | kB), 1)};
    for sd = 1:2
        mm = m & kS{sd}; has(j, sd) = any(mm);
        for b = 1:nsub, S(j, :, sd, b) = sum(X(mm & sub == b, :), 1); end
    end
end
f = @(S) local_statB(S, nsub);
r.d = f(S);
if ~isfinite(r.d), return; end
pd = find(has(:,1) & has(:,2));
if numel(pd) >= 2
    [Sg, exact] = local_flips(numel(pd), nperm); null = nan(size(Sg, 1), 1);
    for k = 1:size(Sg, 1)
        Sw = S; ix = pd(Sg(k, :)); Sw(ix, :, :, :) = S(ix, :, [2 1], :); null(k) = f(Sw);
    end
    r.design = sprintf('paired, %d %ss', numel(pd), lower(unitField));
else
    fr = find(~(has(:,1) & has(:,2))); labB = has(fr, 2); Tot = S(fr, :, 1, :) + S(fr, :, 2, :);
    exact = nchoosek(numel(fr), sum(labB)) <= nperm;
    if exact, C = nchoosek(1:numel(fr), sum(labB)); else, C = zeros(nperm, sum(labB));
        for k = 1:nperm, pp = randperm(numel(fr)); C(k, :) = pp(1:sum(labB)); end
    end
    null = nan(size(C, 1), 1);
    for k = 1:size(C, 1)
        lb = false(numel(fr), 1); lb(C(k, :)) = true;
        Sw = S; Sw(fr, :, :, :) = 0;
        Sw(fr(~lb), :, 1, :) = Tot(~lb, :, 1, :); Sw(fr(lb), :, 2, :) = Tot(lb, :, 1, :);
        null(k) = f(Sw);
    end
    r.design = sprintf('unpaired, %d vs %d %ss', sum(~labB), sum(labB), lower(unitField));
end
r.p = local_p(null, r.d, exact);
[ua, ~, ja] = unique(ani); stA = zeros(numel(ua), 1);
for a = 1:numel(ua), stA(a) = 1*any(has(ja == a, 1)) + 2*any(has(ja == a, 2)); end
bs = nan(nboot, 1);
for b = 1:nboot
    pick = [];
    for st = 1:3
        ca = find(stA == st);
        if ~isempty(ca), pick = [pick; ca(randi(numel(ca), numel(ca), 1))]; end %#ok<AGROW>
    end
    ix = cell2mat(arrayfun(@(a) find(ja == a), pick, 'UniformOutput', false));
    bs(b) = f(S(ix, :, :, :));
end
bs = bs(isfinite(bs));
if ~isempty(bs), r.lo = prctile(bs, 2.5); r.hi = prctile(bs, 97.5); end
end


function v = local_statB(S, nsub)
v = zeros(1, 2);
for sd = 1:2
    if nsub == 1, v(sd) = local_D(S(:, :, sd, 1));
    else, v(sd) = local_D(S(:, :, sd, 2)) - local_D(S(:, :, sd, 1)); end
end
v = v(2) - v(1);
end


%% ---- =CR= helpers, unchanged from compareReactivation.m --------------------------------------
function iv = local_merge(mask)
% union of [start end] rows
iv = zeros(0, 2);
if isempty(mask), return; end
mask = sortrows(mask); iv = mask(1,:);
for k = 2:size(mask, 1)
    if mask(k,1) <= iv(end,2), iv(end,2) = max(iv(end,2), mask(k,2)); else, iv(end+1,:) = mask(k,:); end %#ok<AGROW>
end
end


function [mask, src] = local_runMask(fld, a, b, vmin)
% 1 s windows inside [a b] whose mean speed > vmin; src = '' when there is no tracking file
mask = zeros(0, 2); src = '';
d = dir(fullfile(fld, '*Behavior.mat'));
if isempty(d), return; end
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
% (see compareReactivation.m)
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
% home cage between the end of the last context before 'on' and 'on'
lab = ev.description;
pe = ev.time(startsWith(lab, 'Con') & endsWith(lab, '_end') & ev.time <= on + 1);
prevEnd = 0; if ~isempty(pe), prevEnd = max(pe); end
IV = [max(HC(:,1), prevEnd), min(HC(:,2), on)];
IV = IV(diff(IV, 1, 2) > 0, :);
end


function HC = local_hcBlocks(ev, RP, exclInf, infWin)
% Home-cage intervals = [0, end of recording] minus every ConA*/ConB* interval, minus the infusion
% exclusion window, and minus the marked infusion interval when exclInf.
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
% SALpre | ASV1 | ''
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