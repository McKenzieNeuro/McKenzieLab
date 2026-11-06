%% compareReactivation.m
% Do coactivity patterns from a context exposure recur in the ripples that
% follow it more than in the ripples that preceded it?  Familiar context A
% and novel context B, ASV vs SAL, plus never-infused controls.
%
% SESSION LAYOUT (TransPlasticity, from the .evt.sti markers)
%   HC  ConA1  HC  Infusion  HC  ConB1  HC  ConA2  HC  ConB2  HC
%   The infusion falls between A1 and B1, so A1 is drug-free (within-session
%   baseline) and B1 is encoded with the drug on board.
% CONTROLS (Huszar): familiar and novel [on off] hard-coded below (CONTROLS).
%
% FOR EACH CONTEXT EXPOSURE X
%   template   pairwise Pearson correlations of spike counts in BIN bins during X,
%              restricted to running (1 s windows with mean speed > SPEED_MIN, from
%              <session>.Behavior.mat, used whenever it exists: the controls), else to
%              theta (theta/delta power ratio > THETA_RATIO in 1 s
%              windows, LFP channel from RIPPLES.meta) when an LFP is available,
%              else the whole exposure (thetaSrc says which)
%   PRE / POST pairwise correlations of per-ripple spike counts in the home-cage
%              epoch immediately before / after X, limited to PRE_WIN (1 h) before the exposure
%              starts and POST_WIN (30 min) after it ends; MATCH_N = false uses every ripple
%              in each block, true the same number in both (the last N before, the first N after)
%   rPre, rPost  similarity (Pearson over pairs) of template with PRE / POST
%   EV  = ((rPost - rPre*rPP) / sqrt((1-rPre^2)(1-rPP^2)))^2      (Kudrimoti et al. 1999)
%   REV = ((rPre - rPost*rPP) / sqrt((1-rPost^2)(1-rPP^2)))^2     rPP = corr(PRE, POST)
%   dEV = EV - REV;  dR = rPost - rPre
% TEMPORAL BIAS (Skaggs & McNaughton 1996): for each pair (i,j),
%   bias = (n(0,w] - n[-w,0)) / (n(0,w] + n[-w,0)), n = spikes of j at lag from spikes of i,
%   w = RUN_WIN during the exposure (same running/theta mask) and RIP_WIN in ripples (lags
%   only within the same ripple), PRE / POST ripples as above. Pairs need >= MIN_CC_RUN
%   counts in behaviour and >= MIN_CC_RIP in both PRE and POST.
%   bPre, bPost = Pearson over pairs of behaviour bias with PRE / POST ripple bias;
%   dB = bPost - bPre. Pooled over a group's pairs as well (biasPairs).
% Units: CELLS ('PYR' | 'all'), >= MIN_SPK spikes in the exposure; pairs on
% the same peak channel are dropped (PAIR_EXCLUDE; shared-
% electrode artefacts inflate correlations).
%
% GROUPS (TransPlasticity sessions by their infusion)
%   SAL   SESSION_SET 'all': every saline session (first and repeat, before or after ASV)
%         SESSION_SET 'first': saline sessions before the animal's first ASV
%   ASV   SESSION_SET 'all': every ASV session; 'first': the first-ASV session only
%   (SALpost / ASVrep are no longer produced; their comparisons are skipped)
% Each session also carries Ord (its position among the animal's infusions) and Rep (Ord > 1),
% read from compareInfusionOrder.mat when present, for the order-adjusted ASV vs SAL tests.
% TESTS (per context, on dEV; also printed for dR)
%   within group: dEV > 0, sign-flip over animal means
%   B - A within session (novel vs familiar; A is pre-infusion): sign-flip over animal means
%   between groups: permutation over sessions and LME y ~ Group + (1|Animal)
%   paired: ASV vs SAL and SALpost vs SAL on (B - A), animals with both

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT       = 'R:\TransPlasticity';
CELLS = 'PYR';
BIN = 0.1;                  % s, template bins
THETA_RATIO = 2;            % theta (6-10 Hz) / delta (1-4 Hz) power, 1 s windows
PAIR_EXCLUDE = 'none';          % 'none' (every pair) | 'sameChannel' (cell_metrics.maxWaveformCh1) | 'sameShank' | 'none'
SPEED_MIN = 3;              % cm/s, 1 s windows; from behavior.position.lin (cm) when present, else x/y (metres detected and x100)
MIN_SPK = 20;               % spikes in the exposure (theta part)
MIN_UNITS = 5;
MIN_RIP = 50; MAX_RIP = Inf;
MATCH_N = false;            % false: every ripple in the PRE / POST blocks; true: the same number in both
RUN_WIN = 0.2; RIP_WIN = 0.1;   % s: half-windows for the temporal-bias counts (behaviour / ripples)
MIN_CC_RUN = 20; MIN_CC_RIP = 5; % pair counts within the window needed for a bias
EXCL_INF = false;            % true: the marked Infusion_start..Infusion_end interval is not home cage
INF_EXCL_WIN = [0 1800];     % s [before after] Infusion_start: removed from every PRE/POST block ([0 0] = off).
                             % Fixed window because the marked infusion interval is highly variable. A PRE/POST
                             % block must touch its exposure: if the window reaches the exposure, that block is
                             % empty and the exposure is skipped (ledger), never replaced by a later block.
PRE_B_ALLGAP = true;         % ConB1 PRE = ALL home cage between leaving the previous context (ConA1) and
                             % entering ConB1, minus the infusion exclusion (may be two pieces: before the
                             % infusion and after the window); PRE_WIN does not apply to it. false: the
                             % home-cage block touching ConB1, last PRE_WIN.
PRE_WIN = 3600;              % s: ripples from this long before the exposure (PRE)
POST_WIN = 1800;             % s: ripples from this long after it (POST; 30 min, the drug-timing resolution)
TC_BIN = 3600; TC_MAX = 4*3600;  % s: time course of replay across the home cage after ConB1 (bins from its end)
TC_MIN_RIP = 30;                 % ripples a time-course bin needs; every bin of a session is subsampled to the same N
LATE_START = 2*3600;             % s after ConB1: bins starting here or later form the 'late' replay summary         % s: only ripples within this long before the exposure starts / after it ends
CONTEXTS = {'ConA1','ConB1','ConA2','ConB2'};
NPERM = 1e4;
% Cache mode: reactivationFromCache.m sets REACT_FROM_CACHE (struct) and calls this script; the
% per-pair / per-cell values saved by the last full run (compareReactivation_raw.mat) are re-thresholded
% and every table, test and figure below is rebuilt from them, without loading any session.
FROM_CACHE = exist('REACT_FROM_CACHE', 'var') && isstruct(REACT_FROM_CACHE);
CC = struct(); OUT_TAG = '';
if FROM_CACHE, CC = REACT_FROM_CACHE; end
clear REACT_FROM_CACHE                % so a later direct run recomputes

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
SESSION_SET = 'all';          % 'all' = every SAL and every ASV session | 'first' = SAL before the first ASV + the first ASV
if exist('SHARED', 'var')     % set by runAllAnalyses.m: one source for every shared setting
    CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root;
    if isfield(SHARED, 'sessionSet'), SESSION_SET = SHARED.sessionSet; end
    MIN_SPK = SHARED.minSpk; MIN_UNITS = SHARED.minUnits; PAIR_EXCLUDE = SHARED.pairExclude;
end
opt = struct('runWin', RUN_WIN, 'ripWin', RIP_WIN, 'minCCrun', MIN_CC_RUN, 'minCCrip', MIN_CC_RIP, 'preWin', PRE_WIN, 'postWin', POST_WIN, 'pairExclude', PAIR_EXCLUDE, 'cells', CELLS, 'bin', BIN, 'thetaRatio', THETA_RATIO, 'speedMin', SPEED_MIN, 'minSpk', MIN_SPK, ...
    'minUnits', MIN_UNITS, 'minRip', MIN_RIP, 'maxRip', MAX_RIP, 'matchN', MATCH_N);

if ~FROM_CACHE   % ======== full computation (cache mode: see the else branch) ========
rows = {};
%% ---- controls ------------------------------------------------------------------------
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), warning('react:file', 'missing %s', rf); continue; end
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES;
    nov = CONTROLS{i,2}; fam = CONTROLS{i,3}; rec = RP.meta.recordingSec;
    E = struct('name', {'ConA1','ConB1'}, 'on', {fam(1), nov(1)}, 'off', {fam(2), nov(2)}, 'preIv', {[], []});
    if fam(1) < nov(1)
        E(1).pre = [0 fam(1)];     E(1).post = [fam(2) nov(1)];
        E(2).pre = [fam(2) nov(1)]; E(2).post = [nov(2) rec];
        if PRE_B_ALLGAP, E(2).preIv = [fam(2) nov(1)]; end   % whole gap (no infusion in controls)
    else
        E(2).pre = [0 nov(1)];     E(2).post = [nov(2) fam(1)];
        E(1).pre = [nov(2) fam(1)]; E(1).post = [fam(2) rec];
    end
    for e = E
        y = local_react(RP, fld, e, opt);
        rows(end+1,:) = [{'Control', controlAnimal(rf), rf, e.name}, y, {e}]; %#ok<AGROW>
    end
end

%% ---- TransPlasticity ----------------------------------------------------------------------
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
noCtx = {}; unp = {}; led = cell(0, 4); tim = cell(0, 8);       % ledger: session, group, context, reason not used
for i = 1:numel(fa)
    fld = fileparts(fa{i});
    ev = local_events(fld);
    if ~isempty(ev) && ~isempty(ev.unparsed)
        [~, sn_] = fileparts(fld); unp = [unp; strcat(sn_, ': ', ev.unparsed(:))]; %#ok<AGROW>
    end
    [~, snL] = fileparts(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConA1_start'))
        why_ = 'no event file or annotations'; if ~isempty(ev), why_ = ['no ConA1 in ' ev.source]; end
        noCtx(end+1,:) = {fld, why_}; led(end+1,:) = {snL, '', 'all', why_}; continue %#ok<AGROW>
    end
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    g = local_group(RP, ev, SESSION_SET);
    if isempty(g), led(end+1,:) = {snL, '', 'all', sprintf('not a SAL or ASV session in SESSION_SET ''%s'' (or no infusion time)', SESSION_SET)}; continue; end %#ok<AGROW>
    HC = local_hcBlocks(ev, RP, EXCL_INF, INF_EXCL_WIN);   % home cage = outside ConA/ConB and the infusion exclusion
    for c = 1:numel(CONTEXTS)
        on = ev.time(strcmp(ev.description, [CONTEXTS{c} '_start']));
        off = ev.time(strcmp(ev.description, [CONTEXTS{c} '_end']));
        if isempty(on) || isempty(off)
            if c <= 2, led(end+1,:) = {snL, g, CONTEXTS{c}, 'no start/end marker'}; end %#ok<AGROW>
            continue
        end
        e = struct('name', CONTEXTS{c}, 'on', on(1), 'off', off(1), 'pre', [], 'post', [], 'preIv', []);
        kp = find(abs(HC(:,2) - e.on) <= 1, 1); kq = find(abs(HC(:,1) - e.off) <= 1, 1);   % blocks touching the exposure
        if PRE_B_ALLGAP && strcmp(CONTEXTS{c}, 'ConB1')
            e.preIv = local_gapIv(HC, ev, e.on); kp = [];
            if ~isempty(e.preIv), kp = 0; end
        end
        if isempty(kp) || isempty(kq)
            led(end+1,:) = {snL, g, CONTEXTS{c}, sprintf('no home cage before (%d) / after (%d) the exposure after the context and infusion exclusions', ~isempty(kp), ~isempty(kq))}; %#ok<AGROW>
            continue
        end
        if kp > 0, e.pre = HC(kp,:); else, e.pre = [e.preIv(1,1) e.preIv(end,2)]; end
        e.post = HC(kq,:);
        if c <= 2                                % where PRE/POST sit relative to the infusion (min)
            iOn = ev.time(find(strcmp(ev.description, 'Infusion_start'), 1)); iOff = ev.time(find(strcmp(ev.description, 'Infusion_end'), 1));
            if isempty(iOn), iOn = NaN; end
            if isempty(iOff), iOff = NaN; end
            pr = local_preIv(e, PRE_WIN); po = [e.off, min(e.post(2), e.off + POST_WIN)];
            tim(end+1,:) = {snL, g, CONTEXTS{c}, (iOn - e.off)/60, (e.on - iOn)/60, (iOff - iOn)/60, sum(diff(pr, 1, 2))/60, diff(po)/60}; %#ok<AGROW>
        end
        y = local_react(RP, fld, e, opt);
        rows(end+1,:) = [{g, m.animal, fa{i}, e.name}, y, {e}]; %#ok<AGROW>
    end
end
fprintf('event labels not recognised (or ends inferred), by session:\n'); fprintf('   %s\n', unp{:});
fprintf('\nINFUSION TIMING (min; INF_EXCL_WIN = [-%g +%g] min around Infusion_start, EXCL_INF=%d, PRE_WIN=%g, POST_WIN=%g)\n', ...
    INF_EXCL_WIN/60, EXCL_INF, PRE_WIN/60, POST_WIN/60);
fprintf('   %-34s | %-6s | %-5s | %10s | %10s | %8s | %8s | %8s\n', 'session', 'group', 'ctx', 'infS-ctxEnd', 'ctxOn-infS', 'inf dur', 'PRE', 'POST');
for k = 1:size(tim, 1), fprintf('   %-34s | %-6s | %-5s | %10.1f | %10.1f | %8.1f | %8.1f | %8.1f\n', tim{k,:}); end
fprintf('   (PRE / POST = minutes of home cage actually used, after the context and infusion exclusions)\n');
fprintf('%d TransPlasticity ripple exports without context markers:\n', size(noCtx, 1));
for k = 1:size(noCtx, 1), fprintf('   %-75s %s\n', noCtx{k,1}, noCtx{k,2}); end
else             % ======== cache mode: re-threshold the cached rows ========
    RAW = load(fullfile(CACHE_ROOT, 'compareReactivation_raw.mat'));
    rows = RAW.rows; tcRows = RAW.tcRows; lateRows = RAW.lateRows; gainS = RAW.gainS; opt = RAW.opt;
    PRE_WIN = opt.preWin; POST_WIN = opt.postWin; MIN_UNITS = opt.minUnits;
    TC_BIN = RAW.rawPar.TC_BIN; TC_MAX = RAW.rawPar.TC_MAX; LATE_START = RAW.rawPar.LATE_START;
    EXCL_INF = RAW.rawPar.EXCL_INF; CONTEXTS = RAW.rawPar.CONTEXTS;
    if isfield(RAW.rawPar, 'INF_EXCL_WIN'), INF_EXCL_WIN = RAW.rawPar.INF_EXCL_WIN; else, INF_EXCL_WIN = [NaN NaN]; end
    if ~isfield(CC, 'minSpk'), CC.minSpk = opt.minSpk; end
    if ~isfield(CC, 'minRipSpk'), CC.minRipSpk = 0; end
    if CC.minSpk < opt.minSpk
        error('react:cache', 'cache was built with MIN_SPK = %d; only higher thresholds can be applied', opt.minSpk);
    end
    fprintf('==== FROM CACHE (built %s): MIN_SPK %d -> %d, spikes in PRE and in POST ripples >= %d ====\n', ...
        RAW.rawPar.built, opt.minSpk, CC.minSpk, CC.minRipSpk);
    fprintf('     cached settings: PRE_WIN %g min, POST_WIN %g min, infusion exclusion [-%g +%g] min, EXCL_INF %d, MIN_UNITS %d, CELLS %s, PAIR_EXCLUDE %s\n', ...
        PRE_WIN/60, POST_WIN/60, INF_EXCL_WIN/60, EXCL_INF, MIN_UNITS, opt.cells, opt.pairExclude);
    if isfield(RAW.rawPar, 'PRE_B_ALLGAP'), fprintf('     ConB1 PRE = whole ConA1 -> ConB1 gap: %d\n', RAW.rawPar.PRE_B_ALLGAP); end
    opt.minSpk = CC.minSpk;
    OUT_TAG = sprintf('_spk%d', CC.minSpk); if CC.minRipSpk > 0, OUT_TAG = sprintf('%s_rip%d', OUT_TAG, CC.minRipSpk); end
    % rows: 1-4 Group Animal Session Context | 5 nUnits 6 nPairs 7 nRip 8 thetaFrac 9 thetaSrc 10 rPre 11 rPost
    %       12 rPP 13 EV 14 REV 15 nBias 16 bPre 17 bPost 18 biasPairs 19 pairVec 20 unitN 21 why | 22 E
    for r_ = 1:size(rows, 1)
        if ~isempty(rows{r_,21}), continue; end            % unusable at the original threshold too
        s_ = local_rethresh(rows{r_,19}, rows{r_,20}, rows{r_,18}, CC.minSpk, CC.minRipSpk, MIN_UNITS);
        rows(r_, [5 6 10:17 21]) = {s_.nUnits, s_.nPairs, s_.rPre, s_.rPost, s_.rPP, s_.EV, s_.REV, s_.nBias, s_.bPre, s_.bPost, s_.why};
        rows{r_,18} = {s_.BP}; rows{r_,19} = {s_.PV};   % same 1x1-cell wrapping as a computed row
    end
    % tcRows: 8 EV 9 REV 10 rPre 11 rPost 12 why 13 pairVec 14 unitN;  lateRows: 6 EV 7 REV 8 rPre 9 rPost 10 why 11 pairVec 12 unitN
    for r_ = 1:size(tcRows, 1)
        if ~isempty(tcRows{r_,12}), continue; end
        s_ = local_rethresh(tcRows{r_,13}, tcRows{r_,14}, [], CC.minSpk, CC.minRipSpk, MIN_UNITS);
        tcRows(r_, 8:12) = {s_.EV, s_.REV, s_.rPre, s_.rPost, s_.why}; tcRows{r_,13} = {s_.PV};
    end
    for r_ = 1:size(lateRows, 1)
        if ~isempty(lateRows{r_,10}), continue; end
        s_ = local_rethresh(lateRows{r_,11}, lateRows{r_,12}, [], CC.minSpk, CC.minRipSpk, MIN_UNITS);
        lateRows(r_, 6:10) = {s_.EV, s_.REV, s_.rPre, s_.rPost, s_.why}; lateRows{r_,11} = {s_.PV};
    end
    if isfield(CC, 'gainMinSpk'), GAIN_MIN_SPK_CC = CC.gainMinSpk; end
    led = cell(0, 4);
end              % ======== end of computation / cache branch ========
T = cell2table(rows, 'VariableNames', {'Group','Animal','Session','Context','nUnits','nPairs','nRip', ...
    'thetaFrac','thetaSrc','rPre','rPost','rPP','EV','REV','nBias','bPre','bPost','biasPairs','pairVec','unitN','why','E'});
T.dEV = T.EV - T.REV; T.dR = T.rPost - T.rPre; T.dB = T.bPost - T.bPre;
ok = isfinite(T.dEV);
fprintf('%d context exposures, %d usable\n', height(T), sum(ok));
% ---- ledger of every TransPlasticity export that did not give a usable ConA1 and ConB1 -----------
for r_ = find(~ok)'
    if strcmp(T.Group{r_}, 'Control'), continue; end
    [~, snL] = fileparts(fileparts(T.Session{r_}));
    if isempty(snL), [~, snL] = fileparts(T.Session{r_}); end
    led(end+1,:) = {snL, T.Group{r_}, T.Context{r_}, T.why{r_}}; %#ok<AGROW>
end
fprintf('\nSESSION LEDGER: exports not contributing to the ConA1/ConB1 reactivation (session | group | context | reason)\n');
if isempty(led), fprintf('   none\n'); else
    [~, io] = sortrows(led(:, [2 1 3])); led = led(io, :);
    for r_ = 1:size(led, 1), fprintf('   %-34s | %-6s | %-5s | %s\n', led{r_,:}); end
end
fprintf('usable ConB1: SAL %d, ASV %d; usable ConA1: SAL %d, ASV %d\n', ...
    sum(ok & strcmp(T.Group, 'SAL') & strcmp(T.Context, 'ConB1')), sum(ok & strcmp(T.Group, 'ASV') & strcmp(T.Context, 'ConB1')), ...
    sum(ok & strcmp(T.Group, 'SAL') & strcmp(T.Context, 'ConA1')), sum(ok & strcmp(T.Group, 'ASV') & strcmp(T.Context, 'ConA1')));
if any(~ok)
    w = T.why(~ok); [u, ~, j] = unique(w);
    for k = 1:numel(u), fprintf('   skipped %3d: %s\n', sum(j == k), u{k}); end
end
T = T(ok, :);
G = intersect({'Control','SAL','ASV','SALpost','ASVrep'}, unique(T.Group), 'stable'); G = G(:)';   % row: for-loops iterate columns
% ---- infusion order of each TransPlasticity session, from compareInfusionOrder.mat ----
T.Ord = nan(height(T), 1);
[~, snT] = cellfun(@(f) fileparts(fileparts(f)), T.Session, 'UniformOutput', false);
ordFile = fullfile(CACHE_ROOT, 'compareInfusionOrder.mat');
if exist(ordFile, 'file')
    O_ = load(ordFile, 'R');
    [tf_, loc_] = ismember(snT, O_.R.Session); T.Ord(tf_) = O_.R.Ord(loc_(tf_));
    kI_ = ismember(T.Group, {'SAL', 'ASV'});
    fprintf('infusion order joined for %d of %d infused rows', sum(tf_ & kI_), sum(kI_));
    if any(kI_ & ~tf_), fprintf('; no order for: %s', strjoin(unique(snT(kI_ & ~tf_))', ', ')); end
    fprintf('\n');
else
    fprintf('%s not found: run compareInfusionOrder first for the order-adjusted tests\n', ordFile);
end
T.Rep = double(T.Ord > 1); T.Rep(~isfinite(T.Ord)) = NaN;
fprintf('SESSION_SET ''%s'': %s\n', SESSION_SET, strjoin(cellfun(@(g_) sprintf('%s %d sessions / %d animals', g_, ...
    numel(unique(T.Session(strcmp(T.Group, g_)))), numel(unique(T.Animal(strcmp(T.Group, g_))))), G, 'UniformOutput', false), ', '));
fprintf('theta template source: %s\n', strjoin(cellfun(@(s) sprintf('%s %d', s, sum(strcmp(T.thetaSrc, s))), ...
    unique(T.thetaSrc), 'UniformOutput', false), ', '));

%% ---- per context: group means and within-group tests -------------------------------------------
for mm = {'dEV','dR','dB'}
    nm = mm{1};
    fprintf('\n%s  (mean +- SEM over animal means; [animals]; p = sign-flip over animal means, H0: 0)\n', nm);
    fprintf('%-8s', 'context'); fprintf(' | %-30s', G{:}); fprintf('\n');
    for c = {'ConA1','ConB1','ConA2','ConB2'}
        fprintf('%-8s', c{1});
        for g = G
            k = strcmp(T.Group, g{1}) & strcmp(T.Context, c{1}) & isfinite(T.(nm));
            if ~any(k), fprintf(' | %-30s', '-'); continue; end
            [am, ~] = local_animalMean(T(k,:), nm);
            fprintf(' | %-30s', sprintf('%+.3f +- %.3f [%d] p=%.3f', mean(am), local_sem(am), ...
                numel(am), local_signflip(am)));
        end
        fprintf('\n');
    end
end
fprintf('\nEV / REV means (ConB1):\n');
for g = G
    k = strcmp(T.Group, g{1}) & strcmp(T.Context, 'ConB1');
    if any(k), fprintf('   %-8s EV %.3f  REV %.3f  (n=%d)\n', g{1}, mean(T.EV(k)), mean(T.REV(k)), sum(k)); end
end

%% ---- novel minus familiar within session (A1 is before the infusion) -------------------------
B = local_pairBA(T, 'ConB1', 'ConA1');
fprintf('\nB1 - A1 within session, dEV (A1 = pre-infusion baseline in TransPlasticity)\n');
for g = G
    k = strcmp(B.Group, g{1});
    if ~any(k), continue; end
    am = local_animalMean(B(k,:), 'd');
    fprintf('   %-8s %+.3f +- %.3f (mean +- SEM over animals; %d sessions, %d animals)  sign-flip p = %.3f\n', g{1}, mean(am), ...
        local_sem(am), sum(k), numel(am), local_signflip(am));
end
prs = {'SAL','ASV'; 'SAL','SALpost'; 'Control','SAL'; 'ASV','ASVrep'};
fprintf('\nbetween groups on B1 - A1 (dEV): perm over sessions | LME (1|Animal) | paired sign-flip over animals with both\n');
for j = 1:size(prs, 1)
    if ~all(ismember(prs(j,:), B.Group)), continue; end
    a = B.d(strcmp(B.Group, prs{j,1})); b = B.d(strcmp(B.Group, prs{j,2}));
    fprintf('   %-8s vs %-8s  %+.3f vs %+.3f | perm %.3f | LME %.3f | paired %s\n', prs{j,:}, mean(a), mean(b), ...
        local_perm(a, b, NPERM), local_lme(B, 'd', prs{j,1}, prs{j,2}), local_pairedStr(B, 'd', prs{j,1}, prs{j,2}));
end
fprintf('\nbetween groups on dEV per context (perm | LME | paired)\n');
for c = {'ConA1','ConB1'}
    t = T(strcmp(T.Context, c{1}), :);
    for j = 1:size(prs, 1)
        if ~all(ismember(prs(j,:), t.Group)), continue; end
        a = t.dEV(strcmp(t.Group, prs{j,1})); b = t.dEV(strcmp(t.Group, prs{j,2}));
        fprintf('   %-6s %-8s vs %-8s  %+.3f vs %+.3f | perm %.3f | LME %.3f | paired %s\n', c{1}, prs{j,:}, mean(a), mean(b), ...
            local_perm(a, b, NPERM), local_lme(t, 'dEV', prs{j,1}, prs{j,2}), local_pairedStr(t, 'dEV', prs{j,1}, prs{j,2}));
    end
end
fprintf('\nASV vs SAL with infusion order held fixed: y ~ Group + Ord|Rep + (1|Animal), estimate +- SE (Satterthwaite)\n');
B.Ord = nan(height(B), 1); [tfB, locB] = ismember(B.Session, T.Session); B.Ord(tfB) = T.Ord(locB(tfB));
B.Rep = double(B.Ord > 1); B.Rep(~isfinite(B.Ord)) = NaN;
for cv = {'Ord', 'Rep'}
    fprintf('   B1 - A1 dEV    [%s] %s\n', cv{1}, local_lmeOrder(B, 'd', cv{1}));
    for c = {'ConA1', 'ConB1'}
        fprintf('   %-6s dEV     [%s] %s\n', c{1}, cv{1}, local_lmeOrder(T(strcmp(T.Context, c{1}), :), 'dEV', cv{1}));
    end
end
%% ---- temporal bias (Skaggs & McNaughton 1996) -------------------------------------------------
fprintf('\ntemporal bias, pairs pooled over a group''s sessions: r(behaviour bias, PRE ripple bias) / r(.., POST)\n');
for c = {'ConA1','ConB1','ConA2','ConB2'}
    for g = G
        k = strcmp(T.Group, g{1}) & strcmp(T.Context, c{1});
        if ~any(k), continue; end
        BPg = vertcat(T.biasPairs{k});
        if size(BPg, 1) < 10, continue; end
        [rA, pA] = corr(BPg(:,1), BPg(:,2)); [rB, pB] = corr(BPg(:,1), BPg(:,3));
        fprintf('   %-6s %-8s %4d pairs  PRE r=%+.3f (p=%.3f)  POST r=%+.3f (p=%.3f)\n', c{1}, g{1}, size(BPg, 1), rA, pA, rB, pB);
    end
end
fprintf('(pooled p treats pairs as independent; the per-session dB rows above are the replicate-level test)\n');
okb = isfinite(T.dB);
fprintf('\nbetween groups on dB per context (perm | LME | paired)\n');
for c = {'ConA1','ConB1'}
    t = T(strcmp(T.Context, c{1}) & okb, :);
    for j = 1:size(prs, 1)
        if ~all(ismember(prs(j,:), t.Group)), continue; end
        a = t.dB(strcmp(t.Group, prs{j,1})); b = t.dB(strcmp(t.Group, prs{j,2}));
        fprintf('   %-6s %-8s vs %-8s  %+.3f vs %+.3f | perm %.3f | LME %.3f | paired %s\n', c{1}, prs{j,:}, mean(a), mean(b), ...
            local_perm(a, b, NPERM), local_lme(t, 'dB', prs{j,1}, prs{j,2}), local_pairedStr(t, 'dB', prs{j,1}, prs{j,2}));
    end
end

%% ---- replay time course across the home cage after ConB1 ----------------------------------------
% Templates ConA1 (familiar, encoded before the infusion) and ConB1 (novel, encoded on drug).
% PRE = the last hour before each exposure (as above); POST = successive TC_BIN bins of the home
% cage after ConB1, up to TC_MAX. Each bin: the same number of ripples in PRE and in the bin.
fprintf('\nreplay time course in the home cage after ConB1 (%g min bins)\n', TC_BIN/60);
if ~FROM_CACHE   % ======== full computation of time course, late block, gain cells ========
tcRows = {}; lateRows = {}; gainS = {}; optTC = opt; optTC.minRip = TC_MIN_RIP;
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1}; rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), continue; end
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES; nov = CONTROLS{i,2}; fam = CONTROLS{i,3}; rec = RP.meta.recordingSec;
    if fam(1) > nov(1), continue; end                                 % layout assumed familiar first
    EA = struct('name','ConA1','on',fam(1),'off',fam(2),'pre',[0 fam(1)],'post',[],'postFixed',true,'preIv',[]);
    EB = struct('name','ConB1','on',nov(1),'off',nov(2),'pre',[fam(2) nov(1)],'post',[],'postFixed',true,'preIv',[]);
    if PRE_B_ALLGAP, EB.preIv = [fam(2) nov(1)]; end
    tcRows = [tcRows; local_tc(RP, fld, EA, EB, [nov(2) rec], TC_BIN, TC_MAX, optTC, 'Control', controlAnimal(rf), rf)]; %#ok<AGROW>
    lateRows = [lateRows; local_block(RP, fld, EA, EB, [nov(2) rec], LATE_START, TC_MAX, optTC, 'Control', controlAnimal(rf), rf)]; %#ok<AGROW>
    U_ = local_gainSession(RP, fld, EA, EB, [nov(2) rec], opt, LATE_START, TC_MAX);
    if ~isempty(U_), gainS(end+1,:) = {'Control', controlAnimal(rf), rf, U_}; end %#ok<AGROW>
end
for i = 1:numel(fa)
    fld = fileparts(fa{i}); ev = local_events(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConB1_start')) || ~any(contains(ev.description, 'ConA1_start')), continue; end
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; g = local_group(RP, ev, SESSION_SET);
    if isempty(g), continue; end
    HC = local_hcBlocks(ev, RP, EXCL_INF, INF_EXCL_WIN);   % home cage = outside ConA/ConB and the infusion exclusion
    t = @(lab) ev.time(find(strcmp(ev.description, lab), 1));
    aOn = t('ConA1_start'); aOff = t('ConA1_end'); bOn = t('ConB1_start'); bOff = t('ConB1_end');
    kA = find(abs(HC(:,2) - aOn) <= 1, 1); kB = find(abs(HC(:,2) - bOn) <= 1, 1); kP = find(abs(HC(:,1) - bOff) <= 1, 1);
    pBiv = [];
    if PRE_B_ALLGAP, pBiv = local_gapIv(HC, ev, bOn); kB = []; if ~isempty(pBiv), kB = 0; end; end
    if isempty(kA) || isempty(kB) || isempty(kP), continue; end
    if kB > 0, preB_ = HC(kB,:); else, preB_ = [pBiv(1,1) pBiv(end,2)]; end
    EA = struct('name','ConA1','on',aOn,'off',aOff,'pre',HC(kA,:),'post',[],'postFixed',true,'preIv',[]);
    EB = struct('name','ConB1','on',bOn,'off',bOff,'pre',preB_,'post',[],'postFixed',true,'preIv',pBiv);
    tcRows = [tcRows; local_tc(RP, fld, EA, EB, HC(kP,:), TC_BIN, TC_MAX, optTC, g, RP.meta.animal, fa{i})]; %#ok<AGROW>
    lateRows = [lateRows; local_block(RP, fld, EA, EB, HC(kP,:), LATE_START, TC_MAX, optTC, g, RP.meta.animal, fa{i})]; %#ok<AGROW>
    U_ = local_gainSession(RP, fld, EA, EB, HC(kP,:), opt, LATE_START, TC_MAX);
    if ~isempty(U_), gainS(end+1,:) = {g, RP.meta.animal, fa{i}, U_}; end %#ok<AGROW>
end
% everything reactivationFromCache.m needs to re-threshold without loading sessions
rawPar = struct('TC_BIN', TC_BIN, 'TC_MAX', TC_MAX, 'LATE_START', LATE_START, 'EXCL_INF', EXCL_INF, 'INF_EXCL_WIN', INF_EXCL_WIN, 'PRE_B_ALLGAP', PRE_B_ALLGAP, ...
    'CONTEXTS', {CONTEXTS}, 'built', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
save(fullfile(CACHE_ROOT, 'compareReactivation_raw.mat'), 'rows', 'tcRows', 'lateRows', 'gainS', 'opt', 'rawPar', '-v7.3');
end              % ======== end of computation ========
if isempty(tcRows), error('react:tc', 'no usable ConA1/ConB1 sessions for the time course'); end
TC = cell2table(tcRows, 'VariableNames', {'Group','Animal','Session','Template','Bin','tStart','nRip','EV','REV','rPre','rPost','why','pairVec','unitN'});
TC.dEV = TC.EV - TC.REV; TC.dR = TC.rPost - TC.rPre;
TC = TC(isfinite(TC.dEV), :);
nb = floor(TC_MAX / TC_BIN);
for tp = {'ConA1','ConB1'}
    fprintf('\n template %s: dEV by bin (h after ConB1 ends); mean +- SEM [sessions/animals], p = sign-flip over animals\n', tp{1});
    fprintf('%-9s', 'bin (h)'); fprintf(' | %-30s', G{:}); fprintf('\n');
    for b = 1:nb
        fprintf('%-9s', sprintf('%.1f-%.1f', (b-1)*TC_BIN/3600, b*TC_BIN/3600));
        for g = G
            k = strcmp(TC.Template, tp{1}) & strcmp(TC.Group, g{1}) & TC.Bin == b;
            if ~any(k), fprintf(' | %-30s', '-'); continue; end
            am = local_animalMean(TC(k,:), 'dEV');
            fprintf(' | %-30s', sprintf('%+.3f +- %.3f [%d/%d] p=%.3f', mean(TC.dEV(k)), std(TC.dEV(k))/sqrt(sum(k)), ...
                sum(k), numel(am), local_signflip(am)));
        end
        fprintf('\n');
    end
end
% late replay: mean dEV over bins starting >= LATE_START, per session
LT = varfun(@mean, TC(TC.tStart >= LATE_START, :), 'InputVariables', {'dEV','dR'}, ...
    'GroupingVariables', {'Group','Animal','Session','Template'});
LT.Properties.VariableNames(strcmp(LT.Properties.VariableNames, 'mean_dEV')) = {'dEV'};
LT.Properties.VariableNames(strcmp(LT.Properties.VariableNames, 'mean_dR')) = {'dR'};
fprintf('\nlate replay (bins from %.1f h after ConB1): mean dEV per session\n', LATE_START/3600);
for tp = {'ConA1','ConB1'}
    for g = G
        k = strcmp(LT.Template, tp{1}) & strcmp(LT.Group, g{1});
        if ~any(k), continue; end
        am = local_animalMean(LT(k,:), 'dEV');
        fprintf('   %-6s %-8s %+.3f +- %.3f (%d sessions, %d animals)  animals p=%.3f  sessions p=%.3f\n', tp{1}, g{1}, ...
            mean(LT.dEV(k)), std(LT.dEV(k))/sqrt(sum(k)), sum(k), numel(am), local_signflip(am), local_signflip(LT.dEV(k)));
    end
    t = LT(strcmp(LT.Template, tp{1}), :);
    if all(ismember({'SAL','ASV'}, t.Group))
        a = t.dEV(strcmp(t.Group, 'SAL')); b = t.dEV(strcmp(t.Group, 'ASV'));
        fprintf('   %-6s SAL vs ASV: perm %.3f | LME %.3f | paired %s\n', tp{1}, local_perm(a, b, NPERM), ...
            local_lme(t, 'dEV', 'SAL', 'ASV'), local_pairedStr(t, 'dEV', 'SAL', 'ASV'));
    end
end
fig2 = figure('Color','w','Position',[60 480 1100 380]);
tiledlayout(1, 2, 'TileSpacing','compact');
colTC = containers.Map({'Control','SAL','ASV','SALpost','ASVrep'}, ...
    {[0.4 0.4 0.4], [0.2 0.6 1], [0.85 0.2 0.1], [0 0.25 0.6], [0.55 0.1 0.6]});
for tp = {'ConA1','ConB1'}
    nexttile; hold on
    for g = G
        mu = nan(nb, 1); se = mu;
        for b = 1:nb
            k = strcmp(TC.Template, tp{1}) & strcmp(TC.Group, g{1}) & TC.Bin == b;
            if any(k), mu(b) = mean(TC.dEV(k)); se(b) = std(TC.dEV(k)) / sqrt(sum(k)); end
        end
        x = ((1:nb)' - 0.5) * TC_BIN / 3600;
        errorbar(x, mu, se, '-o', 'Color', colTC(g{1}), 'MarkerFaceColor', colTC(g{1}), 'LineWidth', 1.5, 'CapSize', 0);
    end
    yline(0, 'k:'); xline(LATE_START/3600, 'k--');
    xlabel('h after the novel exposure ends'); ylabel('EV - REV'); title(sprintf('replay of %s template', tp{1}));
    legend(G, 'Location','best', 'Box','off');
end
exportgraphics(fig2, fullfile(CACHE_ROOT, ['compareReactivation_timecourse' OUT_TAG '.pdf']), 'ContentType','vector');

%% ---- diagnostics: where does the spread in one group x context come from? ------------------------
% For each exposure of DIAG_GROUP x DIAG_CONTEXT (and DIAG_REF for comparison):
%   base        the summary value (rPre, rPost, rPP, pairs, units, N)
%   CI          95% bootstrap over pairs (2000 resamples) of dEV
%   LOO         range of dEV when one unit (all its pairs) is left out
%   variants    the same exposure recomputed with one setting changed at a time:
%               N30 / N30even = matched 30 ripples (first 30 / spread over POST_WIN),
%               pre60 / post60 = 60 min windows, spk20 = MIN_SPK 20
DIAG_GROUP = 'SAL'; DIAG_REF = 'ASV'; DIAG_CONTEXT = 'ConB1';
if ~FROM_CACHE   % reloads sessions; skipped in cache mode
fprintf('\n==== diagnostics: %s and %s, %s ====\n', DIAG_GROUP, DIAG_REF, DIAG_CONTEXT);
fprintf('%-8s %-6s %-30s %6s %6s %6s %5s %5s %4s | %7s %17s %17s | %7s %7s %7s %7s %7s\n', 'group', 'animal', 'session', ...
    'rPre', 'rPost', 'rPP', 'pairs', 'units', 'N', 'dEV', 'CI95', 'LOO range', 'N30', 'N30even', 'pre60', 'post60', 'spk20');
vars = {'N30','N30even','pre60','post60','spk20'};
for gg = {DIAG_GROUP, DIAG_REF}
    rr = find(strcmp(T.Group, gg{1}) & strcmp(T.Context, DIAG_CONTEXT));
    for r = rr(:)'
        PV = T.pairVec{r};
        ev = @(V) local_evFromPairs(V);
        d0 = ev(PV);
        bs = arrayfun(@(k) ev(PV(randi(size(PV, 1), size(PV, 1), 1), :)), 1:2000);
        uu = unique(PV(:, 4:5)); lo = arrayfun(@(u) ev(PV(PV(:,4) ~= u & PV(:,5) ~= u, :)), uu);
        L = load(T.Session{r}, 'RIPPLES'); RPd = L.RIPPLES; fldd = fileparts(T.Session{r});
        dv = nan(1, numel(vars));
        for q = 1:numel(vars)
            o2 = opt;
            if iscell(T.E), e2 = T.E{r}; else, e2 = T.E(r); end   % cell2table stores same-field structs as a struct array
            switch vars{q}
                case 'N30', o2.maxRip = 30; o2.minRip = 30; o2.matchN = true;
                case 'N30even'
                    o2.maxRip = 30; o2.minRip = 30; o2.matchN = true;
                    e2.post = [max(e2.post(1), e2.off), min(e2.post(2), e2.off + opt.postWin)]; e2.postFixed = true;
                case 'pre60', o2.preWin = 3600;
                case 'post60', o2.postWin = 3600;
                case 'spk20', o2.minSpk = 20;
            end
            yy = local_react(RPd, fldd, e2, o2);
            dv(q) = yy{9} - yy{10};
        end
        [~, sn] = fileparts(fldd);
        fprintf('%-8s %-6s %-30s %+6.2f %+6.2f %+6.2f %5d %5d %4d | %+7.3f [%+6.3f %+6.3f] [%+6.3f %+6.3f] |', gg{1}, T.Animal{r}, sn, ...
            T.rPre(r), T.rPost(r), T.rPP(r), T.nPairs(r), T.nUnits(r), T.nRip(r), d0, prctile(bs, 2.5), prctile(bs, 97.5), ...
            min(lo), max(lo));
        fprintf(' %+7.3f', dv); fprintf('\n');
    end
end
fprintf(['reading: a wide CI = few / noisy pairs; a wide LOO range = one unit drives it; rPP near 1 = the\n' ...
         'EV/REV denominators amplify small differences; a sign change across variants = the value depends on\n' ...
         'the window / N choice, not on the session.\n']);
end

%% ---- report: early (first POST_WIN after each exposure) and late (LATE_START..TC_MAX after ConB1) ----
% Groups Control, SAL, ASV. Early = the summary analysis (T). Late = one block of every ripple
% from LATE_START to TC_MAX after ConB1 ends, PRE as for the summary. Time course = TC (TC_BIN bins).
REPORT_G = {'Control','SAL','ASV'};
if isempty(lateRows), lateRows = cell(0, 12); end
LB = cell2table(lateRows, 'VariableNames', {'Group','Animal','Session','Template','nRip','EV','REV','rPre','rPost','why','pairVec','unitN'});
LB.dEV = LB.EV - LB.REV; LB = LB(isfinite(LB.dEV), :);
EA_ = T(:, {'Group','Animal','Session','Context','dEV'}); EA_.Properties.VariableNames{'Context'} = 'Template';
blocks = {'early', EA_; 'late', LB(:, {'Group','Animal','Session','Template','dEV'})};
fprintf('\n==== REPORT: EV - REV, groups %s ====\n', strjoin(REPORT_G, ', '));
fprintf('early = first %g min after each exposure; late = %g-%g h after ConB1 (all ripples); PRE = %g min before each exposure\n', ...
    POST_WIN/60, LATE_START/3600, TC_MAX/3600, PRE_WIN/60);
RPT = struct();
for tp = {'ConB1','ConA1'}
    for bk = 1:2
        D = blocks{bk,2}; D = D(strcmp(D.Template, tp{1}) & ismember(D.Group, REPORT_G), :);
        fprintf('\n%s template, %s\n', tp{1}, blocks{bk,1});
        for g = REPORT_G
            x = D.dEV(strcmp(D.Group, g{1}));
            if isempty(x), fprintf('   %-8s none\n', g{1}); continue; end
            am = local_animalMean(D(strcmp(D.Group, g{1}), :), 'dEV');
            fprintf('   %-8s %+.3f +- %.3f  (%d sessions, %d animals)  > 0: animals p=%.3f, sessions p=%.3f\n', g{1}, mean(x), ...
                std(x)/sqrt(numel(x)), numel(x), numel(am), local_signflip(am), local_signflip(x));
        end
        cp = {'SAL','ASV'; 'Control','SAL'; 'Control','ASV'};
        for j = 1:size(cp, 1)
            a = D.dEV(strcmp(D.Group, cp{j,1})); b = D.dEV(strcmp(D.Group, cp{j,2}));
            if numel(a) < 2 || numel(b) < 2, continue; end
            bd = arrayfun(@(k) mean(b(randi(numel(b), numel(b), 1))) - mean(a(randi(numel(a), numel(a), 1))), 1:5000);
            fprintf('   %-8s vs %-8s  diff %+.3f [95%% CI %+.3f %+.3f] | perm %.3f | LME %.3f | paired %s\n', cp{j,:}, ...
                mean(b) - mean(a), prctile(bd, 2.5), prctile(bd, 97.5), local_perm(a, b, NPERM), ...
                local_lme(D, 'dEV', cp{j,1}, cp{j,2}), local_pairedStr(D, 'dEV', cp{j,1}, cp{j,2}));
        end
        RPT.(tp{1}).(blocks{bk,1}) = D;
    end
end
fprintf('(CI: bootstrap over sessions within each group; difference = second minus first)\n');

colR = containers.Map({'Control','SAL','ASV'}, {[0.4 0.4 0.4], [0.2 0.6 1], [0.85 0.2 0.1]});
figR = figure('Color','w','Position',[60 60 1300 700]);
tlR = tiledlayout(2, 3, 'TileSpacing','compact');
nbR = floor(TC_MAX / TC_BIN);
for r = 1:2
    tp = {'ConB1','ConA1'}; tp = tp{r};
    nexttile; hold on
    for g = REPORT_G
        mu = nan(nbR, 1); se = mu;
        for b = 1:nbR
            k = strcmp(TC.Template, tp) & strcmp(TC.Group, g{1}) & TC.Bin == b;
            if any(k), mu(b) = mean(TC.dEV(k)); se(b) = std(TC.dEV(k)) / sqrt(sum(k)); end
        end
        errorbar(((1:nbR)' - 0.5) * TC_BIN / 3600, mu, se, '-o', 'Color', colR(g{1}), 'MarkerFaceColor', colR(g{1}), ...
            'LineWidth', 1.5, 'CapSize', 0);
    end
    yline(0, 'k:'); xline(LATE_START/3600, 'k--'); xlim([0 TC_MAX/3600]);
    xlabel('h after the novel exposure (ConB1) ends'); ylabel('EV - REV');
    title(sprintf('%s template: time course (%g h bins)', tp, TC_BIN/3600), 'FontWeight','normal');
    if r == 1, legend(REPORT_G, 'Location','best', 'Box','off'); end
    ttl = {sprintf('%s: first %g min after %s', tp, POST_WIN/60, tp), sprintf('%s: %g-%g h after ConB1', tp, LATE_START/3600, TC_MAX/3600)};
    for bk = 1:2
        nexttile; hold on
        D = RPT.(tp).(blocks{bk,1});
        for j = 1:numel(REPORT_G)
            x = D.dEV(strcmp(D.Group, REPORT_G{j}));
            plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', colR(REPORT_G{j}), 'MarkerFaceColor', colR(REPORT_G{j}), 'MarkerSize', 5);
            if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', colR(REPORT_G{j}), 'LineWidth', 2); end
        end
        yline(0, 'k:'); xlim([0.5 numel(REPORT_G)+0.5]); xticks(1:numel(REPORT_G)); xticklabels(REPORT_G);
        ylabel('EV - REV'); title(ttl{bk}, 'FontWeight','normal');
    end
end
title(tlR, 'Reactivation of familiar (ConA1) and novel (ConB1) coactivity in ripples');
exportgraphics(figR, fullfile(CACHE_ROOT, ['reactivationReport' OUT_TAG '.pdf']), 'ContentType','vector');

%% ---- gain component of reactivation -------------------------------------------------------------
% EV is built from Pearson correlations, so it is blind to changes in how much each cell fires.
% Here: does a cell's change in ripple firing after ConB1 scale with how active it was in ConB1
% (experience-tied gain), separately from a change shared by all cells (global gain)?
% Per session, over pyramidal cells (GAIN_CELLS), for y = change from BASE to POST block of
%   lrip  log10(in-ripple rate + 0.1)       lout  log10(out-of-ripple rate + 0.01)
%   xp    ripple participation beyond chance at the cell's own out-of-ripple rate in that block
% fit y = a + bB*zB + bA*zA, z = within-session z-score of log10(rate in the exposure + 0.01),
% running/theta part (same mask as the EV template). bB = gain tied to novel-context activity,
% controlling for familiar-context activity (bA); a = change common to all cells.
% BASE: preB (30 min before ConB1, drug on board as in POST) or preA (before ConA1, drug-free).
% POST: early (first POST_WIN after ConB1) or late (LATE_START..TC_MAX after ConB1).
GAIN_CELLS = 'PYR';                 % 'PYR' | 'PYRdeep' | 'PYRsup'
GAIN_MIN_SPK = 10;                  % spikes a cell needs in each of BASE and POST
if exist('GAIN_MIN_SPK_CC', 'var'), GAIN_MIN_SPK = GAIN_MIN_SPK_CC; clear GAIN_MIN_SPK_CC; end
GAIN_MIN_CELLS = 8;
GR = {'Control','SAL','ASV'};   % rows printed for all three; the TESTS are SAL vs ASV only
                                    % (the never-infused arm differs in too much to carry this contrast;
                                    %  it is kept as a check that the decomposition behaves)
ys = {'lrip','lout','xp'}; posts = {'early','late'};
% BASE -> POST pairs. The last is a PLACEBO: the change from preA to preB happens before ConB1, so
% ConB1 activity cannot cause it. Slow excitability drift predicts placebo bB > 0 and bA < 0 (cells
% rising over the session are low in ConA1, high in ConB1); an experience-driven gain predicts bB ~ 0.
% In TransPlasticity the placebo interval contains the infusion (descriptive there).
CMB = {'preB','early'; 'preB','late'; 'preA','early'; 'preA','late'; 'preA','preB'};
gRows = {};
for i = 1:size(gainS, 1)
    U = gainS{i,4};
    sel = U.(GAIN_CELLS);
    for cb = 1:size(CMB, 1)
            nb = [CMB{cb,1} '_']; np = [CMB{cb,2} '_'];
            ok = sel & U.([nb 'n']) >= GAIN_MIN_SPK & U.([np 'n']) >= GAIN_MIN_SPK & isfinite(U.rA) & isfinite(U.rB);
            if sum(ok) < GAIN_MIN_CELLS, continue; end
            zB = zscore(log10(U.rB(ok) + 0.01)); zA = zscore(log10(U.rA(ok) + 0.01));
            for yy = 1:numel(ys)
                y = U.([np ys{yy}])(ok) - U.([nb ys{yy}])(ok);
                kk = isfinite(y);
                if sum(kk) < GAIN_MIN_CELLS, continue; end
                Xd = [ones(sum(kk), 1), zB(kk), zA(kk)];
                cf = Xd \ y(kk);
                gRows(end+1,:) = {gainS{i,1}, gainS{i,2}, gainS{i,3}, CMB{cb,1}, CMB{cb,2}, ys{yy}, sum(kk), cf(1), cf(2), cf(3)}; %#ok<AGROW>
            end
    end
end
if isempty(gRows), gRows = cell(0, 10); end
GN = cell2table(gRows, 'VariableNames', {'Group','Animal','Session','Base','Post','Y','nCells','a','bB','bA'});
fprintf('\n==== gain component (%s; y = POST - BASE per cell; y ~ a + bB*zB + bA*zA) ====\n', GAIN_CELLS);
fprintf('a = shared change; bB = change per SD of ConB1 activity (controlling ConA1); bA = per SD of ConA1 activity\n');
for cb = 1:size(CMB, 1)
        for yy = 1:numel(ys)
            t = GN(strcmp(GN.Base, CMB{cb,1}) & strcmp(GN.Post, CMB{cb,2}) & strcmp(GN.Y, ys{yy}) & ismember(GN.Group, GR), :);
            if isempty(t), continue; end
            plc = ''; if strcmp(CMB{cb,2}, 'preB'), plc = '   [PLACEBO: before ConB1]'; end
            fprintf('\n%s -> %s, %s%s\n', CMB{cb,1}, CMB{cb,2}, ys{yy}, plc);
            for g = GR
                k = strcmp(t.Group, g{1});
                if ~any(k), continue; end
                am = local_animalMean(t(k,:), 'bB');
                fprintf('   %-8s a %+.3f +- %.3f | bB %+.3f +- %.3f (>0: animals p=%.3f, sessions p=%.3f) | bA %+.3f +- %.3f | %d sessions, median %d cells\n', ...
                    g{1}, mean(t.a(k)), std(t.a(k))/sqrt(sum(k)), mean(t.bB(k)), std(t.bB(k))/sqrt(sum(k)), ...
                    local_signflip(am), local_signflip(t.bB(k)), mean(t.bA(k)), std(t.bA(k))/sqrt(sum(k)), sum(k), round(median(t.nCells(k))));
            end
            for cp = {'SAL','ASV'}'
                a1 = t.bB(strcmp(t.Group, cp{1})); b1 = t.bB(strcmp(t.Group, cp{2}));
                a2 = t.a(strcmp(t.Group, cp{1}));  b2 = t.a(strcmp(t.Group, cp{2}));
                if numel(a1) < 2 || numel(b1) < 2, continue; end
                fprintf('   %s vs %s: bB %+.3f vs %+.3f perm %.3f LME %.3f paired %s | a %+.3f vs %+.3f perm %.3f\n', cp{1}, cp{2}, ...
                    mean(a1), mean(b1), local_perm(a1, b1, NPERM), local_lme(t, 'bB', cp{1}, cp{2}), local_pairedStr(t, 'bB', cp{1}, cp{2}), ...
                    mean(a2), mean(b2), local_perm(a2, b2, NPERM));
            end
        end
end
% drift check: experience-tied bB (preA -> early) vs placebo bB (preA -> preB), same session and cells' z
fprintf('\ndrift check, lrip: bB real (preA -> early) vs placebo (preA -> preB); bA alongside; p = sign-flip over sessions of real - placebo\n');
for g = GR
    tr = GN(strcmp(GN.Group, g{1}) & strcmp(GN.Base, 'preA') & strcmp(GN.Post, 'early') & strcmp(GN.Y, 'lrip'), :);
    tp = GN(strcmp(GN.Group, g{1}) & strcmp(GN.Base, 'preA') & strcmp(GN.Post, 'preB') & strcmp(GN.Y, 'lrip'), :);
    tr = tr(:, {'Session','bB','bA'}); tr.Properties.VariableNames = {'Session','bB_left','bA_left'};
    tp = tp(:, {'Session','bB','bA'}); tp.Properties.VariableNames = {'Session','bB_right','bA_right'};
    J = innerjoin(tr, tp, 'Keys', 'Session');
    if height(J) < 2, continue; end
    fprintf('   %-8s real bB %+.3f bA %+.3f | placebo bB %+.3f bA %+.3f | real - placebo bB %+.3f p=%.3f (%d sessions)\n', g{1}, ...
        mean(J.bB_left), mean(J.bA_left), mean(J.bB_right), mean(J.bA_right), mean(J.bB_left - J.bB_right), ...
        local_signflip(J.bB_left - J.bB_right), height(J));
end
% does experience-tied gain go with pattern reactivation in the same session?
fprintf('\ncoupling of gain (bB, lrip, preB -> early / late) with pattern reactivation (dEV ConB1, early / late), over sessions\n');
for pp = 1:2
    t = GN(strcmp(GN.Base, 'preB') & strcmp(GN.Post, posts{pp}) & strcmp(GN.Y, 'lrip'), :);
    if pp == 1, R_ = T(strcmp(T.Context, 'ConB1'), {'Session','dEV'}); else, R_ = LB(strcmp(LB.Template, 'ConB1'), {'Session','dEV'}); end
    J = innerjoin(t, R_, 'Keys', 'Session');
    for gset = {{'SAL'}, {'ASV'}, {'Control'}}
        k = ismember(J.Group, gset{1}) & isfinite(J.dEV) & isfinite(J.bB);
        if sum(k) < 4, continue; end
        [rr, pr] = corr(J.bB(k), J.dEV(k));
        fprintf('   %-5s %-16s r = %+.2f (p = %.3f), n = %d sessions\n', posts{pp}, strjoin(gset{1}, '+'), rr, pr, sum(k));
    end
end

% figure: binned change in ripple firing vs ConB1 activity (pooled cells, sessions weighted equally)
figG = figure('Color','w','Position',[60 60 1300 380]);
tlG = tiledlayout(1, 3, 'TileSpacing','compact');
colG = containers.Map({'Control','SAL','ASV'}, {[0.4 0.4 0.4], [0.2 0.6 1], [0.85 0.2 0.1]});
yl = {'lrip','lout','xp'}; ylab = {'\Delta log_{10} in-ripple rate','\Delta log_{10} out-of-ripple rate','\Delta participation beyond chance'};
for yy = 1:3
    nexttile; hold on
    for g = GR(~strcmp(GR, 'Control'))
        k = find(strcmp(gainS(:,1), g{1}));
        X = []; Y = []; Wt = [];
        for q = k(:)'
            U = gainS{q,4};
            ok = U.(GAIN_CELLS) & U.preB_n >= GAIN_MIN_SPK & U.early_n >= GAIN_MIN_SPK & isfinite(U.rB);
            if sum(ok) < GAIN_MIN_CELLS, continue; end
            z = zscore(log10(U.rB(ok) + 0.01)); y = U.(['early_' yl{yy}])(ok) - U.(['preB_' yl{yy}])(ok);
            X = [X; z]; Y = [Y; y]; Wt = [Wt; repmat(1/sum(ok), sum(ok), 1)]; %#ok<AGROW>
        end
        if isempty(X), continue; end
        ed = [-Inf -1 -0.33 0.33 1 Inf]; [~, ~, bi] = histcounts(X, ed);
        mu = accumarray(bi(bi > 0), Y(bi > 0) .* Wt(bi > 0), [5 1]) ./ accumarray(bi(bi > 0), Wt(bi > 0), [5 1]);
        plot(1:5, mu, '-o', 'Color', colG(g{1}), 'MarkerFaceColor', colG(g{1}), 'LineWidth', 1.5);
    end
    yline(0, 'k:'); xticks(1:5); xticklabels({'<-1','-1..-0.3','mid','0.3..1','>1'});
    xlabel('ConB1 activity (z, log rate)'); ylabel(ylab{yy}, 'Interpreter','tex');
    if yy == 1, legend(GR(~strcmp(GR, 'Control')), 'Location','best', 'Box','off'); end
end
title(tlG, sprintf('Change from the 30 min before ConB1 to the first %g min after, %s', POST_WIN/60, GAIN_CELLS));
exportgraphics(figG, fullfile(CACHE_ROOT, ['reactivationGain' OUT_TAG '.pdf']), 'ContentType','vector');
save(fullfile(CACHE_ROOT, ['compareReactivation' OUT_TAG '.mat']), 'T', 'B', 'TC', 'LT', 'LB', 'RPT', 'GN', 'gainS', 'opt', 'CONTROLS');

%% ---- figure ---------------------------------------------------------------------------------
col = containers.Map({'Control','SAL','ASV','SALpost','ASVrep'}, ...
    {[0.4 0.4 0.4], [0.2 0.6 1], [0.85 0.2 0.1], [0 0.25 0.6], [0.55 0.1 0.6]});
fig = figure('Color','w','Position',[60 60 1300 380]);
tiledlayout(1, 3, 'TileSpacing','compact');
for c = 1:2
    nexttile; hold on
    cn = CONTEXTS{c};
    for j = 1:numel(G)
        k = strcmp(T.Group, G{j}) & strcmp(T.Context, cn);
        x = T.dEV(k);
        plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col(G{j}), 'MarkerFaceColor', col(G{j}), 'MarkerSize', 4);
        if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', col(G{j}), 'LineWidth', 2); end
    end
    yline(0, 'k:'); xlim([0.5 numel(G)+0.5]); xticks(1:numel(G)); xticklabels(G);
    title(sprintf('%s: EV - REV', cn)); ylabel('dEV');
end
nexttile; hold on
for j = 1:numel(G)
    x = B.d(strcmp(B.Group, G{j}));
    plot(j + 0.15*(rand(size(x)) - 0.5), x, 'o', 'Color', col(G{j}), 'MarkerFaceColor', col(G{j}), 'MarkerSize', 4);
    if ~isempty(x), plot(j + [-0.25 0.25], mean(x)*[1 1], '-', 'Color', col(G{j}), 'LineWidth', 2); end
end
yline(0, 'k:'); xlim([0.5 numel(G)+0.5]); xticks(1:numel(G)); xticklabels(G);
title('novel - familiar (B1 - A1), same session'); ylabel('\Delta dEV');
exportgraphics(fig, fullfile(CACHE_ROOT, ['compareReactivation' OUT_TAG '.pdf']), 'ContentType','vector');


%% =========================================================================
function y = local_react(RP, fld, e, o)
% {nUnits, nPairs, nRip, thetaFrac, thetaSrc, rPre, rPost, rPP, EV, REV, nBias, bPre, bPost, biasPairs, pairVec, unitN, why}
% pairVec = [cT cPre cPost i j] per valid pair (i, j = positions in the kept units)
% biasPairs = [bRun bPre bPost i j];  unitN = [spikes in the template, in PRE ripples, in POST ripples] per kept unit
y = {0, 0, 0, NaN, '', NaN, NaN, NaN, NaN, NaN, 0, NaN, NaN, {zeros(0, 5)}, {zeros(0, 5)}, {zeros(0, 3)}, ''};
sf = dir(fullfile(fld, '*spikes.cellinfo.mat')); cf = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(sf) || isempty(cf), y{end} = 'no spikes or cell_metrics'; return; end
S = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); sp = S.spikes;
C = load(fullfile(cf(1).folder, cf(1).name), 'cell_metrics'); cm = C.cell_metrics;
ts = cellfun(@(x) reshape(x(~isnan(x)), [], 1), sp.times(:), 'UniformOutput', false);
nU = numel(ts);
clu = []; if isfield(sp, 'UID') && numel(sp.UID) == nU, clu = sp.UID(:); end
if isfield(sp, 'cluID') && numel(sp.cluID) == nU, clu = sp.cluID(:); end
loc = matchCellMetrics(cm, (1:nU)', clu, nU);
typ = repmat({''}, nU, 1); okL = loc > 0;
ct = cellstr(string(cm.putativeCellType(:))); typ(okL) = ct(loc(okL));
grp = nan(nU, 1);                     % pairs sharing this label are dropped
gf = '';
switch o.pairExclude
    case 'sameChannel', if isfield(cm, 'maxWaveformCh1'), gf = 'maxWaveformCh1'; elseif isfield(cm, 'maxWaveformCh'), gf = 'maxWaveformCh'; end
    case 'sameShank',   if isfield(cm, 'shankID'), gf = 'shankID'; end
end
if ~isempty(gf) && numel(cm.(gf)) == numel(ct), gv = double(cm.(gf)(:)); grp(okL) = gv(loc(okL)); end

% ---- theta part of the exposure ----
[mask, src] = local_runMask(fld, e.on, e.off, o.speedMin);         % 1 s windows: [start end] rows
if isempty(src), [mask, src] = local_thetaMask(RP, e.on, e.off, o.thetaRatio); end
y{4} = sum(diff(mask, 1, 2)) / (e.off - e.on); y{5} = src;
if isempty(mask), y{end} = 'no theta in exposure'; return; end
edges = []; for k = 1:size(mask, 1), edges = [edges; (mask(k,1):o.bin:mask(k,2)-o.bin)']; end %#ok<AGROW>
if numel(edges) < 20, y{end} = 'exposure too short'; return; end
X = zeros(numel(edges), nU);
for u = 1:nU
    t = ts{u}; t = t(t >= edges(1) & t < edges(end) + o.bin);
    if isempty(t), continue; end
    [~, ~, b] = histcounts(t, [edges; edges(end) + o.bin]);
    b = b(:); v = b > 0; v(v) = t(v) < edges(b(v)) + o.bin;   % spikes in gaps between theta windows fall out
    if any(v), X(:, u) = accumarray(b(v), 1, [numel(edges) 1]); end
end
keep = sum(X, 1)' >= o.minSpk;
if strcmp(o.cells, 'PYR'), keep = keep & contains(typ, 'Pyramidal'); end
% ---- ripples: per-ripple counts for the kept units ----
si = RP.units.SpikesIdx(:);
if max(si) > nU
    y{end} = sprintf('RIPPLES indexes unit %d, spikes file has %d (different sort): %s', max(si), nU, fld);
    return
end
uMap = zeros(nU, 1); uMap(si) = 1:numel(si);
keep = keep & uMap > 0;
idx = find(keep);
if numel(idx) < o.minUnits, y{end} = sprintf('< %d units', o.minUnits); return; end
preIv = local_preIv(e, o.preWin);       % PRE_WIN before the exposure, or the given intervals (e.preIv)
if ~isfield(e, 'postFixed') || ~e.postFixed
    e.post = [max(e.post(1), e.off), min(e.post(2), e.off + o.postWin)];   % POST_WIN after it
end
inPre = find(local_ripInIv(RP, preIv));
inPost = find(RP.ts1 >= e.post(1) & RP.ts1Off <= e.post(2));
if ~o.matchN                                   % every ripple in each block
    N = min(numel(inPre), numel(inPost));
    if N < o.minRip, y{end} = sprintf('< %d ripples in pre or post', o.minRip); return; end
    N = numel(inPost);                         % reported: ripples in the POST block
else                                           % the same number in PRE and POST
    N = min([numel(inPre), numel(inPost), o.maxRip]);
    if N < o.minRip, y{end} = sprintf('< %d ripples in pre or post', o.minRip); return; end
    inPre = inPre(end-N+1:end);
    if isfield(e, 'postFixed') && e.postFixed
        inPost = inPost(unique(round(linspace(1, numel(inPost), N))));   % time-course bin: spread over the bin
    else
        inPost = inPost(1:N);
    end
end
Q = local_ripCounts(RP, inPre, uMap(idx)); P = local_ripCounts(RP, inPost, uMap(idx));
cT = corr(X(:, idx)); cA = corr(Q); cB = corr(P);
pairMask = triu(true(numel(idx)), 1);
if any(isfinite(grp(idx))), pairMask = pairMask & (grp(idx) ~= grp(idx)'); end
[pI, pJ] = find(pairMask);
v = [cT(pairMask), cA(pairMask), cB(pairMask)];
okv = all(isfinite(v), 2); v = v(okv, :); pI = pI(okv); pJ = pJ(okv);
if size(v, 1) < 10, y{end} = '< 10 valid pairs'; return; end
r = corr(v);
rPre = r(1,2); rPost = r(1,3); rPP = r(2,3);
EV  = ((rPost - rPre*rPP) / sqrt((1 - rPre^2) * (1 - rPP^2)))^2;
REV = ((rPre - rPost*rPP) / sqrt((1 - rPost^2) * (1 - rPP^2)))^2;
% ---- Skaggs & McNaughton temporal bias ----
iv = local_merge(mask);
tRun = cellfun(@(t) t(local_inIntervals(t, iv)), ts(idx), 'UniformOutput', false);
tA = local_ripTimes(RP, inPre, uMap(idx)); tB = local_ripTimes(RP, inPost, uMap(idx));
[pi_, pj_] = find(pairMask);
BP = nan(numel(pi_), 3); nn = zeros(numel(pi_), 3);
for q = 1:numel(pi_)
    [BP(q,1), nn(q,1)] = local_bias(tRun{pi_(q)}, tRun{pj_(q)}, o.runWin);
    [BP(q,2), nn(q,2)] = local_bias(tA{pi_(q)}, tA{pj_(q)}, o.ripWin);
    [BP(q,3), nn(q,3)] = local_bias(tB{pi_(q)}, tB{pj_(q)}, o.ripWin);
end
okB = nn(:,1) >= o.minCCrun & nn(:,2) >= o.minCCrip & nn(:,3) >= o.minCCrip & all(isfinite(BP), 2);
BP = [BP(okB, :), pi_(okB), pj_(okB)]; bPre = NaN; bPost = NaN;
if size(BP, 1) >= 10, rb = corr(BP(:, 1:3)); bPre = rb(1,2); bPost = rb(1,3); end
UN = [sum(X(:, idx), 1)', sum(Q, 1)', sum(P, 1)'];
y = {numel(idx), size(v, 1), N, y{4}, src, rPre, rPost, rPP, EV, REV, size(BP, 1), bPre, bPost, {BP}, {[v, pI, pJ]}, {UN}, ''};
end


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


function [b, n] = local_bias(a, c, w)
% (spikes of c in (0, w] after spikes of a - in (-w, 0]) / their sum
b = NaN; n = 0;
if isempty(a) || isempty(c), return; end
c = unique(c(:)); a = a(:);
cnt = @(x) local_countLE(x, c);
nA = sum(cnt(a + w) - cnt(a)); nB = sum(cnt(a) - cnt(a - w));
n = nA + nB;
if n > 0, b = (nA - nB) / n; end
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


function g = local_group(RP, ev, sessionSet)
% Group label of one TransPlasticity session: 'SAL', 'ASV' or '' (not used).
%   sessionSet 'all'   : every SAL session is 'SAL' and every ASV session is 'ASV' (first and repeat).
%   sessionSet 'first' : the earlier design -- SAL sessions before the animal's first ASV, and the
%                        first-ASV session; SAL after the first ASV and repeat ASV are ''.
% A session needs an infusion time (RIPPLES.meta, else Infusion_start, else ConA1_end) either way.
if nargin < 3 || isempty(sessionSet), sessionSet = 'all'; end
m = RP.meta; g = '';
inf_ = m.infusionSec;
if ~isfinite(inf_) && nargin > 1 && ~isempty(ev)
    k = find(strcmp(ev.description, 'Infusion_start'), 1);
    if isempty(k), k = find(strcmp(ev.description, 'ConA1_end'), 1); end
    if ~isempty(k), inf_ = ev.time(k); end
end
if ~isfinite(inf_), return; end
if strcmp(sessionSet, 'all')
    if ismember(m.condition, {'SAL', 'ASV'}), g = m.condition; end
    return
end
noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
off = RP.ts2(:) - RP.ts1(:); off = median(off(isfinite(off)));
if noASV || isempty(off) || ~isfinite(off), after = false; else, after = inf_ + off > 0; end
isFirst = ~noASV && strcmp(m.folder, m.firstASV.session);
switch m.condition
    case 'SAL'
        if ~after, g = 'SAL'; end
    case 'ASV'
        if isFirst, g = 'ASV'; end
end
end


function B = local_pairBA(T, cb, ca)
% per session: dEV(context cb) - dEV(context ca)
s = unique(T.Session); rows = {};
for i = 1:numel(s)
    kb = strcmp(T.Session, s{i}) & strcmp(T.Context, cb); ka = strcmp(T.Session, s{i}) & strcmp(T.Context, ca);
    if any(kb) && any(ka)
        r = find(kb, 1);
        rows(end+1,:) = {T.Group{r}, T.Animal{r}, s{i}, T.dEV(r) - T.dEV(find(ka, 1))}; %#ok<AGROW>
    end
end
B = cell2table(rows, 'VariableNames', {'Group','Animal','Session','d'});
if isempty(rows), B = table(cell(0,1), cell(0,1), cell(0,1), zeros(0,1), 'VariableNames', {'Group','Animal','Session','d'}); end
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


function e = local_sem(d)
% standard error of the mean over finite values (NaN for fewer than 2)
d = d(isfinite(d)); e = NaN; if numel(d) > 1, e = std(d) / sqrt(numel(d)); end
end


function s = local_lmeOrder(t, nm, cv)
% ASV - SAL from y ~ Group + cv + (1|Animal) on the SAL and ASV rows with a finite covariate
s = 'not fitted';
k = ismember(t.Group, {'SAL', 'ASV'}) & isfinite(t.(nm)) & isfinite(t.(cv));
if sum(k) < 6 || numel(unique(t.Group(k))) < 2, s = sprintf('too few rows (%d)', sum(k)); return; end
D = table(t.(nm)(k), categorical(t.Group(k), {'SAL', 'ASV'}), t.(cv)(k), categorical(t.Animal(k)), ...
    'VariableNames', {'y', 'Group', 'cov', 'Animal'});
try
    lme = fitlme(D, 'y ~ Group + cov + (1|Animal)');
    [~, ~, C] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
    if isa(C, 'dataset'), C = dataset2table(C); end
    nmC = cellstr(string(C.Name)); kg = find(startsWith(nmC, 'Group'), 1); kc = find(strcmp(nmC, 'cov'), 1);
    s = sprintf('ASV - SAL = %+.3f +- %.3f, p %.3f | %s coef %+.3f +- %.3f, p %.3f | %d sessions, %d animals', ...
        C.Estimate(kg), C.SE(kg), C.pValue(kg), cv, C.Estimate(kc), C.SE(kc), C.pValue(kc), height(D), numel(unique(D.Animal)));
catch ME
    s = ME.message;
end
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


function rows = local_tc(RP, fld, EA, EB, postIv, binW, tmax, o, g, animal, file)
% one row per template x bin of the home cage after ConB1
rows = {};
tEnd = min(postIv(2), postIv(1) + tmax);
edges = postIv(1):binW:tEnd;
if tEnd - edges(end) >= binW / 2, edges(end+1) = tEnd; end    % keep a final partial bin of >= half a bin
for E = [EA, EB]
    % the same number of ripples in every bin (and in PRE): N = the smallest bin that passes o.minRip
    nPre = sum(local_ripInIv(RP, local_preIv(E, o.preWin)));
    nBin = arrayfun(@(b) sum(RP.ts1 >= edges(b) & RP.ts1Off <= edges(b+1)), 1:numel(edges) - 1);
    N = min([nPre, nBin(nBin >= o.minRip)]);
    oo = o; if o.matchN && ~isempty(N), oo.maxRip = N; end
    for b = 1:numel(edges) - 1
        e = E; e.post = edges(b:b+1);
        y = local_react(RP, fld, e, oo);
        rows(end+1,:) = {g, animal, file, E.name, b, edges(b) - postIv(1), y{3}, y{9}, y{10}, y{6}, y{7}, y{end}, y{15}, y{16}}; %#ok<AGROW>
    end
end
end


function s = local_rethresh(PV, UN, BP, minSpk, minRipSpk, minUnits)
% Drop units below the new thresholds and recompute EV / REV (and the temporal-bias correlations)
% from the cached pair values. Exact: each pair's template / PRE / POST correlations depend only on
% its two units, so removing units removes pairs and changes nothing else.
% PV = [cT cPre cPost i j], BP = [bRun bPre bPost i j], UN = [nTemplate nPreRip nPostRip] per unit.
if iscell(PV), PV = PV{1}; end
if iscell(UN), UN = UN{1}; end
if iscell(BP), BP = BP{1}; end
s = struct('nUnits', 0, 'nPairs', 0, 'rPre', NaN, 'rPost', NaN, 'rPP', NaN, 'EV', NaN, 'REV', NaN, ...
    'nBias', 0, 'bPre', NaN, 'bPost', NaN, 'PV', zeros(0, 5), 'BP', zeros(0, 5), 'why', '');
if isempty(UN) || isempty(PV), s.why = 'no cached pairs / unit counts'; return; end
ku = UN(:,1) >= minSpk & UN(:,2) >= minRipSpk & UN(:,3) >= minRipSpk;
s.nUnits = sum(ku);
if s.nUnits < minUnits, s.why = sprintf('< %d units', minUnits); return; end
PV = PV(ku(PV(:,4)) & ku(PV(:,5)), :); s.PV = PV; s.nPairs = size(PV, 1);
if s.nPairs < 10, s.why = '< 10 valid pairs'; return; end
r = corr(PV(:, 1:3)); s.rPre = r(1,2); s.rPost = r(1,3); s.rPP = r(2,3);
s.EV  = ((s.rPost - s.rPre*s.rPP) / sqrt((1 - s.rPre^2) * (1 - s.rPP^2)))^2;
s.REV = ((s.rPre - s.rPost*s.rPP) / sqrt((1 - s.rPost^2) * (1 - s.rPP^2)))^2;
if size(BP, 2) >= 5
    BP = BP(ku(BP(:,4)) & ku(BP(:,5)), :); s.BP = BP; s.nBias = size(BP, 1);
    if s.nBias >= 10, rb = corr(BP(:, 1:3)); s.bPre = rb(1,2); s.bPost = rb(1,3); end
end
end


function d = local_evFromPairs(V)
% EV - REV from pair vectors [cT cPre cPost ...]
d = NaN;
if size(V, 1) < 5, return; end
r = corr(V(:, 1:3)); a = r(1,2); b = r(1,3); c = r(2,3);
EV  = ((b - a*c) / sqrt((1 - a^2) * (1 - c^2)))^2;
REV = ((a - b*c) / sqrt((1 - b^2) * (1 - c^2)))^2;
d = EV - REV;
end


function rows = local_block(RP, fld, EA, EB, postIv, t0, t1, o, g, animal, file)
% one row per template: every ripple from t0 to t1 after the start of the home cage after ConB1
rows = {};
blk = [postIv(1) + t0, min(postIv(2), postIv(1) + t1)];
if diff(blk) <= 0, return; end
for E = [EA, EB]
    e = E; e.post = blk;
    y = local_react(RP, fld, e, o);
    rows(end+1,:) = {g, animal, file, E.name, y{3}, y{9}, y{10}, y{6}, y{7}, y{end}, y{15}, y{16}}; %#ok<AGROW>
end
end


function U = local_gainSession(RP, fld, EA, EB, postIv, o, lateStart, tMax)
% per unit: exposure rates (running/theta part of ConA1, ConB1) and, in blocks preA, preB, early,
% late, in-ripple rate, out-of-ripple rate, ripple participation and its Poisson chance level
U = table();
sf = dir(fullfile(fld, '*spikes.cellinfo.mat')); cf = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(sf) || isempty(cf), return; end
S = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); sp = S.spikes;
C = load(fullfile(cf(1).folder, cf(1).name), 'cell_metrics'); cm = C.cell_metrics;
ts = cellfun(@(x) sort(reshape(x(~isnan(x)), [], 1)), sp.times(:), 'UniformOutput', false);
nU = numel(ts);
clu = []; if isfield(sp, 'UID') && numel(sp.UID) == nU, clu = sp.UID(:); end
if isfield(sp, 'cluID') && numel(sp.cluID) == nU, clu = sp.cluID(:); end
loc = matchCellMetrics(cm, (1:nU)', clu, nU); okL = loc > 0;
ct = repmat({''}, nU, 1); ly = ct;
c0 = cellstr(string(cm.putativeCellType(:))); ct(okL) = c0(loc(okL));
if isfield(cm, 'deepSuperficial') && numel(cm.deepSuperficial) == numel(c0)
    l0 = cellstr(string(cm.deepSuperficial(:))); ly(okL) = l0(loc(okL));
end
isP = contains(ct, 'Pyramidal');
U = table((1:nU)', isP, isP & strcmp(ly, 'Deep'), isP & strcmp(ly, 'Superficial'), 'VariableNames', {'unit','PYR','PYRdeep','PYRsup'});
% exposure rates in the running / theta windows
for E = [EA, EB]
    [mask, src] = local_runMask(fld, E.on, E.off, o.speedMin);
    if isempty(src), mask = local_thetaMask(RP, E.on, E.off, o.thetaRatio); end
    iv = local_merge(mask); L = sum(diff(iv, 1, 2));
    r = nan(nU, 1);
    if L > 0, r = cellfun(@(t) sum(local_inIntervals(t, iv)), ts) / L; end
    if strcmp(E.name, 'ConA1'), U.rA = r; else, U.rB = r; end
end
% blocks
pA = local_preIv(EA, o.preWin); pB = local_preIv(EB, o.preWin);   % n x 2 intervals
st = max(postIv(1), EB.off);
blk = {'preA', pA; 'preB', pB; 'early', [st, min(postIv(2), st + o.postWin)]; ...
       'late', [postIv(1) + lateStart, min(postIv(2), postIv(1) + tMax)]};
for k = 1:size(blk, 1)
    w = blk{k,2}; nm = blk{k,1}; w = w(diff(w, 1, 2) > 0, :); dur = sum(diff(w, 1, 2));
    v = nan(nU, 5);                                   % rRip, rOut, p, pChance, nSpk
    if dur > 0
        in = find(local_ripInIv(RP, w));
        on = RP.ts1(in); of = RP.ts1Off(in);
        [on, o_] = sort(on(:)); of = of(o_);
        kp = of > on & [true; on(2:end) > of(1:end-1)];    % strictly increasing edges
        on = on(kp); of = of(kp); d = of - on; ripT = sum(d);
        if numel(in) >= 10 && ripT > 0
            ed = reshape([on(:)'; of(:)'], [], 1);      % ripples do not overlap: odd bins are inside
            for u = 1:nU
                t = ts{u}; t = t(local_inIntervals(t, w));
                if numel(ed) < 2, continue; end
                cnt = histcounts(t, ed);
                inR = cnt(1:2:end); nIn = sum(inR); nOut = numel(t) - nIn;
                rOut = nOut / max(dur - ripT, eps);
                v(u,:) = [nIn / ripT, rOut, mean(inR > 0), mean(1 - exp(-rOut * d)), numel(t)];
            end
        end
    end
    U.([nm '_lrip']) = log10(v(:,1) + 0.1);
    U.([nm '_lout']) = log10(v(:,2) + 0.01);
    U.([nm '_xp']) = v(:,3) - v(:,4);
    U.([nm '_n']) = v(:,5);
end
end
