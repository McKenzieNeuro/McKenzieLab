%% reactivationCoupling.m
% Does the global change after saline (rates up, shared gain a up) enter EV / REV?
%
% WHY. After saline, novel-context EV - REV is low because REV is high (0.068, against 0.008
% after ASV and 0.006 in controls): the ripples before ConB1 already resemble the ConB1
% template. Those PRE ripples are the home cage between ConA1 and ConB1, i.e. after the infusion,
% when saline rates are rising. A rise in global coupling adds a shared component to every
% pairwise correlation - in the template (theta bins) and in PRE and POST ripples alike - which
% raises rPre and rPP together, and with them REV, whatever the experience did. A rate rise also
% sharpens correlation estimates in the block where it occurs. Either would make the saline value
% unreadable, and the comparison it would distort is the one the design exists for (SAL vs ASV
% within animal).
%
% WHAT. The template / PRE / POST count matrices of compareReactivation.m (local_react, same
% blocks, masks, units and pairs), scored four ways:
%   raw        as published (check: must reproduce compareReactivation)
%   cond       each cell's counts residualised on the leave-one-out population count of the same
%              bin / ripple before correlating: the correlation not carried by population
%              co-fluctuation (population coupling, Okun et al. 2015), in all three blocks
%   thin       in whichever of PRE / POST has the higher mean count per ripple for a cell, its
%              spikes are binomially thinned to the other block's mean: no rate asymmetry between
%              the blocks (N_THIN draws averaged)
%   condOther  as cond, but the population count comes from units OUTSIDE the scored set
%              (interneurons, cells below threshold). With few scored cells the leave-one-out count
%              is mostly the other members of the same assembly, so cond removes the assembly's own
%              covariance and is harsher the lower the yield; condOther cannot contain the
%              scored assembly.
% CtrlMatch: the control sessions re-scored on N_DRAW random subsets of their units, each the
% median iPKR unit count, averaged. Any iPKR-vs-control difference that appears under cond but
% not between Control and CtrlMatch-vs-iPKR is yield, not drug.
% and the coupling itself: rbar = mean pairwise correlation of raw counts over the scored pairs,
% for the template and PRE / POST ripples. dPre = rbar(PRE of ConB1) - rbar(PRE of ConA1), same
% session: drug on board minus drug free, which is the global change in question.
%
% READING (stated before the run)
%  1. dPre SAL - ASV > 0: saline raises coupling in exactly the block that feeds REV.
%  2. cond lowers saline REV toward ASV / control and raises saline dEV: the low saline value was
%     the global change. Then cond dEV SAL = ASV -> reactivation is preserved and orthogonal to
%     the global change; cond dEV SAL > ASV -> ASV reduces pattern reactivation and the global
%     change was masking it (and ASV = control must then be re-checked under cond as well).
%  3. cond leaves saline REV high: REV reflects structure in the PRE ripples other than
%     population coupling, and the global change is not the explanation.
%  4. thin moving dEV: the rate asymmetry between PRE and POST shifts EV on its own.
%  5. (added after run 2, where cond dropped ASV dEV 0.080 -> 0.020 and left controls at 0.085)
%     CtrlMatch dEV falling under cond as ASV's did -> that drop is yield; ASV vs CtrlMatch
%     equal under cond and condOther -> reactivation preserved. ASV below CtrlMatch under
%     condOther -> ASV's raw reactivation is carried by population co-fluctuation.
%  Across infused sessions, corr(dPre, REV) and corr(dPre, dEV) are descriptive at this n.
%
% Needs on the path: getAllExtFiles, matchCellMetrics, controlAnimal, permtest.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT       = 'R:\TransPlasticity';
% ---- must equal compareReactivation.m ----
CELLS = 'PYR'; BIN = 0.1; THETA_RATIO = 2; PAIR_EXCLUDE = 'none'; SPEED_MIN = 3;
MIN_SPK = 100; MIN_UNITS = 10;   % the thresholds behind the published values (Methods)
MIN_RIP = 50; MAX_RIP = Inf; MATCH_N = false;
EXCL_INF = false; INF_EXCL_WIN = [0 1800]; PRE_B_ALLGAP = true; PRE_WIN = 3600; POST_WIN = 1800;
% ---- this script ----
N_THIN = 20;
CTRL_MATCH = true;     % also score controls with their units subsampled to the median iPKR count
N_DRAW = 20;           % random unit subsets per control exposure (CtrlMatch)
CTX = {'ConA1', 'ConB1'};
VAR = {'raw', 'cond', 'condOther', 'thin'};
NPERM = 1e4;
rng(1);
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
opt = struct('preWin', PRE_WIN, 'postWin', POST_WIN, 'pairExclude', PAIR_EXCLUDE, 'cells', CELLS, ...
    'bin', BIN, 'thetaRatio', THETA_RATIO, 'speedMin', SPEED_MIN, 'minSpk', MIN_SPK, ...
    'minUnits', MIN_UNITS, 'minRip', MIN_RIP, 'maxRip', MAX_RIP, 'matchN', MATCH_N, 'nThin', N_THIN, 'nDraw', N_DRAW);

rows = {};
%% ---- TransPlasticity (SAL / ASV; SESSION_SET decides which sessions) ----
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
for i = 1:numel(fa)
    fld = fileparts(fa{i});
    ev = local_events(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConA1_start')), continue; end
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES;
    g = local_group(RP, ev, SESSION_SET);
    if ~ismember(g, {'SAL', 'ASV'}), continue; end
    HC = local_hcBlocks(ev, RP, EXCL_INF, INF_EXCL_WIN);
    for c = 1:numel(CTX)
        on = ev.time(strcmp(ev.description, [CTX{c} '_start']));
        off = ev.time(strcmp(ev.description, [CTX{c} '_end']));
        if isempty(on) || isempty(off), continue; end
        e = struct('name', CTX{c}, 'on', on(1), 'off', off(1), 'pre', [], 'post', [], 'preIv', []);
        kp = find(abs(HC(:,2) - e.on) <= 1, 1); kq = find(abs(HC(:,1) - e.off) <= 1, 1);
        if PRE_B_ALLGAP && strcmp(CTX{c}, 'ConB1')
            e.preIv = local_gapIv(HC, ev, e.on); kp = [];
            if ~isempty(e.preIv), kp = 0; end
        end
        if isempty(kp) || isempty(kq), continue; end
        if kp > 0, e.pre = HC(kp,:); else, e.pre = [e.preIv(1,1) e.preIv(end,2)]; end
        e.post = HC(kq,:);
        s = local_score(RP, fld, e, opt, 0);
        rows(end+1,:) = {g, RP.meta.animal, fa{i}, e.name, s}; %#ok<AGROW>
    end
end
%% ---- controls (after TransPlasticity, so the matched unit count is known) ----
okTP = cellfun(@(x) isempty(x.why), rows(:,5));
mUnits = round(median(cellfun(@(x) x.nUnits, rows(okTP,5))));
fprintf('iPKR exposures: median %d scored units; CtrlMatch subsamples controls to this (%d draws)\n', mUnits, N_DRAW);
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), warning('rc:file', 'missing %s', rf); continue; end
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
        s = local_score(RP, fld, e, opt, 0);
        rows(end+1,:) = {'Control', controlAnimal(rf), rf, e.name, s}; %#ok<AGROW>
        if CTRL_MATCH
            s = local_score(RP, fld, e, opt, mUnits);
            rows(end+1,:) = {'CtrlMatch', controlAnimal(rf), [rf '#matched'], e.name, s}; %#ok<AGROW>
        end
    end
end

%% ---- table ----
ok = cellfun(@(s) isempty(s.why), rows(:,5));
fprintf('%d of %d exposures scored; not scored:\n', sum(ok), numel(ok));
for k = find(~ok)', fprintf('   %-8s %-5s %s\n', rows{k,1}, rows{k,4}, rows{k,5}.why); end
rows = rows(ok, :);
T = cell2table(rows(:,1:4), 'VariableNames', {'Group','Animal','Session','Context'});
S = [rows{:,5}]';
T.nUnits = [S.nUnits]'; T.nOther = [S.nOther]'; T.nPairs = [S.nPairs]'; T.nPre = [S.nPre]'; T.nPost = [S.nPost]';
T.ratePre = [S.ratePre]'; T.ratePost = [S.ratePost]';
T.rbarT = [S.rbarT]'; T.rbarPre = [S.rbarPre]'; T.rbarPost = [S.rbarPost]';
fld5 = {'rPre','rPost','rPP','EV','REV'};
for v = 1:numel(VAR)
    M = vertcat(S.(VAR{v}));
    for f = 1:5, T.([fld5{f} '_' VAR{v}]) = M(:, f); end
    T.(['dEV_' VAR{v}]) = M(:,4) - M(:,5);
end
save(fullfile(CACHE_ROOT, 'reactivationCoupling.mat'), 'T', 'opt');

GR = {'Control', 'CtrlMatch', 'SAL', 'ASV'};
B1 = T(strcmp(T.Context, 'ConB1'), :); A1 = T(strcmp(T.Context, 'ConA1'), :);

%% ---- 1. the global change: coupling and rate in the blocks that feed EV / REV ----
fprintf('\n===== 1. coupling (mean pairwise r of raw counts) and in-ripple rate (spikes / ripple / cell) =====\n');
fprintf('   %-9s %-5s %8s %8s %8s | %8s %8s | %6s %6s | n\n', 'group', 'ctx', 'tmpl', 'PRE', 'POST', 'ratePRE', 'ratePOST', 'units', 'other');
for c = 1:2
    D = T(strcmp(T.Context, CTX{c}), :);
    for g = GR
        k = strcmp(D.Group, g{1}); if ~any(k), continue; end
        fprintf('   %-9s %-5s %8.3f %8.3f %8.3f | %8.3f %8.3f | %6.0f %6.0f | %d\n', g{1}, CTX{c}, mean(D.rbarT(k)), ...
            mean(D.rbarPre(k)), mean(D.rbarPost(k)), mean(D.ratePre(k)), mean(D.ratePost(k)), ...
            median(D.nUnits(k)), median(D.nOther(k)), sum(k));
    end
end
W = innerjoin(B1(:, {'Group','Animal','Session','rbarPre','rbarPost','ratePre','REV_raw','dEV_raw','REV_cond','dEV_cond','dEV_condOther'}), ...
    A1(:, {'Session','rbarPre','ratePre','dEV_raw','dEV_cond','dEV_condOther'}), 'Keys', 'Session');
W.dPre = W.rbarPre_left - W.rbarPre_right;           % PRE of ConB1 - PRE of ConA1 (drug on board - drug free)
W.dRate = log2(W.ratePre_left ./ W.ratePre_right);
fprintf('\nwithin session, PRE of ConB1 minus PRE of ConA1 (infused: drug on board minus drug free)\n');
for nm = {'dPre', 'dRate'}
    for g = GR
        k = strcmp(W.Group, g{1}); if sum(k) < 2, continue; end
        am = local_animalMean(W(k,:), nm{1});
        fprintf('   %-6s %-8s %+.3f +- %.3f  (animals p = %.3f, n = %d/%d)\n', nm{1}, g{1}, mean(W.(nm{1})(k)), ...
            std(W.(nm{1})(k)) / sqrt(sum(k)), local_signflip(am), sum(k), numel(am));
    end
    [md, pp, na] = local_paired(W, nm{1}, 'SAL', 'ASV');
    fprintf('   %-6s ASV - SAL paired %+.3f p = %.3f (%d animals) | LME p = %.3f\n', nm{1}, md, pp, na, local_lme(W, nm{1}, 'SAL', 'ASV'));
end

%% ---- 2. EV / REV for the novel context under each correction ----
fprintf('\n===== 2. ConB1: EV, REV, dEV under each correction =====\n');
for v = 1:numel(VAR)
    fprintf('\n-- %s --\n', VAR{v});
    for g = GR
        k = strcmp(B1.Group, g{1}); if ~any(k), continue; end
        x = B1.(['dEV_' VAR{v}])(k); am = local_animalMean(B1(k,:), ['dEV_' VAR{v}]);
        nf = sum(isfinite(x)); miss = '';
        if nf < sum(k), miss = sprintf('  [%d of %d sessions scored]', nf, sum(k)); end
        fprintf('   %-8s rPre %.3f rPost %.3f rPP %.3f | EV %.3f REV %.3f | dEV %+.3f +- %.3f (>0 animals p = %.3f) n = %d/%d%s\n', g{1}, ...
            local_nmean(B1.(['rPre_' VAR{v}])(k)), local_nmean(B1.(['rPost_' VAR{v}])(k)), local_nmean(B1.(['rPP_' VAR{v}])(k)), ...
            local_nmean(B1.(['EV_' VAR{v}])(k)), local_nmean(B1.(['REV_' VAR{v}])(k)), local_nmean(x), local_nsem(x), ...
            local_signflip(am), sum(k), numel(am), miss);
    end
    nm = ['dEV_' VAR{v}];
    [md, pp, na] = local_paired(B1, nm, 'SAL', 'ASV');
    fprintf('   ASV - SAL: paired %+.3f p = %.3f (%d animals) | LME p = %.3f\n', md, pp, na, local_lme(B1, nm, 'SAL', 'ASV'));
    [md, pp, na] = local_paired(B1, ['REV_' VAR{v}], 'SAL', 'ASV');
    fprintf('   REV ASV - SAL: paired %+.3f p = %.3f (%d animals)\n', md, pp, na);
    for cg = {'Control', 'CtrlMatch'}
        a_ = B1.(nm)(strcmp(B1.Group, 'ASV')); c_ = B1.(nm)(strcmp(B1.Group, cg{1}));
        if numel(c_) < 2, continue; end
        fprintf('   ASV vs %-9s %+.3f vs %+.3f, perm p = %.3f\n', [cg{1} ':'], local_nmean(a_), local_nmean(c_), local_perm(a_, c_, NPERM));
    end
end
% within session: novel (drugged) minus familiar (drug-free), differenced across arms
fprintf('\nwithin session dEV(ConB1) - dEV(ConA1), ASV - SAL paired over animals\n');
for v = {'raw', 'cond', 'condOther'}
    W.(['dd_' v{1}]) = W.(['dEV_' v{1} '_left']) - W.(['dEV_' v{1} '_right']);
    [md, pp, na] = local_paired(W, ['dd_' v{1}], 'SAL', 'ASV');
    fprintf('   %-5s SAL %+.3f ASV %+.3f | ASV - SAL %+.3f p = %.3f (%d animals)\n', v{1}, ...
        local_nmean(W.(['dd_' v{1}])(strcmp(W.Group, 'SAL'))), local_nmean(W.(['dd_' v{1}])(strcmp(W.Group, 'ASV'))), md, pp, na);
    % drug-free familiar exposure as each group's own reference: a deficit present in ConA1 too is
    % cohort or day, not drug; only a deficit in ConB1 beyond ConA1 can be the drug
    for g = GR
        k = strcmp(W.Group, g{1}); if ~any(k), continue; end
        nf = sum(isfinite(W.(['dd_' v{1}])(k)));
        fprintf('        %-9s ConA1 %+.3f  ConB1 %+.3f  B1-A1 %+.3f  (n = %d scored of %d)\n', g{1}, local_nmean(W.(['dEV_' v{1} '_right'])(k)), ...
            local_nmean(W.(['dEV_' v{1} '_left'])(k)), local_nmean(W.(['dd_' v{1}])(k)), nf, sum(k));
    end
    a_ = W.(['dd_' v{1}])(strcmp(W.Group, 'ASV')); c_ = W.(['dd_' v{1}])(strcmp(W.Group, 'CtrlMatch'));
    if numel(c_) >= 2, fprintf('        B1-A1 ASV vs CtrlMatch perm p = %.3f\n', local_perm(a_, c_, NPERM)); end
end

%% ---- 3. does the size of the global change predict REV / dEV across infused sessions? ----
fprintf('\n===== 3. across infused sessions (SAL + ASV pooled, Spearman; descriptive) =====\n');
k = ismember(W.Group, {'SAL', 'ASV'});
for pr = {'dPre','REV_raw'; 'dPre','dEV_raw_left'; 'dPre','REV_cond'; 'dPre','dEV_cond_left'; 'dRate','REV_raw'; 'dRate','dEV_raw_left'}'
    kk = k & isfinite(W.(pr{1})) & isfinite(W.(pr{2}));
    if sum(kk) < 5, continue; end
    [rr, pv] = corr(W.(pr{1})(kk), W.(pr{2})(kk), 'Type', 'Spearman');
    fprintf('   %-6s vs %-14s rho = %+.2f (p = %.3f, n = %d)\n', pr{1}, pr{2}, rr, pv, sum(kk));
end
% shared gain a from compareReactivation (placebo span preA -> preB contains the infusion)
gf = fullfile(CACHE_ROOT, 'compareReactivation.mat');
if exist(gf, 'file')
    G = load(gf, 'GN'); GN = G.GN;
    ga = GN(strcmp(GN.Base, 'preA') & strcmp(GN.Post, 'preB') & strcmp(GN.Y, 'lrip'), {'Session', 'a'});
    J = innerjoin(W(k,:), ga, 'Keys', 'Session');
    for nm = {'REV_raw', 'dEV_raw_left', 'dPre'}
        kk = isfinite(J.a) & isfinite(J.(nm{1}));
        if sum(kk) < 5, continue; end
        [rr, pv] = corr(J.a(kk), J.(nm{1})(kk), 'Type', 'Spearman');
        fprintf('   a(preA->preB) vs %-14s rho = %+.2f (p = %.3f, n = %d)\n', nm{1}, rr, pv, sum(kk));
    end

    %% ---- 4. is what conditioning removes the experience-tied gain? ----
    % bB (per-cell decomposition, lrip) is a rate effect: cells active in ConB1 fire more in the
    % following ripples. In per-ripple counts that is differential coupling to the population,
    % which is what cond / condOther regress out. If the EV lost to conditioning is that gain, the
    % loss (raw - cond, ConB1) should rise with bB across sessions. Read within group: pooling
    % across groups would mix in the cohort difference in coupling.
    fprintf('\n===== 4. EV - REV lost to conditioning vs experience-tied gain bB, across sessions (Spearman) =====\n');
    L4 = B1(:, {'Group', 'Session', 'dEV_raw', 'dEV_cond', 'dEV_condOther'});
    L4.loss_cond = L4.dEV_raw - L4.dEV_cond; L4.loss_other = L4.dEV_raw - L4.dEV_condOther;
    L4.key = erase(L4.Session, '#matched');       % CtrlMatch shares the control session's bB
    for bs = {'preB', 'preA'}
        gb = GN(strcmp(GN.Base, bs{1}) & strcmp(GN.Post, 'early') & strcmp(GN.Y, 'lrip'), {'Session', 'bB'});
        gb.Properties.VariableNames{'Session'} = 'key';
        J4 = innerjoin(L4, gb, 'Keys', 'key');
        fprintf('  bB from %s -> early\n', bs{1});
        for g = {'ASV', 'SAL', 'Control', 'CtrlMatch'}
            for nm = {'loss_cond', 'loss_other'}
                kk = strcmp(J4.Group, g{1}) & isfinite(J4.bB) & isfinite(J4.(nm{1}));
                if sum(kk) < 4, continue; end
                [rr, pv] = corr(J4.bB(kk), J4.(nm{1})(kk), 'Type', 'Spearman');
                fprintf('   %-9s %-10s mean loss %+.3f | rho = %+.2f (p = %.3f, n = %d)\n', g{1}, nm{1}, ...
                    local_nmean(J4.(nm{1})(kk)), rr, pv, sum(kk));
            end
        end
    end
    fprintf('  rho > 0 within ASV: the conditioned-away reactivation is the spared gain component.\n');
    fprintf('  Same sign in controls with a small loss: the cohorts differ in how much reactivation the\n');
    fprintf('  gain carries, not in kind. n = 7-9 sessions per group: report as direction, not test.\n');
else
    fprintf('\n(section 4 skipped: %s not found)\n', gf);
end
fprintf('\nCheck first: raw dEV (ConB1) should reproduce compareReactivation run with the same SESSION_SET (''%s'')\n', SESSION_SET);
fprintf('(the ConB1 dEV row there). If it does not, the blocks here differ and nothing below it is readable.\n');


%% ======================================================================================


function s = local_score(RP, fld, e, o, mUnits)
% mUnits > 0: score o.nDraw random subsets of mUnits of the kept units and average (yield match)
s = struct('why', '', 'nUnits', 0, 'nOther', 0, 'nPairs', 0, 'nPre', 0, 'nPost', 0, 'ratePre', NaN, 'ratePost', NaN, ...
    'rbarT', NaN, 'rbarPre', NaN, 'rbarPost', NaN, 'raw', nan(1,5), 'cond', nan(1,5), 'condOther', nan(1,5), 'thin', nan(1,5));
[X, Q, P, pm, why, O] = local_counts(RP, fld, e, o);
if ~isempty(why), s.why = why; return; end
n = size(X, 2); s.nPre = size(Q, 1); s.nPost = size(P, 1); s.nOther = O.n;
if mUnits > 0 && n > mUnits
    sub = arrayfun(@(d) sort(randperm(n, mUnits)), 1:o.nDraw, 'UniformOutput', false);
else
    sub = {1:n};
end
fl = {'raw', 'cond', 'condOther', 'thin'}; acc = struct();
for f = fl, acc.(f{1}) = nan(numel(sub), 5); end
sc = nan(numel(sub), 6);
for d = 1:numel(sub)
    u = sub{d}; pmd = pm(u, u);
    r = local_variants(X(:,u), Q(:,u), P(:,u), pmd, O, o);
    if isempty(r), continue; end
    for f = fl, acc.(f{1})(d,:) = r.(f{1}); end
    sc(d,:) = [r.nPairs, r.rbarT, r.rbarPre, r.rbarPost, mean(mean(Q(:,u))), mean(mean(P(:,u)))];
end
if all(isnan(sc(:,1))), s.why = '< 10 valid pairs'; return; end
m = mean(sc, 1, 'omitnan');
s.nUnits = numel(sub{1}); s.nPairs = m(1); s.rbarT = m(2); s.rbarPre = m(3); s.rbarPost = m(4);
s.ratePre = m(5); s.ratePost = m(6);
for f = fl, s.(f{1}) = mean(acc.(f{1}), 1, 'omitnan'); end
end


function r = local_variants(X, Q, P, pm, O, o)
% raw | cond (leave-one-out population of the scored units) | condOther (population of the units
% outside the scored set) | thin (rate-matched PRE / POST); [] if too few pairs
r = [];
cT = corr(X); cA = corr(Q); cB = corr(P);
v = [cT(pm), cA(pm), cB(pm)];
pm(pm) = all(isfinite(v), 2);
if sum(pm(:)) < 10, return; end
r.nPairs = sum(pm(:));
r.rbarT = mean(cT(pm)); r.rbarPre = mean(cA(pm)); r.rbarPost = mean(cB(pm));
r.raw = local_ev(cT, cA, cB, pm);
r.cond = local_ev(local_condCorr(X, []), local_condCorr(Q, []), local_condCorr(P, []), pm);
r.condOther = nan(1, 5);
if O.n >= 3, r.condOther = local_ev(local_condCorr(X, O.x), local_condCorr(Q, O.q), local_condCorr(P, O.p), pm); end
th = nan(o.nThin, 5);
for d = 1:o.nThin
    [Qt, Pt] = local_thin(Q, P);
    th(d,:) = local_ev(cT, corr(Qt), corr(Pt), pm);
end
r.thin = mean(th, 1, 'omitnan');
end


function y = local_ev(cT, cA, cB, pm)
% [rPre rPost rPP EV REV] over the pairs in pm (Kudrimoti et al. 1999)
y = nan(1, 5);
v = [cT(pm), cA(pm), cB(pm)]; v = v(all(isfinite(v), 2), :);
if size(v, 1) < 10, return; end
r = corr(v); rPre = r(1,2); rPost = r(1,3); rPP = r(2,3);
EV  = ((rPost - rPre*rPP) / sqrt((1 - rPre^2) * (1 - rPP^2)))^2;
REV = ((rPre - rPost*rPP) / sqrt((1 - rPost^2) * (1 - rPP^2)))^2;
y = [rPre, rPost, rPP, EV, REV];
end


function C = local_condCorr(M, ext)
% correlation of each unit's counts after regressing out a population count of the same bin /
% ripple: ext = [] -> leave-one-out sum of the scored units; else the given external count
s = sum(M, 2); R = zeros(size(M));
for j = 1:size(M, 2)
    if isempty(ext), pop = s - M(:,j); else, pop = ext(:); end
    Xd = [ones(size(s)), pop];
    R(:,j) = M(:,j) - Xd * (Xd \ M(:,j));
end
C = corr(R);
end


function [Qt, Pt] = local_thin(Q, P)
% per unit, thin the block with the higher mean count per ripple to the other's mean
Qt = Q; Pt = P;
mq = mean(Q, 1); mp = mean(P, 1);
for j = 1:size(Q, 2)
    if mq(j) > mp(j) && mq(j) > 0
        Qt(:,j) = binornd(Q(:,j), mp(j) / mq(j));
    elseif mp(j) > mq(j) && mp(j) > 0
        Pt(:,j) = binornd(P(:,j), mq(j) / mp(j));
    end
end
end


function [X, Q, P, pairMask, why, O] = local_counts(RP, fld, e, o)
% template (theta / running bins) and PRE / POST per-ripple counts for the kept units:
% local_react of compareReactivation.m up to the correlations, unchanged
X = []; Q = []; P = []; pairMask = []; why = ''; O = struct('n', 0, 'x', [], 'q', [], 'p', []);
sf = dir(fullfile(fld, '*spikes.cellinfo.mat')); cf = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(sf) || isempty(cf), why = 'no spikes or cell_metrics'; return; end
S = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); sp = S.spikes;
C = load(fullfile(cf(1).folder, cf(1).name), 'cell_metrics'); cm = C.cell_metrics;
ts = cellfun(@(x) reshape(x(~isnan(x)), [], 1), sp.times(:), 'UniformOutput', false);
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
[mask, src] = local_runMask(fld, e.on, e.off, o.speedMin);
if isempty(src), mask = local_thetaMask(RP, e.on, e.off, o.thetaRatio); end
if isempty(mask), why = 'no theta in exposure'; return; end
edges = []; for k = 1:size(mask, 1), edges = [edges; (mask(k,1):o.bin:mask(k,2)-o.bin)']; end %#ok<AGROW>
if numel(edges) < 20, why = 'exposure too short'; return; end
X = zeros(numel(edges), nU);
for u = 1:nU
    t = ts{u}; t = t(t >= edges(1) & t < edges(end) + o.bin);
    if isempty(t), continue; end
    [~, ~, b] = histcounts(t, [edges; edges(end) + o.bin]);
    b = b(:); v = b > 0; v(v) = t(v) < edges(b(v)) + o.bin;
    if any(v), X(:, u) = accumarray(b(v), 1, [numel(edges) 1]); end
end
keep = sum(X, 1)' >= o.minSpk;
if strcmp(o.cells, 'PYR'), keep = keep & contains(typ, 'Pyramidal'); end
si = RP.units.SpikesIdx(:);
if max(si) > nU, why = 'RIPPLES / spikes sort mismatch'; return; end
uMap = zeros(nU, 1); uMap(si) = 1:numel(si);
keep = keep & uMap > 0;
idx = find(keep);
if numel(idx) < o.minUnits, why = sprintf('< %d units', o.minUnits); return; end
preIv = local_preIv(e, o.preWin);
e.post = [max(e.post(1), e.off), min(e.post(2), e.off + o.postWin)];
inPre = find(local_ripInIv(RP, preIv));
inPost = find(RP.ts1 >= e.post(1) & RP.ts1Off <= e.post(2));
if ~o.matchN
    if min(numel(inPre), numel(inPost)) < o.minRip, why = sprintf('< %d ripples in pre or post', o.minRip); return; end
else
    N = min([numel(inPre), numel(inPost), o.maxRip]);
    if N < o.minRip, why = sprintf('< %d ripples in pre or post', o.minRip); return; end
    inPre = inPre(end-N+1:end); inPost = inPost(1:N);
end
Q = local_ripCounts(RP, inPre, uMap(idx)); P = local_ripCounts(RP, inPost, uMap(idx));
oth = find(~keep & uMap > 0 & sum(X, 1)' > 0);     % units in the ripple export outside the scored set
O.n = numel(oth);
if O.n > 0
    O.x = sum(X(:, oth), 2);
    O.q = sum(local_ripCounts(RP, inPre, uMap(oth)), 2); O.p = sum(local_ripCounts(RP, inPost, uMap(oth)), 2);
end
X = X(:, idx);
pairMask = triu(true(numel(idx)), 1);
if any(isfinite(grp(idx))), pairMask = pairMask & (grp(idx) ~= grp(idx)'); end
end


function [md, p, n] = local_paired(T, nm, a, b)
% mean over animals with both arms of (b - a), exact sign-flip
A = T(strcmp(T.Group, a), :); Bt = T(strcmp(T.Group, b), :);
ids = intersect(unique(A.Animal), unique(Bt.Animal)); n = numel(ids); md = NaN; p = NaN;
if n < 2, return; end
d = cellfun(@(x) local_nmean(Bt.(nm)(strcmp(Bt.Animal, x))) - local_nmean(A.(nm)(strcmp(A.Animal, x))), ids);
d = d(isfinite(d)); n = numel(d); md = NaN; p = NaN;
if n < 2, return; end
md = mean(d); p = local_signflip(d);
end


%% ---- copied unchanged from compareReactivation.m ----


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


function [v, ids] = local_animalMean(T, nm)
ids = unique(T.Animal); v = cellfun(@(a) local_nmean(T.(nm)(strcmp(T.Animal, a))), ids);
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


function p = local_perm(x, y, nperm)
% NaN-safe: a session that could not be scored under one correction is dropped, not propagated
p = NaN; x = x(isfinite(x)); y = y(isfinite(y));
if numel(x) >= 2 && numel(y) >= 2, p = permtest(x, y, nperm); end
end


function m = local_nmean(x)
% mean over finite values (NaN when there are none)
x = x(isfinite(x)); m = NaN; if ~isempty(x), m = mean(x); end
end


function e = local_nsem(x)
% standard error over finite values (NaN for fewer than 2)
x = x(isfinite(x)); e = NaN; if numel(x) > 1, e = std(x) / sqrt(numel(x)); end
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
