%% couplingPlasticity.m
% Does experience change each neuron's coupling to the common ripple drive?
%
% WHY. The per-cell decomposition (compareReactivation.m) shows that cells active in ConB1 raise
% their MEAN in-ripple firing afterwards (b_B > 0, spared by ASV). Population conditioning
% (reactivationCoupling.m) removes each cell's slope on the ripple's population count, and in ASV
% sessions the reactivation lost that way tracked b_B across sessions. That association does not
% show that the experience changed the slope itself: cells already strongly coupled before the
% exposure, because they are more excitable, could also be the ones most active during it. This
% script tests the slope directly, cell by cell, before and after the novel exposure.
%
% COUPLING. For unit j in a block of ripples i (PRE or POST), a Poisson GLM with log link
%     log E[x_ij] = alpha_j + gamma_j * log(s_ij + 1)
% x_ij = the unit's spikes in ripple i, s_ij = the population count of that ripple without unit j:
%     'loo'  every other unit in the ripple export
%     'int'  interneurons only (cannot contain the pyramidal assembly being scored)
% gamma_j is an elasticity: rescaling s (a population-wide rate change) moves alpha_j, not gamma_j,
% and a change in the unit's own mean rate is absorbed by alpha_j.
%
% MODELS across the units of one session, dgamma_j = gamma_j(POST) - gamma_j(PRE):
%     M1  dgamma_j = c + d_B z_j^B + d_A z_j^A                       (the requested test)
%     M2  ... + e * dl_j          dl_j = change in log10 in-ripple rate (rate control)
%     M3  ... + e * dl_j + f * gamma_j(PRE)                          (baseline / regression to mean)
% z^X_j = within-session z-score of log10(rate in exposure X + 0.01), running or theta part,
% exactly as in the decomposition. Also reported: Spearman corr(z^B, gamma(PRE)) - are cells
% active in ConB1 ALREADY more coupled before it?
%
% CONTRASTS
%     real     PRE = ripples before ConB1 (whole gap after ConA1 minus the infusion window),
%              POST = first POST_WIN after ConB1
%     placebo  PRE = hour before ConA1, POST = ripples before ConB1. Precedes the novel exposure;
%              in infused sessions it contains the infusion (descriptive), in controls it is clean.
%
% READING (b_B > 0 is established)
%     d_B > 0  experience-active neurons become preferentially coupled to the common ripple drive
%     d_B ~ 0  experience-associated firing rises without a measurable change in coupling
%     d_B < 0  experience-associated firing persists despite weaker coupling
%  A d_B that survives M2 is not a rate-sampling artefact; one present in the placebo too is
%  drift, not experience. corr(z^B, gamma_PRE) > 0 with d_B ~ 0 is the pre-existing-excitability
%  account of the session-level association.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT       = 'R:\TransPlasticity';
THETA_RATIO = 2; SPEED_MIN = 3; CELLS_TYPE = 'Pyramidal';
EXCL_INF = false; INF_EXCL_WIN = [0 1800]; PRE_B_ALLGAP = true; PRE_WIN = 3600; POST_WIN = 1800;
MIN_RIP = 50;             % ripples per block
MIN_RIP_SPK = 30;         % spikes a unit needs in the ripples of each block
MIN_CELLS = 8;            % units per session fit
POPS = {'loo', 'int'};
SESSION_SET = 'all';          % 'all' = every SAL and every ASV session | 'first' = SAL before the first ASV + the first ASV
if exist('SHARED', 'var')     % set by runAllAnalyses.m
    CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root;
    if isfield(SHARED, 'sessionSet'), SESSION_SET = SHARED.sessionSet; end
end
o = struct('thetaRatio', THETA_RATIO, 'speedMin', SPEED_MIN, 'preWin', PRE_WIN, 'postWin', POST_WIN, ...
    'minRip', MIN_RIP, 'minRipSpk', MIN_RIP_SPK, 'minCells', MIN_CELLS, 'cellType', CELLS_TYPE);
o.pops = POPS;
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

rows = {}; skipped = {};
%% ---- controls ----
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), skipped(end+1,:) = {rf, 'missing ripple export'}; continue; end %#ok<AGROW>
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES;
    nov = CONTROLS{i,2}; fam = CONTROLS{i,3}; rec = RP.meta.recordingSec;
    if fam(1) > nov(1), skipped(end+1,:) = {rf, 'novel before familiar'}; continue; end %#ok<AGROW>
    eA = struct('name', 'ConA1', 'on', fam(1), 'off', fam(2), 'pre', [0 fam(1)], 'post', [fam(2) nov(1)], 'preIv', []);
    eB = struct('name', 'ConB1', 'on', nov(1), 'off', nov(2), 'pre', [fam(2) nov(1)], 'post', [nov(2) rec], 'preIv', [fam(2) nov(1)]);
    [R_, why] = local_coupSession(RP, fld, eA, eB, o);
    if ~isempty(why), skipped(end+1,:) = {rf, why}; continue; end %#ok<AGROW>
    for r = 1:numel(R_), rows(end+1,:) = {'Control', controlAnimal(rf), rf, R_(r)}; end %#ok<AGROW>
end
%% ---- TransPlasticity ----
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
for i = 1:numel(fa)
    fld = fileparts(fa{i});
    ev = local_events(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConA1_start')), continue; end
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES;
    g = local_group(RP, ev, SESSION_SET);
    if ~ismember(g, {'SAL', 'ASV'}), continue; end
    HC = local_hcBlocks(ev, RP, EXCL_INF, INF_EXCL_WIN);
    E = struct('name', {}, 'on', {}, 'off', {}, 'pre', {}, 'post', {}, 'preIv', {});
    for c = {'ConA1', 'ConB1'}
        on = ev.time(strcmp(ev.description, [c{1} '_start'])); off = ev.time(strcmp(ev.description, [c{1} '_end']));
        if isempty(on) || isempty(off), break; end
        e = struct('name', c{1}, 'on', on(1), 'off', off(1), 'pre', [], 'post', [], 'preIv', []);
        kp = find(abs(HC(:,2) - e.on) <= 1, 1); kq = find(abs(HC(:,1) - e.off) <= 1, 1);
        if PRE_B_ALLGAP && strcmp(c{1}, 'ConB1')
            e.preIv = local_gapIv(HC, ev, e.on); kp = [];
            if ~isempty(e.preIv), kp = 0; end
        end
        if isempty(kp) || isempty(kq), break; end
        if kp > 0, e.pre = HC(kp,:); else, e.pre = [e.preIv(1,1) e.preIv(end,2)]; end
        e.post = HC(kq,:);
        E(end+1) = e; %#ok<AGROW>
    end
    if numel(E) < 2, skipped(end+1,:) = {fa{i}, 'ConA1 or ConB1 without home cage around it'}; continue; end %#ok<AGROW>
    [R_, why] = local_coupSession(RP, fld, E(1), E(2), o);
    if ~isempty(why), skipped(end+1,:) = {fa{i}, why}; continue; end %#ok<AGROW>
    for r = 1:numel(R_), rows(end+1,:) = {g, RP.meta.animal, fa{i}, R_(r)}; end %#ok<AGROW>
end
fprintf('%d session fits; %d sessions skipped\n', size(rows, 1), size(skipped, 1));
for k = 1:size(skipped, 1), fprintf('   %s: %s\n', skipped{k,:}); end

%% ---- table ----
S = [rows{:,4}]';
T = cell2table(rows(:,1:3), 'VariableNames', {'Group','Animal','Session'});
for f = fieldnames(S)', T.(f{1}) = {S.(f{1})}'; end
for f = fieldnames(S)'
    if all(cellfun(@(x) isnumeric(x) && isscalar(x), T.(f{1}))), T.(f{1}) = cell2mat(T.(f{1})); end
end
T = T(T.nCells >= MIN_CELLS, :);
save(fullfile(CACHE_ROOT, 'couplingPlasticity.mat'), 'T', 'o');

GR = {'Control', 'SAL', 'ASV'};
fprintf('\n===== experience-dependent change in coupling to the common ripple drive =====\n');
fprintf('per session: dgamma_j = c + d_B z_B + d_A z_A [+ e dl_j] [+ f gamma_pre]; mean +- SEM over sessions,\n');
fprintf('p = exact sign-flip over animal means; ASV - SAL paired within animal\n');
for pp = POPS
    for cn = {'real', 'placebo'}
        D = T(strcmp(T.pop, pp{1}) & strcmp(T.contrast, cn{1}), :);
        if isempty(D), continue; end
        fprintf('\n-- population %s | %s --\n', pp{1}, cn{1});
        for g = GR
            k = strcmp(D.Group, g{1}); if ~any(k), continue; end
            Dg = D(k,:);
            fprintf('   %-8s n = %d/%d, median %d cells | gamma PRE %.2f POST %.2f | c %+.3f\n', g{1}, sum(k), ...
                numel(unique(Dg.Animal)), round(median(Dg.nCells)), mean(Dg.gPre), mean(Dg.gPost), mean(Dg.c1));
            for mm = {'dB1', 'M1'; 'dB2', 'M2 (+rate)'; 'dB3', 'M3 (+rate, +baseline)'}'
                x = Dg.(mm{1}); am = local_animalMean(Dg, mm{1});
                fprintf('      %-22s d_B %+.3f +- %.3f  (p = %.3f)\n', mm{2}, mean(x), std(x) / sqrt(numel(x)), local_signflip(am));
            end
            am = local_animalMean(Dg, 'dA1');
            fprintf('      %-22s d_A %+.3f +- %.3f  (p = %.3f)\n', 'M1', mean(Dg.dA1), std(Dg.dA1) / sqrt(height(Dg)), local_signflip(am));
            am = local_animalMean(Dg, 'rhoPre');
            fprintf('      corr(z_B, gamma PRE) rho %+.2f (p = %.3f)   corr(z_B, dl) rho %+.2f\n', mean(Dg.rhoPre, 'omitnan'), ...
                local_signflip(am), mean(Dg.rhoRate, 'omitnan'));
        end
        for mm = {'dB1', 'dB2', 'dB3'}
            [md, pv, na] = local_paired(D, mm{1}, 'SAL', 'ASV');
            fprintf('   ASV - SAL %s: paired %+.3f p = %.3f (%d animals) | LME p = %.3f\n', mm{1}, md, pv, na, local_lme(D, mm{1}, 'SAL', 'ASV'));
        end
    end
    % real minus placebo, same session and cells' z (sign-flip over sessions)
    Rr = T(strcmp(T.pop, pp{1}) & strcmp(T.contrast, 'real'), {'Group', 'Session', 'dB1', 'dB2'});
    Rp = T(strcmp(T.pop, pp{1}) & strcmp(T.contrast, 'placebo'), {'Session', 'dB1', 'dB2'});
    Rr.Properties.VariableNames(3:4) = {'dB1_Rr', 'dB2_Rr'}; Rp.Properties.VariableNames(2:3) = {'dB1_Rp', 'dB2_Rp'};
    J = innerjoin(Rr, Rp, 'Keys', 'Session');
    fprintf('\n   real - placebo d_B (population %s)\n', pp{1});
    for g = GR
        k = strcmp(J.Group, g{1}); if sum(k) < 2, continue; end
        d1 = J.dB1_Rr(k) - J.dB1_Rp(k); d2 = J.dB2_Rr(k) - J.dB2_Rp(k);
        fprintf('      %-8s M1 %+.3f (p = %.3f) | M2 %+.3f (p = %.3f) | %d sessions\n', g{1}, mean(d1), local_signflip(d1), ...
            mean(d2), local_signflip(d2), sum(k));
    end
end

% link to the session-level result: d_B (real, int) against EV - REV lost to conditioning
cf = fullfile(CACHE_ROOT, 'reactivationCoupling.mat');
if exist(cf, 'file')
    C = load(cf, 'T'); RC = C.T(strcmp(C.T.Context, 'ConB1'), {'Session', 'dEV_raw', 'dEV_condOther'});
    RC.loss = RC.dEV_raw - RC.dEV_condOther;
    for pp = POPS
        D = T(strcmp(T.pop, pp{1}) & strcmp(T.contrast, 'real'), {'Group', 'Session', 'dB1'});
        J = innerjoin(D, RC, 'Keys', 'Session');
        for g = GR
            % NaN-safe: a session whose outside-unit conditioning could not be scored is dropped
            k = strcmp(J.Group, g{1}) & isfinite(J.dB1) & isfinite(J.loss);
            if sum(k) < 4, continue; end
            [rr, pv] = corr(J.dB1(k), J.loss(k), 'Type', 'Spearman');
            fprintf('   d_B (%s) vs EV - REV lost to outside-unit conditioning, %-8s rho %+.2f (p = %.3f, n = %d)\n', ...
                pp{1}, g{1}, rr, pv, sum(k));
        end
    end
end


%% ======================================================================================
function [R, why] = local_coupSession(RP, fld, eA, eB, o)
% per-unit coupling gamma in PRE and POST blocks, regressed on z_B, z_A; one struct per (pop, contrast)
R = struct([]); why = '';
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
si = RP.units.SpikesIdx(:);
if max(si) > nU, why = 'RIPPLES / spikes sort mismatch'; return; end
typRP = typ(si); isPyr = contains(typRP, o.cellType); isInt = contains(typRP, 'Interneuron');
% exposure rates (running, else theta part), per RP unit
lam = nan(numel(si), 2); EX = {eA, eB};
for x = 1:2
    [mask, src] = local_runMask(fld, EX{x}.on, EX{x}.off, o.speedMin);
    if isempty(src), mask = local_thetaMask(RP, EX{x}.on, EX{x}.off, o.thetaRatio); end
    if isempty(mask), why = sprintf('no running / theta in %s', EX{x}.name); return; end
    dur = sum(diff(mask, 1, 2));
    for u = 1:numel(si)
        t = ts{si(u)}; n_ = 0;
        for w = 1:size(mask, 1), n_ = n_ + sum(t >= mask(w,1) & t < mask(w,2)); end
        lam(u, x) = n_ / dur;
    end
end
% ripple blocks
blk.preA = find(local_ripInIv(RP, local_preIv(eA, o.preWin)));
blk.preB = find(local_ripInIv(RP, local_preIv(eB, o.preWin)));
pB = [max(eB.post(1), eB.off), min(eB.post(2), eB.off + o.postWin)];
blk.postB = find(RP.ts1 >= pB(1) & RP.ts1Off <= pB(2));
CON = {'real', 'preB', 'postB'; 'placebo', 'preA', 'preB'};
for q = 1:size(CON, 1)
    b0 = blk.(CON{q,2}); b1 = blk.(CON{q,3});
    if numel(b0) < o.minRip || numel(b1) < o.minRip, continue; end
    N0 = local_ripCounts(RP, b0, 1:numel(si)); N1 = local_ripCounts(RP, b1, 1:numel(si));
    d0 = sum(RP.ts1Off(b0) - RP.ts1(b0)); d1 = sum(RP.ts1Off(b1) - RP.ts1(b1));
    cells = find(isPyr & sum(N0, 1)' >= o.minRipSpk & sum(N1, 1)' >= o.minRipSpk & all(isfinite(lam), 2));
    if numel(cells) < o.minCells, continue; end
    zA = zscore(log10(lam(cells, 1) + 0.01)); zB = zscore(log10(lam(cells, 2) + 0.01));
    dl = log10(sum(N1(:, cells), 1)' / d1 + 0.1) - log10(sum(N0(:, cells), 1)' / d0 + 0.1);
    for pp = o.pops
        g0 = local_gamma(N0, cells, isInt, pp{1}); g1 = local_gamma(N1, cells, isInt, pp{1});
        ok = isfinite(g0) & isfinite(g1);
        if sum(ok) < o.minCells, continue; end
        dg = g1(ok) - g0(ok);
        X1 = [ones(sum(ok), 1), zB(ok), zA(ok)];
        X2 = [X1, dl(ok)]; X3 = [X2, g0(ok)];
        b1_ = X1 \ dg; b2_ = X2 \ dg; b3_ = X3 \ dg;
        r = struct('pop', pp{1}, 'contrast', CON{q,1}, 'nCells', sum(ok), 'gPre', mean(g0(ok)), 'gPost', mean(g1(ok)), ...
            'c1', b1_(1), 'dB1', b1_(2), 'dA1', b1_(3), 'dB2', b2_(2), 'eRate', b2_(4), 'dB3', b3_(2), ...
            'rhoPre', corr(zB(ok), g0(ok), 'Type', 'Spearman'), 'rhoRate', corr(zB(ok), dl(ok), 'Type', 'Spearman'));
        if isempty(R), R = r; else, R(end+1) = r; end %#ok<AGROW>
    end
end
if isempty(R), why = sprintf('< %d units with %d ripple spikes per block, or < %d ripples', o.minCells, o.minRipSpk, o.minRip); end
end


function g = local_gamma(N, cells, isInt, pop)
% Poisson GLM slope of each unit's per-ripple count on log(1 + population count without it)
g = nan(numel(cells), 1);
tot = sum(N, 2); totInt = sum(N(:, isInt), 2);
ws = warning('off', 'all'); cl = onCleanup(@() warning(ws));
for c = 1:numel(cells)
    j = cells(c); y = N(:, j);
    switch pop
        case 'loo', s = tot - y;
        case 'int', s = totInt - y .* isInt(j);
    end
    xv = log(s + 1);
    if std(xv) == 0 || sum(y) == 0, continue; end
    try
        b = glmfit(xv, y, 'poisson');
        if all(isfinite(b)) && abs(b(2)) < 20, g(c) = b(2); end
    catch
    end
end
end


%% ---- copied unchanged from reactivationCoupling.m / compareReactivation.m ----
function [md, p, n] = local_paired(T, nm, a, b)
% mean over animals with both arms of (b - a), exact sign-flip
A = T(strcmp(T.Group, a), :); Bt = T(strcmp(T.Group, b), :);
ids = intersect(unique(A.Animal), unique(Bt.Animal)); n = numel(ids); md = NaN; p = NaN;
if n < 2, return; end
d = cellfun(@(x) mean(Bt.(nm)(strcmp(Bt.Animal, x))) - mean(A.(nm)(strcmp(A.Animal, x))), ids);
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
