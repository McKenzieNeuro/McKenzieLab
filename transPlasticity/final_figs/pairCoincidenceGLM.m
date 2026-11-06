%% pairCoincidenceGLM.m
% How long does exploration covariance persist in later ripples? A ripple-level pair GLM.
%
% WHY THIS AND NOT THE SUBSPACE SCORE
% The subspace score (driftInclusionASV) needs time bins and a fixed cell count; at 9 cells it
% could not reproduce its own full-yield answer in the controls. Here every ripple keeps its own
% lag, every scored cell is used, and the null already contains each neuron's rate and its
% coupling to ripple size, so single-neuron changes (gain b_B, coupling d_B) cannot produce the
% effect.
%
% MODEL (per session and exposure; ripples i, scored pyramidal pairs (j,k))
%   1. single-neuron offset, Poisson, fitted on all ripples of the fit range:
%        log mu_ij = log d_i + alpha_j(block of t_i) + gamma_j^{pre|post} log(s_i^{-j} + 1)
%      alpha_j: one level per KNOT_H block on each side (ridge-shrunk), gamma_j: coupling to the
%      ripple's population count, separately before and after the exposure. For the pair the
%      drive excludes both neurons: s_i^{-jk}.
%   2. coincidences y_ijk = x_ij x_ik; under conditional independence E[y] = mu_ij mu_ik.
%   3. log E[y_ijk] = log(mu_ij mu_ik) + u_jk + c_jk [ delta |D_i| + 1(D_i > 0) A exp(-D_i/tau) ]
%      c_jk  = pair partial correlation during exploration (100 ms bins, running/theta part,
%              population count partialled out), z-scored across pairs
%      u_jk  = pair fixed effect: all time-invariant pair structure, profiled out exactly
%              (conditional likelihood of each pair's coincidences given its total)
%      delta = symmetric drift of the alignment with c, the same on both sides
%      A     = extra co-firing of exploration-correlated pairs just after the exposure,
%              log(observed/expected) per SD of c;  tau = how long it lasts
%      D_i   = h after the exposure ends (post), minus h before it starts (pre)
%   tau is profiled on a log grid; at fixed tau the model is a small Newton problem. Sessions
%   share tau within an arm and have their own A and delta, so the arm profile is the sum of the
%   session curves, and the animal bootstrap resamples whole curves.
%
% CONTROLS IN THE SAME RUN
%   HC    c from a duration-matched home-cage block outside the fit range in place of the
%         exploration: if it gives the same A, the effect is state (post-exposure sleep), not
%         experience.
%   joint both c in one model: A_E with the home-cage alignment held fixed.
%   yield never-infused sessions refitted on random subsets of YIELD_SUB scored cells: the
%         subset tau should fall inside the full-yield CI before the infused arms are believed.
%
% LAST TIME WITH EVIDENCE (tau-free)
%   tau is not recoverable at infused yields, so a second fit replaces A exp(-D/tau) by one free
%   excess per post bin (LAST_BIN_H wide), relative to the pre-exposure ripples:
%        ... + c_jk [ delta |D_i| + sum_b A_b 1(D_i in bin b) ]
%   once with delta = 0 ('nodrift', primary) and once with delta free ('drift', identified from
%   the pre side only, so it is extrapolated across the post window). Per animal, bins are pooled
%   by their information. A GREEDY scan then trades resolution for sensitivity: starting at the
%   first bin, a window grows one bin at a time until its animal-level excess is > 0 at sign-flip
%   p < LAST_ALPHA; it is then closed and the next window starts after it. Windows are short where
%   the effect is large and long where it is small. Growing a window until it passes is optional
%   stopping, so each detected window's START a is recalibrated: the statistic is the maximum t
%   over every window [a, e], and its null comes from flipping each unit's sign jointly for all
%   bins >= a (bins before a untouched). The latest a with calibrated p < 0.05 is the last time
%   from which the pairwise signature can be shown to persist (H0: no excess from a onward).
%
% LIMITS, as in driftInclusionASV: ConA1 POST ends at the infusion in infused sessions; ConB1
% PRE is home cage after ConA1 (and drug-free while POST is drugged in ASV). Model SEs are not
% used anywhere: pairs share neurons and ripples share pairs, so all inference is by animal.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT       = 'R:\TransPlasticity';
SESSION_SET = 'all';
if exist('SHARED', 'var')
    CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root;
    if isfield(SHARED, 'sessionSet'), SESSION_SET = SHARED.sessionSet; end
end
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

O = struct();
O.bin        = 0.1;          % s, exploration bins (as the EV template)
O.speedMin   = 3;            % cm/s, running windows
O.thetaRatio = 2;            % theta/delta, when there is no tracking
O.winH       = 4;            % h of home cage used on each side of an exposure
O.infExcl    = [0 1800];     % s [before after] Infusion_start removed from the fit ripples
O.knotH      = 0.5;          % h, blocks for the neuron's rate drift alpha_j
O.ridge      = 1;            % ridge on the block levels (log-rate units^-2)
O.minSpkCtx  = 20;           % spikes in the exploration bins to score a unit
O.minRipSpk  = 10;           % ripple spikes on EACH side to score a unit
O.minUnits   = 5;            % scored units (10 pairs)
O.minRipSide = 50;           % ripples on each side
O.maxPairs   = 1500;         % random cap on pairs (full-yield controls), fixed seed
O.nremOnly   = false;        % fit ripples restricted to NREM where a SleepState file exists
O.tauGrid    = exp(linspace(log(0.25), log(48), 40));   % h. Below ~0.25 h the decay term touches
                             % only the first few ripples and A becomes unconstrained; in simulation
                             % a null session then maximises at the floor, so the floor is not a result
O.tauRef     = 1;            % h, a fixed tau at which A is also tested (model-light amplitude)
O.hcCtrl     = true;
O.lastBinH   = 0.25;         % h, post bins for the tau-free scan
LAST_ALPHA   = 0.05;         % per-window threshold inside the greedy scan (recalibrated afterwards)
LAST_MIN_UNITS = 5;          % units (animals; sessions for never-infused) a window needs
N_FLIP_MAX   = 4096;         % exact flips up to 12 units, else this many random
YIELD_SUB    = [9 20];       % scored cells per draw for the never-infused yield check
N_YIELD_DRAW = 10;
N_BOOT       = 2000;
EQ_RATIO     = 2;            % TOST bound for tau_ASV / tau_SAL
SAVE_OUT     = true;

rng(1);
F = {};        % one struct per session x exposure x cell set
led = cell(0, 3);

%% ---- TransPlasticity ------------------------------------------------------------------------
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
for i = 1:numel(fa)
    fld = fileparts(fa{i}); [~, snL] = fileparts(fld);
    ev = local_events(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConA1_start'))
        led(end+1,:) = {snL, '', 'no event file or no ConA1'}; continue %#ok<AGROW>
    end
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    g = local_group(RP, ev, SESSION_SET);
    if isempty(g), led(end+1,:) = {snL, '', 'not SAL or ASV in this session set'}; continue; end %#ok<AGROW>
    HC = local_hcBlocks(ev, RP, false, [0 0]);
    tOf = @(lab) ev.time(find(strcmp(ev.description, lab), 1));
    infOn = tOf('Infusion_start'); if isempty(infOn), infOn = NaN; end
    allCtx = zeros(0, 2); IV = struct('name', {}, 'iv', {});
    for c = {'ConA1','ConB1','ConA2','ConB2'}
        on = tOf([c{1} '_start']); off = tOf([c{1} '_end']);
        if isempty(on) || isempty(off), continue; end
        allCtx(end+1,:) = [on off]; %#ok<AGROW>
        if any(strcmp(c{1}, {'ConA1','ConB1'})), IV(end+1) = struct('name', c{1}, 'iv', [on off]); end %#ok<AGROW>
    end
    S = local_loadSession(fld, RP);
    if isempty(S), led(end+1,:) = {snL, g, 'no spikes / cell_metrics / PYR units'}; continue; end %#ok<AGROW>
    for w = 1:numel(IV)
        [f, why] = local_fitExposure(S, fld, RP, IV(w), allCtx, HC, infOn, O, []);
        if isempty(f), led(end+1,:) = {snL, g, [IV(w).name ': ' why]}; continue; end %#ok<AGROW>
        f.Group = g; f.Animal = char(string(m.animal)); f.Session = snL; f.Context = IV(w).name; f.CellSet = 'all';
        F{end+1} = local_pack(f); %#ok<AGROW>
        fprintf('%-30s %-4s %s: %d units, %d pairs, %d/%d ripples pre/post\n', snL, g, IV(w).name, f.nUnits, f.nPairs, f.nRipPre, f.nRipPost);
    end
end

%% ---- never-infused ---------------------------------------------------------------------------
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1}; nov = CONTROLS{i,2}; fam = CONTROLS{i,3};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), led(end+1,:) = {local_name(fld), 'Control', 'no rippleEvents_control.mat'}; continue; end %#ok<AGROW>
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES;
    IV = struct('name', {'ConA1','ConB1'}, 'iv', {fam, nov});
    allCtx = sortrows([fam; nov]);
    rec = max([RP.ts1Off(:); 1]) + 1;
    HC = local_complement(allCtx, [0 rec]);
    S = local_loadSession(fld, RP);
    if isempty(S), continue; end
    sn = local_name(fld); mouse = regexprep(sn, '_\d{6}$', '');
    for w = 1:numel(IV)
        [f, why] = local_fitExposure(S, fld, RP, IV(w), allCtx, HC, NaN, O, []);
        if isempty(f), led(end+1,:) = {sn, 'Control', [IV(w).name ': ' why]}; continue; end %#ok<AGROW>
        f.Group = 'Control'; f.Animal = mouse; f.Session = sn; f.Context = IV(w).name; f.CellSet = 'all';
        F{end+1} = local_pack(f); %#ok<AGROW>
        fprintf('%-30s Ctrl %s: %d units, %d pairs, %d/%d ripples pre/post\n', sn, IV(w).name, f.nUnits, f.nPairs, f.nRipPre, f.nRipPost);
        % yield check: the same exposure on random subsets of the scored cells, curves averaged
        for ny = YIELD_SUB
            if f.nUnits <= ny, continue; end
            acc = [];
            for d = 1:N_YIELD_DRAW
                sub = sort(randperm(f.nUnits, ny));
                fy = local_fitExposure(S, fld, RP, IV(w), allCtx, HC, NaN, setfield(O, 'hcCtrl', false), f.unitIdx(sub)); %#ok<SFLD>
                if isempty(fy), continue; end
                if isempty(acc), acc = fy; acc.nDraw = 1;
                else
                    acc.E.L = acc.E.L + fy.E.L; acc.E.th = acc.E.th + fy.E.th; acc.nDraw = acc.nDraw + 1;
                end
            end
            if isempty(acc), continue; end
            acc.E.L = acc.E.L / acc.nDraw; acc.E.th = acc.E.th / acc.nDraw;
            acc.Group = sprintf('Ctrl%d', ny); acc.Animal = mouse; acc.Session = sn; acc.Context = IV(w).name;
            acc.CellSet = sprintf('%d cells x %d draws', ny, acc.nDraw);
            F{end+1} = local_pack(acc); %#ok<AGROW>
        end
    end
end
if isempty(F), error('pcg:none', 'no exposure could be fitted'); end
F = [F{:}];

%% ---- report ----------------------------------------------------------------------------------
R = local_report(F, O, N_BOOT, EQ_RATIO, YIELD_SUB);
R.last = local_lastTime(F, O, LAST_ALPHA, LAST_MIN_UNITS, N_FLIP_MAX);
R.did = local_did(F, O, R);
if ~isempty(led)
    fprintf('\nnot fitted:\n');
    for q = 1:size(led, 1), fprintf('   %-30s %-8s %s\n', led{q,:}); end
end
if SAVE_OUT
    save(fullfile(CACHE_ROOT, 'pairCoincidenceGLM.mat'), 'F', 'R', 'O');
    fprintf('\nsaved pairCoincidenceGLM.mat\n');
end


%% =============================================================================================
%  per exposure
%% =============================================================================================
function [f, why] = local_fitExposure(S, fld, RP, E, allCtx, HC, infOn, O, unitSel)
% unitSel: [] = every eligible unit; else the indices (into S units) to score
f = []; why = '';
on = E.iv(1); off = E.iv(2);
% ---- fit ripples: home cage on each side, stopping at the neighbouring contexts / infusion
tA = on - O.winH * 3600; tB = off + O.winH * 3600;
pv = allCtx(allCtx(:,2) <= on, 2); if ~isempty(pv), tA = max(tA, max(pv)); end
nx = allCtx(allCtx(:,1) >= off, 1); if ~isempty(nx), tB = min(tB, min(nx)); end
if isfinite(infOn) && on < infOn, tB = min(tB, infOn); end        % a pre-infusion exposure stops at it
t = S.tOn;
inHC = local_inIv(t, HC);
keep = inHC & ((t >= tA & t < on) | (t >= off & t < tB));
if isfinite(infOn), keep = keep & ~(t >= infOn - O.infExcl(1) & t < infOn + O.infExcl(2)); end
if O.nremOnly && ~isempty(S.nrem), keep = keep & local_inIv(t, S.nrem); end
ri = find(keep);
D = nan(size(t));
D(t >= off) = (t(t >= off) - off) / 3600; D(t < on) = -(on - t(t < on)) / 3600;
D = D(ri); post = D > 0;
if sum(~post) < O.minRipSide || sum(post) < O.minRipSide
    why = sprintf('%d/%d ripples pre/post', sum(~post), sum(post)); return
end
X = S.X(ri, :); dur = S.dur(ri);
% ---- exploration template: running (or theta) part of the exposure, 100 ms bins
[mask, src] = local_runMask(fld, on, off, O.speedMin);
if isempty(src), mask = local_thetaMask(RP, on, off, O.thetaRatio); end
Z = local_binMask(S.ts, mask, O.bin);
if size(Z, 1) < 100, why = 'too few exploration bins'; return; end
okU = sum(Z, 1)' >= O.minSpkCtx & sum(X(~post,:), 1)' >= O.minRipSpk & sum(X(post,:), 1)' >= O.minRipSpk;
if ~isempty(unitSel), tmp = false(size(okU)); tmp(unitSel) = okU(unitSel); okU = tmp; end
U = find(okU);
if numel(U) < O.minUnits, why = sprintf('%d scored units', numel(U)); return; end
cE = local_partialCorr(Z, U);
% ---- home-cage template: duration-matched block outside the fit range
cH = [];
if O.hcCtrl
    blk = local_hcBlock(HC, [tA tB], off - on, infOn, O.infExcl);
    if ~isempty(blk)
        ZH = local_binMask(S.ts, blk, O.bin);
        if size(ZH, 1) >= 100, cH = local_partialCorr(ZH, U); end
    end
end
% ---- pairs
[jj, kk] = find(triu(true(numel(U)), 1));
P = [jj kk];
okP = isfinite(cE(sub2ind(size(cE), jj, kk)));
P = P(okP, :);
if size(P, 1) > O.maxPairs, P = P(sort(randperm(size(P, 1), O.maxPairs)), :); end
cEp = cE(sub2ind(size(cE), P(:,1), P(:,2)));
cEp = (cEp - mean(cEp)) / max(std(cEp), eps);
% ---- single-neuron offsets (every scored unit, drive = all recorded PYR units)
s = sum(X, 2);
[aIdx, Gam] = local_neuronOffsets(X(:, U), s, dur, D, post, O);
% ---- pair matrices
xj = X(:, U(P(:,1))); xk = X(:, U(P(:,2)));
Y = xj .* xk;
sJK = max(s - xj - xk, 0);
Oo = aIdx(:, P(:,1)) + Gam(:, P(:,1)) .* log(sJK + 1) + aIdx(:, P(:,2)) + Gam(:, P(:,2)) .* log(sJK + 1);
Yt = sum(Y, 1);
kp = Yt > 0;                                % pairs with no coincidence carry no information
Y = Y(:, kp); Oo = Oo(:, kp); cEp = cEp(kp);
if size(Y, 2) < 3, why = 'fewer than 3 pairs with any coincidence'; return; end
absD = abs(D);
f = struct('nUnits', numel(U), 'nPairs', size(Y, 2), 'nRipPre', sum(~post), 'nRipPost', sum(post), ...
    'unitIdx', U, 'tmplSrc', src, 'postRange', [min(D(post)) max(D(post))], 'preRange', [min(-D(~post)) max(-D(~post))], ...
    'nCoinc', sum(Y(:)));
f.E = local_profile(Y, Oo, cEp, [], absD, post, D, O.tauGrid);
f.H = []; f.J = []; f.rEH = NaN;
if ~isempty(cH)
    cHp = cH(sub2ind(size(cH), P(:,1), P(:,2))); cHp = cHp(kp);
    okH = isfinite(cHp);
    if sum(okH) >= 3
        cHp(~okH) = 0;
        cHp = (cHp - mean(cHp(okH))) / max(std(cHp(okH)), eps); cHp(~okH) = 0;
        f.H = local_profile(Y, Oo, cHp, [], absD, post, D, O.tauGrid);
        f.J = local_profile(Y, Oo, cEp, cHp, absD, post, D, O.tauGrid);
        f.rEH = corr(cEp, cHp);
    end
end
% ---- tau-free: one excess per post bin, with and without the drift term
eb = 0:O.lastBinH:O.winH;
f.pwEdges = eb; f.pwA0 = []; f.pwI0 = []; f.pwA1 = []; f.pwI1 = [];
if isempty(unitSel)                          % not for the yield-check subsets
    cHuse = []; if ~isempty(f.J), cHuse = cHp; end
    [f.pwA0, f.pwI0] = local_piecewise(Y, Oo, cEp, cHuse, absD, post, D, eb, false);
    [f.pwA1, f.pwI1] = local_piecewise(Y, Oo, cEp, cHuse, absD, post, D, eb, true);
end
end


function [A, I] = local_piecewise(Y, Oo, c1, c2, absD, post, D, eb, withDrift)
% c_jk [ (delta |D|) + sum_b A_b 1(D in b) ] (+ home-cage nuisance: delta_H |D| + A_H 1(post));
% returns A_b and its model information (used only as a pooling weight, never for inference)
nb = numel(eb) - 1; A = nan(1, nb); I = zeros(1, nb);
bi = discretize(D, eb); bi(~post) = NaN;
Gb = zeros(numel(D), nb);
for b = 1:nb, Gb(:, b) = double(bi == b); end
has = any(Gb, 1) & (sum(Gb' * Y, 2)' > 0);
if ~any(has), return; end
G = Gb(:, has); Cp = repmat(c1, 1, sum(has)); off = 0;
if withDrift, G = [absD, G]; Cp = [c1, Cp]; off = 1; end
if ~isempty(c2), G = [G, absD, double(post)]; Cp = [Cp, c2, c2]; end
[th, ~, H] = local_condFit(Y, Oo, Cp, G, zeros(size(G, 2), 1));
k = off + (1:sum(has));
A(has) = th(k)'; I(has) = max(-diag(H(k, k))', 0);
end


function p = local_pack(f)
% one fixed field set, so exposures, subsets and arms concatenate into one struct array
fn = {'Group','Animal','Session','Context','CellSet','nUnits','nPairs','nRipPre','nRipPost', ...
      'postRange','preRange','nCoinc','tmplSrc','unitIdx','E','H','J','rEH', ...
      'pwA0','pwI0','pwA1','pwI1','pwEdges'};
p = struct();
for q = 1:numel(fn)
    if isfield(f, fn{q}), p.(fn{q}) = f.(fn{q}); else, p.(fn{q}) = []; end
end
end


function C = local_profile(Y, Oo, c1, c2, absD, post, D, grid)
% profile over tau; returns L (concentrated log-lik) and th (grid x M): [delta A] or
% [delta_E A_E delta_H A_H] for the joint model
if isempty(c2), Cp = [c1 c1]; else, Cp = [c1 c1 c2 c2]; end
M = size(Cp, 2);
L = nan(numel(grid), 1); th = nan(numel(grid), M);
t0 = zeros(M, 1);
for q = numel(grid):-1:1                       % from long tau down, warm-starting
    fdec = double(post) .* exp(-max(D, 0) / grid(q));
    if isempty(c2), G = [absD fdec]; else, G = [absD fdec absD fdec]; end
    [t0, L(q)] = local_condFit(Y, Oo, Cp, G, t0);
    th(q, :) = t0';
end
C = struct('L', L, 'th', th, 'M', M);
end


function [th, L, H] = local_condFit(Y, Oo, Cp, G, th)
% maximise sum_p [ sum_i y_ip eta_ip - Y_p log sum_i exp(o_ip + eta_ip) ],
% eta_ip = sum_m G(i,m) Cp(p,m) th(m): the Poisson likelihood with a free intercept per pair,
% profiled out exactly
M = size(Cp, 2); Yt = sum(Y, 1);
Gy = G' * Y;                                    % M x P
fLL = @(t) local_cll(Y, Oo, Cp, G, t, Yt);
[L, W] = fLL(th);
H = zeros(M);
nzG = G ~= 0;
for it = 1:60
    Gw = G' * W;                                % M x P
    grad = sum(Cp' .* (Gy - Yt .* Gw), 2);
    H = zeros(M);
    for a = 1:M
        for b = a:M
            if a == b || any(nzG(:,a) & nzG(:,b))
                cv = ((G(:,a) .* G(:,b))' * W) - Gw(a,:) .* Gw(b,:);
            else
                cv = -Gw(a,:) .* Gw(b,:);              % disjoint indicators: no product term
            end
            H(a,b) = -sum(Yt .* Cp(:,a)' .* Cp(:,b)' .* cv);
            H(b,a) = H(a,b);
        end
    end
    step = -(H - 1e-8 * eye(M)) \ grad;
    if any(~isfinite(step)), break; end
    lam = 1; improved = false;
    for ls = 1:30
        [L1, W1] = fLL(th + lam * step);
        if L1 >= L - 1e-10, improved = true; break; end
        lam = lam / 2;
    end
    if ~improved, break; end
    dth = lam * step; th = th + dth; L = L1; W = W1;
    if max(abs(dth)) < 1e-7, break; end
end
end


function [L, W] = local_cll(Y, Oo, Cp, G, th, Yt)
eta = G * (Cp .* th(:)')';                      % nRip x P
Zz = Oo + eta;
mx = max(Zz, [], 1);
Ex = exp(Zz - mx); Sx = sum(Ex, 1);
W = Ex ./ Sx;
L = sum(Y .* eta, 'all') - sum(Yt .* (log(Sx) + mx));
end


function [Aofs, Gam] = local_neuronOffsets(X, s, dur, D, post, O)
% per neuron: log mu = log d + alpha(block) + gamma_side log(s^{-j}+1); returns, per ripple,
% log d + alpha (Aofs) and gamma_side (Gam), both nRip x nUnits
blk = floor(D / O.knotH); [~, ~, bi] = unique(blk);
nb = max(bi);
B = full(sparse((1:numel(bi))', bi, 1, numel(bi), nb));
nR = size(X, 1); nU = size(X, 2);
Aofs = zeros(nR, nU); Gam = zeros(nR, nU);
for u = 1:nU
    y = X(:, u); su = log(max(s - y, 0) + 1);
    Xd = [ones(nR, 1), B, su .* ~post, su .* post];
    pen = [0, O.ridge * ones(1, nb), 0, 0];
    b = local_poissRidge(Xd, y, log(max(dur, eps)), pen);
    Aofs(:, u) = log(max(dur, eps)) + Xd(:, 1:nb+1) * b(1:nb+1);
    Gam(:, u) = b(nb+2) * ~post + b(nb+3) * post;
end
end


function b = local_poissRidge(Xd, y, off, pen)
b = zeros(size(Xd, 2), 1); b(1) = log(max(sum(y), 0.5) / sum(exp(off)));
Pm = diag(pen);
for it = 1:50
    eta = off + Xd * b; mu = exp(min(eta, 30));
    g = Xd' * (y - mu) - Pm * b;
    H = Xd' * (Xd .* mu) + Pm + 1e-8 * eye(numel(b));
    db = H \ g;
    b = b + db;
    if max(abs(db)) < 1e-8, break; end
end
end


function C = local_partialCorr(Z, U)
% pair correlation of binned counts after removing, from each unit, its regression on the
% summed count of the OTHER units in Z (all recorded PYR units)
s = sum(Z, 2); n = numel(U);
R = zeros(size(Z, 1), n);
for q = 1:n
    z = Z(:, U(q)); x = s - z;
    A = [ones(numel(z), 1), x];
    R(:, q) = z - A * (A \ z);
end
sd = std(R, 0, 1); C = nan(n);
ok = sd > 0;
if sum(ok) >= 2, C(ok, ok) = corr(R(:, ok)); end
end


function Z = local_binMask(ts, mask, binw)
Z = zeros(0, numel(ts));
if isempty(mask), return; end
edges = [];
for k = 1:size(mask, 1), edges = [edges; (mask(k,1):binw:mask(k,2)-binw)']; end %#ok<AGROW>
if numel(edges) < 2, return; end
Z = zeros(numel(edges), numel(ts));
for u = 1:numel(ts)
    tt = ts{u}; tt = tt(tt >= edges(1) & tt < edges(end) + binw);
    if isempty(tt), continue; end
    [~, ~, b] = histcounts(tt, [edges; edges(end) + binw]);
    b = b(:); v = b > 0; v(v) = tt(v) < edges(b(v)) + binw;
    if any(v), Z(:, u) = accumarray(b(v), 1, [numel(edges) 1]); end
end
end


function blk = local_hcBlock(HC, fitRange, dur, infOn, infExcl)
% one contiguous home-cage stretch of length dur, outside the fit range and the infusion
% window; the one farthest from the fit range is taken
blk = [];
H = HC;
H = local_subtract(H, fitRange);
if isfinite(infOn), H = local_subtract(H, [infOn - infExcl(1), infOn + infExcl(2)]); end
H = H(diff(H, 1, 2) >= dur, :);
if isempty(H), return; end
cand = [H(:,1), H(:,1) + dur; H(:,2) - dur, H(:,2)];
dist = min(abs(cand(:,2) - fitRange(1)), abs(cand(:,1) - fitRange(2)));
[~, q] = max(dist);
blk = cand(q, :);
end


function H = local_subtract(H, iv)
out = zeros(0, 2);
for q = 1:size(H, 1)
    a = H(q,1); b = H(q,2);
    if b <= iv(1) || a >= iv(2), out(end+1,:) = [a b]; continue; end %#ok<AGROW>
    if a < iv(1), out(end+1,:) = [a iv(1)]; end %#ok<AGROW>
    if b > iv(2), out(end+1,:) = [iv(2) b]; end %#ok<AGROW>
end
H = out;
end


function H = local_complement(C, rng_)
H = rng_;
for q = 1:size(C, 1), H = local_subtract(H, C(q,:)); end
end


function m = local_inIv(t, iv)
m = false(size(t));
for k = 1:size(iv, 1), m = m | (t >= iv(k,1) & t < iv(k,2)); end
end


function S = local_loadSession(fld, RP)
S = [];
[X, dur, tOn, nU, ts] = local_ripCounts(fld, RP);
if isempty(X) || nU < 2, return; end
S = struct('X', X, 'dur', dur, 'tOn', tOn, 'ts', {ts}, 'nrem', zeros(0, 2));
sf = dir(fullfile(fld, '*.SleepState.states.mat'));
if ~isempty(sf)
    try
        Ls = load(fullfile(sf(1).folder, sf(1).name));
        fn = fieldnames(Ls); st = Ls.(fn{1});
        if isfield(st, 'ints') && isfield(st.ints, 'NREMstate'), S.nrem = st.ints.NREMstate; end
    catch
    end
end
end


%% =============================================================================================
%  report
%% =============================================================================================
function R = local_report(F, O, nBoot, eqRatio, yieldSub)
R = struct();
grid = O.tauGrid(:); gMax = grid(end); gMin = grid(1);
[~, iRef] = min(abs(log(grid) - log(O.tauRef)));
fprintf('\n================ PAIR COINCIDENCE GLM: PERSISTENCE OF EXPLORATION COVARIANCE ================\n');
fprintf('log E[x_j x_k] = log(mu_j mu_k) + u_jk + c_jk [delta |D| + 1(D>0) A exp(-D/tau)]\n');
fprintf('A: log(observed/expected) coincidences per SD of exploration pair correlation, just after the\n');
fprintf('exposure, beyond each pair''s own baseline (u) and symmetric drift (delta). tau shared within an\n');
fprintf('arm; CI by animal bootstrap of the summed profile curves; A tested by sign-flip over animals.\n');
fprintf('Never-infused animals are the three mice; session-level p is printed beside them.\n');
groups = [{'Control'}, arrayfun(@(n) sprintf('Ctrl%d', n), yieldSub, 'UniformOutput', false), {'SAL', 'ASV'}];
groups = groups(ismember(groups, {F.Group}));
models = {'E', 'exploration c'; 'H', 'home-cage c (state control)'; 'J', 'joint: exploration A with home-cage c in the model'};
for c = {'ConA1', 'ConB1'}
    fprintf('\n==================== %s template ====================\n', c{1});
    for g = groups
        k = strcmp({F.Group}, g{1}) & strcmp({F.Context}, c{1});
        if ~any(k), continue; end
        Fg = F(k);
        fprintf('\n  %-8s %d sessions, %d animals | units %s | pairs %s | ripples pre %s, post %s\n', g{1}, numel(Fg), ...
            numel(unique({Fg.Animal})), local_med([Fg.nUnits]), local_med([Fg.nPairs]), local_med([Fg.nRipPre]), local_med([Fg.nRipPost]));
        pr = vertcat(Fg.postRange); pe = vertcat(Fg.preRange);
        fprintf('           lag covered: post %.2f-%.2f h, pre %.2f-%.2f h (medians of each session''s extent)\n', ...
            median(pr(:,1)), median(pr(:,2)), median(pe(:,1)), median(pe(:,2)));
        for mm = 1:size(models, 1)
            md = models{mm, 1};
            okm = arrayfun(@(x) isstruct(x.(md)) && ~isempty(x.(md)), Fg);
            if ~any(okm), continue; end
            Fm = Fg(okm);
            out = local_armFit(Fm, md, grid, nBoot, iRef);
            fprintf('    [%s] %s  (%d sessions)\n', md, models{mm, 2}, numel(Fm));
            local_printArm(out, gMin, gMax, O.tauRef, numel(unique({Fm.Animal})) < 5);
            R.(g{1}).(c{1}).(md) = out;
        end
        if isfield(Fg, 'rEH')
            r = [Fg.rEH]; r = r(isfinite(r));
            if ~isempty(r), fprintf('    corr(exploration c, home-cage c) across pairs: median %.2f\n', median(r)); end
        end
    end
    % ---- ASV vs SAL on the exploration model ----
    kS = strcmp({F.Group}, 'SAL') & strcmp({F.Context}, c{1});
    kA = strcmp({F.Group}, 'ASV') & strcmp({F.Context}, c{1});
    if any(kS) && any(kA)
        fprintf('\n  ASV vs SAL (%s, exploration model)\n', c{1});
        local_compareArms(F(kS), F(kA), grid, nBoot, eqRatio, iRef);
        if isfield(R, 'Control') && isfield(R.Control, c{1})
            mdC = 'E'; if isfield(R.Control.(c{1}), 'J'), mdC = 'J'; end
            tauC = R.Control.(c{1}).(mdC).tau; [~, iC] = min(abs(log(grid) - log(tauC)));
            fprintf('     at the never-infused tau (%s model, %.2f h):\n', mdC, grid(iC));
            for gg = {'SAL', 'ASV'}
                for md = {'E', 'J'}
                    kk = strcmp({F.Group}, gg{1}) & strcmp({F.Context}, c{1}) & ...
                        arrayfun(@(x) isstruct(x.(md{1})) && ~isempty(x.(md{1})), F);
                    if ~any(kk), continue; end
                    Fm = F(kk); [~, ~, ai] = unique({Fm.Animal});
                    a = local_byAnimal(arrayfun(@(x) x.(md{1}).th(iC, 2), Fm(:)), ai);
                    fprintf('        %-3s [%s] A = %+.4f +- %.4f (animals p = %.3f, n = %d)\n', gg{1}, md{1}, ...
                        mean(a), local_sem(a), local_signflip(a), numel(a));
                end
            end
        end
    end
end
fprintf('\n  Yield check: Ctrl9 / Ctrl20 tau should sit inside the Control CI. If it does not, the infused\n');
fprintf('  arms (median ~20 scored cells) cannot be timed with this model either.\n');
end


function out = local_armFit(Fm, md, grid, nBoot, iRef)
% summed profile, bootstrap over animals; A and delta at the arm tau and at the reference tau
Lm = cell2mat(arrayfun(@(x) x.(md).L(:)', Fm(:), 'UniformOutput', false));     % sessions x grid
an = {Fm.Animal}; ua = unique(an);
[~, ai] = ismember(an, ua);
Lsum = sum(Lm, 1, 'omitnan');
[~, qh] = max(Lsum); tau = grid(qh);
tB = nan(nBoot, 1);
for b = 1:nBoot
    pick = randi(numel(ua), numel(ua), 1);
    w = accumarray(pick, 1, [numel(ua) 1]);
    Lb = sum(Lm .* w(ai(:)), 1, 'omitnan');
    [~, q] = max(Lb); tB(b) = grid(q);
end
iA = 2;                                         % [delta A ...]
thH = cell2mat(arrayfun(@(x) x.(md).th(qh, :), Fm(:), 'UniformOutput', false));
thR = cell2mat(arrayfun(@(x) x.(md).th(iRef, :), Fm(:), 'UniformOutput', false));
out = struct('tau', tau, 'tauBoot', tB, 'Lsum', Lsum, ...
    'A', local_byAnimal(thH(:, iA), ai), 'delta', local_byAnimal(thH(:, 1), ai), ...
    'Aref', local_byAnimal(thR(:, iA), ai), 'Ases', thH(:, iA), 'ArefSes', thR(:, iA), ...
    'animals', {ua}, 'ai', ai, 'sessions', {{Fm.Session}});
if size(thH, 2) >= 4
    out.AH = local_byAnimal(thH(:, 4), ai); out.AHref = local_byAnimal(thR(:, 4), ai);
end
end


function x = local_byAnimal(v, ai)
x = accumarray(ai(:), v(:), [], @(z) mean(z, 'omitnan'));
end


function local_printArm(o, gMin, gMax, tauRef, sesLevel)
ci = prctile(o.tauBoot, [2.5 97.5]);
fC = mean(o.tauBoot >= gMax * 0.999); fF = mean(o.tauBoot <= gMin * 1.001);
pA = local_signflip(o.A); pAs = local_signflip(o.Ases);
fprintf('        tau = %.2f h [95%% CI %.2f-%.2f], half-life %.2f h | draws at ceiling %.0f%%, floor %.0f%%\n', ...
    o.tau, ci(1), ci(2), o.tau * log(2), 100 * fC, 100 * fF);
fprintf('        A(tau)   = %+.4f +- %.4f  (animals p = %.3f, n = %d; sessions p = %.3f)\n', ...
    mean(o.A), local_sem(o.A), pA, numel(o.A), pAs);
fprintf('        A(%g h)  = %+.4f +- %.4f  (animals p = %.3f; sessions p = %.3f)\n', tauRef, ...
    mean(o.Aref), local_sem(o.Aref), local_signflip(o.Aref), local_signflip(o.ArefSes));
fprintf('        delta    = %+.4f +- %.4f per h (animals p = %.3f)\n', mean(o.delta), local_sem(o.delta), local_signflip(o.delta));
if isfield(o, 'AH')
    fprintf('        home-cage A in the same model: %+.4f +- %.4f (animals p = %.3f)\n', mean(o.AH), local_sem(o.AH), local_signflip(o.AH));
end
pUse = pA; if sesLevel, pUse = pAs; end
if ~(mean(o.A) > 0 && pUse < 0.05)
    fprintf('        -> A not reliably > 0: nothing to time\n');
elseif fF > 0.025
    fprintf('        -> tau not bounded below: decays faster than the nearest ripples resolve\n');
elseif fC > 0.025
    fprintf('        -> tau not bounded above: no decay detected within the lags covered\n');
else
    fprintf('        -> timeable\n');
end
end


function local_compareArms(FS, FA, grid, nBoot, eqRatio, iRef)
% common tau over both arms, then paired A; tau ratio when both arms are timeable
oS = local_armFit(FS, 'E', grid, nBoot, iRef); oA = local_armFit(FA, 'E', grid, nBoot, iRef);
Lc = oS.Lsum + oA.Lsum; [~, qc] = max(Lc);
aS = arrayfun(@(x) x.E.th(qc, 2), FS(:)); aA = arrayfun(@(x) x.E.th(qc, 2), FA(:));
mS = local_byAnimal(aS, oS.ai); mA = local_byAnimal(aA, oA.ai);
[an, ia, ib] = intersect(oA.animals, oS.animals);
fprintf('     common tau %.2f h: A ASV %+.4f, SAL %+.4f (animal means)\n', grid(qc), mean(mA), mean(mS));
if numel(an) >= 2
    d = mA(ia) - mS(ib);
    fprintf('     paired A, ASV - SAL = %+.4f +- %.4f, p = %.3f (%d animals)\n', mean(d), local_sem(d), local_signflip(d), numel(an));
    dr = local_byAnimal(oA.ArefSes, oA.ai); sr = local_byAnimal(oS.ArefSes, oS.ai);
    d2 = dr(ia) - sr(ib);
    fprintf('     paired A at the reference tau, ASV - SAL = %+.4f +- %.4f, p = %.3f\n', mean(d2), local_sem(d2), local_signflip(d2));
end
if exist('fitlme', 'file') == 2
    try
        T = table([aS; aA], categorical([repmat({'SAL'}, numel(aS), 1); repmat({'ASV'}, numel(aA), 1)], {'SAL','ASV'}), ...
            [{FS.Animal}'; {FA.Animal}'], 'VariableNames', {'A','Group','Animal'});
        lme = fitlme(T, 'A ~ Group + (1|Animal)');
        [b, nm, st] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
        k = strcmp(nm.Name, 'Group_ASV');
        fprintf('     LME A ~ Group + (1|Animal): ASV - SAL = %+.4f, p = %.3f (%d sessions)\n', b(k), st.pValue(k), height(T));
    catch err
        fprintf('     LME not fitted (%s)\n', err.message);
    end
end
okS = local_timeable(oS, grid); okA = local_timeable(oA, grid);
if okS && okA
    an2 = union(oS.animals, oA.animals); r = nan(nBoot, 1);
    LmS = cell2mat(arrayfun(@(x) x.E.L(:)', FS(:), 'UniformOutput', false));
    LmA = cell2mat(arrayfun(@(x) x.E.L(:)', FA(:), 'UniformOutput', false));
    [~, aS2] = ismember({FS.Animal}, an2); [~, aA2] = ismember({FA.Animal}, an2);
    for b = 1:nBoot
        w = accumarray(randi(numel(an2), numel(an2), 1), 1, [numel(an2) 1]);
        [~, q1] = max(sum(LmS .* w(aS2(:)), 1, 'omitnan')); [~, q2] = max(sum(LmA .* w(aA2(:)), 1, 'omitnan'));
        r(b) = grid(q2) / grid(q1);
    end
    ci95 = prctile(r, [2.5 97.5]); ci90 = prctile(r, [5 95]);
    fprintf('     tau ASV / SAL = %.2f [95%% CI %.2f-%.2f; 90%% CI %.2f-%.2f]\n', oA.tau / oS.tau, ci95(1), ci95(2), ci90(1), ci90(2));
    if ci90(1) >= 1 / eqRatio && ci90(2) <= eqRatio
        fprintf('     equivalent within a factor of %g (TOST, alpha 0.05)\n', eqRatio);
    else
        fprintf('     NOT shown equivalent within a factor of %g\n', eqRatio);
    end
else
    fprintf('     tau ratio not estimable: %s\n', strjoin([local_ifc(~okS, {'SAL'}), local_ifc(~okA, {'ASV'})], ' and '));
end
end


function ok = local_timeable(o, grid)
fC = mean(o.tauBoot >= grid(end) * 0.999); fF = mean(o.tauBoot <= grid(1) * 1.001);
ok = mean(o.A) > 0 && local_signflip(o.A) < 0.05 && fC <= 0.025 && fF <= 0.025;
end


function c = local_ifc(t, v)
if t, c = v; else, c = {}; end
end


function s = local_med(x)
s = sprintf('%g (%g-%g)', median(x), min(x), max(x));
end


function e = local_sem(x)
x = x(isfinite(x));
if numel(x) < 2, e = NaN; else, e = std(x) / sqrt(numel(x)); end
end


function p = local_signflip(x)
x = x(isfinite(x)); n = numel(x);
if n < 2, p = NaN; return; end
if n > 20
    s = sign(randn(n, 5000)); s(:,1) = 1;
else
    s = double(dec2bin(0:2^n - 1) - '0')' * 2 - 1;
end
m = abs(x(:)' * s / n);
p = max(mean(m >= abs(mean(x)) - 1e-12), 2 / 2^n);
end


%% =============================================================================================
%  tau-free: last time with evidence
%% =============================================================================================
function out = local_lastTime(F, O, alpha, minU, nFlipMax)
out = struct();
fprintf('\n================ LAST TIME WITH EVIDENCE (tau-free greedy scan) ================\n');
fprintf('Per-bin excess A_b (%.2f h bins) pooled within animal by information. A window grows from\n', O.lastBinH);
fprintf('the first undecided bin until its animal-level excess is > 0 at sign-flip p < %.2f, then closes.\n', alpha);
fprintf('p_cal: max t over all windows starting at that window''s start, against the same maximum with\n');
fprintf('unit signs flipped for every bin from there on (pays for the growing). Never-infused: unit =\n');
fprintf('session (three mice). Primary: no drift term; "drift" refits with delta from the pre side.\n');
for c = {'ConB1', 'ConA1'}
    for g = {'Control', 'SAL', 'ASV'}
        k = strcmp({F.Group}, g{1}) & strcmp({F.Context}, c{1}) & arrayfun(@(x) ~isempty(x.pwA0), F);
        if ~any(k), continue; end
        Fg = F(k); eb = Fg(1).pwEdges;
        if strcmp(g{1}, 'Control'), unit = {Fg.Session}; else, unit = {Fg.Animal}; end
        for v = {'nodrift', 'drift'}
            if strcmp(v{1}, 'nodrift'), Am = vertcat(Fg.pwA0); Im = vertcat(Fg.pwI0);
            else, Am = vertcat(Fg.pwA1); Im = vertcat(Fg.pwI1); end
            [uu, ~, ui] = unique(unit);
            nb = numel(eb) - 1; Au = nan(numel(uu), nb); Iu = zeros(numel(uu), nb);
            for q = 1:numel(uu)
                r = ui == q; a = Am(r, :); w = Im(r, :); w(~isfinite(a)) = 0; a(~isfinite(a)) = 0;
                Iu(q, :) = sum(w, 1);
                Au(q, :) = sum(a .* w, 1) ./ max(Iu(q, :), eps);
                Au(q, Iu(q, :) == 0) = NaN;
            end
            fprintf('\n  %s  %-7s %-7s %d units\n', c{1}, g{1}, v{1}, numel(uu));
            if strcmp(v{1}, 'nodrift')
                fprintf('     per bin (h):');
                for b = 1:nb
                    x = Au(:, b); x = x(isfinite(x));
                    if numel(x) >= 2, fprintf(' %.2f:%+.3f(%d)', eb(b), mean(x), numel(x)); end
                end
                fprintf('\n');
            end
            seg = local_greedy(Au, Iu, 1, alpha, minU);
            lastA = NaN;
            for z = 1:size(seg, 1)
                b1 = seg(z, 1); b2 = seg(z, 2);
                x = local_winEst(Au, Iu, b1, b2);
                if seg(z, 3)
                    pc = local_calib(Au, Iu, b1, minU, nFlipMax);
                    fprintf('     %.2f-%.2f h  A = %+.4f +- %.4f  p = %.3f  p_cal = %.3f  (n = %d)\n', eb(b1), eb(b2+1), ...
                        mean(x), local_sem(x), local_signflip(x), pc, numel(x));
                    if pc < 0.05, lastA = eb(b1); end
                else
                    fprintf('     %.2f-%.2f h  A = %+.4f +- %.4f  p = %.3f  not detected (scan ends)\n', eb(b1), eb(b2+1), ...
                        mean(x), local_sem(x), local_signflip(x));
                end
            end
            if isfinite(lastA)
                fprintf('     -> excess shown to persist into the window starting %.2f h after the exposure\n', lastA);
            else
                fprintf('     -> no window with calibrated evidence\n');
            end
            out.(c{1}).(g{1}).(v{1}) = struct('Au', Au, 'Iu', Iu, 'units', {uu}, 'seg', seg, 'last', lastA, 'edges', eb);
        end
    end
end
end


function seg = local_greedy(Au, Iu, b0, alpha, minU)
% rows [first bin, last bin, detected]
nb = size(Au, 2); seg = zeros(0, 3); s = b0;
while s <= nb
    det = false;
    for e = s:nb
        x = local_winEst(Au, Iu, s, e);
        if numel(x) >= minU && mean(x) > 0 && local_signflip(x) < alpha, det = true; break; end
    end
    seg(end+1, :) = [s e det]; %#ok<AGROW>
    if ~det, break; end
    s = e + 1;
end
end


function x = local_winEst(Au, Iu, b1, b2)
% information-weighted excess over bins b1..b2, one value per unit with any data there
w = Iu(:, b1:b2); a = Au(:, b1:b2); w(~isfinite(a)) = 0; a(~isfinite(a)) = 0;
sw = sum(w, 2); x = sum(a .* w, 2) ./ max(sw, eps); x = x(sw > 0);
end


function pc = local_calib(Au, Iu, b1, minU, nFlipMax)
% max-t over all windows [b1, e]: the family the greedy scan searched, so growing a window until
% it passes is paid for. Signs of every unit are flipped jointly for all bins >= b1.
% (In simulation: a binary "did the flipped scan detect anything" gave p_cal ~ the scan's
% nominal detection rate whatever the data, so the magnitude is what is compared.)
n = size(Au, 1);
if n <= 12, S = double(dec2bin(0:2^n - 1) - '0') * 2 - 1;
else, S = sign(randn(nFlipMax, n)); S(S == 0) = 1; S(1, :) = 1; end
T0 = local_scanT(Au, Iu, b1, minU);
if ~isfinite(T0), pc = NaN; return; end
c = 0;
for q = 1:size(S, 1)
    Af = Au; Af(:, b1:end) = Af(:, b1:end) .* S(q, :)';
    c = c + (local_scanT(Af, Iu, b1, minU) >= T0 - 1e-12);
end
pc = c / size(S, 1);
end


function T = local_scanT(Au, Iu, b1, minU)
T = -Inf;
for e = b1:size(Au, 2)
    x = local_winEst(Au, Iu, b1, e);
    if numel(x) < minU, continue; end
    sd = std(x);
    if sd > 0, T = max(T, mean(x) / (sd / sqrt(numel(x)))); end
end
end


function out = local_did(F, O, R)
% within session: A(ConB1) - A(ConA1), ASV vs SAL, at the 1-h reference tau and at the
% never-infused tau; the drug-minus-day contrast
out = struct();
grid = O.tauGrid(:);
fprintf('\n================ WITHIN-SESSION CONTRAST: A(ConB1) - A(ConA1), ASV vs SAL ================\n');
fprintf('ConA1 is before the infusion in both arms, so this removes session-order and day effects.\n');
taus = O.tauRef; lbl = {sprintf('%g h reference', O.tauRef)};
if isfield(R, 'Control') && isfield(R.Control, 'ConB1') && isfield(R.Control.ConB1, 'J')
    taus(end+1) = R.Control.ConB1.J.tau; lbl{end+1} = 'never-infused ConB1 tau';
end
for md = {'E', 'J'}
    for q = 1:numel(taus)
        [~, it] = min(abs(log(grid) - log(taus(q))));
        rows = {};
        for g = {'SAL', 'ASV'}
            ses = unique({F(strcmp({F.Group}, g{1})).Session});
            for s = 1:numel(ses)
                a = F(strcmp({F.Session}, ses{s}) & strcmp({F.Context}, 'ConA1') & strcmp({F.Group}, g{1}));
                b = F(strcmp({F.Session}, ses{s}) & strcmp({F.Context}, 'ConB1') & strcmp({F.Group}, g{1}));
                if isempty(a) || isempty(b) || isempty(a(1).(md{1})) || isempty(b(1).(md{1})), continue; end
                rows(end+1, :) = {g{1}, b(1).Animal, ses{s}, b(1).(md{1}).th(it, 2) - a(1).(md{1}).th(it, 2)}; %#ok<AGROW>
            end
        end
        if isempty(rows), continue; end
        T = cell2table(rows, 'VariableNames', {'Group','Animal','Session','d'});
        fprintf('\n  [%s] tau = %.2f h (%s)\n', md{1}, grid(it), lbl{q});
        mS = []; mA = []; aS = {}; aA = {};
        for g = {'SAL', 'ASV'}
            k = strcmp(T.Group, g{1}); [ua, ~, ai] = unique(T.Animal(k));
            m = accumarray(ai, T.d(k), [], @mean);
            fprintf('     %-3s B1 - A1 = %+.4f +- %.4f (animals p = %.3f, n = %d)\n', g{1}, mean(m), local_sem(m), local_signflip(m), numel(m));
            if strcmp(g{1}, 'SAL'), mS = m; aS = ua; else, mA = m; aA = ua; end
        end
        [~, ia, ib] = intersect(aA, aS);
        if numel(ia) >= 2
            d = mA(ia) - mS(ib);
            fprintf('     paired ASV - SAL = %+.4f +- %.4f, p = %.3f (%d animals)\n', mean(d), local_sem(d), local_signflip(d), numel(ia));
        end
        if exist('fitlme', 'file') == 2
            try
                T.Group = categorical(T.Group, {'SAL', 'ASV'});
                lme = fitlme(T, 'd ~ Group + (1|Animal)');
                [b, nm, st] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
                k = strcmp(nm.Name, 'Group_ASV');
                fprintf('     LME d ~ Group + (1|Animal): ASV - SAL = %+.4f, p = %.3f (%d sessions)\n', b(k), st.pValue(k), height(T));
            catch err
                fprintf('     LME not fitted (%s)\n', err.message);
            end
        end
        out.(md{1}).(sprintf('tau%d', q)) = T;
    end
end
end


%% =============================================================================================
%  shared with driftInclusionASV.m / compareReactivation.m (copied unchanged)
%% =============================================================================================
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


function s = local_name(fld)
[~, s] = fileparts(fld);
end


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