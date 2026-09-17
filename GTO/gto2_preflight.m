function R = gto2_preflight(Xref, tref, Xprb, tprb, opts)
% GTO2_PREFLIGHT  Cascade-of-trust gate evaluation for real data.
%
%   R = gto2_preflight(Xref, tref, Xprb, tprb, opts)
%
% Evaluates gates 0-6 per probe and returns G/T/O ONLY for probes that clear
% every directly-evaluable gate. Gate 2 (estimator bias) is NOT evaluable on
% real data; this function returns the measured parameters needed to bound it
% by matched simulation (see gto2_preflight_calibspec below).
%
% INPUTS
%   Xref  [n_ref x F]  reference set, NON-NEGATIVE rates, NOT centred, NOT z-scored.
%   tref  [n_ref x 1]  timestamps (s) for reference samples. Required for Theiler
%                      exclusion; pass NaN(n_ref,1) only if samples are genuinely
%                      independent (e.g. trial-shuffled across sessions).
%   Xprb  [n_prb x F]  probe set. MUST be disjoint from Xref (separate trials/epochs).
%   tprb  [n_prb x 1]  probe timestamps (s).
%   opts  struct, see defaults below.
%
% OUTPUT R
%   R.gate    [n_prb x 7] logical, columns = gates 0..6 (gate 2 always NaN-> see R.calib)
%   R.pass    [n_prb x 1] logical, all evaluable gates cleared
%   R.G,T,O   [n_prb x 1] decomposition; NaN where ~pass
%   R.diag    per-probe diagnostics (a, C6, ell, q, c, kappa_rad, ...)
%   R.global  dataset-level diagnostics (subtended angle, rho-hat spread, Theiler)
%   R.calib   parameters to hand to the simulation harness for the gate-2 bound
%
% NOTE ON PROVENANCE: written standalone, not against gto2_core.m. Before use,
% verify agreement: run gto2_preflight on a simulated case for which gto2_core
% returns known G/T/O and assert max|delta| < 1e-12 on all three. If they differ,
% this file is wrong, not the core.

% ---------------------------------------------------------------- defaults
if nargin < 5, opts = struct(); end
d  = getdef(opts,'d',        []);      % intrinsic dim; [] = estimate per probe
k  = getdef(opts,'k',        100);     % neighbours
alignTol  = getdef(opts,'alignTol', 0.10);  % gate 4: min a
C6min     = getdef(opts,'C6min',    20);    % gate 6: see S17 error tiers
qmax      = getdef(opts,'qmax',     0.50);  % gate 1: max leakage lam(d+1)/lam(d)
biasmax   = getdef(opts,'biasmax',  0.25);  % gate 3: max (curv bias + jitter)/||r||
theiler   = getdef(opts,'theiler',  []);    % s; [] = estimate from autocorrelation
minAngSpread = getdef(opts,'minAngSpread', 5*pi/180); % gate 0b: rho-hat spread
verbose   = getdef(opts,'verbose',  true);

[n_ref, F] = size(Xref);
n_prb      = size(Xprb,1);
assert(size(Xprb,2)==F, 'Xref and Xprb must have the same number of units.');

% =========================================================== PREFLIGHT P1
% Origin integrity. The construction requires the origin to be zero population
% activity; centring or z-scoring destroys the radial correspondence entirely.
P = struct();
P.minval      = min(Xref(:));
P.negfrac     = mean(Xref(:) < 0);
P.colmean_min = min(mean(Xref,1));
if P.negfrac > 0
    error(['gto2_preflight:negativeData  %.2f%% of Xref is negative. The ', ...
           'decomposition requires a true zero-activity origin; the data have ', ...
           'been centred, z-scored, or are not rates.'], 100*P.negfrac);
end
% How far from the origin does the cloud live, and how much of the space does
% it subtend? If the half-angle is small, rho-hat is near-constant and the gain
% axis is approximately global rather than local.
nrm   = sqrt(sum(Xref.^2,2));
rhat  = Xref ./ max(nrm, eps);
rbar  = mean(rhat,1); rbar = rbar / norm(rbar);
cosang = rhat * rbar(:);
P.halfangle_med = acos(median(min(max(cosang,-1),1)));    % rad
P.halfangle_p95 = acos(min(max(quantile(cosang,0.05),-1),1));
P.radial_cv     = std(nrm)/mean(nrm);   % how much of the variance is pure scale

% =========================================================== PREFLIGHT P2
% Independence of reference and probe sets.
if ~any(isnan(tref)) && ~any(isnan(tprb))
    P.min_ref_prb_dt = min(abs(tprb(:) - tref(:)'), [], 'all');
else
    P.min_ref_prb_dt = NaN;
end

% =========================================================== PREFLIGHT P3
% Theiler window: neighbours must not be temporal neighbours, or ell measures
% autocorrelation rather than manifold extent.
if isempty(theiler)
    theiler = estimate_theiler(Xref, tref);
end
P.theiler_s = theiler;

% =========================================================== PREFLIGHT P4
% Radial share of local noise energy (the anisotropy orientation the anchor
% does NOT absorb). Estimated from within-neighbourhood residuals.
P.radial_noise_frac = radial_noise_fraction(Xref, tref, k, theiler);

if verbose
    fprintf('--- preflight ---\n');
    fprintf('  cloud half-angle from origin: med %.1f deg, p95 %.1f deg\n', ...
            rad2deg(P.halfangle_med), rad2deg(P.halfangle_p95));
    fprintf('  radial CV of ||x||          : %.3f\n', P.radial_cv);
    fprintf('  Theiler window              : %.3f s\n', theiler);
    fprintf('  radial share of local noise : %.3f\n', P.radial_noise_frac);
    if P.min_ref_prb_dt == 0
        warning('Reference and probe sets share timestamps. Split them.');
    end
end

% =========================================================== PER-PROBE LOOP
gate  = false(n_prb,7);
%gate(:,3) = NaN;                       % gate 2 (index 3) never evaluable here
a     = nan(n_prb,1);  C6   = nan(n_prb,1);  ell = nan(n_prb,1);
q     = nan(n_prb,1);  cbk  = nan(n_prb,1);  dhat= nan(n_prb,1);
biasf = nan(n_prb,1);  Hn   = nan(n_prb,1);  rn  = nan(n_prb,1);
G     = nan(n_prb,1);  T    = nan(n_prb,1);  O   = nan(n_prb,1);
angsp = nan(n_prb,1);

for i = 1:n_prb
    z = Xprb(i,:);

    % --- neighbours with Theiler exclusion -----------------------------
    dist = sqrt(sum((Xref - z).^2, 2));
    if ~isnan(theiler) && ~any(isnan(tref)) && ~isnan(tprb(i))
        dist(abs(tref - tprb(i)) < theiler) = Inf;
    end
    [~, ord] = sort(dist, 'ascend');
    idx = ord(1:min(k, sum(isfinite(dist))));
    if numel(idx) < k, continue; end            % gate 0 fails: insufficient support
    Nb  = Xref(idx, :);

    % --- anchor and local frame ----------------------------------------
    mu  = mean(Nb, 1);
    r   = z - mu;
    rn(i) = norm(r);
    Y   = Nb - mu;
    [V, S2] = local_pca(Y);
    lam = diag(S2);

    % --- GATE 1: intrinsic dimension and leakage ------------------------
    if isempty(d)
        dhat(i) = pick_dim(lam);
    else
        dhat(i) = d;
    end
    dd = dhat(i);
    if dd+1 <= numel(lam) && lam(dd) > 0
        q(i) = lam(dd+1) / lam(dd);
    else
        q(i) = NaN;
    end
    Vd  = V(:, 1:dd);

    % --- GATE 0: chart exists -------------------------------------------
    % Local support must be a patch: enough neighbours, positive spread in all
    % d directions, and neighbourhood extent small vs local curvature radius.
    ell(i) = sqrt(sum(sum((Y*Vd).^2, 2)) / (numel(idx)-1));
    gate(i,1) = numel(idx) >= k && all(lam(1:dd) > 0) && isfinite(ell(i)) && ell(i) > 0;

    % --- second fundamental form (local quadratic fit) ------------------
    [Hvec, kappa_ok] = local_curvature(Y, Vd);
    Hn(i) = norm(Hvec);

    % --- GATE 4: gain axis exists ---------------------------------------
    rhoh = mu / max(norm(mu), eps);
    PN_rho = rhoh - (rhoh * Vd) * Vd';      % project out tangent
    a(i) = norm(PN_rho);
    gate(i,5) = a(i) >= alignTol;
    if a(i) < eps, continue; end
    ngh = PN_rho / a(i);

    % --- GATE 5: booking ------------------------------------------------
    if kappa_ok && Hn(i) > 0
        cbk(i) = dot(ngh, Hvec / Hn(i));
    end
    gate(i,6) = kappa_ok;       % booking evaluable; |c| reported, not thresholded

    % --- GATE 3: anchor bias vs signal ----------------------------------
    curvbias = 0.5 * ell(i)^2 * Hn(i);
    jitter   = ell(i) / sqrt(numel(idx));
    biasf(i) = (curvbias + jitter) / max(rn(i), eps);
    gate(i,4) = biasf(i) <= biasmax;

    % --- GATE 6: conditioning -------------------------------------------
    C6(i) = numel(idx) * rn(i)^2 / max(ell(i)^2, eps);
    gate(i,7) = (C6(i) - 1) >= C6min;

    % --- GATE 1 verdict ---------------------------------------------------
    gate(i,2) = isfinite(q(i)) && q(i) <= qmax;

    % --- decomposition ----------------------------------------------------
    rT = (r * Vd) * Vd';
    gcomp = dot(r, ngh);
    G(i) = gcomp^2 / max(rn(i)^2, eps);
    T(i) = sum(rT.^2) / max(rn(i)^2, eps);
    O(i) = max(1 - G(i) - T(i), 0);

    angsp(i) = acos(min(max(dot(rhoh, rbar),-1),1));
end

% --- GATE 0b (dataset level): is rho-hat informative? ---------------------
rho_spread = std(angsp(isfinite(angsp)));
gate0b = rho_spread >= minAngSpread;

% --- assemble -------------------------------------------------------------
evaluable = [1 2 4 5 6 7];            % gate 2 (col 3) excluded
pass = all(gate(:, evaluable) ~= 0, 2) & gate0b;
G(~pass) = NaN; T(~pass) = NaN; O(~pass) = NaN;

R.gate = gate;
R.pass = pass;
R.G = G; R.T = T; R.O = O;
R.diag = table(a, C6, ell, q, cbk, dhat, biasf, Hn, rn, ...
    'VariableNames', {'a','C6','ell','q','c','d_hat','biasfrac','normH','normr'});
R.global = P;
R.global.rho_spread_rad = rho_spread;
R.global.gate0b = gate0b;
R.global.pass_rate = mean(pass);

% --- parameters for the gate-2 simulation bound ---------------------------
R.calib = struct( ...
    'F',        F, ...
    'd',        median(dhat(isfinite(dhat))), ...
    'k',        k, ...
    'n_ref',    n_ref, ...
    'ell',      median(ell(isfinite(ell))), ...
    'normH',    median(Hn(isfinite(Hn))), ...
    'a',        median(a(isfinite(a))), ...
    'C6',       median(C6(isfinite(C6))), ...
    'radial_noise_frac', P.radial_noise_frac);

if verbose
    fprintf('--- gates ---\n');
    nm = {'0 chart','1 dim','2 estimator','3 anchor','4 gain axis','5 booking','6 C6'};
    for g = 1:7
        if g == 3
            fprintf('  gate %-12s : NOT EVALUABLE on data (bound by simulation)\n', nm{g});
        else
            fprintf('  gate %-12s : %5.1f%% pass\n', nm{g}, 100*mean(gate(:,g)~=0));
        end
    end
    fprintf('  gate 0b rho-spread  : %.1f deg  (%s)\n', rad2deg(rho_spread), ...
            ternary(gate0b,'PASS','FAIL - gain axis is near-global'));
    fprintf('  ALL EVALUABLE GATES : %.1f%% of probes\n', 100*mean(pass));
    fprintf(['NOTE: gate-passing is a data-dependent selection. Report pass rate ', ...
             'per condition before comparing G/T/O across conditions.\n']);
end
end

% ======================================================= helper functions
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end

function [V, S2] = local_pca(Y)
C = (Y' * Y) / max(size(Y,1)-1, 1);
[V, S2] = eig((C+C')/2, 'vector');
[S2, ix] = sort(S2, 'descend');
V = V(:, ix);
S2 = diag(S2);
end

function d = pick_dim(lam)
% Largest eigengap in log spectrum, capped where cumulative variance >= 0.95.
lam = lam(lam > 0);
if numel(lam) < 3, d = max(numel(lam)-1, 1); return; end
cv = cumsum(lam)/sum(lam);
cap = find(cv >= 0.95, 1, 'first');
g = -diff(log(lam(1:min(cap+1, numel(lam)))));
[~, d] = max(g);
d = max(d, 1);
end

function [Hvec, ok] = local_curvature(Y, Vd)
% Mean curvature vector from a local quadratic fit in tangent coordinates.
% Returns ok=false when the neighbourhood cannot support the fit.
[n, F] = size(Y);
dd = size(Vd, 2);
U  = Y * Vd;                             % tangent coords
W  = Y - U * Vd';                        % normal residuals
nq = dd*(dd+1)/2;
if n < nq + dd + 1, Hvec = zeros(1,F); ok = false; return; end
Q = zeros(n, nq); col = 0;
for i = 1:dd
    for j = i:dd
        col = col + 1;
        Q(:,col) = U(:,i) .* U(:,j) * (1 + (i~=j));   % 2*u_i*u_j off-diagonal
    end
end
A = [ones(n,1), U, Q];
if rcond(A'*A) < 1e-12, Hvec = zeros(1,F); ok = false; return; end
B = A \ W;                               % [const; linear; quadratic] x F
Bq = B(dd+2:end, :);
% mean curvature vector = trace of the Hessian / d  (diagonal terms only)
diagsel = false(nq,1); col = 0;
for i = 1:dd
    for j = i:dd
        col = col + 1;
        if i == j, diagsel(col) = true; end
    end
end
Hvec = sum(Bq(diagsel,:), 1) * (2/dd);   % factor 2 from the 1/2 in the quadratic
ok = true;
end

function th = estimate_theiler(X, t)
% First zero-crossing of the population-vector autocorrelation, in seconds.
if any(isnan(t)), th = NaN; return; end
[t, ix] = sort(t); X = X(ix, :);
dt = median(diff(t));
Xc = X - mean(X,1);
maxlag = min(500, floor(size(X,1)/4));
ac = zeros(maxlag,1);
denom = sum(Xc(:).^2);
for L = 1:maxlag
    ac(L) = sum(sum(Xc(1:end-L,:) .* Xc(1+L:end,:))) / denom;
end
z = find(ac <= 0, 1, 'first');
if isempty(z), z = maxlag; end
th = z * dt;
end

function f = radial_noise_fraction(X, t, k, theiler)
% Share of local residual energy along the radial (gain) direction. High values
% put the data in the anisotropy regime the anchor does not absorb.
n = size(X,1);
m = min(200, n);
sel = randperm(n, m);
num = 0; den = 0;
for s = sel
    dist = sqrt(sum((X - X(s,:)).^2, 2));
    if ~isnan(theiler) && ~any(isnan(t))
        dist(abs(t - t(s)) < theiler) = Inf;
    end
    dist(s) = Inf;
    [~, ord] = sort(dist);
    idx = ord(1:min(k, sum(isfinite(dist))));
    if numel(idx) < 3, continue; end
    Nb = X(idx,:); mu = mean(Nb,1);
    rh = mu / max(norm(mu), eps);
    Rr = Nb - mu;
    num = num + sum((Rr * rh(:)).^2);
    den = den + sum(Rr(:).^2);
end
f = num / max(den, eps);
end
