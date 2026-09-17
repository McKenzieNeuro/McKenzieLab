function est = gto2_dim_estimate(Xref, tref, opts)
% GTO2_DIM_ESTIMATE  Independent intrinsic-dimension estimates, to fix
% opts.d in gto2_preflight rather than let the per-probe eigengap heuristic
% vote. Three estimators, deliberately NOT the eigengap method preflight
% uses internally, so agreement across them is real evidence, not circularity.
%
%   est = gto2_dim_estimate(Xref, tref, opts)
%
% INPUTS
%   Xref  [n x F]  reference population rates (non-negative, uncentred), same
%                  matrix you pass as Xref to gto2_preflight.
%   tref  [n x 1]  timestamps (s), for the Theiler exclusion in the local
%                  estimator. Pass NaN(n,1) if genuinely independent samples.
%   opts  struct:
%     .k          neighbours for local participation ratio    (default 100,
%                                                                match preflight's k)
%     .n_sub      subsample size for TwoNN / global PR         (default 3000;
%                                                                O(n_sub^2) memory)
%     .n_local    number of anchor points for local PR         (default 300)
%     .theiler    Theiler window (s); [] = same estimator as preflight
%     .keep_frac  TwoNN: fraction of points kept after discarding the
%                 heaviest tail of mu=r2/r1 (outlier-robust MLE)  (default 0.9)
%
% OUTPUT est
%   .twonn        global ID, Facco et al. 2017 MLE estimator (density-robust,
%                  does not depend on choosing a neighbourhood size)
%   .PR_global    participation ratio of the global covariance spectrum
%                  (sum lam)^2 / sum(lam^2) -- soft, not integer-constrained
%   .PR_local     [n_local x 1] participation ratio computed the SAME way
%                  gate 1 sees it: local k-NN covariance per anchor. This is
%                  the one most directly comparable to gate 1's d_hat/q.
%   .PR_local_med median(PR_local) -- the recommended opts.d if PR_local is
%                  tight; round to nearest integer yourself after inspecting
%                  the spread, do not silently round here.
%   .anchors      indices into Xref used for PR_local (for cross-referencing)
%
% INTERPRETATION
%   - twonn and PR_global agree, both integer-ish, PR_local tight around the
%     same value  -> the manifold has one stable dimension; set opts.d to it.
%   - twonn/PR_global agree but PR_local is spread out across anchors
%     -> locally-varying effective dimension (e.g. curvature, non-stationary
%        density); a single global opts.d is an approximation, and gate 1's
%        q distribution should be read as "distance from flat" rather than a
%        hard fail.
%   - twonn and PR_global disagree substantially -> don't average them; PR is
%     sensitive to a few large eigenvalues (anisotropic noise inflates it),
%     TwoNN is sensitive to density non-uniformity and short-range noise.
%     Report both rather than picking one silently.

if nargin < 3, opts = struct(); end
k         = getdef(opts,'k',         1000);
n_sub     = getdef(opts,'n_sub',     3000);
n_local   = getdef(opts,'n_local',   300);
theiler   = getdef(opts,'theiler',   []);
keep_frac = getdef(opts,'keep_frac', 0.9);

n = size(Xref,1);

% ---------------------------------------------------------- subsample
n_sub = min(n_sub, n);
sub = randperm(n, n_sub);
Xs  = Xref(sub, :);

% ---------------------------------------------------------- TwoNN (global)
D2 = sqdist(Xs, Xs);
D2(1:n_sub+1:end) = Inf;                     % exclude self
[Dsort, ~] = sort(D2, 2, 'ascend');
r1 = sqrt(Dsort(:,1));
r2 = sqrt(Dsort(:,2));
mu = r2 ./ max(r1, eps);
mu = mu(isfinite(mu) & mu > 1);
mu_sorted = sort(mu, 'ascend');
n_keep = max(round(keep_frac * numel(mu_sorted)), 10);
mu_keep = mu_sorted(1:n_keep);
est.twonn = n_keep / sum(log(mu_keep));      % MLE, Facco et al. 2017 eq. 3

% ---------------------------------------------------------- global PR
Xc = Xs - mean(Xs,1);
C  = (Xc' * Xc) / (n_sub - 1);
lam = eig((C+C')/2);
lam = sort(lam(lam > 1e-12*max(lam)), 'descend');
est.PR_global = sum(lam)^2 / sum(lam.^2);

% ---------------------------------------------------------- local PR
if isempty(theiler)
    theiler = local_theiler_like_preflight(Xref, tref);
end
n_local = min(n_local, n);
anchors = randperm(n, n_local);
PR_local = nan(n_local,1);
for ii = 1:n_local
    i = anchors(ii);
    z = Xref(i,:);
    dist = sqrt(sum((Xref - z).^2, 2));
    if ~isnan(theiler) && ~any(isnan(tref))
        dist(abs(tref - tref(i)) < theiler) = Inf;
    end
    dist(i) = Inf;
    [~, ord] = sort(dist);
    idx = ord(1:min(k, sum(isfinite(dist))));
    if numel(idx) < 5, continue; end
    Nb = Xref(idx,:);
    mu_i = mean(Nb,1);
    Y = Nb - mu_i;
    Cl = (Y'*Y)/(numel(idx)-1);
    laml = eig((Cl+Cl')/2);
    laml = sort(laml(laml > 1e-12*max(laml)), 'descend');
    PR_local(ii) = sum(laml)^2 / sum(laml.^2);
end
est.PR_local = PR_local;
est.PR_local_med = median(PR_local, 'omitnan');
est.PR_local_iqr = [quantile(PR_local,0.25), quantile(PR_local,0.75)];
est.anchors = anchors;

fprintf('--- dimension estimate (n_sub=%d, n_local=%d, k=%d) ---\n', n_sub, n_local, k);
fprintf('  TwoNN (global, density-robust)  : %.2f\n', est.twonn);
fprintf('  participation ratio (global)    : %.2f\n', est.PR_global);
fprintf('  participation ratio (local,med) : %.2f  [IQR %.2f-%.2f]\n', ...
    est.PR_local_med, est.PR_local_iqr(1), est.PR_local_iqr(2));
if abs(est.twonn - est.PR_global) > 0.5*mean([est.twonn, est.PR_global])
    fprintf(['  NOTE: TwoNN and global PR disagree by >50%%. Do not average -- ', ...
             'they fail for different reasons (see help text). Report both.\n']);
end
if diff(est.PR_local_iqr) > 0.5*est.PR_local_med
    fprintf(['  NOTE: local PR IQR is wide relative to its median -- effective ', ...
             'dimension may vary across the population, not a single opts.d.\n']);
end
end

% ======================================================= helper functions
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function D2 = sqdist(A, B)
D2 = sum(A.^2,2) + sum(B.^2,2)' - 2*(A*B');
D2 = max(D2, 0);
end

function th = local_theiler_like_preflight(X, t)
if any(isnan(t)), th = NaN; return; end
[t, ix] = sort(t); X = X(ix,:);
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
