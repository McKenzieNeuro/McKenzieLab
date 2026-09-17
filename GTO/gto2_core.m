function D = gto2_core(A, Z, varargin)
% GTO2_CORE  Nested tangent/normal decomposition with a radial gain axis.
%
% D = GTO2_CORE(A, Z) decomposes each test state in Z (nTest-by-F) against a
% sampled reference manifold A (nRef-by-F).
%
% CONSTRUCTION (see Methods sec:bundles, sec:radial)
%   anchor        mu   = centroid of the k nearest reference points to z_B
%   tangent est.  Vd   = leading d right singular vectors of (N_A - mu)
%   normal proj.  P_N  = I - Vd*Vd'
%   radial dir.   rho  = mu/||mu||                    (needs a meaningful origin)
%   gain axis     ng   = P_N*rho / ||P_N*rho||        (radial, tangential part removed)
%   alignment     a    = ||P_N*rho||                  (in [0,1]; gates the gain axis)
%
%   The displacement r = z_B - mu then splits into THREE MUTUALLY ORTHOGONAL
%   pieces, because ng is orthogonal to Vd by construction:
%       r = r_tan  +  (ng'*r)*ng  +  r_novel
%   giving Pythagorean energy fractions that sum to exactly 1:
%       T = ||r_tan||^2 / ||r||^2      on-manifold pattern
%       G = (ng'*r)^2   / ||r||^2      gain-aligned
%       O = ||r_novel||^2 / ||r||^2    off-manifold novelty
%
% WHY ng AND NOT rho DIRECTLY
%   Using the raw radial direction as a gain axis makes the result depend on the
%   order of projection, because rho is not generally orthogonal to Vd. Removing
%   the tangential part first makes the three subspaces mutually orthogonal, so
%   ng'*r and ng'*(P_N*r) are identically equal and no order need be specified.
%   Check C1 in gto2_checks verifies this to machine precision.
%
% DEGENERACY
%   Where a = ||P_N*rho|| -> 0 the radial direction has been absorbed into the
%   estimated tangent space: gain modulation is locally indistinguishable from
%   movement along the manifold, and ng is not merely ill-estimated but
%   undefined. Points with a < 'alignTol' return G = O = NaN (T is still valid,
%   since it does not use ng). This is a LOUD failure by design -- the previous
%   radial-first construction returned a number regardless of alignment.
%
%   THE TOLERANCE IS SCIENTIFIC, NOT NUMERICAL. On a cone -- where the radial
%   direction is tangent everywhere by Euler's relation, so the gain axis is
%   undefined in theory -- the MEASURED alignment is about 0.06, not 1e-16:
%   finite sampling, noise and curvature keep ||P_N*rho|| well off zero even in
%   the fully degenerate case. A machine-epsilon tolerance therefore flags
%   nothing and lets ng be determined by the few percent of rho that happened to
%   survive projection. The default below is set where the gain axis stops
%   carrying usable information rather than where it stops existing; report
%   .align and choose deliberately for a given geometry.
%
% NO MAGNITUDE IS ESTIMATED
%   This returns the composition of a change in direction. It does not estimate
%   alpha, and does not identify the manifold state the change proceeded from.
%
% OPTIONS
%   'k'         40      neighbours defining the local patch
%   'd'         2       assumed tangent dimension
%   'alignTol'  0.1     below this alignment, G and O are returned as NaN
%                       (a scientific threshold, not a numerical one -- see above)
%
% RETURNS a struct with per-test-point columns (nTest-by-1 unless noted):
%   .G .T .O    energy fractions, G+T+O = 1 exactly where defined
%   .align      a = ||P_N*rho||
%   .rnorm      ||r||, the displacement magnitude from the anchor
%   .ell        tangential neighbourhood extent, sqrt(sum ||Vd'(x-mu)||^2/(k-1));
%               this is the ell of Propositions 2-3 and should be used in every
%               prediction (delta = ell^2/2R, tau^2 = ell^2/k)
%   .rnbr       ambient neighbourhood extent, sqrt(sum ||x-mu||^2/(k-1)); includes
%               normal spread and noise in all F dimensions -- NOT for predictions
%   .mu         nTest-by-F, the anchor used for each test point
%   .ng         nTest-by-F, the gain axis used for each test point
%   .degenerate logical, true where alignment fell below alignTol

ip = inputParser;
ip.addParameter('k', 40);
ip.addParameter('d', 2);
ip.addParameter('alignTol', 0.1);
ip.parse(varargin{:});
o = ip.Results;

[M, F] = size(Z);
nRef = size(A, 1);
if size(A,2) ~= F
    error('gto2_core:dim', 'A and Z must have the same number of columns.');
end
if o.d >= F
    error('gto2_core:nogainaxis', ...
        'd = %d leaves no room for a gain axis in F = %d dimensions (need d < F).', o.d, F);
end
k = min(o.k, nRef);

% kNN in chunks, so memory stays bounded at large nRef (a full 4000-by-40000
% distance matrix plus its sort needs several GB, and some pdist2 implementations
% build a chunk-by-nRef-by-F intermediate). Neighbours are identical to sorting
% the full matrix.
knn = zeros(M, k);
chunk = max(1, floor(5e7 / (nRef * F)));
for c0 = 1:chunk:M
    idx = c0:min(M, c0 + chunk - 1);
    [~, ord] = sort(pdist2(Z(idx,:), A), 2);
    knn(idx,:) = ord(:, 1:k);
end

D.G = nan(M,1);  D.T = nan(M,1);  D.O = nan(M,1);
D.align = nan(M,1);  D.rnorm = nan(M,1);  D.rnbr = nan(M,1);
D.mu = nan(M,F);  D.ng = nan(M,F);
D.degenerate = false(M,1);
D.ell = nan(M,1);
for i = 1:M
    loc = A(knn(i,:), :);
    mu  = mean(loc, 1).';
    Xc  = loc - mu.';
    [~, ~, V] = svd(Xc, 'econ');
    dd = min(o.d, size(V,2));
    Vd = V(:, 1:dd);
    % Neighbourhood extent. Both use the unbiased (k-1) normalisation, since the
    % centroid is estimated from the same k points: dividing by k underestimates
    % E||u||^2 by (k-1)/k.
    Xc_tan    = Xc * Vd;                                % k-by-d tangential coordinates
    D.ell(i)  = sqrt(sum(Xc_tan(:).^2) / (k - 1));      % TANGENTIAL extent: the ell of Props 2-3
    D.rnbr(i) = sqrt(sum(Xc(:).^2)     / (k - 1));      % ambient extent (includes normal + noise)
    zB = Z(i,:).';
    r  = zB - mu;
    nr2 = sum(r.^2);

    r_tan = Vd*(Vd.'*r);
    PNr   = r - r_tan;

    % radial direction and its normal-space part
    nmu = norm(mu);
    rho = mu / (nmu + eps);
    PNrho = rho - Vd*(Vd.'*rho);
    a = norm(PNrho);

    D.align(i) = a;
    D.rnorm(i) = sqrt(nr2);
    D.mu(i,:)  = mu.';
    D.T(i)     = sum(r_tan.^2) / (nr2 + eps);

    if a < o.alignTol
        % gain axis undefined: radial direction lies in the estimated tangent
        % space. T remains valid; G and O are not defined here.
        D.degenerate(i) = true;
        continue
    end

    ng = PNrho / a;
    gcomp = ng.' * r;                 % == ng.'*PNr exactly, since ng _|_ Vd
    r_novel = PNr - gcomp*ng;

    D.ng(i,:) = ng.';
    D.G(i)    = gcomp^2 / (nr2 + eps);
    D.O(i)    = sum(r_novel.^2) / (nr2 + eps);
end
end
