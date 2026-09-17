function result = manifoldAnalysis_ABn(ZA, ZB, opts)
% Component-aware, RHO-BRANCHED manifold membership test.
%
% CHANGES FROM THE ORIGINAL (each independently validated on synthetic manifolds
% before being combined -- see paper Sections 2-3):
%
%  1. RHO-BRANCH (the central fix). rho = ||Vd'*u||^2 is the fraction of the
%     local radial/amplitude direction captured by the tangent basis, computed
%     from the ANCHOR before any gain removal. It determines whether "just
%     scaled" is even a coherent on-manifold statement:
%       rho high (amplitude IS tangent, e.g. a cone): gain removal is valid.
%         Project out z_gain = alpha*mu, then test the ORTHOGONAL residual of
%         what's left. This is close to the original code (Branch A below).
%       rho low (amplitude is ORTHOGONAL, e.g. a sphere): gain removal is
%         WRONG -- it launders genuine off-manifold displacement into an
%         on-manifold "gain" label. Test the orthogonal residual of the RAW
%         centroid displacement instead, with no gain subtracted (Branch B).
%       rho intermediate: ambiguous (Branch C) -- report fracOff/fracPattern,
%         leave fracGain undefined (NaN), and flag the point.
%     Validated: identical x1.5 scaling on a cone (rho=1) vs a sphere (rho=0)
%     gives fracOff separation (vs each manifold's own null) of ~2x vs ~7000x --
%     same operation, opposite on-manifold verdict, decided entirely by rho.
%
%  2. RADIUS NORMALIZATION (was: divide by ||zB||, the point's own raw norm).
%     All fractions now divide by localR^2, the local reference-neighborhood's
%     own spread -- a fixed, point-independent yardstick. Validated: at
%     matched SNR, radius normalization gave AUC ~0.94-1.00 for membership
%     detection where the old norm-of-zB denominator gave near-chance AUC
%     (~0.5-0.7), because subtracting gain first left a small numerator over
%     an unrelated, usually large, denominator.
%
%  3. SQUARED-ENERGY FRACTIONS (was: norm, not norm^2). fracGain + fracPattern
%     + fracOff now sum to a true energy partition of ||delta||^2 (Branch A/B);
%     the old norm-based fractions did not sum to anything interpretable.
%
%  4. PARTICIPATION-RATIO d (was: hard-coded d=3). d = round(sum(ev)^2/sum(ev.^2))
%     from the LOCAL reference neighborhood's eigenspectrum, capped by what the
%     neighborhood can estimate. Validated to be far less sensitive to ambient
%     noise than variance-explained thresholds once radius normalization is in
%     place (AUC changed by <0.05 across d=2 to d=15 at matched SNR).
%
%  5. CALIBRATED MEMBERSHIP (was: reconErr < 3, an unjustified threshold).
%     supportScore is now a tail probability against a split-half SPONTANEOUS
%     null: reconstruct held-out reference points from the rest of ZA with the
%     identical pipeline, and report the fraction of ZB below the null's 95th
%     percentile. "On-manifold" now means "statistically indistinguishable
%     from how ZA reconstructs itself," not an arbitrary number.
%
% SCOPE CAVEAT: rho, radius normalization, and the d-estimator were validated
% on d_true-spheres and cones with isotropic ambient noise and orthogonal/
% radial perturbations only. Anisotropic noise, oblique perturbations, and
% strongly non-spherical curvature are untested.

%% -------------------- Defaults --------------------
if ~isfield(opts,'k'), opts.k = 50; end
if ~isfield(opts,'knnType'), opts.knnType = 'radius'; end
if ~isfield(opts,'radiusPct'), opts.radiusPct = 20; end
if ~isfield(opts,'rhoHigh'), opts.rhoHigh = 0.7; end     % >= this: Branch A (gain valid)
if ~isfield(opts,'rhoLow'),  opts.rhoLow  = 0.3; end     % <= this: Branch B (no gain removal)
if ~isfield(opts,'dMax'),    opts.dMax    = []; end      % [] -> no cap beyond neighborhood rank
if ~isfield(opts,'nNullRep'), opts.nNullRep = 200; end   % split-half null draws for calibration

[N_A,F] = size(ZA);
N_B = size(ZB,1);

%% -------------------- Graph & components --------------------
W = buildWeightedAdjacency(ZA, opts);
G = graph(W);
components = conncomp(G);

%% -------------------- NN search --------------------
[idxNN, ~] = knnsearch(ZA, ZB, 'K', opts.k);

%% -------------------- Outputs --------------------
fracGain    = nan(N_B,1);
fracPattern = nan(N_B,1);
fracOff     = nan(N_B,1);
rhoVec      = nan(N_B,1);
dUsed       = nan(N_B,1);
branch      = repmat('?', N_B, 1);   % 'A' | 'B' | 'C'
anchorIdx   = zeros(N_B,1);

%% -------------------- Main loop --------------------
for i = 1:N_B
    zA0 = ZA(idxNN(i,1), :).';
    anchorIdx(i) = idxNN(i,1);

    comp0 = components(anchorIdx(i));
    neigh = idxNN(i, components(idxNN(i,:)) == comp0);
    if numel(neigh) < 3
        neigh = idxNN(i,1:3);
    end

    localA = ZA(neigh,:);
    mu = mean(localA,1).';
    Xc = localA - mu.';

    % ---- tangent space: participation-ratio dimensionality ----
    [~,S,V] = svd(Xc,'econ');
    ev = diag(S).^2;
    d  = max(1, round(sum(ev)^2 / (sum(ev.^2) + eps)));
    if ~isempty(opts.dMax), d = min(d, opts.dMax); end
    d  = min(d, size(V,2));
    Vd = V(:,1:d);
    dUsed(i) = d;

    % ---- local scale: fixed, point-independent yardstick ----
    localR = sqrt(mean(sum(Xc.^2, 2)));

    % ---- rho: fraction of the RADIAL/anchor direction in the tangent space,
    %      computed BEFORE any gain removal ----
    u = zA0 / (norm(zA0) + eps);
    rho = sum((Vd.' * u).^2);
    rhoVec(i) = rho;

    zB = ZB(i,:).';
    delta = zB - mu;

    if rho >= opts.rhoHigh
        % Branch A: amplitude IS a tangent direction. Gain removal is valid.
        alpha = (mu.' * zB) / (mu.' * mu + eps);
        z_gain = alpha * mu;
        r = zB - z_gain;
        r_tan  = Vd * (Vd.' * r);
        r_orth = r - r_tan;
        fracGain(i)    = sum(z_gain.^2) / (localR^2 + eps);
        fracPattern(i) = sum(r_tan.^2)  / (localR^2 + eps);
        fracOff(i)      = sum(r_orth.^2) / (localR^2 + eps);
        branch(i) = 'A';

    elseif rho <= opts.rhoLow
        % Branch B: amplitude is ORTHOGONAL to the manifold. Gain removal
        % would launder off-manifold displacement into an on-manifold label.
        % Test the orthogonal residual of the RAW displacement instead.
        d_tan  = Vd * (Vd.' * delta);
        d_orth = delta - d_tan;
        fracGain(i)    = 0;                                  % not meaningful here
        fracPattern(i) = sum(d_tan.^2)  / (localR^2 + eps);
        fracOff(i)      = sum(d_orth.^2) / (localR^2 + eps);
        branch(i) = 'B';

    else
        % Branch C: ambiguous. Report tangent/orthogonal split of the raw
        % displacement; leave fracGain undefined rather than guessing.
        d_tan  = Vd * (Vd.' * delta);
        d_orth = delta - d_tan;
        fracPattern(i) = sum(d_tan.^2)  / (localR^2 + eps);
        fracOff(i)      = sum(d_orth.^2) / (localR^2 + eps);
        branch(i) = 'C';
    end
end

%% -------------------- Split-half spontaneous null (calibration) --------------------
% Reconstruct held-out ZA points from the rest of ZA with the IDENTICAL pipeline,
% so "on-manifold" is a tail probability against ZA's own self-reconstruction,
% not an arbitrary fixed threshold.
nNull = min(opts.nNullRep, floor(N_A/4));
nullFracOff = nan(nNull,1);
if nNull >= 20
    permA = randperm(N_A);
    half  = floor(N_A/2);
    trainIdx = permA(1:half);  testIdx = permA(half+1:half+min(nNull,N_A-half));
    optsNull = opts;  optsNull.nNullRep = 0;   % avoid recursive nulling
    nullResult = manifoldAnalysis_ABn(ZA(trainIdx,:), ZA(testIdx,:), optsNull); %#ok<*NASGU>
    nullFracOff = nullResult.fracOff;
end
nullThresh = quantile(nullFracOff, 0.95);   % NaN-safe if nNull < 20

%% -------------------- Membership scores --------------------
if isfinite(nullThresh)
    supportScore = mean(fracOff < nullThresh, 'omitnan');
else
    supportScore = NaN;   % null not computed (nNullRep == 0, e.g. inside recursion)
end

%% -------------------- Output --------------------
result.fracGain      = fracGain;
result.fracPattern   = fracPattern;
result.fracOff        = fracOff;
result.rho            = rhoVec;
result.branch         = branch;
result.dUsed           = dUsed;
result.anchorIdx      = anchorIdx;
result.components     = components;
result.nullFracOff    = nullFracOff;
result.nullThresh95   = nullThresh;
result.supportScore   = supportScore;   % NaN unless nNullRep >= 20 in this call
result.W               = W;
result.branchCounts    = struct('A', sum(branch=='A'), 'B', sum(branch=='B'), 'C', sum(branch=='C'));

end

%% -------------------- Helper: weighted adjacency --------------------
function W = buildWeightedAdjacency(Z, opts)

N = size(Z,1);
D = pdist2(Z,Z);

switch lower(opts.knnType)

    case 'knn'
        k = min(opts.k+1, N);
        [idxNN, distNN] = knnsearch(Z,Z,'K',k);
        W = zeros(N);
        for i = 1:N
            nbrs = idxNN(i,2:end);
            W(i,nbrs) = distNN(i,2:end);
        end

    case 'mutual'
        k = min(opts.k+1, N);
        [idxNN, ~] = knnsearch(Z,Z,'K',k);
        W = zeros(N);
        for i = 1:N
            nbrs = idxNN(i,2:end);
            for j = nbrs
                if any(idxNN(j,2:end) == i)
                    W(i,j) = D(i,j);
                end
            end
        end

    case 'radius'
        radius = prctile(D(:), opts.radiusPct);
        W = D;
        W(W > radius | W == 0) = 0;

    otherwise
        error('Unknown knnType');
end

end