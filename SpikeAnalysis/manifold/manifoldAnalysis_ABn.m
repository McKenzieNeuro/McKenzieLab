function result = manifoldAnalysis_ABn(ZA, ZB, opts)
% Component-aware, gain-correct manifold membership test
% Adds fracOnManifold = fraction of vector aligned with manifold (gain + tangent)

%% -------------------- Defaults --------------------
if ~isfield(opts,'k'), opts.k = 50; end
if ~isfield(opts,'knnType'), opts.knnType = 'radius'; end
if ~isfield(opts,'radiusPct'), opts.radiusPct = 20; end

[N_A,F] = size(ZA);
N_B = size(ZB,1);

%% -------------------- Graph & components --------------------
W = buildWeightedAdjacency(ZA, opts);
G = graph(W);
components = conncomp(G);
graphConnected = (max(components) == 1);

%% -------------------- NN search --------------------
[idxNN, ~] = knnsearch(ZA, ZB, 'K', opts.k);

% global scale for normalization
D = pdist2(ZA,ZA);
sigma_self = median(D(:)) + eps;

%% -------------------- Outputs --------------------
fracGain        = zeros(N_B,1);
fracPattern     = zeros(N_B,1);
fracOff         = zeros(N_B,1);
fracOnManifold  = zeros(N_B,1);  % NEW
reconErrB       = zeros(N_B,1);
anchorIdx       = zeros(N_B,1);

%% -------------------- Main loop --------------------
for i = 1:N_B

    % ---- anchor point ----
    zA0 = ZA(idxNN(i,1), :).';
    anchorIdx(i) = idxNN(i,1);

    % ---- component-aware neighborhood ----
    comp0 = components(anchorIdx(i));
    neigh = idxNN(i, components(idxNN(i,:)) == comp0);
    if numel(neigh) < 3
        neigh = idxNN(i,1:3); % fallback
    end

    localA = ZA(neigh,:);
    mu = mean(localA,1).';
    Xc = localA - mu.';

    % ---- tangent space ----
    [~,~,V] = svd(Xc,'econ');
    %d = min(rank(Xc), size(V,2));
    d = 3;
    Vd = V(:,1:d);

    % ---- gain estimation relative to anchor ----
    zB = ZB(i,:).';
    %alpha = (zA0' * zB) / (zA0' * zA0 + eps);
    %z_gain = alpha * zA0;
    % Gain relative to neighborhood centroid
    alpha = (mu' * zB) / (mu' * mu + eps);  
    z_gain = alpha * mu;

    % ---- residual after gain ----
    r = zB - z_gain;

    % ---- tangent / orthogonal split ----
    r_tan  = Vd * (Vd.' * r);
    r_orth = r - r_tan;

    % ---- fractions ----
    nv = norm(zB) + eps;
    fracGain(i)    = norm(z_gain) / nv;
    fracPattern(i) = norm(r_tan)  / nv;
    fracOff(i)     = norm(r_orth)/ nv;

    % ---- NEW: combined on-manifold fraction ----
    % Includes both gain + tangent relative to total vector
    fracOnManifold(i) = (norm(z_gain) + norm(r_tan)) / nv;

    % ---- reconstruction error ----
    reconErrB(i) = norm(r_orth) / sigma_self;
end

%% -------------------- Membership scores --------------------
membershipLocal = exp(-(reconErrB).^2);
membershipScore = mean(membershipLocal);

%% -------------------- Output --------------------
result.reconErrB       = reconErrB;
result.fracGain        = fracGain;
result.fracPattern     = fracPattern;
result.fracOff         = fracOff;
result.fracOnManifold  = fracOnManifold;  % NEW
result.membershipLocal = membershipLocal;
result.membershipScore = membershipScore;
result.supportScore    = mean(reconErrB < 3);
result.anchorIdx       = anchorIdx;
%result.graphConnected  = graphConnected;
result.components      = components;
result.W               = W;

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
