function results = testOffManifoldScaling(ZA, ZB, opts)
% Test whether off-manifold deviations are simple scalings of tangent patterns
%
% ZA : N_A x F   (reference manifold)
% ZB : N_B x F   (test points)
% opts.k        : neighbors (default 50)
% opts.nPerm    : permutations (default 100)
%
% OUTPUT: results struct with quantitative evidence

%% ---------------- Defaults ----------------
if ~isfield(opts,'k'), opts.k = 50; end
if ~isfield(opts,'nPerm'), opts.nPerm = 100; end

[N_B, F] = size(ZB);

%% ---------------- Neighbors ----------------
[idxNN,~] = knnsearch(ZA, ZB, 'K', opts.k);

%% ---------------- Allocate ----------------
cosineVT_VP   = zeros(N_B,1);
R2_scaling    = zeros(N_B,1);
alpha_hat     = zeros(N_B,1);
perm_R2       = zeros(N_B, opts.nPerm);

%% ---------------- Main loop ----------------
for i = 1:N_B

    %% --- Local manifold ---
    localA = ZA(idxNN(i,:), :);
    mu = mean(localA,1).';
    Xc = localA - mu.';
    [~,~,V] = svd(Xc,'econ');
    d = min(rank(Xc), opts.k-1);
    T = V(:,1:d);             % F x d

    %% --- Deviation vector ---
    v = (ZB(i,:) - mu.').';
    v = v(:);                 % enforce column

    %% --- Gain (radial) direction ---
    g = mu / (norm(mu)+eps);
    v_gain = g * (g.' * v);

    %% --- Tangent & off-manifold ---
    v_res = v - v_gain;
    v_tan = T * (T.' * v_res);
    v_off = v_res - v_tan;

    %% --- Skip degenerate cases ---
    if norm(v_tan) < 1e-10 || norm(v_off) < 1e-10
        continue
    end

    %% --- Hypothesis tests ---
    % 1. Cosine similarity
    cosineVT_VP(i) = dot(v_tan, v_off) / ...
        (norm(v_tan)*norm(v_off));

    % 2. Linear scaling regression
    alpha_hat(i) = (v_tan.' * v_off) / (v_tan.' * v_tan);
    v_pred = alpha_hat(i) * v_tan;

    R2_scaling(i) = 1 - ...
        norm(v_off - v_pred)^2 / norm(v_off)^2;

    % 3. Permutation control
    for p = 1:opts.nPerm
        idx = randperm(F);
        v_perm = v_tan(idx);
        a_p = (v_perm.' * v_off) / (v_perm.' * v_perm);
        v_p = a_p * v_perm;
        perm_R2(i,p) = 1 - norm(v_off - v_p)^2 / norm(v_off)^2;
    end
end

%% ---------------- Aggregate evidence ----------------
results.meanCosine     = mean(cosineVT_VP,'omitnan');
results.meanR2         = mean(R2_scaling,'omitnan');
results.meanAlpha      = mean(alpha_hat,'omitnan');

results.permMeanR2     = mean(perm_R2(:),'omitnan');
results.p_value        = mean( ...
    mean(perm_R2,2,'omitnan') >= R2_scaling, ...
    'omitnan');

%% ---------------- Interpretation ----------------
results.supportsScaling = ...
    results.meanR2 > 0.5 && ...
    results.meanCosine > 0.7 && ...
    results.p_value < 0.05;

end
