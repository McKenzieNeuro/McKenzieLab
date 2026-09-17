function S = precompute_sessions(ses, nFold, fold_use, min_sp, min_st)
%  Run Ztensor + UMAP once per session, cache everything the GUI needs

nSes = length(ses);
S    = cell(nSes,1);

for i = 1:nSes
    fprintf('  Session %d/%d...', i, nSes);

    % Find fold
    f = fold_use;
    while f <= nFold && ~isfield(ses(i).tensor.block(f),'tensor_flow')
        f = f + 1;
    end
    if f > nFold, fprintf(' skip (no tensor_flow)\n'); continue; end

    tf       = ses(i).tensor.block(f).tensor_flow;
    stim_ix  = ses(i).tensor.block(f).stim_ix;
    clust_id = ses(i).tensor.block(f).kmean_id(:);
    stim_vec = ismember(1:length(clust_id), stim_ix)';

    X_raw = double(ses(i).binnedPopRipple_tensor);
    A_mat = ses(i).tensorModel.A;
    [~,T,R] = size(X_raw);

    Zfull = zeros(size(A_mat,2),T,R);
    for r = 1:R
        Zfull(:,:,r) = A_mat' * squeeze(X_raw(:,:,r));
    end

    spont_all = ~stim_vec;

    % ── Per-ripple mean firing rate ───────────────────────────────────────
    % fbar_r(r) = mean spike count across all neurons and timebins for ripple r
    % This is the scalar we regress out — each ripple's own activity level.
    fbar_r = squeeze(mean(mean(X_raw, 1), 2));   % R x 1

    % ── Motif identity embedding (Zid) ────────────────────────────────────
    C  = ses(i).tensorModel.C;
    Cz = zscore(C,[],1);
    Cn = Cz ./ (vecnorm(Cz,2,2)+1e-12);          % R x K, unit rows

    % Regress per-ripple firing rate out of Cn (OLS across ripples)
    % beta(k) = regression coefficient of component k on fbar_r
    fbar_c  = fbar_r - mean(fbar_r);              % mean-centre predictor
    beta    = (Cn' * fbar_c) / (fbar_c' * fbar_c + 1e-12);  % K x 1
    Cn_rf   = Cn - fbar_c * beta';                % R x K, rate-regressed
    Zid     = Cn_rf ./ (vecnorm(Cn_rf,2,2)+1e-12);           % R x K, renormed

    % ── Per-timepoint normalization ───────────────────────────────────────
    % Apply the same per-ripple rate regression to each timepoint of Zfull,
    % using the same beta so the space is consistent with Zid / UMAP.
    % Each timepoint is first standardized using Cz column statistics,
    % then unit-normed, then the per-ripple rate contribution is removed
    % using beta scaled by fbar_c(r), then renormed.
    K_dim = size(Zfull,1);
    Cz_mu = mean(Cz,1)';    % K x 1
    Cz_sd = std(Cz,0,1)';   % K x 1
    Znorm = zeros(K_dim,T,R);
    for r = 1:R
        for t = 1:T
            z   = Zfull(:,t,r);
            zz  = (z - Cz_mu) ./ (Cz_sd + 1e-12);   % standardize
            zn  = zz / (norm(zz)+1e-12);              % unit norm
            zrf = zn - fbar_c(r) * beta;              % remove per-ripple rate
            zrn = zrf / (norm(zrf)+1e-12);             % renorm
            Znorm(:,t,r) = zrn;
        end
    end

    % UMAP on Zid (trial-mean, rate-free) — clean cluster separation
    tmpl = sprintf('umap_template_ses%d.mat',i);
    [usp_all,~] = run_umap(Zid(spont_all,:),...
        'n_neighbors',15,'min_dist',0.10,'metric','cosine',...
        'n_components',2,'randomize',false,...
        'save_template_file',tmpl,'verbose','none');
    ust_all = zeros(0,2);
    if any(stim_vec)
        ust_all = run_umap(Zid(stim_vec,:),...
            'template_file',tmpl,'verbose','none');
    end
    umap_all              = zeros(R,2);
    umap_all(spont_all,:) = usp_all;
    umap_all(stim_vec, :) = ust_all;

    % Valid motifs
    motifs_sp = tf.spont.motifID(:);
    rf_sp_all = tf.spont.rotFracTotal(:);
    motifs_st = tf.stim.motifID(:);

    valid = false(length(motifs_sp),1);
    for m = 1:length(motifs_sp)
        cl = motifs_sp(m);
        valid(m) = sum(clust_id==cl & ~stim_vec)>=min_sp && ...
                   sum(clust_id==cl &  stim_vec)>=min_st && ...
                   ismember(cl,motifs_st);
    end
    if sum(valid) < 2, fprintf(' skip (<%d valid motifs)\n',2); continue; end

    mot_valid = motifs_sp(valid);
    rf_valid  = rf_sp_all(valid);
    nV        = length(mot_valid);

    % Sort by rotFrac descending for display
    [rf_valid, sord] = sort(rf_valid,'descend');
    mot_valid = mot_valid(sord);

    % Default top 3 by joint score
    umap_centers = zeros(nV,2);
    umap_spread  = zeros(nV,1);
    for m = 1:nV
        cl   = mot_valid(m);
        kp   = clust_id==cl & ~stim_vec;
        upts = umap_all(kp,:);
        umap_centers(m,:) = mean(upts,1);
        umap_spread(m)    = mean(std(upts,0,1))+1e-12;
    end
    rf_norm = (rf_valid-min(rf_valid))/(max(rf_valid)-min(rf_valid)+1e-12);
    combos  = nchoosek(1:nV,min(3,nV));
    jscores = zeros(size(combos,1),1);
    for ci = 1:size(combos,1)
        sub = combos(ci,:);
        md  = inf;
        for m1=1:length(sub)
            for m2=m1+1:length(sub)
                md = min(md, norm(umap_centers(sub(m1),:)-umap_centers(sub(m2),:)));
            end
        end
        jscores(ci) = 0.20*(md/mean(umap_spread(sub))) + 0.55*mean(rf_norm(sub));
    end
    [~,bci]   = max(jscores);
    top_subset = combos(bci,:);
    top_mot    = mot_valid(top_subset);

    D.clust_id    = clust_id;
    D.stim_vec    = stim_vec;
    D.umap_all    = umap_all;
    D.Znorm       = Znorm;
    D.T           = T;
    D.R           = R;
    D.mot_valid   = mot_valid;
    D.rf_valid    = rf_valid;
    D.top_mot     = top_mot;
    D.umap_spread = umap_spread;
    S{i} = D;

    fprintf(' %d motifs, top=[%s]\n', nV, num2str(top_mot'));
end
end
