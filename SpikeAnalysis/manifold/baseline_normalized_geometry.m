function [out] = baseline_normalized_geometry(ses)
% BASELINE_NORMALIZED_GEOMETRY  Relates baseline manifold geometry to
%   stim-induced changes using a median-split approach.
%
% For each session x fold, clusters are split at the median of the baseline
% predictor (width or sharpness). The effect size is:
%   effect = mean(delta | high predictor) - mean(delta | low predictor)
%
% Delta is regression-to-the-mean corrected:
%   d_spread = (spread_stim - spread_test) / spread_train
%   d_volume = (vol_stim    - vol_test)    / vol_train
%   d_nn     = (nn_stim     - nn_test)     / nn_train
%   d_order  =  order_stim  - order_test
%
% Predicted directions:
%   width → Δspread:   positive (broad basins → more acceleration)
%   sharp → Δspread:   negative (sharp basins → less acceleration)
%   base_spread → Δspread: negative (fast at baseline → smaller normalized gain)
%
% Session-level effects (median across folds) tested by sign test and t-test.

nSes  = length(ses);
nFold = length(ses(1).tensor.block);

fields = {'width_spread','width_volume','width_nn','width_order', ...
          'sharp_spread','sharp_volume','sharp_nn','sharp_order', ...
          'base_spread_d_spread'};

effect_ses = struct();
for f = 1:length(fields)
    effect_ses.(fields{f})         = nan(nSes,1);
    effect_ses.([fields{f} '_hi']) = nan(nSes,1);
    effect_ses.([fields{f} '_lo']) = nan(nSes,1);
end

for i = 1:nSes

    C      = ses(i).tensorModel.C;
    Cz     = zscore(C, [], 1);
    Cn     = Cz ./ (vecnorm(Cz, 2, 2) + 1e-12);
    X_raw  = double(ses(i).binnedPopRipple_tensor);
    uRate  = squeeze(sum(sum(X_raw, 1), 2));
    v_rate = (uRate' * Cn)';
    v_rate = v_rate / (norm(v_rate) + 1e-12);
    Cn_rf  = Cn - (Cn * v_rate) * v_rate';
    Chat   = Cn_rf ./ (vecnorm(Cn_rf, 2, 2) + 1e-12);

    A      = ses(i).tensorModel.A;
    [~, T, R] = size(X_raw);
    nComp  = size(A, 2);
    Cz_mu  = mean(Cz, 1)';
    Cz_sd  = std(Cz, 0, 1)';
    Znorm  = zeros(nComp, T, R);
    for r = 1:R
        for t = 1:T
            z   = A' * double(X_raw(:,t,r));
            zz  = (z - Cz_mu) ./ (Cz_sd + 1e-12);
            zn  = zz / (norm(zz) + 1e-12);
            zrf = zn - (zn' * v_rate) * v_rate;
            zrn = zrf / (norm(zrf) + 1e-12);
            Znorm(:,t,r) = zrn(1:nComp);
        end
    end

    fb = struct();
    for f = 1:length(fields)
        fb.(fields{f})         = nan(nFold,1);
        fb.([fields{f} '_hi']) = nan(nFold,1);
        fb.([fields{f} '_lo']) = nan(nFold,1);
    end
    out.session(i) = fb;

    for j = 1:nFold

        blk      = ses(i).tensor.block(j);
        kmean_id = blk.kmean_id(:);
        train_ix = blk.training_ix(:);
        test_ix  = blk.test_ix(:);
        stim_ix  = blk.stim_ix(:);
        motifs   = blk.tensor_flow.spont.motifID(:);
        nM       = length(motifs);
        K_clust  = blk.optimal_k;

        if nM < 4, continue; end

        width      = nan(nM,1);
        separation = nan(nM,1);
        centroids  = nan(nM, size(Chat,2));
        spread_train= nan(nM,1);
        spread_test = nan(nM,1);
        spread_stim = nan(nM,1);
        vol_test    = nan(nM,1);
        vol_stim    = nan(nM,1);
        nn_train    = nan(nM,1);
        nn_test     = nan(nM,1);
        nn_stim     = nan(nM,1);
        ord_test    = nan(nM,1);
        ord_stim    = nan(nM,1);
        N_train_cl  = zeros(nM,1);
        N_test_cl   = zeros(nM,1);
        N_stim_cl   = zeros(nM,1);

        for m = 1:nM
            cl     = motifs(m);
            tr_idx = train_ix(kmean_id(train_ix) == cl);
            te_idx = test_ix( kmean_id(test_ix)  == cl);
            st_idx = stim_ix( kmean_id(stim_ix)  == cl);
            N_train_cl(m) = length(tr_idx);
            N_test_cl(m)  = length(te_idx);
            N_stim_cl(m)  = length(st_idx);

            if N_train_cl(m) < 5, continue; end

            X   = Chat(tr_idx,:);
            mu  = mean(X,1); mu = mu/(norm(mu)+1e-12);
            centroids(m,:) = mu;
            S   = max(min(X*X',1),-1); D = 1-S;
            ut  = D(triu(true(size(D)),1));
            width(m) = mean(ut);

            [U2,ok] = estimate_jpca_plane(Znorm, tr_idx, T, nComp);
            if ~ok, continue; end
            spread_train(m) = compute_spread(Znorm, tr_idx, U2, T);
            if N_test_cl(m) >= 3
                spread_test(m) = compute_spread(Znorm, te_idx, U2, T);
            end
            if N_stim_cl(m) >= 2
                spread_stim(m) = compute_spread(Znorm, st_idx, U2, T);
            end

            if cl >= 1 && cl <= K_clust
                vol_test(m) = blk.jacobian_volume_spont(cl);
                vol_stim(m) = blk.jacobian_volume_stim(cl);
            end

            if N_train_cl(m) >= 5
                Xtr  = Chat(tr_idx,:);
                D_tr = 1 - max(min(Xtr*Xtr',1),-1);
                D_tr(1:N_train_cl(m)+1:end) = Inf;
                nn_train(m) = mean(min(D_tr,[],2));
            end

            nn_all_test = blk.nearest_neighbor_global_spont;
            nn_all_stim = blk.nearest_neighbor_global_stim;
            te_pos = arrayfun(@(x) find(test_ix==x,1), te_idx,'UniformOutput',true);
            st_pos = arrayfun(@(x) find(stim_ix==x,1), st_idx,'UniformOutput',true);
            if ~isempty(te_pos), nn_test(m) = mean(nn_all_test(te_pos)); end
            if ~isempty(st_pos), nn_stim(m) = mean(nn_all_stim(st_pos)); end

            if cl >= 1 && cl <= K_clust
                ord_test(m) = blk.order_corr_spont(cl);
                ord_stim(m) = blk.order_corr_stim(cl);
            end
        end

        % Sharpness
        valid_c = find(~isnan(width) & ~any(isnan(centroids),2));
        if length(valid_c) < 4, continue; end
        C_valid = centroids(valid_c,:);
        for mi = 1:length(valid_c)
            m      = valid_c(mi);
            others = setdiff(1:length(valid_c),mi);
            dists  = 1 - C_valid(mi,:)*C_valid(others,:)';
            separation(m) = min(dists);
        end
        sharpness = separation ./ (width + 1e-12);

        % Normalized deltas
        d_spread = (spread_stim - spread_test) ./ (spread_train + 1e-12);
        d_volume = (vol_stim    - vol_test)    ./ 1;   % vol_train ≈ 1
        d_nn     = (nn_stim     - nn_test)     ./ (nn_train + 1e-12);
        d_order  =  ord_stim    - ord_test;

        pairs = { ...
            'width_spread',        width,        d_spread; ...
            'width_volume',        width,        d_volume; ...
            'width_nn',            width,        d_nn;     ...
            'width_order',         width,        d_order;  ...
            'sharp_spread',        sharpness,    d_spread; ...
            'sharp_volume',        sharpness,    d_volume; ...
            'sharp_nn',            sharpness,    d_nn;     ...
            'sharp_order',         sharpness,    d_order;  ...
            'base_spread_d_spread',spread_train, d_spread};

        for p = 1:size(pairs,1)
            tag = pairs{p,1};
            [eff, hi_m, lo_m] = msplit_bn(pairs{p,2}, pairs{p,3}, ...
                                           N_train_cl, N_test_cl, N_stim_cl);
            fb.(tag)(j)         = eff;
            fb.([tag '_hi'])(j) = hi_m;
            fb.([tag '_lo'])(j) = lo_m;
        end

    end  % fold loop

    for f = 1:length(fields)
        out.session(i).(fields{f}) = fb.(fields{f});
        out.session(i).([fields{f} '_hi']) = fb.([fields{f} '_hi']);
        out.session(i).([fields{f} '_lo']) = fb.([fields{f} '_lo']);
    end

    for f = 1:length(fields)
        fv = fb.(fields{f});
        fv = fv(~isnan(fv));
        if length(fv) >= 3
            effect_ses.(fields{f})(i) = median(fv);
        end
        fh = fb.([fields{f} '_hi']); fh = fh(~isnan(fh));
        fl = fb.([fields{f} '_lo']); fl = fl(~isnan(fl));
        if length(fh) >= 3, effect_ses.([fields{f} '_hi'])(i) = median(fh); end
        if length(fl) >= 3, effect_ses.([fields{f} '_lo'])(i) = median(fl); end
    end

end  % session loop

% ── Summary ───────────────────────────────────────────────────────────────
function [med_e, iqr_e, p_sign, p_tt] = summarise(ev)
    ev = ev(~isnan(ev));
    if length(ev) < 3
        med_e=nan; iqr_e=nan; p_sign=nan; p_tt=nan; return;
    end
    med_e  = median(ev);
    iqr_e  = iqr(ev);
    n_neg  = sum(ev < 0); n_pos = sum(ev > 0);
    p_sign = 2*binocdf(min(n_neg,n_pos), n_neg+n_pos, 0.5);
    [~,p_tt] = ttest(ev);
end

labels = { ...
    'width_spread',        'Width  → Δspread  (pos=broad→more accel)';   ...
    'width_volume',        'Width  → Δvolume  (pos=broad→more expand)';  ...
    'width_nn',            'Width  → ΔNN      (pos=broad→further)';       ...
    'width_order',         'Width  → Δorder   (neg=broad→more degrad)';  ...
    'sharp_spread',        'Sharp  → Δspread  (neg=sharp→less accel)';   ...
    'sharp_volume',        'Sharp  → Δvolume  (neg=sharp→less expand)';  ...
    'sharp_nn',            'Sharp  → ΔNN      (neg=sharp→less distant)'; ...
    'sharp_order',         'Sharp  → Δorder   (pos=sharp→less degrad)';  ...
    'base_spread_d_spread','Baseline spread → Δspread (neg=fast→less gain)'};

fprintf('\n%s\n', repmat('=',1,72));
fprintf('Baseline geometry median-split — N sessions: %d\n', nSes);
fprintf('delta = (stim - test) / train\n');
fprintf('%s\n', repmat('-',1,72));
fprintf('%-45s  %+8s  %8s  %7s  %7s\n','Analysis','median','IQR','p(sign)','p(ttest)');
fprintf('%s\n', repmat('-',1,72));

for k = 1:length(labels)
    fld = labels{k,1};
    ev  = effect_ses.(fld);
    hi  = effect_ses.([fld '_hi']);
    lo  = effect_ses.([fld '_lo']);
    [med_e,iqr_e,p_sign,p_tt] = summarise(ev);
    hi_m = median(hi(~isnan(hi)));
    lo_m = median(lo(~isnan(lo)));
    fprintf('%-45s  %+8.4f  %8.4f  %7.4f  %7.4f  [hi=%.3f lo=%.3f]\n', ...
        labels{k,2}, med_e, iqr_e, p_sign, p_tt, hi_m, lo_m);
end
fprintf('%s\n', repmat('=',1,72));

out.effect_ses = effect_ses;
out.labels     = labels;
for k = 1:length(labels)
    fld = labels{k,1};
    ev  = effect_ses.(fld);
    hi  = effect_ses.([fld '_hi']);
    lo  = effect_ses.([fld '_lo']);
    [med_e,iqr_e,p_sign,p_tt] = summarise(ev);
    out.(fld) = struct('effect_ses',ev,'median',med_e,'iqr',iqr_e, ...
                       'p_sign',p_sign,'p_ttest',p_tt, ...
                       'hi_ses',hi,'lo_ses',lo, ...
                       'hi_median',median(hi(~isnan(hi))), ...
                       'lo_median',median(lo(~isnan(lo))));
end

end  % main

% =========================================================================
function [U2, ok] = estimate_jpca_plane(Znorm, idx, T, D)
    ok = false; U2 = [];
    if length(idx) < 5, return; end
    Rk = length(idx); Xp = []; Vp = [];
    for r = 1:Rk
        tr = Znorm(:,:,idx(r));
        xp = tr(:,1); vp = tr(:,2)-tr(:,1);
        for t = 2:T-1
            xp = [xp, tr(:,t)];
            vp = [vp, (tr(:,t+1)-tr(:,t-1))/2];
        end
        xp = [xp, tr(:,T)]; vp = [vp, tr(:,T)-tr(:,T-1)];
        Xp = [Xp, xp]; Vp = [Vp, vp];
    end
    M  = (Vp*Xp')/(Xp*Xp'+1e-6*eye(D));
    Ms = (M-M')/2;
    [Ve,Ev] = eig(Ms);
    [~,si]  = sort(abs(imag(diag(Ev))),'descend');
    v1 = real(Ve(:,si(1))); v2 = imag(Ve(:,si(1)));
    v1 = v1/(norm(v1)+1e-12);
    v2 = v2-(v2'*v1)*v1; v2 = v2/(norm(v2)+1e-12);
    U2 = [v1,v2]; ok = true;
end

function phi = compute_spread(Znorm, idx, U2, T)
    if isempty(idx), phi = nan; return; end
    Rk = length(idx); phi_t = nan(T-2,1);
    for t = 2:T-1
        Xt = zeros(Rk,2); Vt = zeros(Rk,2);
        for r = 1:Rk
            tr = Znorm(:,:,idx(r));
            Xt(r,:) = (U2'*tr(:,t))';
            Vt(r,:) = (U2'*(tr(:,t+1)-tr(:,t-1))/2)';
        end
        Xt = Xt-mean(Xt,1); Vt = Vt-mean(Vt,1);
        lambda = max(1e-3*trace(Xt'*Xt)/2, 1e-6);
        B = (Vt'*Xt)/(Xt'*Xt+lambda*eye(2));
        A = (B-B')/2; S = (B+B')/2;
        AX = (A*Xt')'; SX = (S*Xt')';
        phi_t(t-1) = mean(sum(AX.^2,2)+sum(SX.^2,2));
    end
    phi = nanmean(phi_t);
end

function [eff, hi_m, lo_m] = msplit_bn(predictor, outcome, N_train, N_test, N_stim)
% Median split with fixed validity criteria: N_train>=5, N_test>=3, N_stim>=2
    eff = nan; hi_m = nan; lo_m = nan;
    valid = ~isnan(predictor) & ~isnan(outcome) & ...
            N_train >= 5 & N_test >= 3 & N_stim >= 2;
    if sum(valid) < 4, return; end
    p   = predictor(valid);
    y   = outcome(valid);
    med = median(p);
    hi  = p >= med;
    lo  = ~hi;
    if sum(hi) < 1 || sum(lo) < 1, return; end
    hi_m = mean(y(hi));
    lo_m = mean(y(lo));
    eff  = hi_m - lo_m;
end