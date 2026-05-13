function [out] = estimate_basin_geometry(ses)
% ESTIMATE_BASIN_GEOMETRY  Relate basin sharpness to geometric sensitivity
%   and occupancy redistribution using a median-split approach.
%
% For each session x fold, clusters are split at the median of sharpness.
% Effect sizes compare the outcome (dTS, occ_ratio, etc.) between halves:
%   effect = mean(outcome | high sharpness) - mean(outcome | low sharpness)
%
% For dTS and spread: negative effect expected (sharp → less disrupted)
% For occ:           positive effect expected  (sharp → more visited by stim)
% For baseline occ:  near-zero expected        (sharp irrelevant spontaneously)
%
% Session-level effects (median across folds) are tested by sign test
% and t-test across sessions.

nSes  = length(ses);
nFold = length(ses(1).tensor.block);

fields = {'sharp_dTS','sharp_occ','sharp_base', ...
          'occ_dTS',  ...       % median split on sharpness, outcome = dTS
          'sharp_spont_spread','sharp_stim_spread'};

% One value per session
effect_ses = struct();
for f = 1:length(fields)
    effect_ses.(fields{f}) = nan(nSes,1);
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
    Chat   = ses(i).Ztensor;   % use stored rate-free vectors

    % Per-fold buffers
    fb = struct();
    for f = 1:length(fields)
        fb.(fields{f}) = nan(nFold,1);
    end

    for j = 1:nFold

        blk      = ses(i).tensor.block(j);
        tf       = blk.tensor_flow;
        kmean_id = blk.kmean_id(:);
        train_ix = blk.training_ix(:);
        stim_ix  = blk.stim_ix(:);
        test_ix  = blk.test_ix(:);
        motifs   = tf.spont.motifID(:);
        nM       = length(motifs);

        if nM < 4, continue; end   % need ≥2 per half

        % ── Per-cluster basin geometry from training ripples ───────────────
        width      = nan(nM,1);
        centroids  = nan(nM, size(Chat,2));
        N_spont_cl = zeros(nM,1);
        N_stim_cl  = zeros(nM,1);
        N_test_cl  = zeros(nM,1);

        for m = 1:nM
            cl     = motifs(m);
            sp_idx = train_ix(kmean_id(train_ix) == cl);
            st_idx = stim_ix( kmean_id(stim_ix)  == cl);
            te_idx = test_ix( kmean_id(test_ix)  == cl);
            N_spont_cl(m) = length(sp_idx);
            N_stim_cl(m)  = length(st_idx);
            N_test_cl(m)  = length(te_idx);

            if N_spont_cl(m) < 5, continue; end
            X  = Chat(sp_idx,:);
            mu = mean(X,1); mu = mu/(norm(mu)+1e-12);
            centroids(m,:) = mu;
            S  = max(min(X*X',1),-1);
            D  = 1 - S;
            ut = D(triu(true(size(D)),1));
            width(m) = mean(ut);
        end

        % Separation and sharpness
        valid_c = find(~isnan(width) & ~any(isnan(centroids),2));
        if length(valid_c) < 4, continue; end

        separation = nan(nM,1);
        C_valid    = centroids(valid_c,:);
        for mi = 1:length(valid_c)
            m      = valid_c(mi);
            others = setdiff(1:length(valid_c),mi);
            dists  = 1 - C_valid(mi,:)*C_valid(others,:)';
            separation(m) = min(dists);
        end
        sharpness = separation ./ (width + 1e-12);

        % ── Occupancy ratios ───────────────────────────────────────────────
        p_stim  = N_stim_cl  / (sum(N_stim_cl)  + 1e-12);
        p_spont = N_spont_cl / (sum(N_spont_cl) + 1e-12);
        p_test  = N_test_cl  / (sum(N_test_cl)  + 1e-12);
        occ_stim = p_stim ./ (p_spont + 1e-12);
        occ_base = p_test  ./ (p_spont + 1e-12);

        % ── delta_totalSpread ──────────────────────────────────────────────
        dTS          = tf.stim.totalSpread(:) - tf.spont.totalSpread(:);
        spont_spread = tf.spont.totalSpread(:);
        stim_spread  = tf.stim.totalSpread(:);

        % ── Median split: sharpness as predictor ──────────────────────────
        fb.sharp_dTS(j)          = msplit_geo(sharpness, dTS,          N_spont_cl, N_stim_cl, 10, 2);
        fb.sharp_occ(j)          = msplit_geo(sharpness, occ_stim,     N_spont_cl, N_stim_cl, 10, 2);
        fb.sharp_base(j)         = msplit_geo(sharpness, occ_base,     N_spont_cl, N_test_cl,  10, 2);
        fb.sharp_spont_spread(j) = msplit_geo(sharpness, spont_spread, N_spont_cl, N_stim_cl,  5, 0);
        fb.sharp_stim_spread(j)  = msplit_geo(sharpness, stim_spread,  N_spont_cl, N_stim_cl,  5, 2);

        % ── Median split: dTS as predictor, occ as outcome ────────────────
        fb.occ_dTS(j) = msplit_geo(dTS, occ_stim, N_spont_cl, N_stim_cl, 5, 2);

    end  % fold loop

    out.session(i) = fb;

    % Aggregate folds → session median
    for f = 1:length(fields)
        fv = fb.(fields{f});
        fv = fv(~isnan(fv));
        if length(fv) >= 3
            effect_ses.(fields{f})(i) = median(fv);
        end
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
    'sharp_dTS',          'Sharp → Δspread  (neg=sharp less disrupted)'; ...
    'sharp_occ',          'Sharp → stim occ (pos=sharp more visited)';   ...
    'sharp_base',         'Sharp → base occ (neg=AAC-specific if ~0)';   ...
    'occ_dTS',            'High ΔTS → occ   (neg=disrupted avoided)';    ...
    'sharp_spont_spread', 'Sharp → spont spread (neg=sharp slower)';      ...
    'sharp_stim_spread',  'Sharp → stim spread  (neg=sharp slower stim)'};

fprintf('\n%s\n', repmat('=',1,72));
fprintf('Basin geometry median-split — N sessions: %d\n', nSes);
fprintf('Effect = mean(outcome | high predictor) - mean(outcome | low predictor)\n');
fprintf('%s\n', repmat('-',1,72));
fprintf('%-45s  %+8s  %8s  %7s  %7s\n','Analysis','median','IQR','p(sign)','p(ttest)');
fprintf('%s\n', repmat('-',1,72));

for k = 1:length(labels)
    fld = labels{k,1};
    ev  = effect_ses.(fld);
    [med_e,iqr_e,p_sign,p_tt] = summarise(ev);
    fprintf('%-45s  %+8.4f  %8.4f  %7.4f  %7.4f\n', ...
        labels{k,2}, med_e, iqr_e, p_sign, p_tt);
end
fprintf('%s\n', repmat('=',1,72));

out.effect_ses = effect_ses;
for k = 1:length(labels)
    fld = labels{k,1};
    ev  = effect_ses.(fld);
    [med_e,iqr_e,p_sign,p_tt] = summarise(ev);
    out.(fld) = struct('effect_ses',ev,'median',med_e,'iqr',iqr_e, ...
                       'p_sign',p_sign,'p_ttest',p_tt);
end

end

% ── Local helpers ─────────────────────────────────────────────────────────

function eff = msplit_geo(predictor, outcome, N_sp, N_st, min_sp, min_st)
% Median split: effect = mean(outcome|high predictor) - mean(outcome|low)
    eff   = nan;
    valid = ~isnan(predictor) & ~isnan(outcome) & ...
            N_sp >= min_sp & N_st >= min_st;
    if sum(valid) < 4, return; end
    p   = predictor(valid);
    y   = outcome(valid);
    med = median(p);
    hi  = p >= med;
    lo  = ~hi;
    if sum(hi) < 1 || sum(lo) < 1, return; end
    eff = mean(y(hi)) - mean(y(lo));
end