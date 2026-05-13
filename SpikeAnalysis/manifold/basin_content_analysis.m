function [out] = basin_content_analysis(ses)
% BASIN_CONTENT_ANALYSIS  Relates basin sharpness to sequence content
%   and within-cluster variance using a median-split approach.
%
% For each session x fold, clusters are split at the median of sharpness.
% Effect = mean(outcome | high sharpness) - mean(outcome | low sharpness)
%
% Predicted directions:
%   sharp_order_spont:  positive (sharp → more faithful sequences)
%   sharp_order_stim:   positive
%   sharp_delta_order:  positive (sharp → less degraded)
%   sharp_var_spont:    negative (sharp → tighter spont cloud)
%   sharp_var_stim:     negative
%   sharp_delta_var:    negative (sharp → less variance increase)
%
% Session-level effects (median across folds) are tested by sign test
% and t-test across sessions.

nSes  = length(ses);
nFold = length(ses(1).tensor.block);

fields = {'sharp_order_spont','sharp_order_stim','sharp_delta_order', ...
          'sharp_var_spont','sharp_var_stim','sharp_delta_var'};

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

    fb = struct();
    for f = 1:length(fields)
        fb.(fields{f}) = nan(nFold,1);
    end

    for j = 1:nFold

        blk      = ses(i).tensor.block(j);
        kmean_id = blk.kmean_id(:);
        train_ix = blk.training_ix(:);
        stim_ix  = blk.stim_ix(:);
        motifs   = blk.tensor_flow.spont.motifID(:);
        nM       = length(motifs);

        if nM < 4, continue; end

        try
            ord_sp = blk.order_corr_spont(:);
            ord_st = blk.order_corr_stim(:);
        catch
            continue
        end

        % ── Per-cluster quantities ─────────────────────────────────────────
        var_spont = nan(nM,1);
        var_stim  = nan(nM,1);
        ord_sp_m  = nan(nM,1);
        ord_st_m  = nan(nM,1);
        width_m   = nan(nM,1);
        centroids = nan(nM, size(Chat,2));
        N_spont_m = zeros(nM,1);
        N_stim_m  = zeros(nM,1);

        for m = 1:nM
            cl     = motifs(m);
            sp_idx = train_ix(kmean_id(train_ix) == cl);
            st_idx = stim_ix( kmean_id(stim_ix)  == cl);
            N_spont_m(m) = length(sp_idx);
            N_stim_m(m)  = length(st_idx);

            if cl >= 1 && cl <= length(ord_sp)
                ord_sp_m(m) = ord_sp(cl);
                ord_st_m(m) = ord_st(cl);
            end

            if N_spont_m(m) >= 5
                Xsp = Chat(sp_idx,:);
                mu  = mean(Xsp,1); mu = mu/(norm(mu)+1e-12);
                centroids(m,:) = mu;
                S   = max(min(Xsp*Xsp',1),-1);
                D   = 1-S;
                ut  = D(triu(true(size(D)),1));
                var_spont(m) = mean(ut);
                width_m(m)   = mean(ut);
            end

            if N_stim_m(m) >= 3
                Xst = Chat(st_idx,:);
                S   = max(min(Xst*Xst',1),-1);
                D   = 1-S;
                ut  = D(triu(true(size(D)),1));
                var_stim(m) = mean(ut);
            end
        end

        delta_order = ord_st_m - ord_sp_m;
        delta_var   = var_stim - var_spont;

        % ── Basin sharpness ────────────────────────────────────────────────
        valid_c = find(~isnan(width_m) & ~any(isnan(centroids),2));
        if length(valid_c) < 4, continue; end

        separation = nan(nM,1);
        C_valid    = centroids(valid_c,:);
        for mi = 1:length(valid_c)
            m      = valid_c(mi);
            others = setdiff(1:length(valid_c),mi);
            dists  = 1 - C_valid(mi,:)*C_valid(others,:)';
            separation(m) = min(dists);
        end
        sharpness = separation ./ (width_m + 1e-12);

        % ── Median split ──────────────────────────────────────────────────
        fb.sharp_order_spont(j) = msplit_bc(sharpness, ord_sp_m,  N_spont_m, N_stim_m, 10, 0);
        fb.sharp_order_stim(j)  = msplit_bc(sharpness, ord_st_m,  N_spont_m, N_stim_m, 10, 2);
        fb.sharp_delta_order(j) = msplit_bc(sharpness, delta_order,N_spont_m, N_stim_m, 10, 2);
        fb.sharp_var_spont(j)   = msplit_bc(sharpness, var_spont,  N_spont_m, N_stim_m, 10, 0);
        fb.sharp_var_stim(j)    = msplit_bc(sharpness, var_stim,   N_spont_m, N_stim_m, 10, 3);
        fb.sharp_delta_var(j)   = msplit_bc(sharpness, delta_var,  N_spont_m, N_stim_m, 10, 3);

    end  % fold loop

    out.session(i) = fb;

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
    'sharp_order_spont', 'Sharp → spont order   (pos=sharp→faithful)';   ...
    'sharp_order_stim',  'Sharp → stim order    (pos=sharp→faithful)';   ...
    'sharp_delta_order', 'Sharp → Δorder        (pos=sharp→less degrad)';...
    'sharp_var_spont',   'Sharp → spont var     (neg=sharp→tight)';       ...
    'sharp_var_stim',    'Sharp → stim var      (neg=sharp→tight stim)';  ...
    'sharp_delta_var',   'Sharp → Δvar          (neg=sharp→less increase)'};

fprintf('\n%s\n', repmat('=',1,72));
fprintf('Basin content median-split — N sessions: %d\n', nSes);
fprintf('Effect = mean(outcome | high sharpness) - mean(outcome | low sharpness)\n');
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

% ── Local helper ──────────────────────────────────────────────────────────

function eff = msplit_bc(predictor, outcome, N_sp, N_st, min_sp, min_st)
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