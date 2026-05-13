function [out] = correlate_flow_occupancy(ses)
% FLOW_OCCUPANCY_CLUSTER  Cluster-level association between stim occupancy
%   ratio and per-cluster metric deltas (stim - spont).
%
% For each session x fold x embedding, clusters are split at the median of
% each geometry metric delta. The effect size is:
%   effect = mean(occ_ratio | top half) - mean(occ_ratio | bottom half)
%
% Positive = clusters with larger metric delta are more visited.
% Negative = clusters with larger metric delta are less visited (expected).
%
% Per-fold effects are aggregated to session level by median.
% Session-level effects are tested by sign test and t-test across sessions.
%
% Metrics:
%   (1) Rate              rate_stim - rate_spont
%   (2) Order             order_corr_stim - spont
%   (3) NN                nearest_neighbor_global stim - spont  [per-ripple]
%   (4) GeodesicRatio     geodesic_ratio stim - spont
%   (5) AnisotropicExp    jacobian_condition stim - spont
%   (6) Volume            jacobian_volume stim - spont
%   (7) OffManifold       manifold stim - spont                 [per-ripple]
%   (8) FlowSpin          tensor_flow totalSpread stim - spont  [tensor only]

metrics = struct();
metrics(1).name   = 'Rate';
metrics(1).stim   = 'rate_stim';
metrics(1).spont  = 'rate_spont';
metrics(1).perRip = false;

metrics(2).name   = 'Order';
metrics(2).stim   = 'order_corr_stim';
metrics(2).spont  = 'order_corr_spont';
metrics(2).perRip = false;

metrics(3).name   = 'NN';
metrics(3).stim   = 'nearest_neighbor_global_stim';
metrics(3).spont  = 'nearest_neighbor_global_spont';
metrics(3).perRip = true;

metrics(4).name   = 'GeodesicRatio';
metrics(4).stim   = 'geodesic_ratio_stim';
metrics(4).spont  = 'geodesic_ratio_spont';
metrics(4).perRip = false;

metrics(5).name   = 'AnisotropicExp';
metrics(5).stim   = 'jacobian_condition_stim';
metrics(5).spont  = 'jacobian_condition_spont';
metrics(5).perRip = false;

metrics(6).name   = 'Volume';
metrics(6).stim   = 'jacobian_volume_stim';
metrics(6).spont  = 'jacobian_volume_spont';
metrics(6).perRip = false;

metrics(7).name   = 'OffManifold';
metrics(7).stim   = 'manifold_stim';
metrics(7).spont  = 'manifold_spont';
metrics(7).perRip = true;

metrics(8).name   = 'FlowSpin';
metrics(8).stim   = '';
metrics(8).spont  = '';
metrics(8).perRip = false;

nMetrics = length(metrics);
embs     = {'PV','tensor'};
nEmb     = 2;
nSes     = length(ses);
nFold    = length(ses(1).tensor.block);

% Storage: effect_folds{c,e,i} = nFold x 3 [effect, mean_occ_high, mean_occ_low]
effect_folds = cell(nMetrics, nEmb, nSes);
for c = 1:nMetrics; for e = 1:nEmb; for i = 1:nSes
    effect_folds{c,e,i} = nan(nFold,3);
end; end; end

effect_ses = nan(nMetrics, nEmb, nSes, 3);  % [:,:,:,1]=effect, 2=hi, 3=lo

% ── Main loop ─────────────────────────────────────────────────────────────
for i = 1:nSes
    for j = 1:nFold
        for e = 1:nEmb
            embField = embs{e};
            blk_emb  = ses(i).(embField).block(j);
            kmean_id = blk_emb.kmean_id(:);
            stim_ix  = blk_emb.stim_ix(:);
            train_ix = blk_emb.training_ix(:);
            test_ix  = blk_emb.test_ix(:);
            K        = blk_emb.optimal_k;
            R_tot    = length(kmean_id);

            stim_vec  = false(R_tot,1); stim_vec(stim_ix)  = true;
            spont_vec = false(R_tot,1); spont_vec(train_ix) = true;
            test_vec  = false(R_tot,1); test_vec(test_ix)   = true;

            N_stim  = accumarray(kmean_id(stim_vec),  1, [K 1], @sum, 0);
            N_spont = accumarray(kmean_id(spont_vec), 1, [K 1], @sum, 0);

            p_stim    = N_stim  / (sum(N_stim)  + 1e-12);
            p_spont   = N_spont / (sum(N_spont) + 1e-12);
            occ_ratio = p_stim ./ (p_spont + 1e-12);

            for c = 1:nMetrics
                if c == 8 && e == 1, continue; end

                try
                    if c == 8
                        tf      = ses(i).tensor.block(j).tensor_flow;
                        motifs  = tf.spont.motifID(:);
                        nM      = length(motifs);
                        dMetric = tf.stim.totalSpread(:) - tf.spont.totalSpread(:);
                        clust_delta = nan(K,1);
                        for m = 1:nM
                            cl = motifs(m);
                            if cl >= 1 && cl <= K
                                clust_delta(cl) = dMetric(m);
                            end
                        end

                    elseif metrics(c).perRip
                        sv = ses(i).(embField).block(j).(metrics(c).stim)(:);
                        pv = ses(i).(embField).block(j).(metrics(c).spont)(:);
                        mu_stim  = accumarray(kmean_id(stim_vec), sv, [K 1], @nanmean, nan);
                        mu_spont = accumarray(kmean_id(test_vec),  pv, [K 1], @nanmean, nan);
                        clust_delta = mu_stim - mu_spont;

                    else
                        sv = ses(i).(embField).block(j).(metrics(c).stim)(:);
                        pv = ses(i).(embField).block(j).(metrics(c).spont)(:);
                        clust_delta = sv - pv;
                    end

                    % ── Median split on clust_delta ───────────────────────
                    valid = N_stim >= 2 & N_spont >= 5 & ...
                            ~isnan(clust_delta) & ~isnan(occ_ratio);
                    if sum(valid) < 4, continue; end

                    delta_v = clust_delta(valid);
                    occ_v   = occ_ratio(valid);
                    med     = median(delta_v);
                    top     = delta_v >= med;   % high metric delta
                    bot     = ~top;

                    if sum(top) < 1 || sum(bot) < 1, continue; end

                    % Negative effect = disrupted clusters are avoided
                    effect_folds{c,e,i}(j,1) = mean(occ_v(top)) - mean(occ_v(bot));
                    effect_folds{c,e,i}(j,2) = mean(occ_v(top));   % high delta half
                    effect_folds{c,e,i}(j,3) = mean(occ_v(bot));   % low delta half

                catch
                end
            end
        end
    end

    % Aggregate folds → session median (all three columns)
    for c = 1:nMetrics
        for e = 1:nEmb
            if c == 8 && e == 1, continue; end
            fv = effect_folds{c,e,i};   % nFold x 3
            for col = 1:3
                vals = fv(:,col);
                vals = vals(~isnan(vals));
                if length(vals) >= 3
                    effect_ses(c,e,i,col) = median(vals);
                end
            end
        end
    end
end

% ── Summary ───────────────────────────────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,72));
fprintf('Occupancy median-split by metric delta — N sessions: %d\n', nSes);
fprintf('Effect = mean(occ | high delta) - mean(occ | low delta)\n');
fprintf('Negative = disrupted clusters avoided (predicted direction)\n');
fprintf('%s\n', repmat('-',1,72));
fprintf('%-20s  %-8s  %+8s  %8s  %7s  %7s\n', ...
    'Metric','Emb','median','IQR','p(sign)','p(ttest)');
fprintf('%s\n', repmat('-',1,72));

out.results = struct();

for c = 1:nMetrics
    for e = 1:nEmb
        if c == 8 && e == 1, continue; end

        ev    = squeeze(effect_ses(c,e,:,1));
        ev_hi = squeeze(effect_ses(c,e,:,2));
        ev_lo = squeeze(effect_ses(c,e,:,3));
        ev = ev(~isnan(ev));
        if length(ev) < 3, continue; end

        med_e  = median(ev);
        iqr_e  = iqr(ev);
        hi_m   = median(ev_hi(~isnan(ev_hi)));
        lo_m   = median(ev_lo(~isnan(ev_lo)));

        % Two-sided sign test
        n_neg  = sum(ev < 0);
        n_pos  = sum(ev > 0);
        n_tot  = n_neg + n_pos;
        p_sign = 2 * binocdf(min(n_neg,n_pos), n_tot, 0.5);

        [~, p_tt] = ttest(ev);

        fprintf('%-20s  %-8s  %+8.4f  %8.4f  %7.4f  %7.4f  [hi=%.3f lo=%.3f]\n', ...
            metrics(c).name, embs{e}, med_e, iqr_e, p_sign, p_tt, hi_m, lo_m);

        tag = sprintf('%s_%s', metrics(c).name, embs{e});
        out.results.(tag).effect_ses = ev;
        out.results.(tag).median     = med_e;
        out.results.(tag).iqr        = iqr_e;
        out.results.(tag).p_sign     = p_sign;
        out.results.(tag).p_ttest    = p_tt;
        out.results.(tag).n_neg      = n_neg;
        out.results.(tag).n_pos      = n_pos;
        out.results.(tag).occ_hi     = hi_m;   % mean occ ratio, high-delta clusters
        out.results.(tag).occ_lo     = lo_m;   % mean occ ratio, low-delta clusters
    end
end
fprintf('%s\n', repmat('=',1,72));

out.metrics      = metrics;
out.embs         = embs;
out.effect_ses   = effect_ses;
out.effect_folds = effect_folds;

end