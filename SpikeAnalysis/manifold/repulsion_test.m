%% ============================================================
% REPULSION ANALYSIS v2
%
% Core question: does the stim event cause the post-stim ripple
% to move away from the anchor along the anchor→stim axis?
%
% Design:
%   For each stim triplet (anchor=spont, middle=stim, post=spont):
%     1. Define axis = anchor → stim (unit vector)
%     2. Compute projection of (post - anchor) onto this axis
%     3. Test whether mean projection < 0 (repulsion) or > 0 (attraction)
%
%   Control: for matched spontaneous triplets, project
%   (post - anchor) onto the SAME axis from the nearest
%   stim triplet — asking whether spontaneous flow already
%   has a component along this axis
%
%   Null: permute which ripple is j+2 within IRI-matched bins,
%   keeping the anchor→stim axis fixed
% ============================================================

clear repulsion

nPerm     = 1000;
max_iri   = 30;      % seconds
iri_bins  = [0 logspace(log10(0.05), log10(max_iri), 8) Inf];

for i = 1:length(ses)

    fprintf('Session %d / %d\n', i, length(ses))

    % ---- representations ----
    R_pv     = ses(i).binnedPopRipplez_norm;      % nRip x nNeurons
    R_tensor = ses(i).tensorModel.C;
    Rtz      = zscore(R_tensor,[],1);
    R_tensor = Rtz ./ (vecnorm(Rtz,2,2) + 1e-12);

    nRip  = length(ses(i).rip_ts);
    in    = ismember(1:nRip, ses(i).PV.block(1).stim_ix)';
    clust = ses(i).PV.block(1).kmean_id(:);

    dt_all = abs(repmat(ses(i).rip_ts', nRip, 1) - ...
                 repmat(ses(i).rip_ts,  1,    nRip));

    for rep = 1:2

        switch rep
            case 1; R = R_pv;     rlab = 'pv';
            case 2; R = R_tensor; rlab = 'tensor';
        end

        % ----------------------------------------------------------
        % PASS 1: collect stim triplets
        % Store anchor index, stim index, post index, and the
        % anchor→stim axis for each triplet
        % ----------------------------------------------------------

        stim_anchor_idx = [];
        stim_mid_idx    = [];
        stim_post_idx   = [];
        stim_iri        = [];   % IRI from anchor to post (for matching)

        for j = 1:nRip-2
            if in(j),   continue; end   % anchor must be spont
            if ~in(j+1), continue; end  % middle must be stim
            if in(j+2),  continue; end  % post must be spont
            if isnan(clust(j)) || isnan(clust(j+2)), continue; end
            if dt_all(j,j+2) > max_iri, continue; end

            stim_anchor_idx = [stim_anchor_idx; j];
            stim_mid_idx    = [stim_mid_idx;    j+1];
            stim_post_idx   = [stim_post_idx;   j+2];
            stim_iri        = [stim_iri;         dt_all(j,j+2)];
        end

        nStim = length(stim_anchor_idx);

        % ----------------------------------------------------------
        % PASS 2: collect spontaneous triplets
        % ----------------------------------------------------------

        spont_anchor_idx = [];
        spont_mid_idx    = [];
        spont_post_idx   = [];
        spont_iri        = [];

        for j = 1:nRip-2
            if in(j),    continue; end  % anchor must be spont
            if in(j+1),  continue; end  % middle must be spont
            if in(j+2),  continue; end  % post must be spont
            if isnan(clust(j)) || isnan(clust(j+2)), continue; end
            if dt_all(j,j+2) > max_iri, continue; end

            spont_anchor_idx = [spont_anchor_idx; j];
            spont_mid_idx    = [spont_mid_idx;    j+1];
            spont_post_idx   = [spont_post_idx;   j+2];
            spont_iri        = [spont_iri;         dt_all(j,j+2)];
        end

        nSpont = length(spont_anchor_idx);

        if nStim < 5 || nSpont < 5
            repulsion.(rlab).proj_stim(i)       = nan;
            repulsion.(rlab).proj_spont(i)       = nan;
            repulsion.(rlab).perm_p(i)           = nan;
            repulsion.(rlab).basin_cross_stim(i) = nan;
            repulsion.(rlab).basin_cross_spont(i)= nan;
            continue
        end

        % ----------------------------------------------------------
        % COMPUTE PROJECTIONS FOR STIM TRIPLETS
        %
        % For each stim triplet:
        %   axis    = unit(R[stim] - R[anchor])
        %   proj    = dot(R[post] - R[anchor], axis)
        %
        % Positive = post moved toward stim (attraction)
        % Negative = post moved away from stim (repulsion)
        % ----------------------------------------------------------

        proj_stim    = zeros(nStim, 1);
        axis_mat     = zeros(nStim, size(R,2));  % store axes for control

        for k = 1:nStim
            v_anchor  = R(stim_anchor_idx(k), :);
            v_stim    = R(stim_mid_idx(k),    :);
            v_post    = R(stim_post_idx(k),   :);

            ax        = v_stim - v_anchor;
            ax        = ax / (norm(ax) + 1e-12);

            proj_stim(k)   = dot(v_post - v_anchor, ax);
            axis_mat(k,:)  = ax;
        end

        % ----------------------------------------------------------
        % COMPUTE PROJECTIONS FOR SPONT TRIPLETS
        % Project (post - anchor) onto the axis from the nearest
        % stim triplet (matched by IRI)
        % This gives the spontaneous baseline along the SAME axes,
        % asking: does spontaneous flow already move along this axis?
        % ----------------------------------------------------------

        proj_spont = zeros(nSpont, 1);
        iri_bin_st = discretize(stim_iri,  iri_bins);
        iri_bin_sp = discretize(spont_iri, iri_bins);

        for k = 1:nSpont
            % find nearest stim triplet in same IRI bin
            same_bin = find(iri_bin_st == iri_bin_sp(k));
            if isempty(same_bin)
                proj_spont(k) = nan;
                continue
            end
            % pick randomly among matches
            pick = same_bin(randi(length(same_bin)));

            % use that stim triplet's axis
            ax = axis_mat(pick,:);

            v_anchor = R(spont_anchor_idx(k), :);
            v_post   = R(spont_post_idx(k),   :);

            proj_spont(k) = dot(v_post - v_anchor, ax);
        end

        proj_spont = proj_spont(~isnan(proj_spont));

        % ----------------------------------------------------------
        % OBSERVED STATISTICS
        % ----------------------------------------------------------

        mu_stim  = nanmean(proj_stim);
        mu_spont = nanmean(proj_spont);

        repulsion.(rlab).proj_stim(i)  = mu_stim;
        repulsion.(rlab).proj_spont(i) = mu_spont;
        repulsion.(rlab).proj_delta(i) = mu_stim - mu_spont;

        % ----------------------------------------------------------
        % PERMUTATION TEST
        % Null: shuffle which ripple is the post-event within
        % IRI-matched bins, keeping anchor→stim axis fixed
        % ----------------------------------------------------------

        % pool of candidate post-event ripples (spontaneous only)
        spont_pool    = find(~in & ~isnan(clust));
        iri_to_pool   = dt_all(stim_anchor_idx, spont_pool);  % nStim x nPool

        perm_means = zeros(nPerm, 1);

        for pp = 1:nPerm
            proj_perm = zeros(nStim, 1);
            for k = 1:nStim
                % find pool ripples within same IRI bin
                iri_k    = iri_to_pool(k,:);
                bin_k    = discretize(iri_k, iri_bins);
                cands    = spont_pool(bin_k == iri_bin_st(k));
                if isempty(cands)
                    proj_perm(k) = nan;
                    continue
                end
                pick     = cands(randi(length(cands)));
                ax       = axis_mat(k,:);
                v_anchor = R(stim_anchor_idx(k),:);
                v_perm   = R(pick,:);
                proj_perm(k) = dot(v_perm - v_anchor, ax);
            end
            perm_means(pp) = nanmean(proj_perm);
        end

        % two-tailed permutation p-value
        repulsion.(rlab).perm_p(i) = mean(abs(perm_means) >= abs(mu_stim));
        repulsion.(rlab).perm_dist{i} = perm_means;

        % ----------------------------------------------------------
        % BASIN CROSSING
        % Does the post-stim ripple cross a cluster boundary
        % more often than post-spont?
        % ----------------------------------------------------------

        bc_stim  = clust(stim_post_idx)  ~= clust(stim_anchor_idx);
        bc_spont = clust(spont_post_idx) ~= clust(spont_anchor_idx);

        repulsion.(rlab).basin_cross_stim(i)  = nanmean(bc_stim);
        repulsion.(rlab).basin_cross_spont(i) = nanmean(bc_spont);
        repulsion.(rlab).basin_cross_delta(i) = nanmean(bc_stim) - nanmean(bc_spont);

        % ----------------------------------------------------------
        % PERPENDICULAR DISPLACEMENT
        % How much does the post-event ripple move off the
        % anchor→stim axis? If repulsion is directional, perp
        % should not differ between conditions
        % ----------------------------------------------------------

        perp_stim = zeros(nStim,1);
        for k = 1:nStim
            ax       = axis_mat(k,:);
            v_anchor = R(stim_anchor_idx(k),:);
            v_post   = R(stim_post_idx(k),:);
            disp     = v_post - v_anchor;
            perp_stim(k) = norm(disp - dot(disp,ax)*ax);
        end

        perp_spont = zeros(nSpont,1);
        for k = 1:nSpont
            iri_bin_k = iri_bin_sp(k);
            same_bin  = find(iri_bin_st == iri_bin_k);
            if isempty(same_bin)
                perp_spont(k) = nan; continue
            end
            pick     = same_bin(randi(length(same_bin)));
            ax       = axis_mat(pick,:);
            v_anchor = R(spont_anchor_idx(k),:);
            v_post   = R(spont_post_idx(k),:);
            disp     = v_post - v_anchor;
            perp_spont(k) = norm(disp - dot(disp,ax)*ax);
        end

        repulsion.(rlab).perp_stim(i)  = nanmean(perp_stim);
        repulsion.(rlab).perp_spont(i) = nanmean(perp_spont(~isnan(perp_spont)));
        repulsion.(rlab).perp_delta(i) = repulsion.(rlab).perp_stim(i) - ...
                                         repulsion.(rlab).perp_spont(i);

    end % rep
end % session

%% ============================================================
% STATISTICS
% ============================================================

reps   = {'pv','tensor'};
fprintf('\n============ REPULSION RESULTS ============\n')

for r = 1:2
    rl = reps{r};
    fprintf('\n--- %s ---\n', upper(rl))

    % 1. Is mean projection negative (repulsion) or positive (attraction)?
    proj_st = repulsion.(rl).proj_stim(:);
    [~,p1,~,st1] = ttest(proj_st);
    d1 = nanmean(proj_st) / (nanstd(proj_st) + 1e-12);
    fprintf('Projection onto anchor→stim axis (0=no effect):\n')
    fprintf('  mean=%.4f, d=%.3f, t=%.3f, p=%.4f\n', ...
            nanmean(proj_st), d1, st1.tstat, p1)
    if nanmean(proj_st) < 0 && p1 < 0.05
        fprintf('  ** REPULSION: post-stim ripple moves AWAY from stim\n')
    elseif nanmean(proj_st) > 0 && p1 < 0.05
        fprintf('  ** ATTRACTION: post-stim ripple moves TOWARD stim\n')
    else
        fprintf('  No significant directional effect\n')
    end

    % 2. Stim vs matched spontaneous baseline along same axis
    delta = repulsion.(rl).proj_delta(:);
    [~,p2,~,st2] = ttest(delta);
    d2 = nanmean(delta) / (nanstd(delta) + 1e-12);
    fprintf('Stim projection vs spont baseline on same axis:\n')
    fprintf('  mean delta=%.4f, d=%.3f, t=%.3f, p=%.4f\n', ...
            nanmean(delta), d2, st2.tstat, p2)

    % 3. Permutation p-values across sessions
    perm_p = repulsion.(rl).perm_p(:);
    fprintf('Permutation p-values (median across sessions): %.4f\n', ...
            nanmedian(perm_p))
    fprintf('Sessions significant at p<0.05: %d / %d\n', ...
            sum(perm_p < 0.05), sum(~isnan(perm_p)))

    % 4. Basin crossing
    bc_delta = repulsion.(rl).basin_cross_delta(:);
    [~,p3,~,st3] = ttest(bc_delta);
    d3 = nanmean(bc_delta) / (nanstd(bc_delta) + 1e-12);
    fprintf('Basin crossing delta (stim - spont):\n')
    fprintf('  mean=%.4f, d=%.3f, t=%.3f, p=%.4f\n', ...
            nanmean(bc_delta), d3, st3.tstat, p3)

    % 5. Perpendicular displacement
    perp_delta = repulsion.(rl).perp_delta(:);
    [~,p4,~,st4] = ttest(perp_delta);
    d4 = nanmean(perp_delta) / (nanstd(perp_delta) + 1e-12);
    fprintf('Perpendicular displacement delta (should be ~0 if directional):\n')
    fprintf('  mean=%.4f, d=%.3f, t=%.3f, p=%.4f\n', ...
            nanmean(perp_delta), d4, st4.tstat, p4)

end

%% ============================================================
% PLOTS
% ============================================================

figure('Position',[100 100 1400 600]);
cols_cond = [0.5 0.5 0.5; 0.2 0.5 0.85];   % gray=spont, blue=stim

for r = 1:2
    rl  = reps{r};
    row = r;

    % --- Panel 1: projection distributions across sessions ---
    subplot(2,4,(row-1)*4+1)
    dat = [repulsion.(rl).proj_spont(:) repulsion.(rl).proj_stim(:)];
    for k = 1:2
        bar(k, nanmean(dat(:,k)), 0.5, ...
            'FaceColor', cols_cond(k,:), 'EdgeColor','none'); hold on
        errorbar(k, nanmean(dat(:,k)), ...
                 nanstd(dat(:,k))/sqrt(sum(~isnan(dat(:,k)))), ...
                 'k','LineWidth',1.5)
        plot(k + 0.1*randn(size(dat,1),1), dat(:,k), 'o', ...
             'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.8 0.8 0.8], ...
             'MarkerSize',5)
    end
    yline(0,'--k','LineWidth',1)
    [~,p] = ttest(dat(:,2));   % test stim projection against 0
    ymax  = max(abs(nanmean(dat)))*1.4 + 0.01;
    text(2, ymax, sprintf('p=%.3f vs 0',p), ...
         'HorizontalAlignment','center','FontSize',9)
    set(gca,'XTick',1:2,'XTickLabel',{'Spont\n(matched axis)','Stim'}, ...
            'FontSize',10)
    ylabel('Mean projection onto anchor→stim axis')
    title(sprintf('%s: Directional projection\n(- = repulsion)', upper(rl)))
    box off

    % --- Panel 2: permutation distributions (first session shown) ---
    subplot(2,4,(row-1)*4+2)
    ses_show = find(~isnan(repulsion.(rl).perm_p), 1);
    if ~isempty(ses_show)
        pd = repulsion.(rl).perm_dist{ses_show};
        histogram(pd, 30, 'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none')
        hold on
        xline(repulsion.(rl).proj_stim(ses_show), '-b', 'LineWidth', 2)
        xline(0, '--k')
        xlabel('Mean projection (permuted)')
        ylabel('Count')
        title(sprintf('%s: Permutation null\n(session %d)', upper(rl), ses_show))
        legend({'Null','Observed'},'FontSize',8,'Location','northwest')
    end
    box off

    % --- Panel 3: basin crossing ---
    subplot(2,4,(row-1)*4+3)
    dat = [repulsion.(rl).basin_cross_spont(:) ...
           repulsion.(rl).basin_cross_stim(:)];
    for k = 1:2
        bar(k, nanmean(dat(:,k)), 0.5, ...
            'FaceColor', cols_cond(k,:), 'EdgeColor','none'); hold on
        errorbar(k, nanmean(dat(:,k)), ...
                 nanstd(dat(:,k))/sqrt(sum(~isnan(dat(:,k)))), ...
                 'k','LineWidth',1.5)
        plot(k + 0.1*randn(size(dat,1),1), dat(:,k), 'o', ...
             'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.8 0.8 0.8], ...
             'MarkerSize',5)
    end
    [~,p] = ttest(dat(:,1)-dat(:,2));
    ymax  = max(nanmean(dat))*1.2 + 0.01;
    text(1.5, ymax, sprintf('p=%.3f',p), ...
         'HorizontalAlignment','center','FontSize',9)
    set(gca,'XTick',1:2,'XTickLabel',{'Spont','Stim'},'FontSize',10)
    ylabel('Basin crossing probability')
    title(sprintf('%s: Basin crossing', upper(rl)))
    box off

    % --- Panel 4: perpendicular displacement ---
    subplot(2,4,(row-1)*4+4)
    dat = [repulsion.(rl).perp_spont(:) repulsion.(rl).perp_stim(:)];
    for k = 1:2
        bar(k, nanmean(dat(:,k)), 0.5, ...
            'FaceColor', cols_cond(k,:), 'EdgeColor','none'); hold on
        errorbar(k, nanmean(dat(:,k)), ...
                 nanstd(dat(:,k))/sqrt(sum(~isnan(dat(:,k)))), ...
                 'k','LineWidth',1.5)
        plot(k + 0.1*randn(size(dat,1),1), dat(:,k), 'o', ...
             'Color',[0.6 0.6 0.6],'MarkerFaceColor',[0.8 0.8 0.8], ...
             'MarkerSize',5)
    end
    [~,p] = ttest(dat(:,1)-dat(:,2));
    ymax  = max(nanmean(dat))*1.2 + 0.01;
    text(1.5, ymax, sprintf('p=%.3f',p), ...
         'HorizontalAlignment','center','FontSize',9)
    set(gca,'XTick',1:2,'XTickLabel',{'Spont','Stim'},'FontSize',10)
    ylabel('Perpendicular displacement')
    title(sprintf('%s: Off-axis displacement\n(should be ~0 if directional)', upper(rl)))
    box off

end

sgtitle('Repulsion hypothesis: absolute displacement along anchor→stim axis', ...
        'FontSize',13,'FontWeight','bold')