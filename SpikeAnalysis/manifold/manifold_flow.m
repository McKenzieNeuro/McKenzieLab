%% ============================================================
% SESSION MANIFOLD FLOW  v2
%
% For each session:
%   1. Select the 2 tensor dimensions that maximise
%      inter-cluster centre separation (Fisher criterion)
%      among the top n_top candidate motifs
%   2. Plot those motifs in uncentered native space
%   3. Density cloud + flow field both colored by cluster ID
%   4. Stim trials overlaid per-cluster in a lighter tint
% ============================================================

nFold         = 10;
fold_use      = 1;
min_sp        = 5;
min_st        = 2;
n_top         = 3;
n_dims_search = 15;   % search more components for better separation

col_st_alpha = 0.75;

% Qualitative cluster palette — well-separated hues
CPAL = [0.15 0.35 0.85;   % blue
        0.10 0.65 0.35;   % green
        0.65 0.10 0.75];  % purple

% Stim palette: warm/red tints matched to each cluster
% (same hue family but shifted toward red-orange)
CPAL_ST = [0.85 0.30 0.15;   % warm red  (for blue cluster)
           0.90 0.65 0.05;   % amber     (for green cluster)
           0.95 0.20 0.50];  % hot pink  (for purple cluster)

nSes = length(ses);
figure('Position',[20 20 340*nSes 380],'Color','w');

for i = 1:nSes

    %% ---- Find fold ----
    f = fold_use;
    while f <= nFold && ~isfield(ses(i).tensor.block(f),'tensor_flow')
        f = f + 1;
    end
    if f > nFold, continue; end

    tf       = ses(i).tensor.block(f).tensor_flow;
    stim_ix  = ses(i).tensor.block(f).stim_ix;
    clust_id = ses(i).tensor.block(f).kmean_id(:);
    stim_vec = ismember(1:length(clust_id), stim_ix)';

    X_raw = double(ses(i).binnedPopRipple_tensor);
    A_mat = ses(i).tensorModel.A;
    [~,T,R] = size(X_raw);
    nC_avail = min(n_dims_search, size(A_mat,2));

    % Project all ripples into full tensor space
    Zfull = zeros(size(A_mat,2), T, R);
    for r = 1:R
        Zfull(:,:,r) = A_mat' * squeeze(X_raw(:,:,r));
    end

    %% ---- Get UMAP embedding for this fold ----
    % UMAP coords stored in ses(i).PV.block(f).umap_proj — R x 2
    % (all ripples: spont training + held-out + stim projected in)
    umap_all = ses(i).tensor.block(f).umap_proj;   % R x 2

    %% ---- Joint cluster selection in UMAP space ----
    % Separability = min pairwise UMAP centroid distance / mean within-cluster spread
    % Score = 0.7 * norm_separability + 0.3 * norm_mean_rotFrac

    motifs_sp = tf.spont.motifID(:);
    rf_sp_all = tf.spont.rotFracTotal(:);
    motifs_st = tf.stim.motifID(:);

    valid = false(length(motifs_sp),1);
    for m = 1:length(motifs_sp)
        cl  = motifs_sp(m);
        nsp = sum(clust_id==cl & ~stim_vec);
        nst = sum(clust_id==cl &  stim_vec);
        valid(m) = nsp>=min_sp && nst>=min_st && ismember(cl,motifs_st);
    end
    if sum(valid) < 2, continue; end

    mot_valid = motifs_sp(valid);
    rf_valid  = rf_sp_all(valid);
    nV        = length(mot_valid);

    % Per-motif UMAP centroid and spread (spont trials only)
    umap_centers = zeros(nV, 2);
    umap_spread  = zeros(nV, 1);
    for m = 1:nV
        cl   = mot_valid(m);
        kp   = clust_id==cl & ~stim_vec;
        upts = umap_all(kp,:);
        umap_centers(m,:) = mean(upts, 1);
        umap_spread(m)    = mean(std(upts, 0, 1)) + 1e-12;
    end

    % Enumerate all subsets of size n_top
    combos   = nchoosek(1:nV, min(n_top,nV));
    nCombos  = size(combos,1);
    rf_norm  = (rf_valid-min(rf_valid)) / (max(rf_valid)-min(rf_valid)+1e-12);

    sep_scores = zeros(nCombos,1);
    rf_scores  = zeros(nCombos,1);
    for ci = 1:nCombos
        sub = combos(ci,:);
        % Min pairwise centroid distance normalised by mean within-spread
        min_dist = inf;
        for m1 = 1:length(sub)
            for m2 = m1+1:length(sub)
                d = norm(umap_centers(sub(m1),:) - umap_centers(sub(m2),:));
                min_dist = min(min_dist, d);
            end
        end
        sep_scores(ci) = min_dist / mean(umap_spread(sub));
        rf_scores(ci)  = mean(rf_norm(sub));
    end

    % Normalise and combine
    sep_n = sep_scores / (max(sep_scores)+1e-12);
    joint_scores = 0.7*sep_n + 0.3*rf_scores;
    [best_joint, best_ci] = max(joint_scores);
    best_subset = combos(best_ci,:);

    top_mot = mot_valid(best_subset);
    top_rf  = rf_valid(best_subset);
    nM      = length(top_mot);

    fprintf('Session %d: motifs=%s  rF=[%s]  sep=%.2f  joint=%.3f\n',...
        i, num2str(top_mot'), num2str(top_rf','%.2f '), ...
        sep_scores(best_ci), best_joint)

    %% ---- Use UMAP coords as display space ----
    % Flow field is computed in tensor space then mapped to UMAP via
    % local linear approximation (Jacobian of UMAP map)
    % Display coords: umap_all (R x 2), indexed by ripple



    %% ---- Plot ----
    ax = subplot(1, nSes, i);
    hold on

    for m = 1:nM
        cl   = top_mot(m);
        mcol = CPAL(m,:);

        kp_sp = clust_id==cl & ~stim_vec;
        kp_st = clust_id==cl &  stim_vec;

        % UMAP display coordinates (ripple-level, not time-resolved)
        usp = umap_all(kp_sp,:);   % nSp x 2
        ust = umap_all(kp_st,:);   % nSt x 2
        Rsp = size(usp,1);
        Rst = size(ust,1);

        %% -- Density cloud in UMAP space --
        draw_density_cloud2(usp, mcol, 0.9, 0.85)
        hold on

        %% -- Flow field in UMAP space --
        % Compute position-velocity in tensor space, then project
        % velocity into UMAP via local Jacobian:
        %   J(x) = d(umap)/d(tensor) estimated from nearest neighbours

        % Tensor coords at middle timepoint for spont ripples
        t_mid   = ceil(T/2);
        Ztensor = squeeze(Zfull(:, t_mid, kp_sp))';   % Rsp x nDims

        % Position-velocity in tensor space (use all timepoints)
        Xp = []; Vp = []; uXp = [];
        for rr = 1:Rsp
            rip_idx = find(kp_sp);
            tr = squeeze(Zfull(:,:,rip_idx(rr)));   % nDims x T
            u0 = usp(rr,:);                          % UMAP coord for this ripple
            for t = 2:T-1
                Xp  = [Xp,  tr(:,t)];
                Vp  = [Vp,  (tr(:,t+1)-tr(:,t-1))/2];
                uXp = [uXp; u0];   % UMAP position for this sample
            end
            Xp  = [Xp,  tr(:,T)];
            Vp  = [Vp,  tr(:,T)-tr(:,T-1)];
            uXp = [uXp; u0];
        end
        if size(Xp,2) < 4, continue; end

        % Grid in UMAP space
        n_sp_pts = size(usp,1);
        bw_kde   = std(usp) * n_sp_pts^(-1/5) * 1.5;
        bw_kde(bw_kde<1e-6) = 0.05;
        x1r = [min(usp(:,1))-3*bw_kde(1), max(usp(:,1))+3*bw_kde(1)];
        x2r = [min(usp(:,2))-3*bw_kde(2), max(usp(:,2))+3*bw_kde(2)];
        gx  = linspace(x1r(1), x1r(2), 14);
        gy  = linspace(x2r(1), x2r(2), 14);
        [GX,GY] = meshgrid(gx, gy);

        % KDE density mask
        dens_grid = zeros(numel(GX),1);
        for jj = 1:n_sp_pts
            d2v = ((GX(:)-usp(jj,1))/bw_kde(1)).^2 + ...
                  ((GY(:)-usp(jj,2))/bw_kde(2)).^2;
            dens_grid = dens_grid + exp(-0.5*d2v);
        end
        dens_mask = reshape(dens_grid/max(dens_grid), size(GX)) > 0.003;

        % At each grid point: estimate local Jacobian J (2 x nDims)
        % from k nearest UMAP-tensor pairs, then project tensor velocity
        % into UMAP: v_umap = J * v_tensor
        lbw_u = 0.6 * mean([range(x1r) range(x2r)]);
        k_nn  = min(15, Rsp);   % neighbours for Jacobian
        U = zeros(size(GX));
        V = zeros(size(GY));

        for gi = 1:numel(GX)
            if ~dens_mask(gi), continue; end
            gpt_u = [GX(gi); GY(gi)];

            % UMAP-space weights for local affine tensor->velocity fit
            dist2_u = sum((uXp - gpt_u').^2, 2);
            w_u     = exp(-dist2_u / (2*lbw_u^2));
            if sum(w_u) < 1e-6, continue; end

            % Local Jacobian: umap_pos ~ J * tensor_pos + c
            % Fit J from (tensor_pos, umap_pos) pairs weighted by w_u
            Wu     = diag(w_u);
            Xtaug  = [Xp; ones(1,size(Xp,2))];
            lamJ   = 1e-3 * trace(Xtaug*Wu*Xtaug') / (size(Xp,1)+1);
            J_aug  = (uXp'*Wu*Xtaug') / (Xtaug*Wu*Xtaug' + lamJ*eye(size(Xtaug,1)));
            J      = J_aug(:, 1:size(Xp,1));   % 2 x nDims

            % Local affine flow in tensor space at nearest tensor point
            % Find tensor position closest to this UMAP grid point
            [~, nn_idx] = min(dist2_u);
            x_tensor = Xp(:, nn_idx);

            Xtaug2 = [Xp; ones(1,size(Xp,2))];
            dist2_t = sum((Xp - x_tensor).^2, 1)';
            w_t     = exp(-dist2_t / (2*(lbw_u*2)^2));
            Wt      = diag(w_t);
            lamB    = 1e-3 * trace(Xtaug2*Wt*Xtaug2') / (size(Xp,1)+1);
            Baug    = (Vp*Wt*Xtaug2') / (Xtaug2*Wt*Xtaug2' + lamB*eye(size(Xtaug2,1)));
            v_tensor = Baug * [x_tensor; 1];

            % Project into UMAP space
            v_umap  = J * v_tensor;
            U(gi)   = v_umap(1);
            V(gi)   = v_umap(2);
        end

        % Scale arrows
        mag  = sqrt(U.^2+V.^2);
        p95  = prctile(mag(dens_mask), 95);
        ascl = 0.13 * max(range(x1r),range(x2r)) / (p95+1e-12);
        Un   = U*ascl;  Vn = V*ascl;
        mag_n = mag / (p95+1e-12);

        for gi = 1:numel(GX)
            if ~dens_mask(gi), continue; end
            a  = 0.30 + 0.70*min(1, mag_n(gi));
            ac = mcol*a + [1 1 1]*(1-a);
            quiver(GX(gi),GY(gi),Un(gi),Vn(gi),0,...
                'Color',ac,'LineWidth',1.5,'MaxHeadSize',2.5,'AutoScale','off')
        end

        %% -- Stim ripples: scatter in UMAP space (no time axis) --
        stim_col = CPAL_ST(m,:);
        scatter(ust(:,1), ust(:,2), 40, stim_col, 'filled',...
            'MarkerFaceAlpha', col_st_alpha, 'MarkerEdgeColor','w',...
            'LineWidth', 0.5)

        %% -- Cluster label at UMAP centroid --
        cx = mean(usp(:,1));
        cy = mean(usp(:,2));
        text(cx, cy, sprintf('M%d  rF=%.2f', cl, top_rf(m)), ...
            'FontSize',7,'Color',mcol*0.6,...
            'HorizontalAlignment','center','FontWeight','bold',...
            'BackgroundColor','w','Margin',1)
    end

    %% ---- Cosmetics ----
    grid on
    set(ax,'FontSize',8,'GridAlpha',0.10,'Color','w',...
        'XTickLabel','','YTickLabel','')
    xlabel('UMAP 1','FontSize',9)
    ylabel('UMAP 2','FontSize',9)
    title(sprintf('Session %d  (joint=%.2f)', i, best_joint),...
        'FontSize',8,'FontWeight','bold')
    box off
end

sgtitle(['Per-session manifold flow  |  UMAP embedding space  |  '...
    'Cluster selection: 0.7×separability + 0.3×rotFrac  |  '...
    'Color = cluster ID  |  Dots = stim ripples  |  '...
    'Arrow brightness = flow speed'],...
    'FontSize',8)