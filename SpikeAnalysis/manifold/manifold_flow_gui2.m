function manifold_flow_gui(ses)
%% MANIFOLD_FLOW_GUI  Interactive cluster selector for per-session manifold flow
%
%  manifold_flow_gui(ses)
%
%  Left panel  : session list — click to load session
%  Centre panel: cluster checkboxes (all valid clusters, with rotFrac + nSpont/nStim)
%                spinner per cluster to set number of stim examples shown
%  Right panel : live UMAP plot, redraws on any change
%
%  Workflow:
%    1. Select session from list
%    2. Check/uncheck clusters to include
%    3. Adjust stim-example count per cluster (1–5)
%    4. Hit "Export Figure" to write a clean figure to file

% ── Parameters ────────────────────────────────────────────────────────────
nFold    = 10;
fold_use = 1;
min_sp   = 5;
min_st   = 2;
CPAL     = [0.15 0.35 0.85;   % blue
            0.10 0.65 0.35;   % green
            0.65 0.10 0.75;   % purple
            0.85 0.55 0.05;   % amber
            0.05 0.65 0.65;   % teal
            0.80 0.10 0.30;   % crimson
            0.40 0.20 0.80];  % indigo

nSes = length(ses);

% ── Pre-compute per-session UMAP + cluster info (cached in S struct) ──────
fprintf('Pre-computing UMAP for %d sessions...\n', nSes);
S = precompute_sessions(ses, nFold, fold_use, min_sp, min_st);
fprintf('Done.\n');

% ── Build GUI ─────────────────────────────────────────────────────────────
fig = uifigure('Name','Manifold Flow — Cluster Selector',...
    'Position',[40 40 1300 700],'Color',[0.12 0.12 0.14]);

% ── Left: session list ────────────────────────────────────────────────────
uilabel(fig,'Text','Sessions','Position',[10 665 120 22],...
    'FontColor','w','FontWeight','bold','FontSize',11,...
    'BackgroundColor','none');
ses_list = uilistbox(fig,...
    'Position',[10 80 130 580],...
    'BackgroundColor',[0.18 0.18 0.20],...
    'FontColor','w',...
    'FontSize',10);
ses_items = cell(nSes,1);
for k = 1:nSes
    if isempty(S{k})
        ses_items{k} = sprintf('S%d  [no data]', k);
    else
        ses_items{k} = sprintf('S%d  (%d clusters)', k, length(S{k}.mot_valid));
    end
end
ses_list.Items = ses_items;
ses_list.Value = ses_items{1};

% ── Centre: cluster controls panel ───────────────────────────────────────
ctrl_panel = uipanel(fig,'Position',[150 10 320 680],...
    'BackgroundColor',[0.15 0.15 0.17],...
    'BorderType','none');
uilabel(ctrl_panel,'Text','Clusters','Position',[10 650 280 22],...
    'FontColor','w','FontWeight','bold','FontSize',11,...
    'BackgroundColor','none');
uilabel(ctrl_panel,'Text','☑ Include   Motif   rF   nSp  nSt   #',...
    'Position',[10 628 300 18],'FontColor',[0.6 0.6 0.6],'FontSize',9,...
    'BackgroundColor','none');

% Scrollable container for cluster rows
MAX_CLUST  = 20;
row_h      = 34;
scroll_h   = 570;   % visible height of scroll area
inner_h    = MAX_CLUST * row_h + 10;   % total content height

scroll_panel = uipanel(ctrl_panel, ...
    'Position',     [0 55 320 scroll_h], ...
    'BackgroundColor', [0.15 0.15 0.17], ...
    'BorderType',   'none', ...
    'Scrollable',   'on');

cluster_rows = struct();
for r = 1:MAX_CLUST
    % Stack rows top-to-bottom inside the scrollable panel
    y = inner_h - r*row_h;
    cluster_rows(r).panel = uipanel(scroll_panel,...
        'Position',[4 y 304 row_h-2],...
        'BackgroundColor',[0.20 0.20 0.22],...
        'BorderType','line','BorderWidth',1,...
        'Visible','off');
    cluster_rows(r).check = uicheckbox(cluster_rows(r).panel,...
        'Position',[8 7 20 20],...
        'Text','','Value',false,...
        'FontColor','w');
    cluster_rows(r).swatch = uipanel(cluster_rows(r).panel,...
        'Position',[32 8 16 18],...
        'BackgroundColor',[0.5 0.5 0.5],'BorderType','none');
    cluster_rows(r).label = uilabel(cluster_rows(r).panel,...
        'Position',[54 5 178 24],...
        'Text','','FontColor','w','FontSize',9,...
        'BackgroundColor','none');
    cluster_rows(r).spinner = uispinner(cluster_rows(r).panel,...
        'Position',[238 6 58 22],...
        'Limits',[1 20],'Value',3,...
        'Step',1,'FontSize',9,...
        'BackgroundColor',[0.25 0.25 0.27],'FontColor','w');
end

% Export / Redraw buttons sit below the scroll area
export_btn = uibutton(ctrl_panel,'push',...
    'Text','Export Figure',...
    'Position',[10 10 140 34],...
    'BackgroundColor',[0.20 0.45 0.80],...
    'FontColor','w','FontWeight','bold','FontSize',10);
redraw_btn = uibutton(ctrl_panel,'push',...
    'Text','↺  Redraw',...
    'Position',[160 10 140 34],...
    'BackgroundColor',[0.25 0.25 0.27],...
    'FontColor','w','FontWeight','bold','FontSize',10);

% ── Right: axes ───────────────────────────────────────────────────────────
ax = uiaxes(fig,'Position',[480 20 800 650],...
    'Color',[0.97 0.97 0.98],...
    'XColor',[0.4 0.4 0.4],'YColor',[0.4 0.4 0.4],...
    'FontSize',9,'GridAlpha',0.12,'Box','off');
xlabel(ax,'UMAP 1'); ylabel(ax,'UMAP 2');
title(ax,'Select a session','FontSize',10,'FontWeight','bold','Color',[0.2 0.2 0.2]);

% ── State ─────────────────────────────────────────────────────────────────
state.ses_idx      = 1;
state.checked      = false(MAX_CLUST,1);
state.n_stim       = 3*ones(MAX_CLUST,1);
state.n_clust      = 0;

% ── Callbacks ─────────────────────────────────────────────────────────────
ses_list.ValueChangedFcn     = @(src,evt) on_session_change(src);
export_btn.ButtonPushedFcn   = @(~,~)    on_export();
redraw_btn.ButtonPushedFcn   = @(~,~)    do_redraw();

% Initial load
on_session_change(ses_list);

% ══════════════════════════════════════════════════════════════════════════
    function on_session_change(src)
        idx = find(strcmp(src.Items, src.Value));
        if isempty(idx), return; end
        state.ses_idx = idx;
        D = S{idx};
        if isempty(D)
            title(ax, sprintf('Session %d — no valid data',idx),...
                'FontSize',10,'FontWeight','bold');
            hide_all_rows();
            return
        end

        nC = length(D.mot_valid);
        state.n_clust = nC;

        % Resize scrollable content to fit actual number of rows
        actual_h = max(scroll_h, nC * row_h + 10);
        scroll_panel.Position(4) = scroll_h;   % visible area stays fixed
        % reposition rows to fill from top
        for r = 1:MAX_CLUST
            y = actual_h - r*row_h;
            cluster_rows(r).panel.Position(2) = y;
        end

        % Default: top 3 by joint score checked
        top3 = D.top_mot;

        for r = 1:MAX_CLUST
            if r <= nC
                cl  = D.mot_valid(r);
                rf  = D.rf_valid(r);
                nsp = sum(D.clust_id==cl & ~D.stim_vec);
                nst = sum(D.clust_id==cl &  D.stim_vec);
                col = CPAL(mod(r-1,size(CPAL,1))+1,:);

                cluster_rows(r).panel.Visible        = 'on';
                cluster_rows(r).swatch.BackgroundColor = col;
                cluster_rows(r).label.Text           = ...
                    sprintf('M%-3d  rF=%.2f  %3dsp %2dst', cl, rf, nsp, nst);
                cluster_rows(r).check.Value           = ismember(cl, top3);
                cluster_rows(r).spinner.Value         = 3;

                % Callbacks — need to capture r in closure
                rr = r;
                cluster_rows(r).check.ValueChangedFcn   = @(~,~) on_ctrl_change(rr);
                cluster_rows(r).spinner.ValueChangedFcn = @(~,~) on_ctrl_change(rr);
            else
                cluster_rows(r).panel.Visible = 'off';
            end
        end

        do_redraw();
    end

    function on_ctrl_change(~)
        do_redraw();
    end

    function hide_all_rows()
        for r = 1:MAX_CLUST
            cluster_rows(r).panel.Visible = 'off';
        end
    end

    function do_redraw()
        rng('shuffle');   % advance RNG so randperm gives new selection every call
        idx = state.ses_idx;
        D   = S{idx};
        if isempty(D), return; end

        % Collect checked clusters + their stim counts
        nC_cur = length(D.mot_valid);
        sel_cl    = [];
        sel_nstim = [];
        for r = 1:min(nC_cur, MAX_CLUST)
            if cluster_rows(r).check.Value
                sel_cl(end+1)    = D.mot_valid(r);  %#ok
                sel_nstim(end+1) = cluster_rows(r).spinner.Value; %#ok
            end
        end

        cla(ax);
        if isempty(sel_cl)
            title(ax, sprintf('Session %d — no clusters selected',idx),...
                'FontSize',10,'FontWeight','bold','Color',[0.2 0.2 0.2]);
            return
        end

        hold(ax,'on');
        ax.XTickLabel = ''; ax.YTickLabel = '';
        grid(ax,'on');

        nM = length(sel_cl);
        for m = 1:nM
            cl    = sel_cl(m);
            n_st  = sel_nstim(m);
            mcol  = CPAL(mod(m-1,size(CPAL,1))+1,:);

            kp_sp = D.clust_id==cl & ~D.stim_vec;
            kp_st = D.clust_id==cl &  D.stim_vec;
            usp   = D.umap_all(kp_sp,:);
            ust   = D.umap_all(kp_st,:);
            Rsp   = size(usp,1);
            Rst   = size(ust,1);

            if Rsp < 3, continue; end

            % ── Density cloud (imagesc, drawn first) ─────────────────────
            if Rsp >= 3
                bw_cl = std(usp) * Rsp^(-1/5) * 1.5;
                bw_cl(bw_cl<1e-6) = 0.05;
                gx_cl = linspace(min(usp(:,1))-2*bw_cl(1), max(usp(:,1))+2*bw_cl(1), 80);
                gy_cl = linspace(min(usp(:,2))-2*bw_cl(2), max(usp(:,2))+2*bw_cl(2), 80);
                [GXc,GYc] = meshgrid(gx_cl, gy_cl);
                kde = zeros(size(GXc));
                for jj = 1:Rsp
                    d2  = ((GXc-usp(jj,1))/bw_cl(1)).^2 + ...
                          ((GYc-usp(jj,2))/bw_cl(2)).^2;
                    kde = kde + exp(-0.5*d2);
                end
                kde = kde / max(kde(:));
                img = ones(size(kde,1), size(kde,2), 3);
                for ch = 1:3
                    img(:,:,ch) = 1 - kde * (1 - mcol(ch)) * 0.55;
                end
                img_mask = repmat(kde < 0.05, [1 1 3]);
                img(img_mask) = 1;
                ih = image(ax, gx_cl, gy_cl, img);
                ih.AlphaData = min(kde * 0.7, 0.7);
                uistack(ih, 'bottom');
            end

            % ── Flow field: K-space velocity, K-space weighted Jacobian ──
            % uXp = trial-mean UMAP positions (one per ripple, repeated for T)
            % Jacobian weighted by K-space proximity — avoids embedding gradient artifact
            rip_sp_idx = find(kp_sp);
            Xp = []; Vp = []; uXp = [];
            for rr = 1:Rsp
                tr = squeeze(D.Zfull(:,:,rip_sp_idx(rr)));  % K x T
                u0 = usp(rr,:);
                for t = 2:D.T-1
                    Xp  = [Xp,  tr(:,t)];
                    Vp  = [Vp,  (tr(:,t+1)-tr(:,t-1))/2];
                    uXp = [uXp; u0];
                end
                Xp  = [Xp,  tr(:,D.T)];
                Vp  = [Vp,  tr(:,D.T)-tr(:,D.T-1)];
                uXp = [uXp; u0];
            end
            if size(Xp,2) < 4, continue; end

            % Grid defined over Ztensor UMAP positions (clean cluster display)
            bw_kde = std(usp) * Rsp^(-1/5) * 1.5;
            bw_kde(bw_kde<1e-6) = 0.05;
            x1r = [min(usp(:,1))-3*bw_kde(1), max(usp(:,1))+3*bw_kde(1)];
            x2r = [min(usp(:,2))-3*bw_kde(2), max(usp(:,2))+3*bw_kde(2)];
            gx  = linspace(x1r(1),x1r(2),20);
            gy  = linspace(x2r(1),x2r(2),20);
            [GX,GY] = meshgrid(gx,gy);

            dens_grid = zeros(numel(GX),1);
            for jj = 1:Rsp
                d2v = ((GX(:)-usp(jj,1))/bw_kde(1)).^2 + ...
                      ((GY(:)-usp(jj,2))/bw_kde(2)).^2;
                dens_grid = dens_grid + exp(-0.5*d2v);
            end
            dens_norm = reshape(dens_grid/max(dens_grid), size(GX));
            dens_mask = dens_norm > 0.08;

            % Jacobian bandwidth in K-space
            lbw_k = median(std(Xp',0,1)) * 2.0;
            U = zeros(size(GX)); V_field = zeros(size(GY));

            for gi = 1:numel(GX)
                if ~dens_mask(gi), continue; end
                gpt_u   = [GX(gi); GY(gi)];
                % Find nearest K-space point to this grid location via UMAP proxy
                dist2_u = sum((uXp - gpt_u').^2, 2);
                [~, nn] = min(dist2_u);
                x_ref   = Xp(:, nn);   % nearest K-space state
                % Weight by K-space proximity — avoids embedding gradient artifact
                dist2_k = sum((Xp - x_ref).^2, 1)';
                w_k     = exp(-dist2_k / (2*lbw_k^2));
                if sum(w_k) < 1e-6, continue; end
                Wk     = diag(w_k);
                Xtaug  = [Xp; ones(1,size(Xp,2))];
                % Jacobian J: maps K-space -> UMAP via locally weighted regression
                lamJ   = 1e-3*trace(Xtaug*Wk*Xtaug')/(size(Xp,1)+1);
                J_aug  = (uXp'*Wk*Xtaug')/(Xtaug*Wk*Xtaug'+lamJ*eye(size(Xtaug,1)));
                J      = J_aug(:,1:size(Xp,1));
                % K-space velocity at nearest neighbor
                dist2_k2 = sum((Xp - x_ref).^2, 1)';
                w_v      = exp(-dist2_k2 / (2*(lbw_k*2)^2));
                Wv       = diag(w_v);
                Xtaug2   = [Xp; ones(1,size(Xp,2))];
                lamB     = 1e-3*trace(Xtaug2*Wv*Xtaug2')/(size(Xp,1)+1);
                Baug     = (Vp*Wv*Xtaug2')/(Xtaug2*Wv*Xtaug2'+lamB*eye(size(Xtaug2,1)));
                v_k      = Baug * [x_ref; 1];
                % Project K-space velocity into UMAP via J
                v_u      = J * v_k;
                U(gi)        = v_u(1);
                V_field(gi)  = v_u(2);
            end

            mag   = sqrt(U.^2+V_field.^2);
            p95   = prctile(mag(dens_mask),95)+1e-12;
            gs    = mean([range(x1r)/(size(GX,2)-1), range(x2r)/(size(GX,1)-1)]);
            ascl  = 4.0*gs/p95;
            Un    = U*ascl; Vn = V_field*ascl;
            mag_n = mag/p95;

            for gi = 1:numel(GX)
                if ~dens_mask(gi), continue; end
                d_n = dens_norm(gi);
                s_n = min(1,mag_n(gi));
                a   = min(0.95, d_n^0.5*(0.35+0.45*s_n));
                ac  = mcol*a + [1 1 1]*(1-a);
                quiver(ax,GX(gi),GY(gi),Un(gi),Vn(gi),0,...
                    'Color',ac,'LineWidth',2.0,'MaxHeadSize',3,'AutoScale','off')
            end

        end  % for m = 1:nM (cloud + flow field)

        % ── Stim trajectories: J*(z_t - z_mean) ─────────────────────────
        all_traj_data = {};
        for m2 = 1:nM
            cl2    = sel_cl(m2);
            n_st2  = sel_nstim(m2);   % correct spinner value for this cluster
            kp_sp2 = D.clust_id==cl2 & ~D.stim_vec;
            kp_st2 = D.clust_id==cl2 &  D.stim_vec;
            usp2   = D.umap_all(kp_sp2,:);
            ust2   = D.umap_all(kp_st2,:);
            Rst2   = size(ust2,1);
            Rsp2   = size(usp2,1);
            if Rst2==0 || Rsp2<4, continue; end

            rip_sp_idx2 = find(kp_sp2);
            rip_st_idx2 = find(kp_st2);

            % Build Xp2/uXp2 for this cluster
            Xp2 = []; uXp2 = [];
            for rr2 = 1:Rsp2
                tr2 = squeeze(D.Zfull(:,:,rip_sp_idx2(rr2)));
                u0  = usp2(rr2,:);
                for t2 = 2:D.T-1
                    Xp2  = [Xp2,  tr2(:,t2)];
                    uXp2 = [uXp2; u0];
                end
                Xp2  = [Xp2,  tr2(:,D.T)];
                uXp2 = [uXp2; u0];
            end
            if size(Xp2,2) < 4, continue; end

            lbw_k2 = median(std(Xp2',0,1)) * 2.0;

            tvar2 = zeros(Rst2,1);
            for rr = 1:Rst2
                tr = squeeze(D.Zfull(:,:,rip_st_idx2(rr)))';
                tvar2(rr) = sum(var(tr,0,1));
            end
            % Variance-weighted random sampling across ALL stim ripples
            % Higher variance = more likely to be selected, but all are eligible
            tvar2_w  = tvar2 / (sum(tvar2) + 1e-12);   % normalize to probabilities
            cum_w    = cumsum(tvar2_w);
            n_pick   = min(n_st2, Rst2);
            show_rr2 = zeros(n_pick, 1);
            available = true(Rst2, 1);
            for pp = 1:n_pick
                % Resample from remaining ripples
                w_rem = tvar2_w .* available;
                w_rem = w_rem / (sum(w_rem) + 1e-12);
                r_pick = find(rand() <= cumsum(w_rem), 1, 'first');
                if isempty(r_pick), r_pick = find(available, 1); end
                show_rr2(pp) = r_pick;
                available(r_pick) = false;
            end

            for rr = show_rr2'
                tr    = squeeze(D.Zfull(:,:,rip_st_idx2(rr)))';
                tr_z  = (tr-mean(D.Cn,1))./(std(D.Cn,0,1)+1e-12);
                tr_n  = tr_z./(vecnorm(tr_z,2,2)+1e-12);
                tr_rf = tr_n-(tr_n*D.uRate)*D.uRate';
                tr_rn = tr_rf./(vecnorm(tr_rf,2,2)+1e-12);
                u_anc  = ust2(rr,:);
                z_mean = mean(tr_rn,1);

                % Jacobian weighted by K-space proximity to trial mean
                dist2_k = sum((Xp2 - z_mean').^2, 1)';
                w_k     = exp(-dist2_k/(2*lbw_k2^2));
                if sum(w_k)<1e-6, continue; end
                Wk    = diag(w_k);
                Xtaug = [Xp2; ones(1,size(Xp2,2))];
                lamJ  = 1e-3*trace(Xtaug*Wk*Xtaug')/(size(Xp2,1)+1);
                J_aug = (uXp2'*Wk*Xtaug')/(Xtaug*Wk*Xtaug'+lamJ*eye(size(Xtaug,1)));
                J     = J_aug(:,1:size(Xp2,1));

                dZ         = tr_rn - z_mean;
                u_traj_raw = u_anc + dZ*J';
                all_traj_data{end+1} = struct('u_traj', u_traj_raw);
            end
        end

        if ~isempty(all_traj_data)
            for ti = 1:length(all_traj_data)
                td     = all_traj_data{ti};
                u_traj = td.u_traj;
                t_orig = 1:D.T;
                t_fine = linspace(1,D.T,30);
                xf = interp1(t_orig,u_traj(:,1),t_fine,'pchip');
                yf = interp1(t_orig,u_traj(:,2),t_fine,'pchip');
                ph = plot(ax,xf,yf,'-','Color',[0 0 0],'LineWidth',2.0);
                ph.Color(4) = 0.85;
                % Open circle at start
                plot(ax,u_traj(1,1), u_traj(1,2), 'o','Color',[0 0 0],...
                    'MarkerSize',5,'MarkerFaceColor','w','LineWidth',1.5)
                % Arrowhead at end: tiny quiver along final trajectory direction
                dx = xf(end) - xf(end-2);
                dy = yf(end) - yf(end-2);
                len = sqrt(dx^2+dy^2)+1e-12;
                % Normalise to a fixed display length then let MaxHeadSize=1 fill it
                arrow_len = mean([range(xlim(ax)) range(ylim(ax))]) * 0.03;
                dx_n = dx/len * arrow_len;
                dy_n = dy/len * arrow_len;
                ah = quiver(ax, xf(end)-dx_n, yf(end)-dy_n, dx_n, dy_n, 0, ...
                    'Color',[0 0 0],'LineWidth',2.0,...
                    'MaxHeadSize', 1.5, 'AutoScale','off');
                ah.Color(4) = 0.85;
            end
        end

        hold(ax,'off');
        rF_str = '';
        for m = 1:nM
            idx2 = find(D.mot_valid==sel_cl(m));
            if ~isempty(idx2)
                rF_str = [rF_str sprintf('M%d:%.2f  ',sel_cl(m),D.rf_valid(idx2))]; %#ok
            end
        end
        title(ax, sprintf('Session %d  |  %s', idx, strtrim(rF_str)),...
            'FontSize',9,'FontWeight','bold','Color',[0.15 0.15 0.15]);
    end

    function on_export()
        idx = state.ses_idx;
        efig = figure('Color','w','Position',[100 100 700 650]);
        eax  = axes(efig,'Color',[0.97 0.97 0.98],'FontSize',9,...
            'GridAlpha',0.12,'Box','off');
        xlabel(eax,'UMAP 1'); ylabel(eax,'UMAP 2');

        % Copy axes content
        kids = ax.Children;
        copyobj(kids, eax);
        axis(eax,'tight');
        title(eax, ax.Title.String,'FontSize',10,'FontWeight','bold');

        fname = sprintf('manifold_flow_S%d_%s.pdf', idx, datestr(now,'yyyymmdd_HHMM'));
        exportgraphics(efig, fname, 'ContentType','vector');
        fprintf('Exported: %s\n', fname);
    end

end % main function


%% ══════════════════════════════════════════════════════════════════════════
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

    % Ztensor: rate-free trial weight embedding (R x K)
    C     = ses(i).tensorModel.C;
    Cz    = zscore(C,[],1);
    Cn    = Cz ./ (vecnorm(Cz,2,2)+1e-12);
    uRate = mean(Cn(spont_all,:),1)';
    uRate = uRate/(norm(uRate)+1e-12);
    Cn_rf = Cn - (Cn*uRate)*uRate';
    Ztensor = Cn_rf ./ (vecnorm(Cn_rf,2,2)+1e-12);   % R x K

    % UMAP on trial-mean vectors — clean cluster separation
    tmpl = sprintf('umap_template_ses%d.mat',i);
    [usp_all,~] = run_umap(Ztensor(spont_all,:),...
        'n_neighbors',15,'min_dist',0.10,'metric','cosine',...
        'n_components',2,'randomize',false,...
        'save_template_file',tmpl,'verbose','none');
    ust_all = zeros(0,2);
    if any(stim_vec)
        ust_all = run_umap(Ztensor(stim_vec,:),...
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
    D.Zfull       = Zfull;
    D.Cn          = Cn;
    D.uRate       = uRate;
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