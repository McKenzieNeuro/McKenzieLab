%% ============================================================
% INTUITION VISUALIZATION FOR PER-MOTIF jPCA GEOMETRY
%
% Four panels:
%   1. Scatter: omegaMean vs alphaMean per motif (spont)
%      with rotFracTotal as color — are they consistent?
%
%   2. Per-motif jPCA plane trajectories for 4 example motifs
%      (hi-rot, lo-rot, hi-exp, mixed) — what does the plane show?
%
%   3. rotFrac across time (t=2,3,4) per motif — is it stable
%      or does it vary? Stim vs spont overlaid.
%
%   4. PCA plane vs jPCA plane omega comparison (scatter)
%      — how much did the fix matter?
%
% Requires: ses struct with tensor.block(f).tensor_flow
%           Run after compute_timeResolved_geometry reruns.
% ============================================================

ses_show = 1;   % change to session with most motifs
fold_show = 1;

tf = ses(ses_show).tensor.block(fold_show).tensor_flow;

%% ---- Pull spont and stim metrics ----

omega_sp  = tf.spont.omegaMean;
alpha_sp  = tf.spont.alphaMean;
rot_sp    = tf.spont.rotFracTotal;
motif_sp  = tf.spont.motifID;

omega_st  = tf.stim.omegaMean;
alpha_st  = tf.stim.alphaMean;
rot_st    = tf.stim.rotFracTotal;
motif_st  = tf.stim.motifID;

% Matched pairs (motifs present in both)
[common, ia, ib] = intersect(motif_sp, motif_st);
omega_sp_m = omega_sp(ia); alpha_sp_m = alpha_sp(ia); rot_sp_m = rot_sp(ia);
omega_st_m = omega_st(ib); alpha_st_m = alpha_st(ib); rot_st_m = rot_st(ib);

fprintf('Session %d, fold %d: %d spont motifs, %d stim motifs, %d matched\n',...
    ses_show, fold_show, length(motif_sp), length(motif_st), length(common))

%% ---- Recompute jPCA planes for example motifs ----
% (to show actual trajectories in the plane)

nComp = 6;
X  = double(ses(ses_show).binnedPopRipple_tensor);
A  = ses(ses_show).tensorModel.A;
[~,T,R] = size(X);
Z  = zeros(size(A,2),T,R);
for r=1:R, Z(:,:,r) = A'*squeeze(X(:,:,r)); end
Z  = Z(1:nComp,:,:);

stim_vec = ismember(1:R, ses(ses_show).PV.block(fold_show).stim_ix)';
clust_id = ses(ses_show).PV.block(fold_show).kmean_id(:);

% Pick 4 example motifs spanning the omega/rot space
% hi-rot, lo-rot, hi-alpha, lo-both
if length(common) >= 4
    [~,i_hirot] = max(rot_sp_m);
    [~,i_lorot] = min(rot_sp_m);
    [~,i_hialp] = max(alpha_sp_m);
    [~,i_loalp] = min(alpha_sp_m);
    ex_idx = unique([i_hirot i_lorot i_hialp i_loalp]);
    ex_idx = ex_idx(1:min(4,end));
    ex_labels = {'Hi rotFrac','Lo rotFrac','Hi alpha','Lo alpha'};
    ex_labels = ex_labels(1:length(ex_idx));
else
    ex_idx = 1:min(4,length(common));
    ex_labels = arrayfun(@(x) sprintf('Motif %d',x), ex_idx,'uni',0);
end

%% ============================================================
% FIGURE 1: omega vs alpha scatter, colored by rotFrac
%% ============================================================

figure('Position',[50 50 1400 350],'Color','w');

subplot(1,4,1); hold on
scatter(omega_sp, alpha_sp, 60, rot_sp, 'filled', 'MarkerFaceAlpha',0.8)
colormap(gca, hot(256)); cb=colorbar; cb.Label.String='rotFracTotal';
clim([0 1])
xlabel('\omega_{mean}','FontSize',11)
ylabel('\alpha_{mean}','FontSize',11)
title('Spont: \omega vs \alpha','FontSize',11,'FontWeight','bold')
% Add correlation
kp = ~isnan(omega_sp)&~isnan(alpha_sp)&~isnan(rot_sp);
[r1,p1]=corr(omega_sp(kp),rot_sp(kp),'type','spearman');
[r2,p2]=corr(alpha_sp(kp),rot_sp(kp),'type','spearman');
text(0.05,0.95,sprintf('r(\\omega,rot)=%.2f\nr(\\alpha,rot)=%.2f',r1,r2),...
    'Units','normalized','FontSize',9,'VerticalAlignment','top')
box off

%% ============================================================
% FIGURE 1 panel 2: rotFrac time profile (spont vs stim)
%% ============================================================

subplot(1,4,2); hold on

% Pool across all matched motifs
rotFrac_sp_t = tf.spont.rotFrac(ia,:);   % nMatched x T
rotFrac_st_t = tf.stim.rotFrac(ib,:);

t_ax = 1:T;
mu_sp = nanmean(rotFrac_sp_t,1);
se_sp = nanstd(rotFrac_sp_t,[],1)./sqrt(sum(~isnan(rotFrac_sp_t)));
mu_st = nanmean(rotFrac_st_t,1);
se_st = nanstd(rotFrac_st_t,[],1)./sqrt(sum(~isnan(rotFrac_st_t)));

fill([t_ax fliplr(t_ax)],[mu_sp+se_sp fliplr(mu_sp-se_sp)],...
    [0.2 0.4 0.85],'FaceAlpha',0.2,'EdgeColor','none')
plot(t_ax,mu_sp,'-','Color',[0.2 0.4 0.85],'LineWidth',2)

fill([t_ax fliplr(t_ax)],[mu_st+se_st fliplr(mu_st-se_st)],...
    [0.85 0.25 0.1],'FaceAlpha',0.2,'EdgeColor','none')
plot(t_ax,mu_st,'-','Color',[0.85 0.25 0.1],'LineWidth',2)

xlabel('Time bin','FontSize',11)
ylabel('rotFrac','FontSize',11)
title('rotFrac over time','FontSize',11,'FontWeight','bold')
legend({'Spont','','Stim',''},'Location','best','FontSize',9)
xlim([1 T]); box off

%% ============================================================
% FIGURE 1 panel 3: omega stim vs spont per matched motif
%% ============================================================

subplot(1,4,3); hold on

% Scatter matched pairs
scatter(omega_sp_m, omega_st_m, 50, rot_sp_m, 'filled','MarkerFaceAlpha',0.8)
colormap(gca,hot(256)); clim([0 1])
lims=[min([omega_sp_m;omega_st_m]) max([omega_sp_m;omega_st_m])];
plot(lims,lims,'--k','LineWidth',1)
xlabel('\omega_{spont}','FontSize',11)
ylabel('\omega_{stim}','FontSize',11)
title('\omega: stim vs spont','FontSize',11,'FontWeight','bold')
[r3,p3]=corr(omega_sp_m(~isnan(omega_sp_m)&~isnan(omega_st_m)),...
             omega_st_m(~isnan(omega_sp_m)&~isnan(omega_st_m)),'type','spearman');
d_omega = nanmean(omega_st_m-omega_sp_m)/nanstd(omega_st_m-omega_sp_m);
text(0.05,0.95,sprintf('r=%.2f (p=%.3f)\nd=%.2f',r3,p3,d_omega),...
    'Units','normalized','FontSize',9,'VerticalAlignment','top')
box off

%% ============================================================
% FIGURE 1 panel 4: example motif trajectories in jPCA plane
%% ============================================================

subplot(1,4,4); hold on
cols_ex = lines(length(ex_idx));
leg_str = {};

for ei = 1:length(ex_idx)

    m_idx = common(ex_idx(ei));
    kp_sp = clust_id==m_idx & ~stim_vec;
    Zk    = Z(:,:,kp_sp);
    Rk    = size(Zk,3);
    if Rk < 3, continue; end

    % Fit jPCA plane for this motif
    Xp=[]; Vp=[];
    for r=1:Rk
        tr=Zk(:,:,r);
        for t=2:T-1
            Xp=[Xp,tr(:,t)];
            Vp=[Vp,(tr(:,t+1)-tr(:,t-1))/2];
        end
        Xp=[Xp,tr(:,T)];
        Vp=[Vp,tr(:,T)-tr(:,T-1)];
    end
    M  = (Vp*Xp')/(Xp*Xp'+1e-6*eye(nComp));
    Ms = (M-M')/2;
    [Ve,~]=eig(Ms);
    [~,si]=sort(abs(imag(diag(eig(Ms)))),'descend');
    v1=real(Ve(:,si(1))); v2=imag(Ve(:,si(1)));
    v1=v1/(norm(v1)+1e-12);
    v2=v2-(v2'*v1)*v1; v2=v2/(norm(v2)+1e-12);
    U2=[v1 v2];

    % Project and plot each trial
    for r=1:Rk
        p = (U2'*Zk(:,:,r))';   % T x 2
        p = p-mean(p,1);
        plot(p(:,1),p(:,2),'-','Color',[cols_ex(ei,:) 0.3],'LineWidth',0.8)
    end
    % Mean trajectory
    mu_tr = mean(Zk,3);
    pm = (U2'*mu_tr)'; pm=pm-mean(pm,1);
    plot(pm(:,1),pm(:,2),'-','Color',cols_ex(ei,:),'LineWidth',2.5)
    plot(pm(1,1),pm(1,2),'o','Color',cols_ex(ei,:),'MarkerSize',7,...
        'MarkerFaceColor',cols_ex(ei,:))
    plot(pm(end,1),pm(end,2),'s','Color',cols_ex(ei,:),'MarkerSize',7,...
        'MarkerFaceColor',cols_ex(ei,:))

    leg_str{end+1} = sprintf('%s (\\omega=%.2f, rot=%.2f)',...
        ex_labels{ei}, omega_sp_m(ex_idx(ei)), rot_sp_m(ex_idx(ei)));
end

xlabel('jPC_1','FontSize',11); ylabel('jPC_2','FontSize',11)
title('Example motifs in jPCA plane','FontSize',11,'FontWeight','bold')
legend(leg_str,'Location','best','FontSize',7)
axis equal; box off; xline(0,'Color',[0.85 0.85 0.85]); yline(0,'Color',[0.85 0.85 0.85])

sgtitle(sprintf('Per-motif jPCA geometry — Session %d Fold %d',ses_show,fold_show),...
    'FontSize',12,'FontWeight','bold')

%% ============================================================
% FIGURE 2: rotFrac heatmap — motifs x time, sorted by omega
%% ============================================================

figure('Position',[50 450 900 400],'Color','w');

[~,sort_idx] = sort(omega_sp,'descend');
kp_valid = ~all(isnan(tf.spont.rotFrac),2);
sort_idx  = sort_idx(kp_valid(sort_idx));

subplot(1,2,1)
imagesc(tf.spont.rotFrac(sort_idx,:))
colormap(hot(256)); colorbar; clim([0 1])
xlabel('Time bin','FontSize',11)
ylabel('Motifs (sorted by \omega)','FontSize',11)
title('rotFrac — Spont','FontSize',11,'FontWeight','bold')
set(gca,'YTick',[],'FontSize',10)

subplot(1,2,2)
% delta rotFrac (stim - spont) for matched motifs
delta_rf = tf.stim.rotFrac(ib,:) - tf.spont.rotFrac(ia,:);
[~,sort_idx2] = sort(omega_sp_m,'descend');
imagesc(delta_rf(sort_idx2,:))
colormap(gca,redblue(256)); colorbar
sym_lim = max(abs(abs(delta_rf(:))));
if sym_lim>0, clim([-sym_lim sym_lim]); end
xlabel('Time bin','FontSize',11)
ylabel('Matched motifs (sorted by \omega_{spont})','FontSize',11)
title('\Delta rotFrac (stim - spont)','FontSize',11,'FontWeight','bold')
set(gca,'YTick',[],'FontSize',10)

sgtitle('rotFrac heatmaps — motifs sorted by \omega','FontSize',12,'FontWeight','bold')

%% helper — redblue colormap
