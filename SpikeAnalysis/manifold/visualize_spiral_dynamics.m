function visualize_spiral_dynamics(Z, replayID, stim, stats, motifID, condition)
% condition: 'spont' or 'stim'

[D,T,R] = size(Z);

% --- Select condition ---
if strcmp(condition,'spont')
    kp = ~stim;
elseif strcmp(condition,'stim')
    kp = stim;
else
    error('condition must be spont or stim');
end

% --- Extract motif ---
idx = replayID==motifID & kp;
Zk = Z(:,:,idx);
Rk = size(Zk,3);

if Rk < 5
    error('Not enough trials');
end

alpha_t = stats.(condition).alpha_t(stats.(condition).motifID==motifID,:);
omega_t = stats.(condition).omega_t(stats.(condition).motifID==motifID,:);

spiralMask = (alpha_t>0) & (omega_t>0);

%% ============================================================
% 1) Timecourse diagnostics
% ============================================================

figure;
subplot(2,1,1)
plot(alpha_t,'LineWidth',2); hold on
yline(0,'k--');
ylabel('\alpha (expansion rate)')
title(['Motif ' num2str(motifID) ' - ' condition])

subplot(2,1,2)
plot(omega_t,'LineWidth',2); hold on
ylabel('\omega (rotation magnitude)')
xlabel('Time')

%% ============================================================
% 2) Dominant 2D plane projection
% ============================================================

% Stack all timepoints across trials
Xall = reshape(Zk, D, T*Rk)';
Xall = Xall - mean(Xall,1);

% PCA
[coeff,~,~] = pca(Xall);
W = coeff(:,1:2);

% Project mean trajectory
mu = mean(Zk,3);
mu2 = (mu' - mean(Xall,1)) * W;

figure;
plot(mu2(:,1), mu2(:,2),'k','LineWidth',2); hold on
scatter(mu2(spiralMask,1), mu2(spiralMask,2),50,'r','filled')

xlabel('PC1'); ylabel('PC2');
title(['Dominant plane - Motif ' num2str(motifID)])
axis equal

%% ============================================================
% 3) Local flow field overlay (optional sparse sampling)
% ============================================================

step = 3;
for t = 2:step:T-1
    
    Xt_prev = squeeze(Zk(:,t-1,:))';
    Xt      = squeeze(Zk(:,t,:))';
    Xt_next = squeeze(Zk(:,t+1,:))';
    
    Xc = Xt - mean(Xt,1);
    Vc = (Xt_next - Xt_prev)/2;
    Vc = Vc - mean(Vc,1);
    
    lambda = 1e-3 * trace(Xc'*Xc)/D;
    J = (Vc' * Xc) / (Xc' * Xc + lambda*eye(D));
    
    z = mu(:,t)';
    v = (mu(:,t+1)-mu(:,t-1))'/2;
    
    z2 = (z - mean(Xall,1)) * W;
    v2 = (v * W);
    
    quiver(z2(1), z2(2), v2(1), v2(2), ...
        'Color',[0 0.4 1],'MaxHeadSize',1,'AutoScale','off');
end

legend({'Mean trajectory','Expanding spiral points','Flow vectors'})