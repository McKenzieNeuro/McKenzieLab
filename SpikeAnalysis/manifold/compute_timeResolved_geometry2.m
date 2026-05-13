function [stats] = compute_timeResolved_geometry2(Z, replayID, stim)

% Z         : K x T x R  per-timepoint normalized latent trajectories
%             (rate axis removed, unit norm — output of Znorm pipeline)
% replayID  : R x 1  motif cluster labels
% stim      : R x 1  logical, true = stimulated ripple
%
% Dynamics are quantified as mean squared velocity in full K-space.
% No dimensionality reduction or regression is performed — this avoids
% the ill-conditioning of fitting a K×K dynamics matrix from O(R_m)
% observations, and avoids the arbitrary choice of projection plane.
%
% Primary summary statistic:
%   totalSpread(m) = mean over timebins of mean squared velocity
%                  = mean kinetic energy of the trajectory in K-space

[D,T,R] = size(Z);
motifs   = unique(replayID);
nM       = length(motifs);

for s = 1:2

    if s==1, kp = ~stim; outfield = 'spont';
    else,    kp =  stim; outfield = 'stim';
    end

    % Time-resolved
    flowEnergy  = nan(nM,T);   % mean squared velocity at each timebin
    pathLength  = nan(nM,1);   % mean total path length across ripples

    % Motif-level summaries
    totalSpread = nan(nM,1);   % mean flowEnergy across timebins
    effRank     = nan(nM,1);   % effective dimensionality in K-space

    for m = 1:nM

        idx = replayID == motifs(m) & kp;
        Zk  = Z(:,:,idx);
        Rk  = size(Zk,3);
        if Rk < 5, continue; end

        % ── Per-timebin mean squared velocity ─────────────────────────────
        for t = 2:T-1
            V = zeros(D,Rk);
            for r = 1:Rk
                V(:,r) = (Zk(:,t+1,r) - Zk(:,t-1,r)) / 2;
            end
            flowEnergy(m,t) = nanmean(sum(V.^2, 1));
        end

        % ── Mean total path length ─────────────────────────────────────────
        % Sum of step distances along each ripple trajectory, averaged across ripples
        pl = zeros(Rk,1);
        for r = 1:Rk
            for t = 1:T-1
                pl(r) = pl(r) + norm(Zk(:,t+1,r) - Zk(:,t,r));
            end
        end
        pathLength(m) = mean(pl);

        % ── Motif-level summary ────────────────────────────────────────────
        totalSpread(m) = nanmean(flowEnergy(m,:));

        % ── Effective rank (full K-space) ──────────────────────────────────
        Xfull = reshape(Zk,D,[]);
        Xfull = Xfull - nanmean(Xfull,2);
        evC   = eig(cov(Xfull'));
        evC   = evC(evC>0);
        pC    = evC/nansum(evC);
        effRank(m) = exp(-sum(pC.*log(pC+1e-12)));

    end

    % ── Store outputs ──────────────────────────────────────────────────────
    stats.(outfield).flowEnergy   = flowEnergy;
    stats.(outfield).totalSpread  = totalSpread;
    stats.(outfield).pathLength   = pathLength;
    stats.(outfield).effectiveRank= effRank;
    stats.(outfield).motifID      = motifs;

end

% ── Delta metrics (stim - spont) ──────────────────────────────────────────
stats.delta.motifID           = stats.spont.motifID;
stats.delta.delta_totalSpread = stats.stim.totalSpread(:) - stats.spont.totalSpread(:);
stats.delta.delta_pathLength  = stats.stim.pathLength(:)  - stats.spont.pathLength(:);
stats.delta.delta_effRank     = stats.stim.effectiveRank(:) - stats.spont.effectiveRank(:);

end