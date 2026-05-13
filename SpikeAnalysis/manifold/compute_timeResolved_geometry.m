function [stats] = compute_timeResolved_geometry(Z, replayID, stim)

% Z : D x T x R latent trajectories (fixed representational space)
% replayID : motif labels
% stim : logical vector (true = stim condition)

[D,T,R] = size(Z);
motifs = unique(replayID);
nM = length(motifs);

for s = 1:2
    
    if s==1
        kp = ~stim; outfield='spont';
    else
        kp = stim;  outfield='stim';
    end
    
    % Time-resolved
    alpha_t      = nan(nM,T);
    omega_t      = nan(nM,T);
    spiralIndex  = nan(nM,T);
    rotFrac      = nan(nM,T);
    transVarAll  = nan(nM,T);
    
    % Energy metrics
    rotEnergy    = nan(nM,T);
    symEnergy    = nan(nM,T);
    radGrowth    = nan(nM,T);
    spiralAmp    = nan(nM,T);
    
    % Means
    alphaMean = nan(nM,1);
    omegaMean = nan(nM,1);
    spiralScore = nan(nM,1);
    
    rotFracTotal = nan(nM,1);
    totalSpread  = nan(nM,1);
    effRank      = nan(nM,1);
    
    for m = 1:nM
        
        idx = replayID==motifs(m) & kp;
        Zk = Z(:,:,idx);
        Rk = size(Zk,3);
        if Rk < 5, continue; end
        
        % ---------------------------------
        % Local tangent plane (rank-2)
        % ---------------------------------
        Zall = reshape(Zk,D,[]);
        Zall = Zall - mean(Zall,2);
        
      %  [U,~,~] = svd(Zall,'econ');
      %  U2 = U(:,1:2);   % local 2D plane

      % PROPOSED — maximum rotational variance plane (jPCA)
      Xp =[];Vp =[];
      for r=1:Rk
          tr = Zk(:,:,r);
          Xp = [tr(:,1)];
          Vp = tr(:,2) - tr(:,1);




          for t = 2:T-1
              Xp = [Xp, tr(:,t)];
              Vp = [Vp, (tr(:,t+1) - tr(:,t-1))/2];   % centered
          end
          % Add endpoint with forward difference from t-1
          Xp = [Xp, tr(:,T)];
          Vp = [Vp, tr(:,T) - tr(:,T-1)];



      end
      M = (Vp*Xp')/(Xp*Xp'+1e-6*eye(D));
      Ms = (M-M')/2;
      [Ve,~] = eig(Ms);
      [~,si] = sort(abs(imag(diag(eig(Ms)))),'descend');
      v1=real(Ve(:,si(1))); v2=imag(Ve(:,si(1)));
      v1=v1/(norm(v1)+1e-12);
      v2=v2-(v2'*v1)*v1; v2=v2/(norm(v2)+1e-12);
      U2=[v1 v2];

        
        % Project into 2D
        Z2 = zeros(2,T,Rk);
        for r = 1:Rk
            Z2(:,:,r) = U2' * Zk(:,:,r);
        end
        
        transVar = nan(T,1);
        
        for t = 2:T-1
            
            Xt = squeeze(Z2(:,t,:))';        % Rk x 2
            Vt = squeeze((Z2(:,t+1,:) - Z2(:,t-1,:))/2)';  % Rk x 2
            
            Xt = Xt - mean(Xt,1);
            Vt = Vt - mean(Vt,1);
            msv(t) = mean(sum(Vt.^2, 2));
            % Ridge regression (2x2)
            lambda = 1e-3 * trace(Xt'*Xt)/2;
            B = (Vt' * Xt) / (Xt'*Xt + lambda*eye(2));
            
            S = (B + B')/2;
            A = (B - B')/2;
            
            % Rotation parameter (exact in 2D)
            omega = A(2,1);
            
            % Dominant symmetric eigenvalue
            evS = eig(S);
            alpha = max(real(evS));
            
            alpha_t(m,t) = alpha;
            omega_t(m,t) = omega;
            
            spiralIndex(m,t) = abs(omega)/(abs(alpha)+1e-8);
            rotFrac(m,t) = norm(A,'fro')/(norm(B,'fro')+1e-8);
            
            % Energies
            AX = (A * Xt')';
            SX = (S * Xt')';
            
            Erot = mean(sum(AX.^2,2));
            Esym = mean(sum(SX.^2,2));
            
            rotEnergy(m,t) = Erot;
            symEnergy(m,t) = Esym;
            
            % Radial growth proxy
            g = mean(sum(Xt .* Vt,2));
            radGrowth(m,t) = g;
            
            spiralAmp(m,t) = omega^2 * max(0,g);
            
            % Transverse variance
            C = cov(Xt);
            transVar(t) = trace(C);
            transVarAll(m,t) = transVar(t);
            
        end
        
        % ---------------------------------
        % Motif-level summaries
        % ---------------------------------
        alphaMean(m) = nanmean(alpha_t(m,:));
        omegaMean(m) = nanmean(abs(omega_t(m,:)));
        
        totalSpread(m) = nanmean(msv(2:T-1));
        totalSpread(m) = nanmean(rotEnergy(m,:)) + nanmean(symEnergy(m,:));
        rotFracTotal(m) = ...
            nanmean(rotEnergy(m,:)) / (totalSpread(m)+1e-12);
        
        spiralScore(m) = ...
            nanmean(spiralAmp(m,:));
        
        % Effective rank (full latent space)
        Xfull = reshape(Zk,D,[]);
        Xfull = Xfull - mean(Xfull,2);
        Cfull = cov(Xfull');
        evC = eig(Cfull);
        evC = evC(evC>0);
        pC = evC/sum(evC);
        H = -sum(pC.*log(pC+1e-12));
        effRank(m) = exp(H);
        
    end
    
    % ---------------------------------
    % Store outputs (compatible structure)
    % ---------------------------------
    stats.(outfield).alpha_t = alpha_t;
    stats.(outfield).omega_t = omega_t;
    stats.(outfield).spiralIndex = spiralIndex;
    stats.(outfield).rotFrac = rotFrac;
    stats.(outfield).transverseVar = transVarAll;
    
    stats.(outfield).alphaMean = alphaMean;
    stats.(outfield).omegaMean = omegaMean;
    stats.(outfield).spiralScore = spiralScore;
    stats.(outfield).motifID = motifs;
    
    stats.(outfield).rotEnergy = rotEnergy;
    stats.(outfield).symEnergy = symEnergy;
    stats.(outfield).radialGrowth = radGrowth;
    stats.(outfield).spiralAmp_t = spiralAmp;
    
    stats.(outfield).totalSpread = totalSpread;
    stats.(outfield).rotFracTotal = rotFracTotal;
    stats.(outfield).effectiveRank = effRank;
    
end

  % Ensure vectors are column-aligned
    rot_spont   = stats.spont.rotFracTotal(:);
    rot_stim    = stats.stim.rotFracTotal(:);
    
    alpha_spont = stats.spont.alphaMean(:);
    alpha_stim  = stats.stim.alphaMean(:);
    
    omega_spont = stats.spont.omegaMean(:);
    omega_stim  = stats.stim.omegaMean(:);
    
    spiral_spont = stats.spont.spiralScore(:);
    spiral_stim  = stats.stim.spiralScore(:);
    
    % ---- Δ rotation dominance ----
    delta_rotFracTotal = rot_stim - rot_spont;
    
    % ---- ω/α ratios (protect division) ----
    ratio_spont = omega_spont ./ (alpha_spont + 1e-12);
    ratio_stim  = omega_stim  ./ (alpha_stim  + 1e-12);
    
    delta_omegaAlphaRatio = ratio_stim - ratio_spont;
    
    % ---- Δ spiral amplification ----
    delta_spiralScore = spiral_stim - spiral_spont;
    
    % ---- Store ----
    stats.delta.motifID = stats.spont.motifID;
    stats.delta.delta_rotFracTotal = delta_rotFracTotal;
    stats.delta.delta_omegaAlphaRatio = delta_omegaAlphaRatio;
    stats.delta.delta_spiralScore = delta_spiralScore;
    
    % Optional: store raw ratios for inspection
    stats.delta.omegaAlphaRatio_spont = ratio_spont;
    stats.delta.omegaAlphaRatio_stim  = ratio_stim;

    

end