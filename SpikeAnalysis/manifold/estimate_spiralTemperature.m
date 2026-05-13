function [temp] = estimate_spiralTemperature(Z, replayID, stim, varargin)
% ESTIMATE_SPIRALTEMPERATURE  One temperature scalar per ripple.
%
%   temp = estimate_spiralTemperature(Z, replayID, stim)
%
% Temperature is the mean squared speed of a ripple's trajectory in
% the complex jPCA plane:  z(t) = x(t) + i*y(t),  w = dz/dt.
%
%   r = |w|^2  averaged over time, normalised so mean ripple = 1.
%
% The complex velocity also yields the instantaneous spiral parameters:
%   w/z = alpha + i*omega   (complex growth rate)
%
% INPUTS
%   Z         : D x T x R  latent trajectories
%   replayID  : R x 1  motif labels per trajectory
%   stim      : R x 1  logical (true = stim condition)
%
% OPTIONAL
%   'min_trials'   : minimum trials per motif  (default: 5)
%
% OUTPUTS  (all per-ripple vectors are R x 1, original indexing)
%   .r             temperature: normalised |w|^2
%   .logr          log(r)
%   .speed_raw     un-normalised mean |w|^2
%   .alpha         mean Re(w/z)  — expansion rate
%   .omega         mean Im(w/z)  — rotation rate
%   .pitch         alpha / |omega|
%   .rotFrac       |tangential v|^2 / |w|^2
%   .r_polar_mean  mean |z| in jPCA plane
%   .motifLabel    motif ID
%   .isStim        logical
%   .valid         logical
%   .motif         per-motif summary struct

% -- Parse -----------------------------------------------------------------
p = inputParser;
addParameter(p, 'min_trials', 5);
parse(p, varargin{:});
opts = p.Results;

[D, T, R] = size(Z);
motifs     = unique(replayID);
nM         = length(motifs);
stim       = stim(:);
replayID   = replayID(:);

% =========================================================================
% PASS 1 :  Per-motif geometry (centers, jPCA planes)
% =========================================================================
motif_mu_t  = cell(nM, 1);
motif_U2    = cell(nM, 1);
motif_valid = false(nM, 1);
motif_idx   = cell(nM, 1);

for m = 1:nM
    idx_all = find(replayID == motifs(m));
    R_all   = length(idx_all);
    if R_all < opts.min_trials, continue; end
    
    motif_valid(m) = true;
    motif_idx{m}   = idx_all;
    Zk = Z(:, :, idx_all);
    
    % Time-resolved center
    motif_mu_t{m} = mean(Zk, 3);   % D x T
    
    % jPCA plane (maximal rotational variance)
    Xp = []; Vp = [];
    for ri = 1:R_all
        tr = Zk(:, :, ri);
        Xp = [Xp, tr(:,1)]; %#ok
        Vp = [Vp, tr(:,2) - tr(:,1)]; %#ok
        for t = 2:T-1
            Xp = [Xp, tr(:,t)]; %#ok
            Vp = [Vp, (tr(:,t+1) - tr(:,t-1))/2]; %#ok
        end
        Xp = [Xp, tr(:,T)]; %#ok
        Vp = [Vp, tr(:,T) - tr(:,T-1)]; %#ok
    end
    M_dyn = (Vp * Xp') / (Xp * Xp' + 1e-6*eye(D));
    Ms    = (M_dyn - M_dyn') / 2;
    [Ve, ~] = eig(Ms);
    eigvals = diag(eig(Ms));
    [~, si] = sort(abs(imag(eigvals)), 'descend');
    v1 = real(Ve(:, si(1)));
    v2 = imag(Ve(:, si(1)));
    v1 = v1 / (norm(v1) + 1e-12);
    v2 = v2 - (v2'*v1)*v1;
    v2 = v2 / (norm(v2) + 1e-12);
    motif_U2{m} = [v1, v2];
end

% =========================================================================
% PASS 2 :  Complex velocity per ripple
% =========================================================================
speed_raw     = nan(R, 1);
alpha_out     = nan(R, 1);
omega_out     = nan(R, 1);
pitch_out     = nan(R, 1);
rotFrac_out   = nan(R, 1);
r_polar_out   = nan(R, 1);
valid_out     = false(R, 1);

for m = 1:nM
    if ~motif_valid(m), continue; end
    
    idx_all = motif_idx{m};
    mu_t    = motif_mu_t{m};
    U2      = motif_U2{m};
    
    for ri = 1:length(idx_all)
        idx_r  = idx_all(ri);
        z_full = Z(:, :, idx_r);   % D x T
        
        % == Project into complex jPCA plane ===============================
        %   z(t) = x(t) + i*y(t),  centred on motif mean trajectory
        proj = U2' * z_full;            % 2 x T
        cent = U2' * mu_t;              % 2 x T
        dx   = proj - cent;             % 2 x T,  deviation from center
        
        zc = complex(dx(1,:), dx(2,:)); % 1 x T,  complex position
        
        % == Complex velocity  w = dz/dt ===================================
        w = nan(1, T, 'like', zc);
        
        w(1) = zc(2) - zc(1);                      % forward
        for t = 2:T-1
            w(t) = (zc(t+1) - zc(t-1)) / 2;       % centered
        end
        w(T) = zc(T) - zc(T-1);                    % backward
        
        % == Temperature:  mean |w|^2 ======================================
        speed2 = abs(w).^2;                          % |w|^2 per timebin
        speed_raw(idx_r) = mean(speed2);
        
        % == Complex growth rate:  w/z = alpha + i*omega ====================
        %   Only at timebins where |z| is large enough
        growth = nan(1, T);
        v_rad  = nan(1, T);    % radial speed^2
        v_tan  = nan(1, T);    % tangential speed^2
        
        for t = 1:T
            rr = abs(zc(t));
            if rr < 1e-6 || isnan(w(t)), continue; end
            
            g = w(t) / zc(t);              % complex growth rate
            growth(t) = g;
            
            % Decompose |w|^2 = v_rad^2 + v_tan^2
            %   v_rad = Re(w * conj(z)) / |z|
            %   v_tan = Im(w * conj(z)) / |z|
            wz = w(t) * conj(zc(t));
            v_rad(t) = (real(wz) / rr)^2;
            v_tan(t) = (imag(wz) / rr)^2;
        end
        
        alpha_out(idx_r)   = nanmean(real(growth));
        omega_out(idx_r)   = nanmean(imag(growth));
        pitch_out(idx_r)   = nanmean(real(growth)) / ...
                              (abs(nanmean(imag(growth))) + 1e-8);
        
        % Rotation fraction: fraction of KE in tangential direction
        tot = nansum(v_rad + v_tan);
        if tot > 0
            rotFrac_out(idx_r) = nansum(v_tan) / tot;
        end
        
        r_polar_out(idx_r) = mean(abs(zc));
        valid_out(idx_r)   = true;
    end
end

% =========================================================================
% Normalise temperature:  r = speed_raw / pooled_mean
% =========================================================================
ref_speed = nanmean(speed_raw(valid_out));

temp.r            = speed_raw / ref_speed;
temp.logr         = log(max(temp.r, 1e-12));
temp.speed_raw    = speed_raw;
temp.alpha        = alpha_out;
temp.omega        = omega_out;
temp.pitch        = pitch_out;
temp.rotFrac      = rotFrac_out;
temp.r_polar_mean = r_polar_out;
temp.motifLabel   = replayID;
temp.isStim       = stim;
temp.valid        = valid_out;
temp.ref_speed    = ref_speed;

% =========================================================================
% Per-motif summaries
% =========================================================================
temp.motif.motifID        = motifs(:);
temp.motif.mu_t           = motif_mu_t;
temp.motif.U2             = motif_U2;
temp.motif.r_mean_spont   = nan(nM, 1);
temp.motif.r_mean_stim    = nan(nM, 1);
temp.motif.r_median_spont = nan(nM, 1);
temp.motif.r_median_stim  = nan(nM, 1);

for m = 1:nM
    msk_sp = replayID == motifs(m) & ~stim & valid_out;
    msk_st = replayID == motifs(m) &  stim & valid_out;
    if any(msk_sp)
        temp.motif.r_mean_spont(m)   = nanmean(temp.r(msk_sp));
        temp.motif.r_median_spont(m) = nanmedian(temp.r(msk_sp));
    end
    if any(msk_st)
        temp.motif.r_mean_stim(m)   = nanmean(temp.r(msk_st));
        temp.motif.r_median_stim(m) = nanmedian(temp.r(msk_st));
    end
end

temp.motif.delta_r = temp.motif.r_mean_stim - temp.motif.r_mean_spont;
temp.motif.temp_modulation = ...
    (temp.motif.r_mean_stim - temp.motif.r_mean_spont) ./ ...
    (temp.motif.r_mean_stim + temp.motif.r_mean_spont + 1e-12);

end