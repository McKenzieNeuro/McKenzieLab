 clc; rng(2);

%% ============================================================
% PARAMETERS
%% ============================================================
N_CA3 = 200; 
N_motifs = 20; 
r = 10; % manifold dimension
nRipples = 1200; 
T_ripple = 80; 
dt = 0.5; 
steps = T_ripple/dt; 
g_CA3 = 1.2; 
tau_CA3 = 10; 
tau_CA1 = 10; 
tau_I = 20; 
noise_level = 0.002; 
alpha_list = linspace(1.5,0,8); % inhibition sweep

%% ============================================================
% 1. BUILD CA3 LOW-RANK STORED MOTIFS
%% ============================================================
U = randn(N_CA3,r); [U,~] = qr(U,0); 
coeff = randn(r,N_motifs); 
mu = U * coeff; 
W_CA3 = mu * mu' / N_CA3; 
W_CA3 = W_CA3 / max(abs(eig(W_CA3)));

%% ============================================================
% 2. CA1 STRUCTURE
%% ============================================================
r_CA1 = r; 
W_CA1 = U'; % feedforward projection
R = randn(r_CA1); [Q,~] = qr(R); 
J = 1.02 * Q; % slightly expansive rotation

%% ============================================================
% 3. BASELINE HETEROGENEOUS INHIBITION
%% ============================================================
delta_inh = 0.05*randn(r_CA1,1); % latent baseline heterogeneity

%% ============================================================
% 4. DEFINE EXPLICIT OPTO PERTURBATIONS
%% ============================================================
u_opsin_none      = zeros(r_CA1,1);                   % baseline
u_opsin_silence   = -1.0*ones(r_CA1,1);              % hyperpolarize all inhibitory axes
u_opsin_activate  = 1.0 + 0.5*randn(r_CA1,1);        % depolarize heterogeneous

%% ============================================================
% STORAGE FOR METRICS
%% ============================================================
mean_rate = zeros(length(alpha_list),3); 
total_variance = zeros(length(alpha_list),3); 
procrustes_volume = zeros(length(alpha_list),3); 
participation_ratio = zeros(length(alpha_list),3); 
isotropy_ratio = zeros(length(alpha_list),3);

opto_conditions = {'baseline','silence','activate'};

%% ============================================================
% INHIBITION SWEEP WITH EXPLICIT OPTO
%% ============================================================
for b = 1:3
    switch b
        case 1, u_opsin = u_opsin_none;
        case 2, u_opsin = u_opsin_silence;
        case 3, u_opsin = u_opsin_activate;
    end
    
    for a = 1:length(alpha_list)
        alpha_inh = alpha_list(a); 
        fprintf('Condition %s | alpha = %.2f\n',opto_conditions{b},alpha_inh); 
        
        Pop_CA1 = zeros(r_CA1,nRipples); 
        firing_rate = 0; 
        
        for k = 1:nRipples
            m = randi(N_motifs); 
            x = mu(:,m) + 0.1*randn(N_CA3,1); 
            z = zeros(r_CA1,1); 
            y = zeros(r_CA1,1); % vector inhibition
            z_sum = zeros(r_CA1,1); 
            
            for t = 1:steps
                %% ---- CA3 dynamics ----
                dx = (-x + g_CA3 * W_CA3 * tanh(x)) / tau_CA3; 
                x = x + dt*dx + noise_level*randn(N_CA3,1);
                
                %% ---- CA1 feedforward ----
                u = W_CA1 * tanh(x);
                
                %% ---- Vector inhibition dynamics ----
                dy = (-y + z.^2) / tau_I; 
                y = y + dt*dy;
                
                %% ---- CA1 update with heterogeneous inhibition + explicit opsin ----
                dz = (-z + J*z + u - alpha_inh.*(y + delta_inh + u_opsin).*z) / tau_CA1; 
                z = z + dt*dz;
                
                z_sum = z_sum + z; 
                firing_rate = firing_rate + mean(abs(z));
            end
            Pop_CA1(:,k) = z_sum; 
        end
        
        %% ---- Compute mean firing rate ----
        firing_rate = firing_rate / (nRipples*steps); 
        mean_rate(a,b) = firing_rate; 
        
        %% ---- Geometry Metrics ----
        C = cov(Pop_CA1'); 
        eigvals = eig(C); 
        total_variance(a,b) = sum(eigvals); 
        
        % participation ratio
        participation_ratio(a,b) = (sum(eigvals)^2) / sum(eigvals.^2); 
        
        % isotropy (condition number)
        isotropy_ratio(a,b) = max(eigvals) / min(eigvals); 
        
        % Procrustes volume (sqrt of det covariance)
        procrustes_volume(a,b) = sqrt(det(C));
    end
end