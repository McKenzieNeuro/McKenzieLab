clear; clc; rng(2);

%% ============================================================
% PARAMETERS
%% ============================================================

N_CA3 = 200;
N_motifs = 20;
r = 10;                      % manifold dimension
nRipples = 1200;
T_ripple = 80;
dt = 0.5;
steps = T_ripple/dt;

g_CA3 = 1.2;
tau_CA3 = 10;

tau_CA1 = 10;
tau_I = 20;

noise_level = 0.002;

alpha_list = linspace(1.5,0,8);   % inhibition sweep

%% ============================================================
% 1. BUILD CA3 LOW-RANK STORED MOTIFS
%% ============================================================

U = randn(N_CA3,r);
[U,~] = qr(U,0);

coeff = randn(r,N_motifs);
mu = U * coeff;

W_CA3 = mu * mu' / N_CA3;
W_CA3 = W_CA3 / max(abs(eig(W_CA3)));

%% ============================================================
% 2. CA1 STRUCTURE
%% ============================================================

r_CA1 = r;
W_CA1 = U';

% near-orthogonal slightly expansive rotation
R = randn(r_CA1);
[Q,~] = qr(R);
J = 1.02 * Q;

%% ============================================================
% STORAGE FOR METRICS
%% ============================================================

mean_rate = zeros(size(alpha_list));
total_variance = zeros(size(alpha_list));
procrustes_volume = zeros(size(alpha_list));
participation_ratio = zeros(size(alpha_list));
isotropy_ratio = zeros(size(alpha_list));

%% ============================================================
% INHIBITION SWEEP
%% ============================================================

for a = 1:length(alpha_list)
    
    alpha_inh = alpha_list(a);
    fprintf('Running alpha = %.2f\n',alpha_inh);
    
    Pop_CA1 = zeros(r_CA1,nRipples);
    firing_rate = 0;
    
    for k = 1:nRipples
        
        m = randi(N_motifs);
        x = mu(:,m) + 0.1*randn(N_CA3,1);
        z = zeros(r_CA1,1);
        y = 0;  % scalar inhibition
        
        z_sum = zeros(r_CA1,1);
        
        for t = 1:steps
            
            % ---- CA3 ----
            dx = (-x + g_CA3 * W_CA3 * tanh(x)) / tau_CA3;
            x = x + dt*dx + noise_level*randn(N_CA3,1);
            
            % ---- CA1 feedforward ----
            u = W_CA1 * tanh(x);
            
            % ---- Scalar radial inhibition ----
            dy = (-y + (z'*z)) / tau_I;
            y = y + dt*dy;
            
            dz = (-z + J*z + u - alpha_inh*y*z) / tau_CA1;
            z = z + dt*dz;
            
            z_sum = z_sum + z;
            firing_rate = firing_rate + mean(abs(z));
        end
        
        Pop_CA1(:,k) = z_sum;
    end
    
    firing_rate = firing_rate / (nRipples*steps);
    mean_rate(a) = firing_rate;
    
    % ===== Geometry Metrics =====
    
    C = cov(Pop_CA1');
    eigvals = eig(C);
    
    total_variance(a) = sum(eigvals);
    
    % participation ratio
    participation_ratio(a) = (sum(eigvals)^2) / sum(eigvals.^2);
    
    % isotropy (condition number)
    isotropy_ratio(a) = max(eigvals) / min(eigvals);
    
    % Procrustes volume (sqrt of det covariance)
    procrustes_volume(a) = sqrt(det(C));
    
end

%% ============================================================
% PLOTS
%% ============================================================

figure;

subplot(2,3,1)
plot(alpha_list,mean_rate,'-o')
xlabel('Inhibition strength')
ylabel('Mean CA1 rate')
title('Rate vs inhibition')
set(gca,'XDir','reverse')

subplot(2,3,2)
plot(alpha_list,total_variance,'-o')
xlabel('Inhibition strength')
ylabel('Total variance')
title('Variance')
set(gca,'XDir','reverse')

subplot(2,3,3)
plot(alpha_list,procrustes_volume,'-o')
xlabel('Inhibition strength')
ylabel('Procrustes volume')
title('Volume')
set(gca,'XDir','reverse')

subplot(2,3,4)
plot(alpha_list,participation_ratio,'-o')
xlabel('Inhibition strength')
ylabel('Participation ratio')
title('Dimensionality')
set(gca,'XDir','reverse')

subplot(2,3,5)
plot(alpha_list,isotropy_ratio,'-o')
xlabel('Inhibition strength')
ylabel('Condition number')
title('Isotropy')
set(gca,'XDir','reverse')