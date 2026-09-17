function out = gto2_paramsweep(mode, varargin)
% GTO2_PARAMSWEEP  Free-parameter dependence of the anchor bias, the tangent-rotation
%                  term of Gate 2, and the Eq. 24 recovery prediction.
%
%   gto2_paramsweep('anchor')    % ell, delta, jitter vs nRef / R / cap / gain / noise   (Table S9)
%   gto2_paramsweep('rotation')  % ||Vd' n||^2 vs nRef, k and R, over seeds              (Table S4)
%   gto2_paramsweep('cor3')      % Eq. 24 vs the exact form vs a Props 2-3 anchor model (Table S16)
%   gto2_paramsweep('all')
%
% ELL CONVENTION (used throughout the suite)
%   ell^2 is the TANGENTIAL second moment of the k neighbours about their centroid,
%   normalised by (k-1). Dividing by k underestimates E||u||^2 by (k-1)/k, which on
%   its own produced a spurious k/(k-1) excess (1.026 at k = 40) in every
%   delta/(ell^2/2R) ratio of the earlier tables.
%
% ANCHOR MODE. For a 2-sphere of radius R uniformly sampled with nRef points,
%       ell^2 = 2*k*R^2/nRef,  delta = ell^2/(2R) = k*R/nRef,  jitter = ell/sqrt(k).
%   On a closed sphere at fixed nRef and k the whole configuration scales with R, so
%   delta/R = k/nRef is invariant; delta grows with R only relative to a FIXED
%   absolute probe magnitude or noise. With ell measured tangentially, the exact
%   sphere geometry adds a known next-order term,
%       delta = (ell^2/2R) * (1 + ell^2/(3R^2)),
%   (uniform-disc neighbourhood: E rho^4 = (4/3) ell^4). 'ratio2' uses this form.
%   The gain and noise blocks test two channels that are absent in this
%   construction: radial probes on an origin-centred sphere leave the kNN ordering
%   invariant, and reference noise enters as fluctuation and through the measured ell.
%
% ROTATION MODE. The PCA plane of a sampled spherical patch is the least-squares
%   fit, whose tilt is set by the fluctuation of Cov(u, |u|^2), not by the centroid
%   offset alone. For a uniform disc of radius a this gives, per tangent direction,
%   Var(slope) = a^2/(12 k R^2), hence
%       E ||Vd'n||^2 = ell^2/(3 k R^2) = 2/(3 nRef),
%   independent of k and R, and (slopes ~ Gaussian, d = 2) an exponential law with
%   median = ln2 * mean = 0.462/nRef. Every cell is reseeded, so the R block uses
%   identical draws at every R: its rows must agree to rounding (the construction is
%   exactly scale-covariant), and any difference would be an absolute-scale leak.
%
% COR3 MODE. Eq. 24 (first order) is compared against (i) the exact Prop. 1 form at
%   the predicted delta and tau^2 and (ii) an anchor MODEL: for each probe the anchor
%   is placed at (R - delta)*e + ubar, ubar ~ N(0, tau^2/2 I_2) in the tangent plane
%   (Props 2-3 and nothing else), and O is computed with the exact sphere frame at
%   that anchor. Comparisons are mean-to-mean: per-probe tau^2 is exponentially
%   distributed, so the median of O is not the formula at the mean tau^2.
%
% RETURNS out.<sweep>.rows with columns in out.<sweep>.colnames.

if nargin < 1, mode = 'all'; end
out = struct();
switch lower(mode)
    case 'anchor',   out.anchor   = sweep_anchor(varargin{:});
    case 'rotation', out.rotation = sweep_rotation(varargin{:});
    case 'cor3',     out.cor3     = sweep_cor3(varargin{:});
    case 'all'
        out.anchor   = sweep_anchor(varargin{:});
        out.rotation = sweep_rotation(varargin{:});
        out.cor3     = sweep_cor3(varargin{:});
    otherwise, error('mode must be anchor, rotation, cor3 or all');
end
end

% =====================================================================
function s = sweep_anchor(varargin)
ip = inputParser;
ip.addParameter('F', 30);
ip.addParameter('R', 5.0);
ip.addParameter('nRef', 4000);
ip.addParameter('k', 40);
ip.addParameter('nTrial', 250);
ip.addParameter('seed', 3);
ip.parse(varargin{:});
o = ip.Results;

s.rows = [];
s.colnames = {'sweep','value','ell','delta','deltaPred','ratio','jitter', ...
              'deltaPred2','ratio2','ratio2SE','jitterPred'};
fprintf('\n=== Anchor bias vs the free parameters (2-sphere) ===\n');
fprintf('  defaults: F = %d, R = %.1f, nRef = %d, k = %d, %d trials/cell; every cell reseeded\n', ...
    o.F, o.R, o.nRef, o.k, o.nTrial);
fprintf('  ratio = delta/(ell^2/2R);  ratio2 = delta/[(ell^2/2R)(1 + ell^2/3R^2)];  SE of ratio2\n');

blocks = { ...
  'nRef',  [500 1000 2000 4000 8000 16000]; ...
  'R',     [3 5 8 12 20]; ...
  'cap',   [0.3 0.6 1.0 1.5 pi]; ...
  'gain',  [-1.0 -0.5 -0.2 0 0.2 0.5 1.0 2.0]; ...
  'noise', [0 0.01 0.05 0.1 0.2] };

for bi = 1:size(blocks,1)
    name = blocks{bi,1};  vals = blocks{bi,2};
    fprintf('\n  -- %s --\n', name);
    fprintf('  %8s %9s %10s %10s %7s %7s %7s %9s %9s\n', ...
        name, 'ell', 'delta', 'ell^2/2R', 'ratio', 'ratio2', 'SE', 'jitter', 'ell/sqrtk');
    for v = vals
        p = o;  cap = pi;  gain = 0;  noise = 0;
        switch name
            case 'nRef',  p.nRef = v;
            case 'R',     p.R    = v;
            case 'cap',   cap    = v;
            case 'gain',  gain   = v;
            case 'noise', noise  = v;
        end
        rng(o.seed);
        [ell, delta, jit, dSE] = anchor_cell(p, cap, gain, noise);
        pred  = ell^2 / (2*p.R);
        pred2 = pred * (1 + ell^2/(3*p.R^2));
        jp    = ell / sqrt(p.k);
        s.rows(end+1,:) = [bi v ell delta pred delta/pred jit pred2 delta/pred2 dSE/pred2 jp]; %#ok<AGROW>
        fprintf('  %8.3f %9.4f %10.5f %10.5f %7.3f %7.3f %7.3f %9.5f %9.5f\n', ...
            v, ell, delta, pred, delta/pred, delta/pred2, dSE/pred2, jit, jp);
    end
end
fprintf(['\n  ratio2 removes the tangent-projection term; residual departures should be\n' ...
         '  within a few SE. In the R block delta/R is invariant (scale covariance).\n\n']);
end

% ---------------------------------------------------------------------
function [ell, delta, jit, dSE] = anchor_cell(o, capang, gain, noise)
ell2 = zeros(o.nTrial,1); rad = zeros(o.nTrial,1); jit2 = zeros(o.nTrial,1);
u = [0 0 1].';
for t = 1:o.nTrial
    if capang >= pi
        v = randn(o.nRef,3); v = v ./ vecnorm(v,2,2);
    else
        ca = cos(capang);
        cz = ca + (1-ca)*rand(o.nRef,1);
        st = sqrt(1 - cz.^2);
        ph = 2*pi*rand(o.nRef,1);
        v  = [st.*cos(ph), st.*sin(ph), cz];
    end
    X = [o.R*v, zeros(o.nRef, o.F-3)];
    if noise > 0, X = X + noise*randn(size(X)); end
    z = zeros(1,o.F); z(1:3) = (o.R + gain)*u.';
    [~, ord] = sort(sum((X - z).^2, 2));
    loc = X(ord(1:o.k), :);
    mu  = mean(loc, 1);
    m   = zeros(1,o.F); m(1:3) = o.R*u.';
    tan_ = loc(:,1:2) - mu(1:2);
    ell2(t) = sum(sum(tan_.^2, 2)) / (o.k - 1);     % tangential, unbiased
    rad(t)  = -(mu(1:3) - m(1:3)) * u;               % inward positive
    jit2(t) = sum(mu(1:2).^2);
end
ell = sqrt(mean(ell2));  delta = mean(rad);  jit = sqrt(mean(jit2));
dSE = std(rad) / sqrt(o.nTrial);
end

% =====================================================================
function s = sweep_rotation(varargin)
ip = inputParser;
ip.addParameter('F', 30);
ip.addParameter('nTest', 1000);
ip.addParameter('seeds', 0:4);
ip.parse(varargin{:});
o = ip.Results;
nS = numel(o.seeds);

cells = [1  1000  5 40; 1  2000  5 40; 1  4000  5 40; 1  8000  5 40; 1 16000  5 40; ...
         2  4000  5 20; 2  4000  5 40; 2  4000  5 80; 2  4000  5 160; ...
         3  4000  3 40; 3  4000  5 40; 3  4000 12 40; 3  4000 20 40];
bname = {'nRef (R = 5, k = 40)', 'k (nRef = 4000, R = 5)', 'R (nRef = 4000, k = 40; identical draws at every R)'};

s.rows = [];
s.colnames = {'block','nRef','R','k','leak','leak_times_nRef', ...
              'leakMean','medxN_sd','meanxN','meanxN_sd'};
fprintf('\n=== Gate 2: tangent rotation, ||Vd''n||^2 (%d probes/cell, seeds %s) ===\n', ...
    o.nTest, mat2str(o.seeds));
fprintf('  columns are x nRef, mean (sd) over seeds. Prediction: mean 2/3 = %.3f, median 2ln2/3 = %.3f\n', ...
    2/3, 2*log(2)/3);
lastB = 0;
for c = 1:size(cells,1)
    blk = cells(c,1); nRef = cells(c,2); R = cells(c,3); k = cells(c,4);
    if blk ~= lastB
        fprintf('\n  -- %s --\n', bname{blk});
        fprintf('  %6s %6s %5s %12s %8s %8s %8s %8s\n', 'nRef','R','k','median','(sd)','mean','(sd)','mean/pred');
        lastB = blk;
    end
    med = zeros(nS,1);  mn = zeros(nS,1);
    for si = 1:nS
        rng(o.seeds(si));
        Lv = rotation_cell(nRef, R, k, o);
        med(si) = median(Lv);  mn(si) = mean(Lv);
    end
    s.rows(end+1,:) = [blk nRef R k mean(med) mean(med)*nRef mean(mn) std(med*nRef) ...
                       mean(mn)*nRef std(mn*nRef)]; %#ok<AGROW>
    fprintf('  %6d %6.1f %5d %12.3f %8.3f %8.3f %8.3f %8.3f\n', nRef, R, k, ...
        mean(med)*nRef, std(med*nRef), mean(mn)*nRef, std(mn*nRef), mean(mn)*nRef/(2/3));
end
fprintf('\n');
end

% ---------------------------------------------------------------------
function Lv = rotation_cell(nRef, R, k, o)
[Qemb,~] = qr(randn(o.F));
S = Qemb(:,1:3);
v = randn(nRef,3); v = v ./ vecnorm(v,2,2);
A = [R*v, zeros(nRef, o.F-3)] * Qemb.';

b = randn(o.nTest,3); b = b ./ vecnorm(b,2,2);
Z = zeros(o.nTest, o.F);  N = zeros(o.nTest, o.F);
for i = 1:o.nTest
    m = Qemb * [R*b(i,:).'; zeros(o.F-3,1)];
    N(i,:) = (m / norm(m)).';                  % true normal at the base point
    ov = randn(o.F,1);  ov = ov - S*(S.'*ov);  ov = ov / norm(ov);
    Z(i,:) = (m + 0.5*ov).';                   % pure off-manifold probe (kNN-invariant)
end

Lv = zeros(o.nTest,1);
[~, ord] = sort(pdist2(Z, A), 2);
for i = 1:o.nTest
    loc = A(ord(i,1:k), :);
    Xc  = loc - mean(loc,1);
    [~,~,V] = svd(Xc, 'econ');
    Vd = V(:,1:2);
    Lv(i) = sum((Vd.' * N(i,:).').^2);
end
end

% =====================================================================
function s = sweep_cor3(varargin)
ip = inputParser;
ip.addParameter('F', 30);   ip.addParameter('R', 5.0);
ip.addParameter('nu', 1.0); ip.addParameter('Ostar', 0.5);
ip.addParameter('nTest', 400);
ip.addParameter('cells', [4000 40; 8000 40; 16000 40; 16000 20]);
ip.addParameter('seeds', 0:4);
ip.parse(varargin{:});
o = ip.Results;
nS = numel(o.seeds);
g0 = o.nu*sqrt(1-o.Ostar);  b0 = o.nu*sqrt(o.Ostar);

s.rows = [];
s.colnames = {'nRef','k','delta','tau2','dO','predDelta','ratioDelta','predFull','ratioFull', ...
              'dOmean','predExact','dOmodelMean','ratioMeanFull','ratioMeanExact', ...
              'ratioMeanModel','ratioMeanModel_sd','dOmodelMed','ratioMedModel','ratioMeanDelta'};
fprintf('\n=== Eq. 24 against the exact form and a Props 2-3 anchor model ===\n');
fprintf('  sphere, R = %.1f, F = %d, noiseless, nu = %.2f, O* = %.2f, %d probes, seeds %s\n', ...
    o.R, o.F, o.nu, o.Ostar, o.nTest, mat2str(o.seeds));
fprintf('  delta = ell^2/2R and tau^2 = ell^2/k from D.ell; the exact form and the model use\n');
fprintf('  delta*(1 + ell^2/3R^2). Ratios are measured/predicted, averaged over seeds.\n\n');
fprintf('  %6s %4s %8s %8s | %10s %10s | %8s %8s %8s %7s %7s | %8s\n', ...
    'nRef','k','delta','tau^2','dO med','dO mean','1st:dlt','1st:full','exact', ...
    'model','(sd)','med/mdl');
for c = 1:size(o.cells,1)
    nRef = o.cells(c,1);  k = o.cells(c,2);
    acc = zeros(nS, 12);
    for si = 1:nS
        rng(o.seeds(si));
        [Q,~] = qr(randn(o.F));  S = Q(:,1:3);
        v = randn(nRef,3); v = v ./ vecnorm(v,2,2);
        A = [o.R*v, zeros(nRef,o.F-3)] * Q.';
        b = randn(o.nTest,3); b = b ./ vecnorm(b,2,2);
        Z = zeros(o.nTest,o.F);
        for j = 1:o.nTest
            m   = Q*[o.R*b(j,:).'; zeros(o.F-3,1)];
            rad = m/norm(m);
            ov  = randn(o.F,1); ov = ov - S*(S.'*ov); ov = ov/norm(ov);
            Z(j,:) = (m + g0*rad + b0*ov).';
        end
        D = gto2_core(A, Z, 'k', k, 'd', 2);
        keep = ~D.degenerate;
        ell    = median(D.ell(keep));
        delta  = ell^2/(2*o.R);
        tau2   = ell^2/k;
        delta2 = delta*(1 + ell^2/(3*o.R^2));
        dOmed  = median(D.O(keep)) - o.Ostar;
        dOmean = mean(D.O(keep))   - o.Ostar;
        pd = -2*o.Ostar*sqrt(1-o.Ostar)*delta/o.nu;
        pf = pd - o.Ostar*tau2/o.nu^2;
        pe = b0^2/((g0 + delta2)^2 + b0^2 + tau2) - o.Ostar;
        Om = anchor_model_O(Z, S, o.R, delta2, tau2, 0);
        dMmean = mean(Om(keep)) - o.Ostar;
        dMmed  = median(Om(keep)) - o.Ostar;
        acc(si,:) = [delta tau2 dOmed pd pf dOmean pe dMmean dMmed ...
                     dOmean/pf dOmean/pe dOmean/dMmean];
    end
    a = mean(acc,1);
    rMod = acc(:,12);
    s.rows(end+1,:) = [nRef k a(1) a(2) a(3) a(4) mean(acc(:,3)./acc(:,4)) a(5) ...
                       mean(acc(:,3)./acc(:,5)) a(6) a(7) a(8) a(10) a(11) mean(rMod) ...
                       std(rMod) a(9) mean(acc(:,3)./acc(:,9)) mean(acc(:,6)./acc(:,4))]; %#ok<AGROW>
    fprintf('  %6d %4d %8.5f %8.5f | %10.6f %10.6f | %8.3f %8.3f %8.3f %7.3f %7.3f | %8.3f\n', ...
        nRef, k, a(1), a(2), a(3), a(6), mean(acc(:,6)./acc(:,4)), a(10), a(11), ...
        mean(rMod), std(rMod), mean(acc(:,3)./acc(:,9)));
end
fprintf(['\n  1st:dlt and 1st:full are mean(dO) against Eq. 24 without and with the tangential\n' ...
         '  term; exact is the Prop. 1 form at the predicted delta, tau^2; model propagates\n' ...
         '  Props 2-3 through exact geometry per probe. med/mdl compares medians like-for-like.\n\n']);
end

% ---------------------------------------------------------------------
function Om = anchor_model_O(Z, S, R, delta, tau2, sigMu)
% O under Propositions 2-3 alone. For each probe, p = projection of z_B's sphere
% component onto the sphere; the model anchor is (R - delta)*e + ubar (+ isotropic
% anchor noise sigMu in all F dims), ubar ~ N(0, tau2/2 * I_2) in the tangent plane
% at p. O is then computed with the EXACT sphere frame at the model anchor.
[n, F] = size(Z);
Om = nan(n,1);
for i = 1:n
    z  = Z(i,:).';
    c  = S.'*z;  e3 = c / norm(c);
    T3 = null(e3.');
    mu3 = (R - delta)*e3 + T3*(sqrt(tau2/2)*randn(2,1));
    mu  = S*mu3;
    if sigMu > 0, mu = mu + sigMu*randn(F,1); end
    Om(i) = exact_sphere_O(mu, z - mu, S);
end
end

% ---------------------------------------------------------------------
function Ofrac = exact_sphere_O(mu, r, S)
p   = S*(S.'*mu);
rad = p / norm(p);
B   = S - rad*(rad.'*S);
[Tb,~] = qr(B, 0);  Tb = Tb(:,1:2);
r_nov = r - Tb*(Tb.'*r) - (rad.'*r)*rad;
Ofrac = sum(r_nov.^2) / (sum(r.^2) + eps);
end
