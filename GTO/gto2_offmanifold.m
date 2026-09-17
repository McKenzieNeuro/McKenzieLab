function out = gto2_offmanifold(varargin)
% GTO2_OFFMANIFOLD  Recovery of off-manifold energy against a gain-modulated anchor.
%
% A test that imposes a displacement at a manifold point m and scores the recovery
% of its (G,T,O) split cannot work, and not because the estimator is weak: the
% anchor mu is the centroid of the neighbours of z_B, so it TRACKS any on-manifold
% component of the move. A tangential displacement relocates mu with it and is
% annihilated by construction -- exactly the tier-three statement of Methods
% sec:identifiable. Imposing tangential energy and scoring its recovery measures
% anchor drift and reports it as estimator error.
%
% THE QUESTION THAT IS WELL POSED
%   Given a state that has been gain-modulated AND pushed off the manifold, does
%   the off-manifold fraction O recover the off-manifold share of the energy?
%   Here the anchor is itself gain-modulated in the sense that matters: z_B sits
%   off the sampled manifold along the radial direction, and mu is drawn to the
%   manifold points beneath it.
%
% CONSTRUCTION
%   Reference manifold: 2-sphere of radius R centred at the ORIGIN, embedded in
%   F dimensions. Origin-centred so the radial direction IS the normal and the
%   alignment a = ||P_N rho|| is 1: the gain axis is well defined, which is the
%   precondition for asking the question at all (contrast gto2_checks C2).
%
%   Each probe takes a base point m on the sphere and adds
%       nu * ( sgn*sqrt(1-O)*rad  +  sqrt(O)*odir )
%   where rad = m/||m|| (radial = gain), odir is an ambient direction orthogonal
%   to the 3-space the sphere occupies (genuinely novel), and sgn = +1 for
%   amplifying and -1 for suppressive gain. NO TANGENTIAL COMPONENT IS IMPOSED,
%   so nothing in the imposed split can be absorbed by anchor tracking and the
%   ground truth O is well posed at the anchor.
%
% WHAT LIMITS RECOVERY -- ONE MECHANISM, IN CLOSED FORM
%   mu is a kNN centroid inside a convex manifold, so it sits a sagitta
%       delta = (rho^2 / 2R) * (1 + rho^2/3R^2)
%   BELOW the surface, where rho is the tangential neighbourhood extent D.ell
%   (unbiased, (k-1) normalisation) and the bracket is the exact-sphere next order. On an
%   origin-centred sphere that displacement is exactly RADIAL, i.e. collinear
%   with the gain axis. It therefore adds to the gain leg rather than to the
%   novel leg, and O is biased through its denominator. With
%       g = sgn*nu*sqrt(1-O)   (imposed gain leg)
%       b = nu*sqrt(O)         (imposed novel leg)
%       j ~ rho/sqrt(k)        (finite-sample anchor jitter, tangential)
%   the predicted recovery is
%       Opred = ( b^2 + (F-3)*sigma^2 ) / ( (g+delta)^2 + b^2 + (F-3)*sigma^2 + j^2 + sigma^2 )
%   The bias is SIGNED BY THE GAIN DIRECTION: delta adds to |g| under
%   amplification (O under-reported) and cancels part of it under suppression
%   (O over-reported). O is not conservative in general.
%
%   Opred is the formula at the MEAN jitter. Per-probe jitter is exponentially
%   distributed, so it is not the prediction for a median; where g + delta ~ 0 the
%   denominator is dominated by that fluctuation (Gate 6) and the two diverge.
%
% MODEL ARM
%   For each probe the anchor is also placed by Propositions 2-3 alone: at
%   (R - delta)*e + ubar + sigma/sqrt(k)*noise, ubar ~ N(0, rho^2/(2k) I_2) in the
%   tangent plane at the probe's projection p, and O is computed with the exact
%   sphere frame there. O_model has the full per-probe distribution, so its median
%   is compared with the median of O_rec and its mean with the mean.
%
% CONDITIONING
%   C6 = median( k*||r||^2 / ell^2 ) is the observable Gate-6 index: ||r||^2 is
%   approximately (g+delta)^2 + b^2 + ell^2/k, so C6 - 1 >> 1 is Eq. 25.
%
%   For a d = 2 manifold rho grows as sqrt(k), so j = rho/sqrt(k) is
%   k-INVARIANT and every k-dependence in the table is sagitta.
%
% ORACLE ARM
%   For each probe the displacement r = z_B - mu is ALSO decomposed against the
%   exact sphere geometry at mu (true radial, true tangent plane, true ambient
%   complement) instead of the estimated Vd, rho, ng. Comparing the two isolates
%   estimator error from anchor bias: if O_rec = O_oracle the estimator is exact
%   and everything left is the geometry of mu.
%
% OPTIONS
%   'F'       30       ambient dimension
%   'R'       5.0      sphere radius
%   'nRef'    4000     reference points
%   'nTest'   300      probes per cell
%   'kgrid'   [40 160] neighbourhood sizes
%   'd'       2        assumed tangent dimension
%   'noise'   0.01     reference/observation noise
%   'nugrid'  [0.3 1.0] total displacement magnitudes
%   'Ogrid'   [0 0.25 0.5 0.75 1.0]  imposed off-manifold fractions
%   'signs'   [1 -1]   amplifying / suppressive gain
%   'seed'    0
%
% RETURNS out.rows (columns in out.colnames) and prints a table per (sgn,k,nu).

ip = inputParser;
ip.addParameter('F', 30);
ip.addParameter('R', 5.0);
ip.addParameter('nRef', 4000);
ip.addParameter('nTest', 300);
ip.addParameter('kgrid', [40 160]);
ip.addParameter('d', 2);
ip.addParameter('noise', 0.01);
ip.addParameter('nugrid', [0.3 1.0]);
ip.addParameter('Ogrid', [0 0.25 0.5 0.75 1.0]);
ip.addParameter('signs', [1 -1]);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;
rng(o.seed);

% fixed embedding: the sphere occupies the first three coordinates, rotated
[Qemb,~] = qr(randn(o.F));
S = Qemb(:,1:3);                                   % the sphere's 3-space
A = [sample_sphere(o.nRef, o.R), zeros(o.nRef, o.F-3)] * Qemb.';
if o.noise > 0, A = A + o.noise*randn(size(A)); end   % noise in ALL F dims:
% confining reference noise to the sphere's own 3 coordinates would put Vd
% exactly inside that subspace and make tangent estimation error identically zero


out.rows = [];
fprintf('\n=== Off-manifold energy recovery, gain-modulated anchor (F=%d, R=%.1f, d=%d) ===\n', ...
    o.F, o.R, o.d);
for sgn = o.signs(:).'
    for k = o.kgrid(:).'
        for nu = o.nugrid(:).'
            fprintf('\n  gain %s   k = %d   nu = %.2f\n', ...
                ternary(sgn > 0, 'amplifying', 'suppressive'), k, nu);
            fprintf('  %7s %7s %7s %7s %7s | %7s %7s | %7s %7s %7s %7s\n', ...
                'O_true', 'O_rec', 'O_orac', 'O_pred', 'O_mdl', 'G_rec', 'T_rec', ...
                '|r|', 'rho', 'delta', 'C6');
            for Otrue = o.Ogrid(:).'
                [Z, ~] = make_probes(A, Qemb, S, sgn, nu, Otrue, o);
                D = gto2_core(A, Z, 'k', k, 'd', o.d);
                keep = ~D.degenerate;

                % oracle: same displacement, exact geometry at mu
                Oorc = nan(size(Z,1),1);
                for i = 1:size(Z,1)
                    Oorc(i) = oracle_O(D.mu(i,:).', Z(i,:).' - D.mu(i,:).', S);
                end

                rho   = median(D.ell(keep));
                delta = rho^2 / (2*o.R) * (1 + rho^2/(3*o.R^2));
                jit   = rho / sqrt(k);
                g = sgn * nu * sqrt(1 - Otrue);
                b = nu * sqrt(Otrue);
                bnov2 = b^2 + (o.F - 3) * o.noise^2;
                Opred = bnov2 / ((g + delta)^2 + bnov2 + jit^2 + o.noise^2);

                st = rng;                % model draws must not shift the probe stream
                Omdl = anchor_model_O(Z, S, o.R, delta, jit^2, o.noise/sqrt(k));
                rng(st);

                mO = median(D.O(keep));  mOo = median(Oorc(keep));
                mG = median(D.G(keep));  mT = median(D.T(keep));
                mr = median(D.rnorm(keep));  ma = median(D.align(keep));
                C6 = median(k * D.rnorm(keep).^2 ./ D.ell(keep).^2);

                out.rows(end+1,:) = [sgn k nu Otrue mO mOo Opred mG mT mr rho delta jit ma ...
                                     mean(D.O(keep)) median(Omdl(keep)) mean(Omdl(keep)) C6]; %#ok<AGROW>
                fprintf('  %7.2f %7.3f %7.3f %7.3f %7.3f | %7.3f %7.3f | %7.3f %7.3f %7.3f %7.1f\n', ...
                    Otrue, mO, mOo, Opred, median(Omdl(keep)), mG, mT, mr, rho, delta, C6);
            end
        end
    end
end
out.colnames = {'sgn','k','nu','O_true','O_rec','O_oracle','O_pred', ...
                'G_rec','T_rec','rnorm','rho','delta','jitter','align', ...
                'O_rec_mean','O_model_med','O_model_mean','C6'};

% --- summary ---------------------------------------------------------
X = out.rows;
estErr  = abs(X(:,5) - X(:,6));       % estimator vs exact geometry
predErr = abs(X(:,5) - X(:,7));       % median recovered vs formula at mean jitter
mdlErr  = abs(X(:,5) - X(:,16));      % median recovered vs median of the model
mdlErrM = abs(X(:,15) - X(:,17));     % mean recovered vs mean of the model
rawErr  = X(:,5) - X(:,4);            % recovered vs imposed (signed)
amp  = X(:,1) > 0;
sgnd = X(:,4) < 1;                    % at O_true = 1 the gain leg is zero: sign undefined
fprintf('\n  max |O_rec - O_oracle|          = %.4f   (estimator error: Vd, rho, ng)\n', max(estErr));
fprintf('  max |O_rec - O_pred|            = %.4f   (median vs formula at mean jitter)\n', max(predErr));
fprintf('  max |O_rec - O_model|, medians  = %.4f\n', max(mdlErr));
fprintf('  max |O_rec - O_model|, means    = %.4f\n', max(mdlErrM));
if any(amp & sgnd) && any(~amp & sgnd)
    fprintf('  O_true < 1 only -- mean signed (O_rec - O_true): amplifying %+.4f, suppressive %+.4f\n', ...
        mean(rawErr(amp & sgnd)), mean(rawErr(~amp & sgnd)));
    fprintf('                  -- mean |O_rec - O_true|:      amplifying %.4f, suppressive %.4f\n', ...
        mean(abs(rawErr(amp & sgnd))), mean(abs(rawErr(~amp & sgnd))));
end
fprintf('  mean |O_rec - O_true| (all cells) = %.4f  (anchor bias, not estimator error)\n', mean(abs(rawErr)));
edges = [0 2 5 20 Inf];
fprintf('  by conditioning C6:   bin        n   mean|O_rec-O_true|  mean|O_rec-O_pred|  mean|O_rec-O_model|\n');
for b = 1:numel(edges)-1
    sel = X(:,18) >= edges(b) & X(:,18) < edges(b+1);
    if ~any(sel), continue, end
    fprintf('                   [%4g,%4g)  %3d   %18.4f  %18.4f  %19.4f\n', edges(b), edges(b+1), ...
        sum(sel), mean(abs(rawErr(sel))), mean(predErr(sel)), mean(mdlErr(sel)));
end
fprintf('\n');
end

% ---------------------------------------------------------------------
function Om = anchor_model_O(Z, S, R, delta, tau2, sigMu)
% O under Propositions 2-3 alone (see header, MODEL ARM).
[n, F] = size(Z);
Om = nan(n,1);
for i = 1:n
    z  = Z(i,:).';
    c  = S.'*z;  e3 = c / norm(c);
    T3 = null(e3.');
    mu3 = (R - delta)*e3 + T3*(sqrt(tau2/2)*randn(2,1));
    mu  = S*mu3;
    if sigMu > 0, mu = mu + sigMu*randn(F,1); end
    Om(i) = oracle_O(mu, z - mu, S);
end
end

% =====================================================================
function [Z, Mtrue] = make_probes(A, Qemb, S, sgn, nu, Otrue, o) %#ok<INUSL>
% Probes with NO imposed tangential component: a radial (gain) leg and an
% ambient novel leg, in the prescribed energy ratio.
n = o.nTest;
base = sample_sphere(n, o.R);
Z = zeros(n, o.F);  Mtrue = zeros(n, o.F);
for i = 1:n
    m   = Qemb * [base(i,:).'; zeros(o.F-3,1)];
    rad = m / norm(m);
    ov  = randn(o.F,1);
    ov  = ov - S*(S.'*ov);                 % orthogonal to the sphere's 3-space
    ov  = ov / norm(ov);
    Mtrue(i,:) = m.';
    Z(i,:) = (m + nu*(sgn*sqrt(1-Otrue)*rad + sqrt(Otrue)*ov)).';
end
if o.noise > 0, Z = Z + o.noise*randn(size(Z)); end
end

% ---------------------------------------------------------------------
function Ofrac = oracle_O(mu, r, S)
% Off-manifold fraction of r at anchor mu using the EXACT sphere geometry:
% radial = the mu direction inside S, tangent plane = its complement in S,
% novel = everything orthogonal to S. No estimated quantities.
p   = S*(S.'*mu);
rad = p / norm(p);
B   = S - rad*(rad.'*S);
[Tb,~] = qr(B, 0);  Tb = Tb(:,1:2);
r_nov = r - Tb*(Tb.'*r) - (rad.'*r)*rad;
Ofrac = sum(r_nov.^2) / (sum(r.^2) + eps);
end

% ---------------------------------------------------------------------
function X = sample_sphere(n, R)
v = randn(n,3);  v = v ./ vecnorm(v,2,2);
X = R * v;
end

% ---------------------------------------------------------------------
function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
