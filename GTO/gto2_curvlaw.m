function out = gto2_curvlaw(varargin)
% GTO2_CURVLAW  Direct test of the anchor-displacement law (Methods sec:anchorbias).
%
% Proposition 1 states that the kNN centroid used as an anchor is displaced from
% the manifold by
%       E[mu] - m  =  (ell^2 / 2) * Hvec   + O(ell^3),
% Hvec the MEAN CURVATURE VECTOR and ell^2 = E||u||^2 the squared tangential
% extent of the neighbourhood. Proposition 2 states that the tangential
% fluctuation of the anchor has rms ell/sqrt(k).
%
% This file tests both directly, on the anchor itself. It does NOT call
% gto2_core: the law is a property of mu, and testing it through the full
% decomposition would confound it with tangent-space estimation. The consequences
% for G and O are tested separately in gto2_offmanifold and gto2_anchorbias.
%
% SURFACE FAMILY
%   A quadratic patch through the origin with prescribed principal curvatures,
%   embedded in F dimensions by a fixed random rotation:
%
%     same-normal:  x(u) = ( u1, u2, (k1*u1^2 + k2*u2^2)/2, 0, ... )
%                   both curvatures share one normal direction, so
%                   Hvec = ((k1+k2)/2) * n3   and k1 = -k2 gives a MINIMAL
%                   surface with Hvec = 0 exactly.
%
%     biaxial:      x(u) = ( u1, u2, k1*u1^2/2, k2*u2^2/2, 0, ... )
%                   the curvatures occupy TWO orthogonal normal directions, so
%                   Hvec = (k1/2)*n3 + (k2/2)*n4 and no choice of k1,k2 other
%                   than both zero makes Hvec vanish. This is the case that
%                   matters for O: a normal direction carrying curvature that is
%                   not the gain axis puts anchor displacement into the residual
%                   normal space, i.e. into apparent novelty.
%
%   With d = 2 the prediction is componentwise:
%       same-normal  ->  <mu, n3> = ell^2 * (k1 + k2) / 4
%       biaxial      ->  <mu, n3> = ell^2 * k1 / 4 ,  <mu, n4> = ell^2 * k2 / 4
%
% OPTIONS
%   'F'       30      ambient dimension
%   'nRef'    4000    reference points per trial
%   'a'       1.0     chart radius (points uniform in the disc of radius a)
%   'kgrid'   [20 40 80 160 320]
%   'nTrial'  400     independent reference draws per cell
%   'gain'    0       normal offset of the probe along n3 (tests the coupling
%                     between the probe's own offset and neighbour selection)
%   'noise'   0       isotropic ambient noise
%   'seed'    0
%
% RETURNS out.rows / out.colnames and prints two tables: a curvature sweep at
% fixed k, and a k sweep at fixed curvature.

% Below this magnitude a predicted value is display-zero (rounds to
% 0.000000 at the precision this function prints), so a ratio against it
% is not a measurement of anything -- print '---' instead of a number.
PRED_FLOOR = 1e-9;

ip = inputParser;
ip.addParameter('F', 30);
ip.addParameter('nRef', 4000);
ip.addParameter('a', 1.0);
ip.addParameter('kgrid', [20 40 80 160 320]);
ip.addParameter('nTrial', 400);
ip.addParameter('gain', 0);
ip.addParameter('noise', 0);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;
rng(o.seed);

[Q,~] = qr(randn(o.F));            % fixed embedding; the law is rotation-equivariant
E = Q(:,1:2);                      % tangent directions at the origin
n3 = Q(:,3);  n4 = Q(:,4);         % the two normal directions used below

% ---- curvature sweep at k = 40 --------------------------------------
cases = { %  k1     k2    sameNormal   label
    0.20,  0.20,  true,  'isotropic  k1=k2=0.2, one normal'
    0.20, -0.20,  true,  'MINIMAL    k1=-k2=0.2, one normal  (Hvec = 0)'
    0.40,  0.00,  true,  'anisotropic k1=0.4, k2=0, one normal'
    0.50,  0.50,  true,  'isotropic  k1=k2=0.5, one normal'
    0.20, -0.20,  false, 'biaxial    k1=0.2 on n3, k2=-0.2 on n4'
    0.20,  0.20,  false, 'biaxial    k1=k2=0.2, two normals'
};

k0 = 40;
out.rows = [];
fprintf('\n=== Proposition 1: anchor displacement = (ell^2/2)*Hvec ===\n');
fprintf('  F = %d, nRef = %d, chart radius = %.2f, k = %d, %d trials/cell, probe gain = %.2f\n\n', ...
    o.F, o.nRef, o.a, k0, o.nTrial, o.gain);
fprintf('  %-44s %9s %11s %11s %7s %11s %11s\n', ...
    'surface', 'ell^2', '<mu,n3>', 'predicted', 'ratio', '<mu,n4>', 'predicted');
for c = 1:size(cases,1)
    [k1, k2, sameN, lab] = cases{c,:};
    [m3, m4, l2, ~] = run_cell(k1, k2, sameN, k0, Q, E, n3, n4, o);
    if sameN, p3 = mean(l2)*(k1+k2)/4;  p4 = 0;
    else,     p3 = mean(l2)*k1/4;       p4 = mean(l2)*k2/4;
    end
    if abs(p3) < PRED_FLOOR
        rat_str = '    ---';
    else
        rat_str = sprintf('%7.3f', mean(m3) / p3);
    end
    out.rows(end+1,:) = [k1 k2 sameN k0 mean(l2) mean(m3) p3 mean(m4) p4]; %#ok<AGROW>
    fprintf('  %-44s %9.5f %11.6f %11.6f %7s %11.6f %11.6f\n', ...
        lab, mean(l2), mean(m3), p3, rat_str, mean(m4), p4);
end

% ---- k sweep at fixed isotropic curvature ---------------------------
fprintf('\n=== Proposition 2: rms tangential jitter = ell/sqrt(k) ===\n');
fprintf('  isotropic surface, k1 = k2 = 0.20, one normal\n\n');
fprintf('  %6s %10s %11s %11s %7s %11s %12s %7s\n', ...
    'k', 'ell^2', '<mu,n3>', 'predicted', 'ratio', 'rms jitter', 'ell/sqrt(k)', 'ratio');
for k = o.kgrid(:).'
    [m3, ~, l2, j2] = run_cell(0.20, 0.20, true, k, Q, E, n3, n4, o);
    L  = mean(l2);
    p3 = L * 0.40 / 4;
    rms_jit  = sqrt(mean(j2));
    pred_jit = sqrt(L/k);
    fprintf('  %6d %10.5f %11.6f %11.6f %7.3f %11.5f %12.5f %7.3f\n', ...
        k, L, mean(m3), p3, mean(m3)/p3, rms_jit, pred_jit, rms_jit/pred_jit);
end
out.colnames = {'k1','k2','sameNormal','k','ell2','mu_n3','pred_n3','mu_n4','pred_n4'};
fprintf(['\n  ell^2 scales linearly in k (ell ~ k^(1/d), d = 2), so the bias grows as k\n' ...
         '  while the jitter is k-invariant. Neither is reduced by a larger neighbourhood.\n\n']);
end

% =====================================================================
function [m3, m4, l2, j2] = run_cell(k1, k2, sameN, k, Q, E, n3, n4, o)
% One cell: nTrial independent reference draws, each giving one anchor.
m3 = zeros(o.nTrial,1);  m4 = zeros(o.nTrial,1);
l2 = zeros(o.nTrial,1);  j2 = zeros(o.nTrial,1);
z  = o.gain * n3;                                  % probe, offset along n3
for t = 1:o.nTrial
    % uniform in the disc of radius a
    rr = o.a * sqrt(rand(o.nRef,1));
    th = 2*pi*rand(o.nRef,1);
    u  = [rr.*cos(th), rr.*sin(th)];

    if sameN
        h = [0.5*(k1*u(:,1).^2 + k2*u(:,2).^2), zeros(o.nRef,1)];
    else
        h = [0.5*k1*u(:,1).^2, 0.5*k2*u(:,2).^2];
    end
    X = u*E.' + h(:,1)*n3.' + h(:,2)*n4.';         % embedded patch
    if o.noise > 0, X = X + o.noise*randn(size(X)); end

    dz = X - z.';
    [~, ord] = sort(sum(dz.^2, 2));
    loc = X(ord(1:k), :);
    mu  = mean(loc, 1).';

    uk  = loc * E;                                  % chart coords of the neighbours
    l2(t) = mean(sum(uk.^2, 2));                    % ell^2 = E||u||^2
    m3(t) = n3.' * mu;
    m4(t) = n4.' * mu;
    j2(t) = sum((E.' * mu).^2);                     % squared tangential jitter
end
end