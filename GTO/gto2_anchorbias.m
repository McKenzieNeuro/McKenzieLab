function out = gto2_anchorbias(varargin)
% GTO2_ANCHORBIAS  Where the anchor's own offset goes.
%
% WHY THE ABSOLUTE NOVEL ENERGY, NOT THE FRACTION
%   It is tempting to measure the off-manifold FRACTION O produced by a pure
%   on-manifold move and read it as a false-positive rate for novelty. That is a
%   mistake, and this file is arranged to make the mistake visible. Reporting the
%   ABSOLUTE novel energy instead shows ||r_novel|| flat at sqrt(F-3)*sigma across
%   every radius and every neighbourhood size -- it is the ambient observation
%   noise and nothing else. There is no curvature-driven novelty false positive on
%   a sphere to tabulate, and any R- or k-dependence seen in the FRACTION is the
%   denominator moving, not the numerator.
%
%   What curvature does produce is REAL, and it is radial. The anchor mu is a kNN
%   centroid inside a convex manifold, so it sits a sagitta
%       delta ~ ell^2 / (2R)
%   below the surface, where ell is the tangential neighbourhood extent D.ell
%   (unbiased, (k-1) normalisation; called rho below only as a local variable name).
%   With ell measured tangentially the exact sphere geometry adds a known next-order
%   term, delta = (ell^2/2R)(1 + ell^2/3R^2); 'ratio2' uses it, with the Monte Carlo
%   standard error of the median. On an origin-centred sphere that
%   offset is exactly collinear with the radial gain axis, so the construction
%   absorbs it into G, leaving O uncontaminated. That is a positive result for the
%   nested construction and a stated limit on G, not a failure mode of O. This
%   sweep therefore measures the anchor offset directly and checks it against
%   ell^2/(2R), and reports the novel energy in absolute units so the two cannot
%   be confused.
%
%   The result is specific to |c| = 1 (H parallel to the gain axis), which holds
%   identically for any 2-manifold whose span is 3-dimensional. See gto2_manifolds
%   for the |c| < 1 case, where the offset does reach O.
%
% CONSTRUCTION
%   Every probe is a PURE on-manifold move: a geodesic step of fixed arc length
%   along a sphere of radius R, so the imposed change has no radial and no novel
%   component. Any radial energy in r = z_B - mu is the anchor's, not the probe's.
%
% OPTIONS
%   'F'      30
%   'Rgrid'  [3 5 8 12 20]     sphere radii (curvature 1/R)
%   'kgrid'  [20 40 80 160]    neighbourhood sizes
%   'nRef'   6000
%   'nTest'  300
%   'd'      2
%   'noise'  0.01
%   'scale'  0.3               on-manifold arc length
%   'seed'   0
%
% RETURNS out.rows (columns in out.colnames) and prints a table.

ip = inputParser;
ip.addParameter('F', 30);
ip.addParameter('Rgrid', [3 5 8 12 20]);
ip.addParameter('kgrid', [20 40 80 160]);
ip.addParameter('nRef', 6000);
ip.addParameter('nTest', 300);
ip.addParameter('d', 2);
ip.addParameter('noise', 0.01);
ip.addParameter('scale', 0.3);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;
rng(o.seed);

% Below this magnitude a predicted value is display-zero, so a ratio against
% it is not a measurement -- print '---' instead (mirrors gto2_curvlaw).
PRED_FLOOR = 1e-9;

[Qemb,~] = qr(randn(o.F));
ambFloor = sqrt(o.F - 3) * o.noise;    % expected ||r_novel|| from noise alone

out.rows = [];
fprintf('\n=== Anchor offset and where it lands (pure on-manifold probes) ===\n');
fprintf('  scale = %.2f, F = %d, d = %d; expected ||r_novel|| from noise = %.4f\n', ...
    o.scale, o.F, o.d, ambFloor);
fprintf('\n  %6s %6s %8s %10s %10s %10s %7s %7s %6s %10s %7s\n', ...
    'R', 'k', 'ell', '|r|', '|r_rad|', 'ell^2/2R', 'ratio', 'ratio2', 'SE', '|r_novel|', 'ratio');
for R = o.Rgrid(:).'
    A = sample_sphere(o.nRef, R, o.F, o.noise) * Qemb.';
    for k = o.kgrid(:).'
        Z = geodesic_probes(Qemb, R, o);
        D = gto2_core(A, Z, 'k', k, 'd', o.d);
        keep = ~D.degenerate;

        rho    = median(D.ell(keep));
        radAbs = median(sqrt(D.G(keep)) .* D.rnorm(keep));   % absolute radial energy
        novAbs = median(sqrt(D.O(keep)) .* D.rnorm(keep));   % absolute novel energy
        rn     = median(D.rnorm(keep));
        sag    = rho^2 / (2*R);
        sag2   = sag * (1 + rho^2/(3*R^2));
        radVec = sqrt(D.G(keep)) .* D.rnorm(keep);
        radSE  = 1.2533 * std(radVec) / sqrt(numel(radVec));  % SE of a median

        if abs(sag) < PRED_FLOOR
            rad_rat_str = '   ---';
        else
            rad_rat_str = sprintf('%6.3f', radAbs / sag);
        end
        if abs(ambFloor) < PRED_FLOOR
            nov_rat_str = '   ---';
        else
            nov_rat_str = sprintf('%6.3f', novAbs / ambFloor);
        end

        out.rows(end+1,:) = [R 1/R k rho rn radAbs sag novAbs median(D.O(keep)) ...
                             sag2 radAbs/sag2 radSE/sag2]; %#ok<AGROW>
        fprintf('  %6.1f %6d %8.3f %10.4f %10.4f %10.4f %7s %7.3f %6.3f %10.4f %7s\n', ...
            R, k, rho, rn, radAbs, sag, rad_rat_str, radAbs/sag2, radSE/sag2, novAbs, nov_rat_str);
    end
end
out.colnames = {'R','kappa','k','rho','rnorm','radAbs','sagittaPred','novAbs','medO', ...
                'sagittaPred2','ratio2','ratio2SE'};

ratio = out.rows(:,11);
zdev  = (ratio - 1) ./ out.rows(:,12);
fprintf(['\n  ratio2 = |r_rad| / [(ell^2/2R)(1+ell^2/3R^2)]: mean %.3f, range %.3f-%.3f;\n' ...
         '  largest |deviation|/SE = %.1f. The anchor offset is radial, so the gain axis absorbs it.\n'], ...
         mean(ratio), min(ratio), max(ratio), max(abs(zdev)));
fprintf(['  ||r_novel||: mean %.4f, range %.4f-%.4f against a noise floor of %.4f -- flat in\n' ...
         '  BOTH R and k. Curvature contributes no off-manifold energy; the O fraction\n' ...
         '  moves only because ||r|| does.\n\n'], ...
         mean(out.rows(:,8)), min(out.rows(:,8)), max(out.rows(:,8)), ambFloor);
end

% =====================================================================
function Z = geodesic_probes(Qemb, R, o)
% Pure on-manifold probes: a geodesic step of arc length o.scale, so the probe
% stays ON the sphere and the imposed (G,T,O) is (0,1,0) exactly. A straight-line
% tangential step would NOT do this -- it leaves the sphere by scale^2/(2R),
% which is genuine radial displacement and would be correctly reported as gain.
n = o.nTest;
u = randn(n,3); u = u ./ vecnorm(u,2,2);
Z = zeros(n, o.F);
ang = o.scale / R;
for i = 1:n
    ui  = u(i,:).';
    t   = null(ui.');
    dir = t * randn(2,1);  dir = dir / norm(dir);
    m3  = R * (cos(ang)*ui + sin(ang)*dir);
    Z(i,:) = (Qemb * [m3; zeros(o.F-3,1)]).';
end
if o.noise > 0, Z = Z + o.noise*randn(size(Z)); end
end

% ---------------------------------------------------------------------
function X = sample_sphere(n, R, F, noise)
v = randn(n,3);  v = v ./ vecnorm(v,2,2);
X = [R*v, zeros(n, F-3)];
if noise > 0, X = X + noise*randn(size(X)); end
end