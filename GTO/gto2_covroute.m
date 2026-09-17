function out = gto2_covroute(mode, varargin)
% GTO2_COVROUTE  The covariance route: how Vd is tilted toward the normal.
%
%   gto2_covroute('gradient')   % curvature gradient, controlled cubic patch
%   gto2_covroute('boundary')   % probe approaching the edge of a sampled region
%   gto2_covroute('all')
%
% TWO MECHANISMS, NOT THREE ROUTES
%   Vd is the top-d eigenspace of the local second-moment matrix. Splitting the
%   local coordinates into tangential u and normal w,
%       M = [ E[u u'] , E[u w'] ;  E[w u'] , E[w w'] ],
%   there are exactly two ways for the leading subspace to rotate toward the
%   normal:
%     VARIANCE   E[w w'] grows until it competes with E[u u'].  Dimension excess
%                and gain-aligned anisotropic noise both act this way.
%     COVARIANCE E[u w'] is nonzero, tilting the subspace even while E[w w']
%                stays small.  Tilt ~ E[u w'] / (E[u u'] - E[w w']).
%
%   On a SYMMETRIC patch of a HOMOGENEOUS manifold the covariance block vanishes
%   by parity: the normal displacement (1/2)*II(u,u) is EVEN in u, so it cannot
%   correlate with u. That is exactly why curvature alone produces a sagitta and
%   no tilt. Every entry below works by breaking that parity.
%
% GRADIENT MODE. Patch  w = (1/2)*kappa*u1^2 + (1/6)*g*u1^3  in one normal
%   direction, uniform disc sampling. The cubic term is ODD in u1, so it
%   correlates with u1 and tilts the fit. With E[u1^4]/E[u1^2] = a^2/2 for a
%   uniform disc of radius a, and ell^2 = a^2/2,
%       theta  =  E[u1 w] / E[u1^2]  =  (g/6) * E[u1^4]/E[u1^2],
%   and for a uniformly filled disc of radius a, E[u1^4]/E[u1^2] = a^2/2 = ell^2,
%   giving theta = g*ell^2/6. THE COEFFICIENT IS A PROPERTY OF THE SAMPLING
%   DISTRIBUTION, not of the geometry: it would differ for a Gaussian
%   neighbourhood. The ell^4 scaling is the robust part.
%       ||Vd'n||^2  ~  (g*ell^2/6)^2   ~  ell^4  ~  k^2.
%   This is the only FRAME term that worsens with neighbourhood size, and it does so
%   QUADRATICALLY (anchor bias also grows, as delta ~ k). On an inhomogeneous
%   manifold the frame cost can outgrow the bias.
%
%   The finite-sample rotation floor adds in expectation (the random tilt has zero
%   mean), and for a uniform disc its mean is known in closed form,
%       E||Vd'n||^2 (g = 0) = ell^2 * kappa^2 / (3k),
%   at leading order in 1/k. At finite k the floor is larger (about 1 + 6/k in Monte
%   Carlo on paraboloid, one-direction patch and sphere alike). Because every cell
%   reuses the same draws, the table subtracts the MEASURED g = 0 rotation at the
%   same k (a paired floor, no added variance) and reports floor_meas/floor_pred as
%   the finite-k factor. ell uses the (k-1) normalisation.
%
% BOUNDARY MODE. A probe near the edge of a sampled region has its neighbourhood
%   truncated on one side, so E[u] ~= 0 and the same parity argument fails. The
%   expected scaling is theta ~ kappa*ell/2, ||Vd'n||^2 ~ ell^2 ~ k, growing as the
%   probe comes within a neighbourhood extent of the edge. Each (k, edge distance)
%   cell uses nProbe probes at random azimuth (the cap is rotationally symmetric)
%   over nSeed reference draws, and reports median, mean, the SE of the median and
%   the interior floor ell^2/(3kR^2) for reference. Sweeping k tests the ell^2 law.
%
% NOT ISOLATED. A density gradient across the neighbourhood belongs to the same
%   class and is NOT separated here. A torus of revolution shows rotation well
%   above the finite-sample floor but with a k-scaling that is neither ell^2 nor
%   ell^4, consistent with a curvature gradient and a density gradient (its
%   density goes as Rc + r*cos(theta)) superposing at different orders. A clean
%   test needs uniform curvature with an imposed density gradient.
%
% CONSEQUENCE FOR LEAKAGE. A tilted basis direction consumes the orthogonal
%   complement in the same way an excess dimension does, so the tilt enters the
%   survival law of gate 1 through the same denominator. Frame tilt and dimension
%   misspecification are the same currency.
%
% RETURNS out.<mode>.rows with columns in out.<mode>.colnames.

if nargin < 1, mode = 'all'; end
out = struct();
switch lower(mode)
    case 'gradient', out.gradient = run_gradient(varargin{:});
    case 'boundary', out.boundary = run_boundary(varargin{:});
    case 'all'
        out.gradient = run_gradient(varargin{:});
        out.boundary = run_boundary(varargin{:});
    otherwise, error('mode must be gradient, boundary or all');
end
end

% =====================================================================
function s = run_gradient(varargin)
ip = inputParser;
ip.addParameter('F', 30);      ip.addParameter('kappa', 0.2);
ip.addParameter('a', 1.0);     ip.addParameter('nRef', 6000);
ip.addParameter('ggrid', [0 0.5 1 2 4]);
ip.addParameter('kgrid', [20 40 80 160 320]);
ip.addParameter('nTrial', 2000); ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;

s.rows = [];
s.colnames = {'block','g','k','ell','rot','pred','ratio','floorPred','ratioSub','floorMeas','floorRatio'};
fprintf('\n=== Covariance route: curvature gradient tilts Vd ===\n');
fprintf('  patch w = kappa*u1^2/2 + g*u1^3/6, uniform disc, F = %d, nRef = %d, %d trials/cell\n', ...
    o.F, o.nRef, o.nTrial);
fprintf(['  Every cell uses the same draws, so the g = 0 rotation at the same k is a PAIRED floor:\n' ...
         '  sub = (rot - floor_meas)/pred. floor_meas/floor_pred measures the finite-k factor of the\n' ...
         '  leading-order floor ell^2 kappa^2/(3k).\n']);

hdr = '  %8s %9s %12s %12s %12s %10s %8s %8s\n';
fprintf('\n  -- g sweep at k = 40 (rotation should scale as g^2) --\n');
fprintf(hdr, 'g', 'ell', '||Vd''n||^2', '(g ell^2/6)^2', 'floor meas', 'meas/pred', 'ratio', 'sub');
[rot0, ~] = cov_cell(0, 40, o);
for g = o.ggrid(:).'
    [rot, ell] = cov_cell(g, 40, o);
    pred = (g*ell^2/6)^2;  flr = ell^2*o.kappa^2/(3*40);
    [rs, ss] = ratio_strings(rot, pred, rot0);
    s.rows(end+1,:) = [1 g 40 ell rot pred rot/max(pred,eps) flr (rot-rot0)/max(pred,eps) rot0 rot0/flr]; %#ok<AGROW>
    fprintf('  %8.2f %9.4f %12.3e %12.3e %12.3e %10.3f %8s %8s\n', g, ell, rot, pred, rot0, rot0/flr, rs, ss);
end

fprintf('\n  -- k sweep at g = 2 (rotation should scale as ell^4 ~ k^2) --\n');
fprintf(hdr, 'k', 'ell', '||Vd''n||^2', '(g ell^2/6)^2', 'floor meas', 'meas/pred', 'ratio', 'sub');
for k = o.kgrid(:).'
    [rot0k, ~]  = cov_cell(0, k, o);
    [rot, ell] = cov_cell(2.0, k, o);
    pred = (2.0*ell^2/6)^2;  flr = ell^2*o.kappa^2/(3*k);
    [rs, ss] = ratio_strings(rot, pred, rot0k);
    s.rows(end+1,:) = [2 2.0 k ell rot pred rot/pred flr (rot-rot0k)/pred rot0k rot0k/flr]; %#ok<AGROW>
    fprintf('  %8d %9.4f %12.3e %12.3e %12.3e %10.3f %8s %8s\n', k, ell, rot, pred, rot0k, rot0k/flr, rs, ss);
end
fprintf('\n');
end

% ---------------------------------------------------------------------
function [rot, ell] = cov_cell(g, k, o)
rng(o.seed);
rots = zeros(o.nTrial,1); ells = zeros(o.nTrial,1);
n = zeros(o.F,1); n(3) = 1;
for t = 1:o.nTrial
    r  = o.a*sqrt(rand(o.nRef,1));
    th = 2*pi*rand(o.nRef,1);
    u  = [r.*cos(th), r.*sin(th)];
    X = zeros(o.nRef, o.F);
    X(:,1:2) = u;
    X(:,3) = 0.5*o.kappa*u(:,1).^2 + (g/6)*u(:,1).^3;
    [~, ord] = sort(sum(X.^2, 2));          % probe at the origin of the chart
    loc = X(ord(1:k), :);
    Xc  = loc - mean(loc,1);
    [~,~,V] = svd(Xc, 'econ');
    rots(t) = sum((V(:,1:2).' * n).^2);
    ells(t) = sum(sum(Xc(:,1:2).^2, 2)) / (k - 1);   % tangential, unbiased
end
rot = mean(rots);  ell = sqrt(mean(ells));
end

% ---------------------------------------------------------------------
function [rs, ss] = ratio_strings(rot, pred, flr)
if pred <= 0
    rs = '---';  ss = '---';
else
    rs = sprintf('%.3f', rot/pred);  ss = sprintf('%.3f', (rot - flr)/pred);
end
end

% =====================================================================
function s = run_boundary(varargin)
ip = inputParser;
ip.addParameter('F', 30);    ip.addParameter('R', 5.0);
ip.addParameter('cap', 0.35);
ip.addParameter('nRef', 8000);
ip.addParameter('kgrid', 40);
ip.addParameter('fracs', [0 0.3 0.6 0.8 0.9 0.95 0.99]);
ip.addParameter('nProbe', 200);
ip.addParameter('nSeed', 3);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;

s.rows = []; s.colnames = {'frac','edgeDist','edgeOverEll','rot','ell','k','rotMean','rotSE','floorPred'};
fprintf('\n=== Covariance route: a boundary tilts Vd ===\n');
fprintf(['  spherical cap, half-angle %.2f, R = %.1f, nRef = %d; noiseless, d = d_true;\n' ...
         '  %d probes at random azimuth x %d reference draws per cell\n'], ...
    o.cap, o.R, o.nRef, o.nProbe, o.nSeed);
for k = o.kgrid(:).'
    fprintf('\n  -- k = %d --\n', k);
    fprintf('  %8s %10s %9s %12s %12s %10s %12s\n', 'ang/cap', 'edge dist', 'edge/ell', ...
        'median', 'mean', 'SE(med)', 'interior');
    for f = o.fracs(:).'
        aa  = o.cap*f;
        rot = zeros(o.nProbe*o.nSeed,1);  ell = rot;  c = 0;
        for sd = 1:o.nSeed
            rng(o.seed + sd - 1);
            a  = o.cap*sqrt(rand(o.nRef,1));
            ph = 2*pi*rand(o.nRef,1);
            u  = [sin(a).*cos(ph), sin(a).*sin(ph), cos(a)];
            A  = [o.R*u, zeros(o.nRef, o.F-3)];
            for ip_ = 1:o.nProbe
                az = 2*pi*rand;
                z  = zeros(1,o.F);  z(1:3) = o.R*[sin(aa)*cos(az), sin(aa)*sin(az), cos(aa)];
                nrm = (z/norm(z)).';
                [~, ord] = sort(sum((A - z).^2, 2));
                loc = A(ord(1:k), :);
                Xc  = loc - mean(loc,1);
                [~,~,V] = svd(Xc, 'econ');
                c = c + 1;
                rot(c) = sum((V(:,1:2).' * nrm).^2);
                ell(c) = sqrt(sum(sum((Xc*V(:,1:2)).^2, 2)) / (k - 1));
            end
        end
        edge = o.R*(o.cap - aa);
        el   = median(ell);
        mr   = median(rot);
        se   = 1.2533*std(rot)/sqrt(numel(rot));
        flr  = el^2/(3*k*o.R^2);
        s.rows(end+1,:) = [f edge edge/el mr el k mean(rot) se flr]; %#ok<AGROW>
        fprintf('  %8.2f %10.4f %9.2f %12.3e %12.3e %10.1e %12.3e\n', f, edge, edge/el, mr, mean(rot), se, flr);
    end
end
fprintf(['\n  Interior cells should sit near the floor (mean). Rotation rises as the probe comes\n' ...
         '  within a neighbourhood extent of the edge; compare across k at matched edge/ell.\n\n']);
end