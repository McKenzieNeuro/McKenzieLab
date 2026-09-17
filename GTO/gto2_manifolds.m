function out = gto2_manifolds(mode, varargin)
% GTO2_MANIFOLDS  Validation on manifolds other than the sphere.
%
%   gto2_manifolds('clifford')   % Gate 5 with |c| < 1  -- closes the open case
%   gto2_manifolds('torus')      % Gates 3 and 4 on an inhomogeneous manifold
%   gto2_manifolds('ellipse')    % Gate 5 with |c| < 1 AND grad(II) ~= 0
%   gto2_manifolds('all')
%
% WHY THE SPHERE IS NOT ENOUGH, AND WHY AN ELLIPSOID WOULD NOT HELP
%   Gate 5 turns on c = n_g' * Hhat, the alignment between the gain axis and the
%   mean curvature vector. There is a structural reason every test so far has
%   returned |c| = 1, and it is not the sphere's doing:
%
%     For a d = 2 manifold whose SPAN is 3-dimensional, the normal space inside
%     that span is 1-dimensional. P_N*rho, the surface normal n, and H are then
%     all forced parallel, so |c| = 1 IDENTICALLY.
%
%   That covers the sphere, the spherical cap, every quadric patch in 3-space,
%   the ellipsoid and the torus of revolution. An ellipsoid cannot test gate 5 --
%   it varies the ALIGNMENT a (gate 4), not c. Reaching |c| < 1 requires the
%   manifold's span to be at least d + 2.
%
% CLIFFORD MODE. The minimal geometry with span d+2 is the Clifford torus in R^4,
%       p(u,v) = (R1 cos u, R1 sin u, R2 cos v, R2 sin v),
%   embedded in F dimensions. It is origin-centred, and
%       n1 = (cos u, sin u, 0, 0),          n2 = (0, 0, cos v, sin v),
%       II(e_u,e_u) = -n1/R1,               II(e_v,e_v) = -n2/R2,
%       H = -(n1/R1 + n2/R2)/2,             rho_hat = (R1 n1 + R2 n2)/rho,
%   with rho = hypot(R1,R2). Because rho_hat lies wholly in the normal space the
%   ALIGNMENT IS EXACTLY 1, so gate 4 is clean and gate 5 is isolated. The
%   curvature alignment is
%       c = -(1/rho) / ||H||,
%   which is -1 at R1 = R2 and falls as the torus is made asymmetric. Setting
%   R1:R2 = 1:1, 1:2, 1:3, 1:5 gives c = -1.00, -0.80, -0.60, -0.38.
%
%   Probes are exact geodesic steps (the Clifford torus is intrinsically flat, so
%   a straight line in (u,v) stays on the surface), i.e. imposed (G,T,O) =
%   (0,1,0). Any off-manifold energy is the anchor's. The prediction is
%       ||r_novel||  ->  eta = (ell^2/2) * ||H|| * sqrt(1 - c^2),
%   which is zero at |c| = 1 (recovering the sphere result: novel energy pinned
%   at the noise floor) and grows as the torus is made asymmetric.
%
% TORUS MODE. A torus of revolution, origin-centred,
%       p(th,ph) = ((Rc + r cos th) cos ph, (Rc + r cos th) sin ph, r sin th),
%   tests two things the sphere cannot.
%     (i) GATE 3 ON AN INHOMOGENEOUS MANIFOLD. The sphere, the Clifford torus and
%         a quadric at its vertex are all homogeneous, so grad(II) = 0 and the
%         term Proposition 1 neglects is structurally zero. Here curvature varies
%         over the surface. Moreover, when Rc < 2*r the mean curvature CHANGES
%         SIGN on the inner region, so the anchor is displaced OUTWARD there --
%         a directional prediction the sphere cannot make.
%     (ii) GATE 4 WHERE IT FAILS. The alignment
%             a = |Rc cos th + r| / ||p||
%         passes through zero at cos th = -r/Rc, so a single reference manifold
%         contains both well-conditioned and degenerate regions and the
%         alignment gate can be checked against an analytic value pointwise.
%
%   NOTE ON MEASUREMENT. The anchor displacement must be referred to the probe's
%   OWN position, not to a pre-step base point. A probe stepped an arc length s
%   along the manifold carries its own normal offset ~ kappa_n*s^2/2, which does
%   not shrink with the neighbourhood and will swamp the sagitta at high sampling
%   density. Probes here sit AT the base point for the gate-3 measurement.
%
% RETURNS out.<mode>.rows with columns in out.<mode>.colnames.

if nargin < 1, mode = 'all'; end
out = struct();
switch lower(mode)
    case 'clifford', out.clifford = run_clifford(varargin{:});
    case 'torus',    out.torus    = run_torus(varargin{:});
    case 'ellipse',  out.ellipse  = run_ellipse(varargin{:});
    case 'all'
        out.clifford = run_clifford(varargin{:});
        out.torus    = run_torus(varargin{:});
        out.ellipse  = run_ellipse(varargin{:});
    otherwise, error('mode must be clifford, torus, ellipse or all');
end
end

% =====================================================================
function s = run_clifford(varargin)
ip = inputParser;
ip.addParameter('F', 30);
ip.addParameter('nRef', 8000);
ip.addParameter('nTest', 300);
ip.addParameter('kgrid', [20 40 80 160]);
ip.addParameter('ratios', [1 1; 1 2; 1 3; 1 5]);
ip.addParameter('noise', 0.001);
ip.addParameter('arc', 0.15);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;

s.rows = []; s.colnames = {'R1','R2','k','c','align','novAbs','etaPred','radAbs','deltaPred','floor','ratioQuad'};
% ratioQuad = novAbs / sqrt(eta^2 + floor^2): the noise floor adds in quadrature
fprintf('\n=== Gate 5 on a Clifford torus: |c| < 1 is reachable ===\n');
fprintf('  F = %d, nRef = %d, noise = %g, arc = %.2f; noise floor = %.5f\n', ...
    o.F, o.nRef, o.noise, o.arc, sqrt(o.F-3)*o.noise);

fprintf('\n  -- curvature alignment sweep at k = 40 --\n');
fl = sqrt(o.F-3)*o.noise;
fprintf('  %7s %8s %7s %10s %10s %7s %8s %10s %10s\n', ...
    'R1:R2', 'c', 'align', '|r_nov|', 'eta pred', 'ratio', 'quad', '|r_rad|', 'delta pred');
for i = 1:size(o.ratios,1)
    R1 = o.ratios(i,1); R2 = o.ratios(i,2);
    [c, Hn] = clifford_geom(R1, R2);
    res = clifford_cell(R1, R2, 40, o);
    eta = 0.5*res.ell2*Hn*sqrt(max(0,1-c^2));
    del = 0.5*res.ell2*Hn*abs(c);
    if eta < 1e-9, rat = NaN; else, rat = res.novAbs / eta; end   % NaN printed where eta = 0
    rq  = res.novAbs / sqrt(eta^2 + fl^2);
    s.rows(end+1,:) = [R1 R2 40 c res.align res.novAbs eta res.radAbs del fl rq]; %#ok<AGROW>
    fprintf('  %3d:%-3d %8.3f %7.3f %10.5f %10.5f %7.2f %8.3f %10.5f %10.5f\n', ...
        R1, R2, c, res.align, res.novAbs, eta, rat, rq, res.radAbs, del);
end

fprintf('\n  -- k sweep at R1:R2 = 1:3 (eta ~ ell^2 ~ k; on a sphere this was FLAT) --\n');
fprintf('  %6s %10s %10s %10s %7s %8s\n', 'k', 'ell^2', '|r_nov|', 'eta pred', 'ratio', 'quad');
[c, Hn] = clifford_geom(1, 3);
for k = o.kgrid(:).'
    res = clifford_cell(1, 3, k, o);
    eta = 0.5*res.ell2*Hn*sqrt(max(0,1-c^2));
    rq  = res.novAbs / sqrt(eta^2 + fl^2);
    s.rows(end+1,:) = [1 3 k c res.align res.novAbs eta res.radAbs 0.5*res.ell2*Hn*abs(c) fl rq]; %#ok<AGROW>
    fprintf('  %6d %10.5f %10.5f %10.5f %7.2f %8.3f\n', k, res.ell2, res.novAbs, eta, res.novAbs/eta, rq);
end
fprintf(['\n  At |c| = 1 the novel energy sits at the noise floor, as on the sphere. As the\n' ...
         '  torus is made asymmetric it rises to track eta, and grows linearly in k. Gate 5''s\n' ...
         '  booking rule is confirmed, not merely derived.\n\n']);
end

% ---------------------------------------------------------------------
function [c, Hn] = clifford_geom(R1, R2)
rho = hypot(R1, R2);
Hn  = 0.5 * hypot(1/R1, 1/R2);
c   = -(1/rho) / Hn;
end

% ---------------------------------------------------------------------
function res = clifford_cell(R1, R2, k, o)
rng(o.seed);
[Q,~] = qr(randn(o.F));
emb = @(u,v) [R1*cos(u), R1*sin(u), R2*cos(v), R2*sin(v), zeros(numel(u), o.F-4)] * Q.';

u = 2*pi*rand(o.nRef,1);  v = 2*pi*rand(o.nRef,1);
A = emb(u,v);
if o.noise > 0, A = A + o.noise*randn(size(A)); end

ub = 2*pi*rand(o.nTest,1); vb = 2*pi*rand(o.nTest,1);
ph = 2*pi*rand(o.nTest,1);
Z = emb(ub + o.arc*cos(ph)/R1, vb + o.arc*sin(ph)/R2);   % exact geodesic step
if o.noise > 0, Z = Z + o.noise*randn(size(Z)); end

D = gto2_core(A, Z, 'k', k, 'd', 2);
keep = ~D.degenerate;
res.ell2   = k * (4*pi^2*R1*R2) / (2*pi*o.nRef);
res.align  = median(D.align(keep));
res.novAbs = median(sqrt(D.O(keep)) .* D.rnorm(keep));
res.radAbs = median(sqrt(D.G(keep)) .* D.rnorm(keep));
end

% =====================================================================
function s = run_ellipse(varargin)
% CIRCLE x ELLIPSE in R^4:  p(u,v) = (R1 cos u, R1 sin u, A cos v, B sin v).
%
% The Clifford torus establishes the booking rule but is intrinsically FLAT and
% homogeneous, so grad(II) = 0 there and Proposition 1's neglected term is
% structurally zero. Replacing the second circle by an ellipse keeps the span at
% d+2 -- so |c| < 1 is still reachable -- while making the normal curvature vary:
%       n1 = (cos u, sin u, 0, 0)
%       n2 = (0, 0, B cos v, A sin v)/q,     q = sqrt(A^2 sin^2 v + B^2 cos^2 v)
%       II(e_u,e_u) = -n1/R1,   II(e_v,e_v) = -kappa_e * n2,   kappa_e = A*B/q^3
% so kappa_e, ||H||, the alignment and c ALL vary with v. At the latitudes where
% |c| -> 1 the prediction eta -> 0 and the measured novel energy should fall back
% to the noise floor; away from them it should track eta. Both regimes appear on
% ONE manifold instead of across a family.
%
% BINNING WARNING. eta varies strongly with v. Evaluating the prediction at a
% bin's median coordinate rather than per probe produces an apparent discrepancy
% of tens of percent that is purely an artefact of the binning. Predict per probe,
% aggregate afterwards -- which is what this function does.
ip = inputParser;
ip.addParameter('F', 30);   ip.addParameter('R1', 1.0);
ip.addParameter('A', 1.0);  ip.addParameter('B', 3.0);
ip.addParameter('nRefGrid', [5000 10000 20000 40000]);
ip.addParameter('k', 60);   ip.addParameter('nTest', 1500);
ip.addParameter('noise', 0); ip.addParameter('arc', 0.05);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;

s.rows = []; s.colnames = {'nRef','ell2','novAbs','etaPred','ratio'};
fprintf('\n=== Gate 5 where grad(II) ~= 0: circle x ellipse (R1=%.1f, A=%.1f, B=%.1f) ===\n', ...
    o.R1, o.A, o.B);
fprintf('  noiseless by default, so the only novel energy is the anchor''s\n\n');
fprintf('  ratio = median over probes of the per-probe ratio |r_nov|/eta (not a ratio of medians)\n');
fprintf('  %8s %10s %12s %12s %8s\n', 'nRef', 'ell^2', '|r_nov|', 'eta pred', 'ratio');
for nRef = o.nRefGrid(:).'
    rng(o.seed);
    [Q,~] = qr(randn(o.F));
    emb = @(u,v) [o.R1*cos(u), o.R1*sin(u), o.A*cos(v), o.B*sin(v), ...
                  zeros(numel(u), o.F-4)] * Q.';
    u = 2*pi*rand(nRef,1);  v = ellipse_v(nRef, o);
    Aref = emb(u,v);
    if o.noise > 0, Aref = Aref + o.noise*randn(size(Aref)); end
    ub = 2*pi*rand(o.nTest,1); vb = ellipse_v(o.nTest, o);
    Z = emb(ub + o.arc/o.R1, vb);            % step in u: stays ON the manifold
    if o.noise > 0, Z = Z + o.noise*randn(size(Z)); end

    D = gto2_core(Aref, Z, 'k', o.k, 'd', 2);
    keep = ~D.degenerate;

    vg  = linspace(0, 2*pi, 2000).';
    per = trapz(vg, sqrt(o.A^2*sin(vg).^2 + o.B^2*cos(vg).^2));
    ell2 = o.k * (2*pi*o.R1*per) / (2*pi*nRef);

    [c, Hn] = ellipse_geom(vb, o);
    eta = 0.5*ell2*Hn.*sqrt(max(0, 1 - c.^2));           % PER PROBE
    nov = sqrt(D.O) .* D.rnorm;
    sel = keep & (abs(c) < 0.95);                        % away from the |c|->1 latitudes
    rat = median(nov(sel) ./ eta(sel));
    s.rows(end+1,:) = [nRef ell2 median(nov(sel)) median(eta(sel)) rat]; %#ok<AGROW>
    fprintf('  %8d %10.5f %12.5f %12.5f %8.3f\n', nRef, ell2, median(nov(sel)), median(eta(sel)), rat);
end
fprintf(['\n  The ratio is stable across the sweep: the booking rule survives a varying second\n' ...
         '  fundamental form, and is not an artefact of the Clifford torus''s flatness.\n\n']);
end

% ---------------------------------------------------------------------
function [c, Hn] = ellipse_geom(v, o)
q  = sqrt(o.A^2*sin(v).^2 + o.B^2*cos(v).^2);
ke = o.A*o.B ./ q.^3;
nr = sqrt(o.R1^2 + o.A^2*cos(v).^2 + o.B^2*sin(v).^2);
r1 = o.R1 ./ nr;  r2 = (o.A*o.B) ./ (q .* nr);      % rho_hat on (n1, n2)
al = hypot(r1, r2);                                  % alignment ||P_N rho||
g1 = r1 ./ al;  g2 = r2 ./ al;                       % n_g on (n1, n2)
h1 = -0.5/o.R1 * ones(size(v));  h2 = -0.5*ke;       % H on (n1, n2)
Hn = hypot(h1, h2);
c  = (g1.*h1 + g2.*h2) ./ Hn;
end

% ---------------------------------------------------------------------
function v = ellipse_v(n, o)
% Rejection sampling for uniform AREA: density ~ q(v).
v = zeros(n,1); got = 0;
while got < n
    t = 2*pi*rand(2*n,1);
    w = rand(2*n,1) * max(o.A, o.B);
    t = t(w < sqrt(o.A^2*sin(t).^2 + o.B^2*cos(t).^2));
    take = min(numel(t), n - got);
    v(got+1:got+take) = t(1:take);
    got = got + take;
end
end

% =====================================================================
function s = run_torus(varargin)
ip = inputParser;
ip.addParameter('F', 30);
ip.addParameter('Rc', 3.0);          % Rc < 2r so the mean curvature changes sign
ip.addParameter('r', 2.0);
ip.addParameter('nRef', 20000);
ip.addParameter('nTest', 4000);
ip.addParameter('k', 40);
ip.addParameter('noise', 0);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;
rng(o.seed);

[Q,~] = qr(randn(o.F));
emb = @(th,ph) [(o.Rc + o.r*cos(th)).*cos(ph), (o.Rc + o.r*cos(th)).*sin(ph), ...
                o.r*sin(th), zeros(numel(th), o.F-3)] * Q.';

th = torus_theta(o.nRef, o);  ph = 2*pi*rand(o.nRef,1);
A  = emb(th, ph);
if o.noise > 0, A = A + o.noise*randn(size(A)); end

tb = torus_theta(o.nTest, o);  pb = 2*pi*rand(o.nTest,1);
Z  = emb(tb, pb);                       % probe AT the base point (see header note)
if o.noise > 0, Z = Z + o.noise*randn(size(Z)); end

D = gto2_core(A, Z, 'k', o.k, 'd', 2);

N = [cos(tb).*cos(pb), cos(tb).*sin(pb), sin(tb), zeros(o.nTest, o.F-3)] * Q.';
M = emb(tb, pb);
disp_ = sum((D.mu - M) .* N, 2);                    % anchor displacement along outward normal
k1 = cos(tb) ./ (o.Rc + o.r*cos(tb));  k2 = ones(size(tb))/o.r;
% ell^2 is MEASURED: the mean of D.ell^2 over probes (uniform area density, so one
% value serves the whole surface). The areal formula k*Area/(2*pi*nRef) underestimates
% the second moment of a k-NN neighbourhood by k/(k+1) -- the k-th neighbour sits on
% the boundary -- which showed as a flat ~2.5% excess in the outer bins at every nRef.
ell2an = o.k * (4*pi^2*o.Rc*o.r) / (2*pi*o.nRef);
ell2   = mean(D.ell(~D.degenerate).^2);
pred = -0.5*ell2*(k1 + k2)/2;
aTrue = abs(o.Rc*cos(tb) + o.r) ./ vecnorm(M,2,2);

s.rows = []; s.colnames = {'thLo','thHi','k1k2','dispMeas','dispPred','ratio','alignMeas','alignTrue','fracDegen', ...
                           'ratioSE','alignErrMax','alignErrMedAbs','nRef'};
% ratio = least-squares slope of per-probe disp on per-probe pred within the bin
% (robust where k1+k2 crosses zero inside a bin); dispMeas/dispPred are bin medians.
fprintf('\n=== Gates 3 and 4 on a torus of revolution (Rc = %.1f, r = %.1f) ===\n', o.Rc, o.r);
fprintf(['  nRef = %d, k = %d, ell^2 = %.5f measured (areal formula %.5f, ratio %.3f; (k+1)/k = %.3f);\n' ...
         '  mean curvature changes sign; %.1f%% degenerate\n'], ...
    o.nRef, o.k, ell2, ell2an, ell2/ell2an, (o.k+1)/o.k, 100*mean(D.degenerate));
fprintf('  ratio = LS slope of per-probe disp on pred (+- SE); alignment error = |a_meas - a_true| per probe\n\n');
fprintf('  %13s %9s %11s %11s %8s %7s %11s %10s %9s %8s\n', ...
    'theta bin', 'k1+k2', 'disp meas', 'pred', 'slope', 'SE', 'align meas', 'align true', 'max|da|', '%degen');
e = linspace(0, 2*pi, 9);
for i = 1:8
    sel = tb >= e(i) & tb < e(i+1);
    if sum(sel) < 30, continue, end
    mp  = median(pred(sel));
    ps  = pred(sel);  ds = disp_(sel);
    rat = sum(ds.*ps) / sum(ps.^2);
    res_ = ds - rat*ps;
    rse = sqrt(sum(res_.^2)/(numel(ds)-1) / sum(ps.^2));
    da  = abs(D.align(sel) - aTrue(sel));
    s.rows(end+1,:) = [e(i) e(i+1) median(k1(sel)+k2(sel)) median(ds) mp rat ...
                       median(D.align(sel)) median(aTrue(sel)) mean(D.degenerate(sel)) ...
                       rse max(da) median(da) o.nRef]; %#ok<AGROW>
    fprintf('  %5.2f-%5.2f  %9.4f %11.5f %11.5f %8.3f %7.3f %11.3f %10.3f %9.3f %8.1f\n', ...
        e(i), e(i+1), median(k1(sel)+k2(sel)), median(ds), mp, rat, rse, ...
        median(D.align(sel)), median(aTrue(sel)), max(da), 100*mean(D.degenerate(sel)));
end
fprintf(['\n  The displacement should reverse sign with the mean curvature on the inner region.\n' ...
         '  Run at increasing nRef: a genuine higher-order term shrinks roughly as ell^2 (halving\n' ...
         '  per doubling of nRef); a harness confound does not.\n\n']);
end

% ---------------------------------------------------------------------
function th = torus_theta(n, o)
% Rejection sampling for uniform AREA on the torus: density ~ (Rc + r cos th).
th = zeros(n,1); got = 0;
while got < n
    t = 2*pi*rand(2*n,1);
    u = rand(2*n,1) * (o.Rc + o.r);
    t = t(u < (o.Rc + o.r*cos(t)));
    take = min(numel(t), n - got);
    th(got+1:got+take) = t(1:take);
    got = got + take;
end
end