function out = gto2_checks(mode)
% GTO2_CHECKS  Invariant checks for the nested decomposition in gto2_core.
%
%   gto2_checks           % runs all
%   gto2_checks('c1')     % order-independence of the gain projection
%   gto2_checks('c2')     % loud failure as alignment -> 0
%   gto2_checks('c3')     % G+T+O = 1 pointwise
%   gto2_checks('c4')     % pure-mechanism attribution on a sphere
%   gto2_checks('c5')     % Proposition 1, exact identity in the anchor frame (Table S1)
%
% EXACT checks are algebraic identities: a FAIL is a bug, not a tolerance
% question. C2 is a behavioural check -- it asserts that the construction
% refuses to answer where the gain axis is undefined, rather than returning a
% number. That distinguishes this construction from the previous radial-first
% one, which returned a gain estimate regardless of alignment.

if nargin < 1, mode = 'all'; end
out = struct();

switch lower(mode)
    case 'c1',  out.c1 = check_order();
    case 'c2',  out.c2 = check_degenerate();
    case 'c3',  out.c3 = check_sum();
    case 'c4',  out.c4 = check_attribution();
    case 'c5',  out.c5 = check_prop1();
    case 'all'
        out.c1 = check_order();
        out.c2 = check_degenerate();
        out.c3 = check_sum();
        out.c4 = check_attribution();
        out.c5 = check_prop1();
        report(out);
    otherwise
        error('mode must be c1, c2, c3, c4, c5 or all');
end
end

% =====================================================================
function s = check_order()
% C1 (EXACT). ng is orthogonal to Vd by construction, so projecting the RAW
% displacement onto ng must equal projecting the tangent-removed displacement
% onto it. If these differ, the gain axis is not being built inside the normal
% space and the decomposition has become order-dependent.
rng(0);
F = 30; d = 2;
A = sphere_cap(3000, F, 5.0, 0.6, 0.01);
Z = sphere_cap(200,  F, 5.0, 0.6, 0.01) * 1.4;      % gain-modulated probes

D = gto2_core(A, Z, 'k', 40, 'd', d);

% recompute both ways directly
worst = 0;
for i = 1:size(Z,1)
    mu = D.mu(i,:).';  ng = D.ng(i,:).';
    if any(isnan(ng)), continue, end
    r = Z(i,:).' - mu;
    % need Vd again to form PNr
    [~, ord] = sort(pdist2(Z(i,:), A), 2);
    loc = A(ord(1:40), :);  muL = mean(loc,1).';
    [~,~,V] = svd(loc - muL.', 'econ');  Vd = V(:,1:d);
    PNr = r - Vd*(Vd.'*r);
    worst = max(worst, abs(ng.'*r - ng.'*PNr));
end
s.name = 'C1 gain projection is order-independent';
s.kind = 'EXACT';
s.pass = worst < 1e-10;
s.detail = sprintf('max |ng''r - ng''(P_N r)| = %.2e (must be ~0 by construction)', worst);
end

% =====================================================================
function s = check_degenerate()
% C2 (BEHAVIOURAL). On a manifold whose tangent space contains the radial
% direction, the gain axis is undefined. The construction must REFUSE (NaN),
% not return a number.
%
% A patch of a cone through the origin has exactly this property: the position
% vector of any point is itself tangent to the cone (Euler's relation for a
% degree-1 homogeneous surface), so P_N*rho -> 0.
rng(1);
F = 20;
A = cone_patch(4000, F, 0.01);
Z = cone_patch(150,  F, 0.01) * 1.3;                % pure radial scaling

D = gto2_core(A, Z, 'k', 40, 'd', 2);

fracDegen = mean(D.degenerate);
medAlign  = median(D.align);
gAllNaN   = all(isnan(D.G(D.degenerate)));
tStillOk  = all(~isnan(D.T));

s.name = 'C2 gain axis fails loudly where alignment vanishes';
s.kind = 'BEHAVIOURAL';
s.pass = medAlign < 0.2 && fracDegen > 0.5 && gAllNaN && tStillOk;
s.detail = sprintf(['median alignment = %.4f on a cone (NOT ~0: sampling and noise ' ...
    'keep it off zero even where the gain axis is undefined in theory); %.0f%% of ' ...
    'points flagged degenerate at the default tolerance; G is NaN wherever flagged ' ...
    '(%d); T still defined everywhere (%d)'], ...
    medAlign, 100*fracDegen, gAllNaN, tStillOk);
end

% =====================================================================
function s = check_sum()
% C3 (EXACT). The three subspaces are mutually orthogonal, so the energy
% fractions must sum to one pointwise.
rng(2);
F = 30;
A = sphere_cap(3000, F, 5.0, 0.6, 0.01);
Zg = sphere_cap(100, F, 5.0, 0.6, 0.01) * 1.4;                       % gain
Zo = sphere_cap(100, F, 5.0, 0.6, 0.01);                             % + novel
[Q,~] = qr(randn(F)); Zo = Zo + 2.0 * repmat(Q(:,end).', 100, 1);
Z = [Zg; Zo];

D = gto2_core(A, Z, 'k', 40, 'd', 2);
ok = ~D.degenerate;
worst = max(abs(D.G(ok) + D.T(ok) + D.O(ok) - 1));

s.name = 'C3 G+T+O = 1 pointwise';
s.kind = 'EXACT';
s.pass = worst < 1e-12;
s.detail = sprintf('max |G+T+O-1| = %.2e over %d points', worst, sum(ok));
end

% =====================================================================
function s = check_attribution()
% C4 (APPROX). On an origin-centred sphere the radial direction IS the normal,
% so alignment is 1 and gain should separate cleanly. A pure radial scaling
% should land in G; a push along an ambient direction outside the local
% (tangent + radial) span should land in O.
rng(3);
F = 30; R = 5.0;
A  = sphere_cap(4000, F, R, 0.6, 0.01);
base = sphere_cap(300, F, R, 0.6, 0.01);

Zg = base * 1.4;                                      % pure gain

% novel direction: orthogonal to the 3-dim subspace the cap occupies
[~,~,V] = svd(A - mean(A,1), 'econ');
novel = V(:, end).';
Zo = base + 2.0 * repmat(novel, 300, 1);              % pure off-manifold

Dg = gto2_core(A, Zg, 'k', 40, 'd', 2);
Do = gto2_core(A, Zo, 'k', 40, 'd', 2);

medAlign = median(Dg.align);
gG = median(Dg.G, 'omitnan');
oO = median(Do.O, 'omitnan');

s.name = 'C4 pure mechanisms attribute to their own component';
s.kind = 'APPROX';
s.pass = medAlign > 0.95 && gG > 0.9 && oO > 0.9;
s.detail = sprintf(['sphere alignment %.3f; pure gain -> median G = %.3f; ' ...
    'pure off-manifold -> median O = %.3f'], medAlign, gG, oO);
end

% =====================================================================
function s = check_prop1()
% C5 (EXACT). Proposition 1 in the anchor frame. Every quantity is computed from
% the REALISED mu, Vd and ng, never from the construction's intent:
%   gamma = ng'(z_B - m),  delta = ng'(m - mu),
%   b = residual-normal part of (z_B - m),  eta = residual-normal part of (m - mu),
%   t = tangential part of (z_B - m) + (m - mu).
% Then ||r||^2 = (gamma+delta)^2 + ||b + eta||^2 + ||t||^2 and G, O, T follow
% exactly. Prop. 1 as printed assumes the imposed change has no tangential part in
% the anchor frame; the share it actually has is reported, not assumed away.
F = 30; R = 5.0; nRef = 4000; k = 40; sig = 0.01; nTest = 300; seeds = 0:2;
tab = zeros(numel(seeds), 6);
for si = 1:numel(seeds)
    rng(seeds(si));
    [Q,~] = qr(randn(F));  S = Q(:,1:3);
    v = randn(nRef,3); v = v ./ vecnorm(v,2,2);
    A = [R*v, zeros(nRef,F-3)]*Q.' + sig*randn(nRef,F);
    b = randn(nTest,3); b = b ./ vecnorm(b,2,2);
    M = [R*b, zeros(nTest,F-3)]*Q.';
    Z = zeros(nTest,F);
    for i = 1:nTest
        m = M(i,:).';  rad = m/R;
        ov = randn(F,1); ov = ov - S*(S.'*ov); ov = ov/norm(ov);
        Os = rand;  nu = 0.3 + 0.7*rand;  sg = sign(randn);
        Z(i,:) = (m + nu*(sg*sqrt(1-Os)*rad + sqrt(Os)*ov)).';
    end
    Z = Z + sig*randn(size(Z));
    D = gto2_core(A, Z, 'k', k, 'd', 2);
    [~, ord] = sort(pdist2(Z, A), 2);
    worst = zeros(1,5);  tshare = nan(nTest,1);
    for i = 1:nTest
        if D.degenerate(i), continue, end
        loc = A(ord(i,1:k), :);  muL = mean(loc,1).';
        [~,~,V] = svd(loc - muL.', 'econ');  Vd = V(:,1:2);
        mu = D.mu(i,:).';  ng = D.ng(i,:).';
        dz = Z(i,:).' - M(i,:).';  dm = M(i,:).' - mu;
        Pt = @(x) Vd*(Vd.'*x);
        gam = ng.'*dz;  del = ng.'*dm;
        bv  = dz - Pt(dz) - gam*ng;
        ev  = dm - Pt(dm) - del*ng;
        tv  = Pt(dz) + Pt(dm);
        r   = Z(i,:).' - mu;  r2 = sum(r.^2);
        r2p = (gam+del)^2 + sum((bv+ev).^2) + sum(tv.^2);
        Gp = (gam+del)^2/r2p;  Op = sum((bv+ev).^2)/r2p;  Tp = sum(tv.^2)/r2p;
        worst = max(worst, [abs(r2p-r2), abs(Gp-D.G(i)), abs(Tp-D.T(i)), abs(Op-D.O(i)), ...
                            abs(D.G(i)+D.T(i)+D.O(i)-1)]);
        tshare(i) = sum(Pt(dz).^2) / r2;
    end
    tab(si,:) = [worst, median(tshare, 'omitnan')];
end
fprintf('\n  C5 / Table S1: maximum absolute departure over %d probes per seed\n', nTest);
fprintf('  %4s %12s %12s %12s %12s %12s %14s\n', 'seed', '|dr^2|', '|dG|', '|dT|', '|dO|', ...
    '|G+T+O-1|', 'imposed T share');
for si = 1:numel(seeds)
    fprintf('  %4d %12.1e %12.1e %12.1e %12.1e %12.1e %14.2e\n', seeds(si), tab(si,:));
end
s.name = 'C5 Proposition 1 is exact in the anchor frame';
s.kind = 'EXACT';
s.table = tab;
s.pass = max(max(tab(:,1:5))) < 1e-10;
s.detail = sprintf(['max departure %.1e over %d seeds; median share of ||r||^2 that the imposed ' ...
    'change places in the ANCHOR tangent space = %.1e (Prop. 1 as printed assumes 0)'], ...
    max(max(tab(:,1:5))), numel(seeds), median(tab(:,6)));
end

% =====================================================================
function X = sphere_cap(n, F, R, capang, noise)
% Uniform-in-angle cap of a sphere of radius R centred at the ORIGIN, embedded
% in F dimensions via a fixed random rotation. Origin-centred by design: the
% radial direction is exactly the normal, so alignment is 1.
st = rng; rng(4242); [Q,~] = qr(randn(F)); rng(st);
ang = capang * sqrt(rand(n,1));                 % ~uniform on the cap for d=2
dir2 = randn(n,2); dir2 = dir2 ./ vecnorm(dir2,2,2);
u = [sin(ang).*dir2, cos(ang), zeros(n, F-3)];
X = R * u * Q.';
if noise > 0, X = X + noise*randn(size(X)); end
end

% ---------------------------------------------------------------------
function X = cone_patch(n, F, noise)
% Patch of a circular cone with apex at the ORIGIN: p(s,theta) =
% s*(cos t, sin t, 1)/sqrt(2). The position vector is tangent to the surface
% everywhere (d p/d s = p/s), so the radial direction lies in the tangent
% space and the gain axis is undefined.
st = rng; rng(777); [Q,~] = qr(randn(F)); rng(st);
s = 3 + 2*rand(n,1);
t = 2*pi*rand(n,1);
u = [s.*cos(t), s.*sin(t), s, zeros(n, F-3)] / sqrt(2);
X = u * Q.';
if noise > 0, X = X + noise*randn(size(X)); end
end

% ---------------------------------------------------------------------
function report(out)
f = fieldnames(out);
fprintf('\n%s\n', repmat('=', 1, 78));
fprintf('  GTO2 INVARIANT CHECKS\n');
fprintf('%s\n', repmat('=', 1, 78));
np = 0;
for i = 1:numel(f)
    s = out.(f{i});
    tag = 'FAIL'; if s.pass, tag = 'PASS'; np = np + 1; end
    fprintf('  [%s] %-6s %-48s\n         %s\n', tag, s.kind, s.name, s.detail);
end
fprintf('%s\n  %d of %d passed\n\n', repmat('-', 1, 78), np, numel(f));
end
