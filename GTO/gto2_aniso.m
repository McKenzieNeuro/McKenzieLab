function out = gto2_aniso(mode, varargin)
% GTO2_ANISO  Anisotropic observation noise: which directions are absorbed.
%
%   gto2_aniso('fractions')  % effect on G, T, O and the alignment  (gate 4)
%   gto2_aniso('anchor')     % effect on the anchor displacement    (gate 3, Prop 1)
%   gto2_aniso('frame')      % effect on the frame and on tau^2     (gates 2 and 6)
%   gto2_aniso('all')
%
% Every other file in this suite assumes isotropic noise, Sigma = sigma^2 * I.
% Real recordings are not isotropic, and this file establishes what that costs.
%
% THE OBVIOUS CORRECTION IS WRONG
%   Corollary 2 carries noise as a scalar times a dimension count in each leg:
%   sigma^2 in the gain leg, (F-d-1)*sigma^2 in the novel leg, d*sigma^2 in the
%   tangent leg. The natural generalisation replaces each by a directional
%   quadratic form -- ng'*Sigma*ng, trace(P_res*Sigma*P_res), trace(Vd'*Sigma*Vd).
%   That substitution is algebraically correct for a FIXED anchor and a FIXED
%   frame, and it predicts effects that DO NOT OCCUR. It is not the right model,
%   and this sweep is arranged to show why.
%
%   The anchor and the tangent basis are estimated from the SAME noisy reference
%   cloud that the probe is drawn from. Noise shared between reference and probe
%   is therefore partly absorbed rather than measured: the neighbourhood is
%   thickened along the noisy direction, and mu and Vd move with it. Whether a
%   given anisotropy matters depends entirely on WHICH direction it occupies.
%
% WHAT THE SWEEP SHOWS
%   Anisotropy is imposed as Sigma = sigma^2 (I + kappa*w*w') for a unit w placed
%   in one of three geometric roles, on an origin-centred spherical cap small
%   enough that the local frame is near-constant:
%
%     w TANGENT  -- no effect on the fractions out to sigma_w/ell > 2. The anchor
%                   tracks tangential displacement, which is what makes T
%                   unidentifiable; the same mechanism makes the decomposition
%                   immune to tangential noise. The cost and the protection are
%                   the same fact.
%
%     w AMBIENT (in the residual normal space) -- almost no effect. The noise
%                   enters the numerator and the denominator of O together, and
%                   the reference cloud is thickened along w as well, so mu
%                   partly follows.
%
%     w ALONG THE GAIN AXIS -- NOT absorbed. The alignment a = ||P_N rho||
%                   degrades and energy moves out of G and into O, monotonically
%                   in sigma_w/ell and with no threshold to hide behind. Gain
%                   modulation is read as novelty.
%
%   The third case is the one that matters for population recordings, because a
%   shared multiplicative fluctuation -- common-mode gain drift across the
%   population -- IS a radial direction, and on an origin-centred manifold the
%   radial direction IS the gain axis. The worst available alignment for this
%   failure is the one the data is most likely to supply.
%
% THE ALIGNMENT TOLERANCE DOES NOT PROTECT
%   The damage is done at alignments far above any sensible alignTol: by the time
%   a would trip a 0.1 threshold, G has already been largely converted to O. The
%   alignment gate detects a gain axis that does not EXIST; it does not detect a
%   gain axis that exists but has been contaminated. Report a alongside G and O
%   and treat a drop from unity as quantitative, not as a pass/fail.
%
% OPTIONS
%   'F' 30, 'R' 5, 'cap' 0.30, 'nRef' 1500, 'k' 100, 'nTest' 300
%   'noise'  0.01      isotropic scale sigma
%   'kgrid'  [0 100 400 900 1600 6400]     anisotropy factors kappa
%   'roles'  {'gain','tangent','ambient'}
%   'nu' 1.0, 'Ostar' 0.5                  imposed displacement and off-manifold share
%   'kAnchor' 40, 'nTrialAnchor' 2000       neighbourhood size and trials for 'anchor' (Table S10)
%   'kFrame'  40                            neighbourhood size for 'frame' (Table S12)
%   'seed' 0                                every cell is reseeded (paired comparisons)
%
% ISOTROPY IS LOAD-BEARING IN MORE THAN ONE PLACE
%   Proposition 1's claim that noise is not a bias channel for the anchor, the
%   gate-2 rotation law, and the gate-6 tangential term are all derived or
%   measured under isotropy. The 'anchor' and 'frame' modes test them directly.
%   Two of the three do not survive:
%     - PROPOSITION 1 FAILS. Radial anisotropy ERASES the sagitta (kNN selection
%       pulls the neighbours' mean radius toward the probe, cancelling the inward
%       displacement); tangential anisotropy INFLATES it (tangentially displaced
%       points from further around the manifold become selectable, so the selected
%       set spans a larger arc than the measured ell reports).
%     - THE GATE-2 ROTATION LAW FAILS, and badly. Radial anisotropy drives
%       ||Vd'n||^2 up by orders of magnitude: the normal direction is captured
%       into the estimated tangent space outright.
%     - GATE 6 REMAINS CHECKABLE. tau^2 does move with the anisotropy direction,
%       but it tracks ell^2/k measured on the SAME neighbourhood, so the observable
%       conditioning check survives. The frame mode prints tau^2/(ell^2/k).
%
% TWO MECHANISMS
%   Vd absorbing the radial direction happens through the VARIANCE route
%   (over-estimated d, gain-aligned anisotropic noise: normal variance competes
%   for a singular-vector slot) or the COVARIANCE route (finite-sample rotation,
%   curvature gradients, boundaries). All present identically -- the alignment
%   falls and energy moves from G to O. Finite-sample rotation is a floor that
%   scales as 1/nRef; the others have no such floor.
%
% RETURNS out.rows with columns in out.colnames.

if nargin < 1, mode = 'fractions'; end
ip = inputParser;
ip.addParameter('F', 30);      ip.addParameter('R', 5.0);
ip.addParameter('cap', 0.30);  ip.addParameter('nRef', 1500);
ip.addParameter('k', 100);     ip.addParameter('nTest', 300);
ip.addParameter('noise', 0.01);
ip.addParameter('kgrid', [0 100 400 900 1600 6400]);
ip.addParameter('roles', {'gain','tangent','ambient'});
ip.addParameter('nu', 1.0);    ip.addParameter('Ostar', 0.5);
ip.addParameter('kAnchor', 40);   ip.addParameter('nTrialAnchor', 2000);
ip.addParameter('kFrame', 40);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results;

switch lower(mode)
    case 'fractions'
    case 'anchor',  out = aniso_anchor(o);  return
    case 'frame',   out = aniso_frame(o);   return
    case 'all'
        out.anchor = aniso_anchor(o);
        out.frame  = aniso_frame(o);
    otherwise, error('mode must be fractions, anchor, frame or all');
end

out.rows = [];
out.colnames = {'role','kappa','sigma_w','ell','sw_over_ell','G','T','O','align','fracDegen'};
fprintf('\n=== Anisotropic noise: Sigma = sigma^2 (I + kappa w w'') ===\n');
fprintf('  cap = %.2f, R = %.1f, F = %d, nRef = %d, k = %d, sigma = %g\n', ...
    o.cap, o.R, o.F, o.nRef, o.k, o.noise);
fprintf('  imposed (G,T,O) = (%.2f, 0, %.2f) at displacement %.2f\n', ...
    1-o.Ostar, o.Ostar, o.nu);

for ri = 1:numel(o.roles)
    role = o.roles{ri};
    fprintf('\n  -- w along the %s direction --\n', upper(role));
    fprintf('  %8s %10s %10s %8s %8s %8s %8s %8s\n', ...
        'kappa', 'sigma_w', 'sw/ell', 'G', 'T', 'O', 'align', '%degen');
    for kap = o.kgrid(:).'
        r = aniso_cell(kap, role, o);
        out.rows(end+1,:) = [ri kap r.sw r.ell r.sw/r.ell r.G r.T r.O r.align r.degen]; %#ok<AGROW>
        fprintf('  %8d %10.4f %10.3f %8.4f %8.4f %8.4f %8.4f %8.1f\n', ...
            kap, r.sw, r.sw/r.ell, r.G, r.T, r.O, r.align, 100*r.degen);
    end
end
fprintf(['\n  Tangential and ambient anisotropy are absorbed; gain-aligned anisotropy is not,\n' ...
         '  and converts G into O with no threshold. Note the alignment is still well above\n' ...
         '  any usable alignTol while the damage is already severe.\n\n']);
end

% =====================================================================
function r = aniso_cell(kappa, role, o)
rng(o.seed);
[Q,~] = qr(randn(o.F));
S    = Q(:,1:3);
axis_ = Q(:,3);  t1 = Q(:,1);  t2 = Q(:,2);  amb = Q(:,8);
switch role
    case 'gain',    w = axis_;      % cap axis: radial, hence the gain axis
    case 'tangent', w = t1;         % in the tangent plane at the cap centre
    case 'ambient', w = amb;        % in the residual normal space
end
w = w / norm(w);
nz = @(n) o.noise * (randn(n, o.F) + sqrt(kappa) * (randn(n,1) * w.'));

cappts = @(n) capsample(n, o, Q);
A = cappts(o.nRef) + nz(o.nRef);

base = cappts(o.nTest);
Z = zeros(o.nTest, o.F);
for i = 1:o.nTest
    m   = base(i,:).';
    rad = m / norm(m);
    ov  = randn(o.F,1);  ov = ov - S*(S.'*ov);  ov = ov / norm(ov);
    Z(i,:) = (m + o.nu*(sqrt(1-o.Ostar)*rad + sqrt(o.Ostar)*ov)).';
end
Z = Z + nz(o.nTest);

D = gto2_core(A, Z, 'k', o.k, 'd', 2);
keep = ~D.degenerate;
r.ell   = median(D.ell(keep));
r.sw    = o.noise * sqrt(1 + kappa);
r.G     = median(D.G(keep));
r.O     = median(D.O(keep));
r.T     = median(D.T(keep));
r.align = median(D.align);
r.degen = mean(D.degenerate);
end

% ---------------------------------------------------------------------
function X = capsample(n, o, Q)
a  = o.cap * sqrt(rand(n,1));
ph = 2*pi*rand(n,1);
u  = [sin(a).*cos(ph), sin(a).*sin(ph), cos(a)];
X  = [o.R*u, zeros(n, o.F-3)] * Q.';
end

% =====================================================================
function s = aniso_anchor(o)
% GATE 3 / PROPOSITION 1 under anisotropy. Measures the anchor displacement
% directly at a single base point, so 'radial' and 'tangential' are well defined.
% Proposition 1 predicts delta = ell^2/(2R) and (verified isotropically) is
% insensitive to the noise scale. Neither survives anisotropy.
kg = [0 25 100 400];  sig = 0.05;  nTrial = o.nTrialAnchor;  k = o.kAnchor;  nRef = 4000;
s.rows = []; s.colnames = {'role','kappa','sigma_w','ell','delta','deltaPred','ratio', ...
                           'deltaPred2','ratio2','ratio2SE'};
fprintf('\n=== Gate 3: does anisotropic noise bias the ANCHOR? ===\n');
fprintf('  R = %.1f, F = %d, nRef = %d, k = %d, sigma = %g, %d trials/cell\n\n', ...
    o.R, o.F, nRef, k, sig, nTrial);
fprintf('  ratio2 = delta/[(ell^2/2R)(1 + ell^2/3R^2)]; SE is the Monte Carlo error of ratio2\n\n');
fprintf('  %9s %8s %9s %8s %10s %10s %8s %8s %7s\n', ...
    'w role', 'kappa', 'sigma_w', 'ell', 'delta', 'ell^2/2R', 'ratio', 'ratio2', 'SE');
roles = {'radial','tangent','ambient'};
for ri = 1:3
    for kap = kg
        [dlt, ell, dse] = anchor_cell_aniso(kap, roles{ri}, sig, nRef, k, nTrial, o);
        pred  = ell^2/(2*o.R);
        pred2 = pred*(1 + ell^2/(3*o.R^2));
        s.rows(end+1,:) = [ri kap sig*sqrt(1+kap) ell dlt pred dlt/pred pred2 dlt/pred2 dse/pred2]; %#ok<AGROW>
        fprintf('  %9s %8d %9.3f %8.4f %10.5f %10.5f %8.3f %8.3f %7.3f\n', ...
            roles{ri}, kap, sig*sqrt(1+kap), ell, dlt, pred, dlt/pred, dlt/pred2, dse/pred2);
    end
end
fprintf(['\n  Radial anisotropy ERASES the sagitta; tangential anisotropy INFLATES it;\n' ...
         '  ambient is inert. Both departures are in the measured displacement itself,\n' ...
         '  not in the convention for ell -- delta falls while ell rises, and vice versa.\n\n']);
end

% ---------------------------------------------------------------------
function [dlt, ell, dse] = anchor_cell_aniso(kappa, role, sig, nRef, k, nTrial, o)
rng(o.seed);                                % paired draws across roles and kappa
u = [0 0 1].';
rad = zeros(o.F,1); rad(1:3) = u;
tan_ = zeros(o.F,1); tan_(1) = 1;
amb = zeros(o.F,1); amb(8) = 1;
switch role, case 'radial', w = rad; case 'tangent', w = tan_; case 'ambient', w = amb; end
dl = zeros(nTrial,1); el = zeros(nTrial,1);
for t = 1:nTrial
    v = randn(nRef,3); v = v ./ vecnorm(v,2,2);
    X = [o.R*v, zeros(nRef,o.F-3)];
    X = X + sig*(randn(nRef,o.F) + sqrt(kappa)*(randn(nRef,1)*w.'));
    z = zeros(1,o.F); z(1:3) = o.R*u.';
    [~, ord] = sort(sum((X - z).^2, 2));
    loc = X(ord(1:k), :);  mu = mean(loc,1);
    m = zeros(1,o.F); m(1:3) = o.R*u.';
    dl(t) = -(mu - m) * rad;
    el(t) = sum(sum((loc(:,1:2) - mu(1:2)).^2, 2)) / (k - 1);   % tangential, unbiased
end
dlt = mean(dl);  ell = sqrt(mean(el));  dse = std(dl)/sqrt(nTrial);
end

% =====================================================================
function s = aniso_frame(o)
% GATES 2 and 6 under anisotropy, on a cap small enough that the local frame is
% near-constant. Reports the rotation ||Vd'n||^2 (gate 2's error term) and the
% tangential energy tau^2 (gate 6's condition). The first is direction-sensitive
% by orders of magnitude; the second is not.
kg = [0 25 100 400];  sig = 0.02;  nRef = 4000; k = o.kFrame; nTest = 200; cap = 0.25;
s.rows = []; s.colnames = {'role','kappa','sigma_w','rot','tau2','ell','tau2_over_ell2k'};
fprintf('\n=== Gates 2 and 6 under anisotropy ===\n');
fprintf('  cap %.2f, R = %.1f, F = %d, nRef = %d, k = %d, sigma = %g, %d probes; medians\n\n', ...
    cap, o.R, o.F, nRef, k, sig, nTest);
fprintf('  %9s %8s %9s %12s %11s %9s %12s\n', ...
    'w role', 'kappa', 'sigma_w', '||Vd''n||^2', 'tau^2', 'ell^2/k', 'tau^2/(l^2/k)');
roles = {'radial','tangent','ambient'};
for ri = 1:3
    for kap = kg
        [rot, tau2, ell] = frame_cell(kap, roles{ri}, sig, nRef, k, nTest, cap, o);
        s.rows(end+1,:) = [ri kap sig*sqrt(1+kap) rot tau2 ell tau2/(ell^2/k)]; %#ok<AGROW>
        fprintf('  %9s %8d %9.3f %12.3e %11.5f %9.5f %12.3f\n', ...
            roles{ri}, kap, sig*sqrt(1+kap), rot, tau2, ell^2/k, tau2/(ell^2/k));
    end
end
fprintf(['\n  Gate 2''s rotation term is NOT direction-blind: radial anisotropy raises it by\n' ...
         '  orders of magnitude, capturing the normal into Vd outright. tau^2 moves with the\n' ...
         '  anisotropy but should track ell^2/k measured on the same neighbourhood.\n\n']);
end

% ---------------------------------------------------------------------
function [rot, tau2, ell] = frame_cell(kappa, role, sig, nRef, k, nTest, cap, o)
rng(o.seed);                                % paired draws across roles and kappa
[Q,~] = qr(randn(o.F));  S = Q(:,1:3);
rad = Q*[0;0;1;zeros(o.F-3,1)];  tan_ = Q*[1;0;0;zeros(o.F-3,1)];  amb = Q(:,8);
switch role, case 'radial', w = rad; case 'tangent', w = tan_; case 'ambient', w = amb; end
w = w/norm(w);
nz = @(n) sig*(randn(n,o.F) + sqrt(kappa)*(randn(n,1)*w.'));
pts = @(n) capsample2(n, cap, o, Q);
A = pts(nRef) + nz(nRef);
base = pts(nTest);
N = base ./ vecnorm(base,2,2);
ov = randn(nTest,o.F);  ov = ov - (ov*S)*S.';  ov = ov ./ vecnorm(ov,2,2);
Z = base + 0.5*ov + nz(nTest);
D = gto2_core(A, Z, 'k', k, 'd', 2);
keep = ~D.degenerate;
[~, ord] = sort(pdist2(Z, A), 2);
rv = zeros(nTest,1);
for i = 1:nTest
    loc = A(ord(i,1:k), :);
    [~,~,V] = svd(loc - mean(loc,1), 'econ');
    rv(i) = sum((V(:,1:2).' * N(i,:).').^2);
end
rot  = median(rv(keep));
tau2 = median(D.T(keep) .* D.rnorm(keep).^2);
ell  = median(D.ell(keep));
end

% ---------------------------------------------------------------------
function X = capsample2(n, cap, o, Q)
a = cap*sqrt(rand(n,1)); ph = 2*pi*rand(n,1);
v = [sin(a).*cos(ph), sin(a).*sin(ph), cos(a)];
X = [o.R*v, zeros(n, o.F-3)] * Q.';
end
