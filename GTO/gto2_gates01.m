function out = gto2_gates01(mode, varargin)
% GTO2_GATES01  The two gates at the top of the cascade.
%
%   gto2_gates01('gate1')   % dimension excess: leakage, survival, alignment
%   gto2_gates01('gate0')   % sheet proximity: where a chart failure is booked
%   gto2_gates01('all')
%
% GATE 1. For an imposed purely off-manifold direction v, the leakage
%       L(v) = ||Vd' vhat||^2
% is the fraction absorbed by the estimated tangent basis. If the q = d - d_true
% excess directions were RANDOM in the ambient complement, the surviving fraction
% would follow the Beta survival law 1 - q/(F - d_true). They are not always
% random, and the difference is the point of this test.
%
%   On a curved manifold the local patch has a genuine third principal direction:
%   the NORMAL, along which the patch bends by the sagitta. Whenever the variance
%   that bending contributes along the normal (sd ell^2/(2*sqrt(3)*R) for a uniform
%   disc) beats the largest sampled noise variance, the first excess singular vector
%   is not a random ambient direction but the normal itself ('capture' measures
%   this). Two consequences, both measured here:
%     (i)  O is BETTER preserved than the survival law predicts at q = 1, because
%          the absorbed direction is not one that carried novelty.
%     (ii) The gain axis is DESTROYED, because ng is built from rho inside the
%          normal space and rho has just been absorbed into Vd. The alignment
%          does not degrade; it COLLAPSES at q = 1 and is flat thereafter.
%   In the opposite regime, noise-dominated, the excess directions are random,
%   the survival law holds, and the alignment degrades smoothly. Over-estimating d
%   is therefore dangerous in a way that depends on the curvature-to-noise ratio,
%   and the damage lands on gate 4 rather than where the survival law points.
%
% GATE 0. Two concentric spherical shells separated by h, both origin-centred, so
% the gain axis is defined. Probes sit on the inner shell and are purely
% on-manifold: any structure in r is an artefact of the chart failing. As h falls
% through the neighbourhood extent ell the kNN set bridges the sheets.
%
%   The result is a booking rule directly parallel to gate 5's. A chart failure is
%   apportioned according to the DIRECTION IN WHICH THE SHEETS ARE STACKED:
%     - stacked RADIALLY (along ng): the novel energy stays pinned at the noise
%       floor at every separation, while the radial energy spikes. Bridging
%       masquerades as GAIN.
%     - stacked in an AMBIENT direction: the reverse. Radial energy stays at or
%       below control while novel energy spikes. Bridging masquerades as NOVELTY.
%   Both peak near h ~ 0.8*ell and subside on either side: at h >> ell there is no
%   bridging, and at h << ell the two sheets are indistinguishable from one
%   slightly thicker sheet.
%
%   PRACTICAL WARNING. The fraction of neighbours drawn from the wrong sheet is
%   NOT a usable diagnostic. At h << ell that fraction reaches 50% while the
%   decomposition is essentially unharmed; at h ~ ell it is around 40% and the
%   artefact is at its worst. Cross-sheet counting identifies mixing, not damage.
%
% OPTIONS (both modes): 'F' 30, 'R' 5, 'nRef', 'k' 40, 'noise', 'nTest', 'seed'.
%
% RETURNS out.<mode>.rows with columns in out.<mode>.colnames.

if nargin < 1, mode = 'all'; end
out = struct();
switch lower(mode)
    case 'gate1', out.gate1 = run_gate1(varargin{:});
    case 'gate0', out.gate0 = run_gate0(varargin{:});
    case 'all'
        out.gate1 = run_gate1(varargin{:});
        out.gate0 = run_gate0(varargin{:});
    otherwise, error('mode must be gate0, gate1 or all');
end
end

% =====================================================================
function s = run_gate1(varargin)
% Survival is compared with FOUR predictions, because the right comparator is the
% question at issue (Table S3):
%   predMean     1 - q/(F - d_true)                      mean of the Haar Beta law
%   predMedian   1 - betaincinv(0.5, q/2, (F-d)/2)        median of the same law
%   predMedFloor, predMeanFloor: the same, with the q = 0 leakage floor of THIS
%                seed folded in, S = (1 - L0)(1 - B), L0 drawn from the per-probe
%                q = 0 leakage and B ~ Beta(q/2,(F-d)/2) (independence assumed).
% 'capture' is the median squared overlap of the (d_true+1)-th singular vector with
% the true normal at the base point: ~1 when curvature wins that slot, ~0 when noise
% does. Sweeping sigma locates the transition; the variance estimate for where it
% sits is sigma ~ ell^2/(2*sqrt(3)*R), up to an O(1) extreme-value factor.
ip = inputParser;
ip.addParameter('F', 30);   ip.addParameter('R', 5.0);
ip.addParameter('nRef', 4000); ip.addParameter('k', 40);
ip.addParameter('nTest', 250);
ip.addParameter('dgrid', 2:6);
ip.addParameter('noiseGrid', [0.01 0.30]);
ip.addParameter('seeds', 0);
ip.addParameter('nMC', 1e5);
ip.parse(varargin{:});
o = ip.Results;
nS = numel(o.seeds);  nd = numel(o.dgrid);  dmax = max(o.dgrid);

s.rows = [];
s.colnames = {'sigma','d','leak','survival','survPred','align','fracDegen','medO', ...
              'survivalSD','leakMean','survMean','predMedian','predMedFloor','predMeanFloor','capture'};
ell2 = 2*o.k*o.R^2/o.nRef;
fprintf('\n=== Gate 1: dimension excess, leakage, and the alignment interaction ===\n');
fprintf('  R = %.1f, F = %d, nRef = %d, k = %d, %d probes x %d seeds\n', o.R, o.F, o.nRef, o.k, o.nTest, nS);
fprintf('  ell^2/2R = %.4f; sagitta sd along the normal ell^2/(2 sqrt3 R) = %.4f\n', ...
    ell2/(2*o.R), ell2/(2*sqrt(3)*o.R));
for sig = o.noiseGrid(:).'
    acc = nan(nS, nd, 12);
    for si = 1:nS
        rng(o.seeds(si));
        [Q,~] = qr(randn(o.F));  S = Q(:,1:3);
        v = randn(o.nRef,3); v = v ./ vecnorm(v,2,2);
        A = [o.R*v, zeros(o.nRef, o.F-3)] * Q.' + sig*randn(o.nRef, o.F);

        b = randn(o.nTest,3); b = b ./ vecnorm(b,2,2);
        Z = zeros(o.nTest,o.F);  V = zeros(o.nTest,o.F);  N = zeros(o.nTest,o.F);
        for i = 1:o.nTest
            m  = Q * [o.R*b(i,:).'; zeros(o.F-3,1)];
            ov = randn(o.F,1); ov = ov - S*(S.'*ov); ov = ov/norm(ov);  % truly ambient
            V(i,:) = ov.';  N(i,:) = (m/norm(m)).';  Z(i,:) = (m + 0.5*ov).';
        end
        Z = Z + sig*randn(size(Z));

        % one SVD per probe serves every d
        [~, ord] = sort(pdist2(Z, A), 2);
        Lcum = zeros(o.nTest, dmax);  cap3 = zeros(o.nTest,1);
        for i = 1:o.nTest
            loc = A(ord(i,1:o.k), :);
            [~,~,Vv] = svd(loc - mean(loc,1), 'econ');
            Lcum(i,:) = cumsum(((Vv(:,1:dmax).' * V(i,:).').^2).');
            cap3(i)   = (Vv(:,3).' * N(i,:).')^2;
        end
        L0 = Lcum(:,2);                          % q = 0 leakage floor, per probe
        rs = rng;  rng(1000 + o.seeds(si));      % separate stream for the Monte Carlo
        L0draw = L0(randi(numel(L0), o.nMC, 1));
        for di = 1:nd
            d = o.dgrid(di);  q = d - 2;
            D = gto2_core(A, Z, 'k', o.k, 'd', d);
            L = Lcum(:,d);
            pMean = 1 - q/(o.F-2);
            if q == 0
                pMed = 1;  B = zeros(o.nMC,1);
            else
                pMed = 1 - betaincinv(0.5, q/2, (o.F-d)/2);
                B = betarnd(q/2, (o.F-d)/2, o.nMC, 1);
            end
            Sfl = (1 - L0draw) .* (1 - B);
            acc(si,di,:) = [median(L), 1-median(L), pMean, median(D.align), mean(D.degenerate), ...
                            median(D.O(~D.degenerate)), mean(L), 1-mean(L), pMed, ...
                            median(Sfl), mean(Sfl), median(cap3)];
        end
        rng(rs);
    end
    fprintf('\n  sigma = %.3f\n', sig);
    fprintf('  %3s %9s %7s %9s | %9s %9s %9s %9s | %8s %7s %8s %8s\n', 'd', 'surv med', '(sd)', 'surv mean', ...
        'pred mean', 'pred med', 'med*flr', 'mean*flr', 'align', '%degen', 'capture', 'median O');
    for di = 1:nd
        a  = squeeze(mean(acc(:,di,:),1)).';
        sd = std(acc(:,di,2));
        if nS == 1, sd = NaN; end
        s.rows(end+1,:) = [sig o.dgrid(di) a(1) a(2) a(3) a(4) a(5) a(6) sd a(7) a(8) a(9) a(10) a(11) a(12)]; %#ok<AGROW>
        fprintf('  %3d %9.4f %7.4f %9.4f | %9.4f %9.4f %9.4f %9.4f | %8.4f %7.1f %8.3f %8.4f\n', ...
            o.dgrid(di), a(2), sd, a(8), a(3), a(9), a(10), a(11), a(4), 100*a(5), a(12), a(6));
    end
end
fprintf(['\n  Compare median survival with the median predictions and mean with mean. capture\n' ...
         '  near 1 means the first excess direction is the normal (alignment collapses at q = 1).\n\n']);
end

% =====================================================================
function s = run_gate0(varargin)
ip = inputParser;
ip.addParameter('F', 30);   ip.addParameter('R', 5.0);
ip.addParameter('nRef', 6000); ip.addParameter('k', 40);
ip.addParameter('noise', 0.01); ip.addParameter('nTest', 300);
ip.addParameter('hgrid', [4.0 1.0 0.5 0.25 0.10 0.05]);
ip.addParameter('nSeed', 5);
ip.parse(varargin{:});
o = ip.Results;

s.rows = []; s.colnames = {'stack','h','hOverEll','fracWrong','radAbs','novAbs','rnorm'};
fprintf('\n=== Gate 0: two sheets, separation h; probes PURELY on-manifold ===\n');
for st = 1:2
    stack = ternary(st==1, 'radial', 'ambient');
    fprintf('\n  -- sheets stacked %s --\n', upper(stack));
    fprintf('  %8s %8s %10s %10s %10s %9s\n', 'h', 'h/ell', 'wrong', '|r_rad|', '|r_nov|', '|r|');
    for h = [o.hgrid, Inf]
        acc = zeros(o.nSeed,5);
        for sd = 1:o.nSeed
            acc(sd,:) = gate0_cell(h, stack, sd, o);
        end
        m = mean(acc,1);
        lab = ternary(isinf(h), '  ctrl', sprintf('%8.2f', h));
        s.rows(end+1,:) = [st h m(1) m(2) m(3) m(4) m(5)]; %#ok<AGROW>
        fprintf('  %8s %8.2f %9.1f%% %10.5f %10.5f %9.5f\n', lab, m(1), 100*m(2), m(3), m(4), m(5));
    end
end
fprintf(['\n  Novel energy is pinned at the noise floor when the sheets are stacked radially,\n' ...
         '  and spikes when they are stacked ambiently; the radial energy does the opposite.\n' ...
         '  A chart failure is booked by the STACKING DIRECTION. Note the wrong-sheet\n' ...
         '  fraction is not diagnostic: it is largest where the damage is smallest.\n\n']);
end

% ---------------------------------------------------------------------
function v = gate0_cell(h, stack, seed, o)
rng(seed);
[Q,~] = qr(randn(o.F));
sh = @(n,rad,off) [rad*(randn(n,3)./vecnorm(randn(n,3),2,2)), zeros(n,o.F-3)]*Q.' + off;
n2 = floor(o.nRef/2);
u1 = randn(n2,3); u1 = u1./vecnorm(u1,2,2);
u2 = randn(n2,3); u2 = u2./vecnorm(u2,2,2);
if strcmp(stack,'radial')
    r2 = o.R + h; off = zeros(1,o.F);
    if isinf(h), r2 = 1e9; end
else
    r2 = o.R;  off = (Q(:,5)*h).';
    if isinf(h), off = (Q(:,5)*1e9).'; end
end
A = [ [o.R*u1, zeros(n2,o.F-3)]*Q.' ; [r2*u2, zeros(n2,o.F-3)]*Q.' + off ];
lab = [zeros(n2,1); ones(n2,1)];
A = A + o.noise*randn(size(A));

b = randn(o.nTest,3); b = b./vecnorm(b,2,2);
Z = [o.R*b, zeros(o.nTest,o.F-3)]*Q.' + o.noise*randn(o.nTest,o.F);

D = gto2_core(A, Z, 'k', o.k, 'd', 2);
[~, ord] = sort(pdist2(Z, A), 2);
fracWrong = mean(mean(lab(ord(:,1:o.k)), 2));
ell = median(D.ell);          % tangential extent (unbiased)
ok = ~D.degenerate;
v = [h/ell, fracWrong, ...
     median(sqrt(D.G(ok)).*D.rnorm(ok)), ...
     median(sqrt(D.O(ok)).*D.rnorm(ok)), ...
     median(D.rnorm(ok))];
end

% ---------------------------------------------------------------------
function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
