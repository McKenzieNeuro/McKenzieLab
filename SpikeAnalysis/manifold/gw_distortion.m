function out = gw_distortion(ZA, ZB, varargin)
% GW_DISTORTION  Do stimulated events PRESERVE or DISTORT the intrinsic geodesic
% geometry of the spontaneous neural manifold? (Gromov-Wasserstein.)
%
%   ZA : spontaneous reference [N_A x F]
%   ZB : stimulation           [N_B x F]
%
% Geodesics are computed WITHIN each condition on its own kNN graph, so GW compares
% the two as abstract metric spaces -- no shared embedding or shared latent needed,
% and relabeling / occupancy ordering cannot masquerade as distortion. Each geodesic
% matrix is normalized by its mean distance before GW, so the GW value reflects SHAPE
% distortion only; global gain is reported separately as the mean-distance ratio.
%
% Inference is against a within-spontaneous split-half null (GW between two disjoint
% spontaneous subsamples) -> z-score / percentile. A within-stim split-half GW is
% also returned; if stim is merely noisy its self-reliability will already be poor.
% Supports are equalized by farthest-point (MaxMin) subsampling to a common n.
%
% NOTE: the GW value is from an ENTROPIC solver (Peyre et al. 2016), so it is
% regularized -- interpret via z / pct against the null (same epsilon throughout),
% not as an absolute number. Solver validated against Python POT before transcription;
% sanity-check in MATLAB by passing two halves of ZA as (ZA,ZB): z should be ~0.
%
% Domain note: GW on cosine geodesics tests DIRECTIONAL/PATTERN geometry and is
% blind to amplitude (cos(c*x,y)=cos(x,y)). Use scaled_version_test for the gain /
% on-manifold question; the two are complementary. Run the on-manifold gate first --
% GW shape distortion is only meaningful for points that are on-manifold.
%
% Validated (Python, swiss roll): preserved z~0, gain x1.5 z~0, arc-length warp z~+4.

ip = inputParser;
ip.addParameter('k', 10);
ip.addParameter('metric', 'cosine');
ip.addParameter('nSub', 300);        % GW problem size (cost ~ O(nSub^2) per Sinkhorn step)
ip.addParameter('nNull', 200);
ip.addParameter('epsilon', 0.05);    % entropic regularization (matrices are mean-normalized)
ip.addParameter('gwIters', 100);
ip.addParameter('skIters', 100);
ip.addParameter('seed', 0);
ip.parse(varargin{:});
o = ip.Results; rng(o.seed);

% within-condition geodesics on a common, coverage-matched support
GA_full = geodesics(ZA, o.k, o.metric);
GB_full = geodesics(ZB, o.k, o.metric);
m  = min([o.nSub, size(GA_full,1), size(GB_full,1)]);
sa = maxmin(GA_full, m);  sb = maxmin(GB_full, m);
GA = GA_full(sa, sa);     GB = GB_full(sb, sb);

[shape_gw, gain, T] = shapeGW(GA, GB, o.epsilon, o.gwIters, o.skIters);

% within-spontaneous split-half null
nullv = zeros(o.nNull, 1);  nA = size(GA_full, 1);  h = floor(nA/2);
for kk = 1:o.nNull
    pr = randperm(nA);  s1 = pr(1:h);  s2 = pr(h+1:end);
    i1 = maxmin(GA_full(s1, s1), m);  i2 = maxmin(GA_full(s2, s2), m);
    nullv(kk) = shapeGW(GA_full(s1(i1), s1(i1)), GA_full(s2(i2), s2(i2)), ...
                        o.epsilon, o.gwIters, o.skIters);
end
mu = mean(nullv);  sd = std(nullv);

% within-stim reliability
nB = size(GB_full, 1);  pB = randperm(nB);  hb = floor(nB/2);
b1 = pB(1:hb);  b2 = pB(hb+1:end);
jb1 = maxmin(GB_full(b1, b1), m);  jb2 = maxmin(GB_full(b2, b2), m);
stim_reliab = shapeGW(GB_full(b1(jb1), b1(jb1)), GB_full(b2(jb2), b2(jb2)), ...
                      o.epsilon, o.gwIters, o.skIters);

% per-spontaneous-point distortion (square-loss factorization; localizes failures)
C1 = GA / meanOffDiag(GA);  C2 = GB / meanOffDiag(GB);
Trow = sum(T, 2);  Tcol = sum(T, 1).';
L = ((C1.^2) * Trow) * ones(1, m) + ones(m, 1) * ((C2.^2) * Tcol).' - 2 * (C1 * T * C2);
local = sum(T .* L, 2);

out.shape_gw         = shape_gw;
out.gain             = gain;                 % mean-distance ratio (isometry-up-to-scale)
out.z                = (shape_gw - mu) / sd; % PRIMARY verdict vs spontaneous noise floor
out.pct              = mean(nullv >= shape_gw);
out.null_mean        = mu;
out.null_sd          = sd;
out.stim_reliability = stim_reliab;          % guardrail: compare to shape_gw
out.T                = T;                    % optimal coupling (rows = support_idx)
out.local_distortion = local;
out.support_idx      = sa;
end

% ---------------------------------------------------------------------------
function G = geodesics(X, k, metric)
% within-condition geodesics: symmetrized union-kNN graph, largest component.
n  = size(X, 1);
D  = squareform(pdist(X, metric));
[~, nb] = sort(D, 2);
src = repmat((1:n).', k, 1);
dst = reshape(nb(:, 2:k+1), [], 1);
wt  = D(sub2ind([n n], src, dst));
A   = sparse(src, dst, wt, n, n);
A   = max(A, A.');                       % undirected union of kNN edges
Gr  = graph(A);
comp = conncomp(Gr);
keep = find(comp == mode(comp));
if numel(keep) < n
    warning('gw:disconnected', 'graph disconnected; keeping largest component (%d/%d). Raise k.', ...
            numel(keep), n);
end
G = distances(subgraph(Gr, keep));       % all-pairs shortest path
end

% ---------------------------------------------------------------------------
function idx = maxmin(G, m)
% farthest-point (MaxMin) subsample of m indices on the geodesic metric G.
n = size(G, 1);
if m >= n, idx = (1:n).'; return; end
idx = zeros(m, 1);  idx(1) = randi(n);  dmin = G(idx(1), :);
for i = 2:m
    [~, j] = max(dmin);  idx(i) = j;  dmin = min(dmin, G(j, :));
end
end

% ---------------------------------------------------------------------------
function [gw, gain, T] = shapeGW(Ga, Gb, epsilon, gwIters, skIters)
% entropic Gromov-Wasserstein (Peyre et al. 2016, square loss) on mean-normalized
% geodesic matrices -> shape distortion. gain = ratio of mean distances.
na = meanOffDiag(Ga);  nb = meanOffDiag(Gb);
C1 = Ga / na;  C2 = Gb / nb;
n = size(C1, 1);  m = size(C2, 1);
p = ones(n, 1) / n;  q = ones(m, 1) / m;
constC = ((C1.^2) * p) * ones(1, m) + ones(n, 1) * (q.' * (C2.^2));
T = p * q.';
for it = 1:gwIters
    Ggrad = 2 * (constC - 2 * (C1 * T * C2));       % square-loss gradient (hC2 = 2*C2)
    K = exp(-(Ggrad - min(Ggrad(:))) / epsilon);    % stabilized kernel
    u = ones(n, 1);
    for s = 1:skIters
        v = q ./ (K.' * u + 1e-300);
        u = p ./ (K  * v + 1e-300);
    end
    T = (u .* K) .* v.';                            % diag(u) K diag(v)
end
gw   = sqrt(max(sum(sum((constC - 2 * (C1 * T * C2)) .* T)), 0));
gain = nb / na;
end

% ---------------------------------------------------------------------------
function s = meanOffDiag(G)
n = size(G, 1);
s = (sum(G(:)) - sum(diag(G))) / (n * (n - 1));
end
