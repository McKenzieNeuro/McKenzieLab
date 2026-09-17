function out = gw_distortion_components(ZA, ZB, varargin)
% GW_DISTORTION_COMPONENTS  Two-level distortion analysis for a DISCONNECTED manifold.
% Self-contained: no dependency on gw_distortion.m (per-component GW is the local
% function runGW below).
%
%   ZA : spontaneous reference [N_A x F]     ZB : stimulation [N_B x F]  (shared ambient space)
%
% Geodesic distance exists only within a connected component, so GW is run PER COMPONENT,
% never globally (ill-posed) and never on tiny neighborhoods (no relational structure).
%   COARSE (report first): component counts/sizes, between-centroid distances, spont<->stim
%   matching. Merges/splits/novel components / components pushed apart are reorganization of
%   coarse structure -- a larger effect than within-component shape distortion.
%   FINE: within-component GW (runGW) for each cleanly-matched pair above minComp.
%
% Matching is point-level: each stim point -> nearest spont point (ambient) -> its spont
% component. >1 spont comp touched = MERGE; one spont comp matched by >1 stim comp = SPLIT;
% points beyond spont support (gapFrac high) = NOVEL. Between-centroid distances (comparable
% since ZA,ZB share coordinates) carry the coarse layout GW cannot see.

ip = inputParser;
ip.addParameter('k', 10);           ip.addParameter('metric', 'cosine');
ip.addParameter('nSub', 300);       ip.addParameter('nNull', 2);
ip.addParameter('epsilon', 0.05);   ip.addParameter('gwIters', 100);
ip.addParameter('skIters', 100);    ip.addParameter('seed', 0);
ip.addParameter('minComp', 150);    ip.addParameter('gapMult', 3);
ip.addParameter('touchFrac', 0.10); ip.addParameter('minPurity', 0.80);
ip.parse(varargin{:});
o = ip.Results; rng(o.seed);

CA = condComponents(ZA, o.k, o.metric);
CB = condComponents(ZB, o.k, o.metric);

% ---- drop sub-threshold fragments before anything else ----
% Singletons and tiny pieces arise from isolated points losing the kNN lottery;
% they carry no manifold structure and should not drive the reorganized flag.
CA.substantial = find(CA.size >= o.minComp);
CB.substantial = find(CB.size >= o.minComp);
nAS = numel(CA.substantial);  nBS = numel(CB.substantial);
if CA.ncomp > nAS
    fprintf('  spont: ignoring %d sub-threshold fragment(s) (sizes: %s)\n', ...
        CA.ncomp - nAS, num2str(CA.size(CA.size < o.minComp)'));
end
if CB.ncomp > nBS
    fprintf('  stim:  ignoring %d sub-threshold fragment(s) (sizes: %s)\n', ...
        CB.ncomp - nBS, num2str(CB.size(CB.size < o.minComp)'));
end

% ---- point-level matching (substantial components only) ----
spontKeep = ismember(CA.label, CA.substantial);
ZA_sub = ZA(spontKeep, :);  labA_sub = CA.label(spontKeep);
Dcross = pdist2(ZB, ZA_sub, o.metric);
[nnd, nn] = min(Dcross, [], 2);
labOfB = labA_sub(nn);
gapThr = o.gapMult * prctile(CA.kthDist(spontKeep), 99);

cont = zeros(nBS, nAS);
gapFrac = zeros(nBS, 1);  purity = zeros(nBS, 1);
matchSpont = zeros(nBS, 1);  nTouch = zeros(nBS, 1);
for ti = 1:nBS
    tc = CB.substantial(ti);
    m = CB.label == tc;  nb = nnz(m);
    for si = 1:nAS
        cont(ti, si) = nnz(labOfB(m) == CA.substantial(si));
    end
    gapFrac(ti) = mean(nnd(m) > gapThr);
    [mx, mi] = max(cont(ti, :));
    matchSpont(ti) = CA.substantial(mi);  purity(ti) = mx / nb;
    nTouch(ti) = nnz(cont(ti, :) > o.touchFrac * nb);
end

isMerge = nTouch > 1;
isNovel = gapFrac > 0.5;
isClean = ~isMerge & ~isNovel & purity >= o.minPurity;
sp = matchSpont(isClean);
splitSpont = unique(sp(arrayfun(@(x) sum(sp == x) > 1, sp)));

% ---- coarse layout (substantial components only) ----
coarse.ncompSpont    = nAS;  coarse.ncompStim    = nBS;
coarse.ncompSpontRaw = CA.ncomp;  coarse.ncompStimRaw = CB.ncomp;
coarse.sizeSpont = CA.size(CA.substantial);
coarse.sizeStim  = CB.size(CB.substantial);
centA = CA.centroid(CA.substantial, :);
centB = CB.centroid(CB.substantial, :);
coarse.centDistSpont = safeSquareform(pdist(centA, o.metric));
coarse.centDistStim  = safeSquareform(pdist(centB, o.metric));
coarse.centDistCross = pdist2(centA, centB, o.metric);

% layout ratio: for each pair of cleanly-matched substantial comps mapping to distinct spont comps
ci = find(isClean); layout = [];
for a = 1:numel(ci)
    for b = a+1:numel(ci)
        ti1 = ci(a); ti2 = ci(b);
        s1 = matchSpont(ti1); s2 = matchSpont(ti2);
        if s1 ~= s2
            si1 = find(CA.substantial == s1); si2 = find(CA.substantial == s2);
            layout(end+1, :) = [s1, s2, CB.substantial(ti1), CB.substantial(ti2), ...
                coarse.centDistStim(ti1, ti2) / coarse.centDistSpont(si1, si2)]; %#ok<AGROW>
        end
    end
end
coarse.layoutRatio = layout;   % [spont1 spont2 stim1 stim2  stimDist/spontDist]

% ---- fine: per-matched-component GW ----
comps = struct('spontComp', {}, 'stimComp', {}, 'nSpont', {}, 'nStim', {}, 'gw', {});
for ti = 1:nBS
    if ~isClean(ti), continue; end
    sc = matchSpont(ti);                    % original spont component label
    tc = CB.substantial(ti);               % original stim component label
    r = runGW(ZA(CA.label == sc, :), ZB(CB.label == tc, :), o);
    comps(end+1) = struct('spontComp', sc, 'stimComp', tc, ...
        'nSpont', CA.size(sc), 'nStim', CB.size(tc), 'gw', r); %#ok<AGROW>
end

match.contingency = cont; match.matchSpont = matchSpont; match.purity = purity;
match.gapFrac = gapFrac;  match.nTouch = nTouch;
match.isMerge = find(isMerge); match.isNovel = find(isNovel);
match.isSplit = splitSpont;    match.isClean = find(isClean);

out.spont = CA; out.stim = CB; out.match = match; out.coarse = coarse; out.components = comps;
% reorganized flag uses SUBSTANTIAL component counts only; raw count includes graph artifacts
out.reorganized = ~isempty(match.isMerge) || ~isempty(match.isNovel) || ...
                  ~isempty(splitSpont) || (nAS ~= nBS);

fprintf('COARSE: spont %d comps, stim %d comps (raw: %d / %d, fragments ignored) | merges=%d novel=%d splits=%d | GW on %d pairs\n', ...
    nAS, nBS, CA.ncomp, CB.ncomp, numel(match.isMerge), numel(match.isNovel), numel(splitSpont), numel(comps));
if out.reorganized
    fprintf('  -> coarse structure REORGANIZED; report this before within-component GW.\n');
end
end

% ===========================================================================
function res = runGW(ZA, ZB, o)
% within-condition GW distortion with split-half spontaneous null (was gw_distortion.m)
GA_full = geodesics(ZA, o.k, o.metric);
GB_full = geodesics(ZB, o.k, o.metric);
m  = min([o.nSub, size(GA_full, 1), size(GB_full, 1)]);
sa = maxmin(GA_full, m);  sb = maxmin(GB_full, m);
GA = GA_full(sa, sa);     GB = GB_full(sb, sb);
[shape_gw, gain, T] = shapeGW(GA, GB, o.epsilon, o.gwIters, o.skIters);

% NULL: random subsampling (NOT MaxMin) so each draw genuinely differs.
% MaxMin on a large manifold is nearly deterministic -- always picks the same
% extremal points -- giving null_sd ~ 0 and meaningless z-scores.
% The actual X-vs-Y/Z comparison keeps MaxMin for coverage; only the null uses random.
nullv = zeros(o.nNull, 1);  nA = size(GA_full, 1);
for kk = 1:o.nNull
    pr = randperm(nA, min(2*m, nA));
    s1 = pr(1:m);
    if numel(pr) >= 2*m
        s2 = pr(m+1:2*m);
    else
        s2 = randperm(nA, m);
    end
    nullv(kk) = shapeGW(GA_full(s1, s1), GA_full(s2, s2), ...
                        o.epsilon, o.gwIters, o.skIters);
end
mu = mean(nullv);  sd = std(nullv);

% within-stim reliability: two random halves of stim
nB = size(GB_full, 1);  hb = min(m, floor(nB/2));
pB = randperm(nB, min(2*hb, nB));
b1 = pB(1:hb);
if numel(pB) >= 2*hb
    b2 = pB(hb+1:2*hb);
else
    b2 = randperm(nB, hb);
end
stim_reliab = shapeGW(GB_full(b1, b1), GB_full(b2, b2), ...
                      o.epsilon, o.gwIters, o.skIters);

nm1 = size(GA,1); C1 = GA / ((sum(GA(:))-sum(diag(GA)))/(nm1*(nm1-1)));
nm2 = size(GB,1); C2 = GB / ((sum(GB(:))-sum(diag(GB)))/(nm2*(nm2-1)));
Trow = sum(T, 2);  Tcol = sum(T, 1).';
L = ((C1.^2) * Trow) * ones(1, m) + ones(m, 1) * ((C2.^2) * Tcol).' - 2 * (C1 * T * C2);
local = sum(T .* L, 2);

% geodesic quantile profiles on the MaxMin-subsampled supports
% use full matrices so quantiles reflect the whole geodesic distribution
qp = [0.1 0.25 0.5 0.75 0.9];
ga_vec = GA_full(triu(true(size(GA_full)),1));
gb_vec = GB_full(triu(true(size(GB_full)),1));
gq.q          = qp;
gq.spont      = quantile(ga_vec, qp);
gq.stim       = quantile(gb_vec, qp);
gq.ratio      = gq.stim ./ gq.spont;   % >1: stim geodesics longer at this percentile
gq.spont_mean = mean(ga_vec);
gq.stim_mean  = mean(gb_vec);

res.shape_gw = shape_gw;  res.gain = gain;  res.z = (shape_gw - mu) / sd;
res.pct = mean(nullv >= shape_gw);  res.null_mean = mu;  res.null_sd = sd;
res.stim_reliability = stim_reliab;  res.T = T;  res.local_distortion = local;
res.support_idx = sa;  res.geodesic_quantiles = gq;
end

% ---------------------------------------------------------------------------
function G = geodesics(X, k, metric)
n = size(X, 1);
D = squareform(pdist(X, metric));
[~, nb] = sort(D, 2);
src = repmat((1:n).', k, 1);
dst = reshape(nb(:, 2:k+1), [], 1);
wt  = D(sub2ind([n n], src, dst));
A   = max(sparse(src, dst, wt, n, n), sparse(dst, src, wt, n, n));
Gr  = graph(A);
comp = conncomp(Gr);
keep = find(comp == mode(comp));
G = distances(subgraph(Gr, keep));
end

% ---------------------------------------------------------------------------
function idx = maxmin(G, m)
n = size(G, 1);
if m >= n, idx = (1:n).'; return; end
idx = zeros(m, 1);  idx(1) = randi(n);  dmin = G(idx(1), :);
for i = 2:m
    [~, j] = max(dmin);  idx(i) = j;  dmin = min(dmin, G(j, :));
end
end

% ---------------------------------------------------------------------------
function [gw, gain, T] = shapeGW(Ga, Gb, epsilon, gwIters, skIters)
n_ = size(Ga,1); na = (sum(Ga(:)) - sum(diag(Ga))) / (n_*(n_-1));
m_ = size(Gb,1); nb = (sum(Gb(:)) - sum(diag(Gb))) / (m_*(m_-1));
C1 = Ga / na;  C2 = Gb / nb;
n = size(C1, 1);  m = size(C2, 1);
p = ones(n, 1) / n;  q = ones(m, 1) / m;
constC = ((C1.^2) * p) * ones(1, m) + ones(n, 1) * (q.' * (C2.^2));
T = p * q.';
for it = 1:gwIters
    Ggrad = 2 * (constC - 2 * (C1 * T * C2));
    K = exp(-(Ggrad - min(Ggrad(:))) / epsilon);
    u = ones(n, 1);
    for s = 1:skIters
        v = q ./ (K.' * u + 1e-300);
        u = p ./ (K  * v + 1e-300);
    end
    T = (u .* K) .* v.';
end
gw   = sqrt(max(sum(sum((constC - 2 * (C1 * T * C2)) .* T)), 0));
gain = nb / na;
end

% ---------------------------------------------------------------------------
function C = condComponents(X, k, metric)
n = size(X, 1);
D = squareform(pdist(X, metric));
[Ds, nb] = sort(D, 2);
src = repmat((1:n).', k, 1);
dst = reshape(nb(:, 2:k+1), [], 1);
wt  = D(sub2ind([n n], src, dst));
A   = max(sparse(src, dst, wt, n, n), sparse(dst, src, wt, n, n));
lab = conncomp(graph(A)).';
nc  = max(lab);
sz  = zeros(nc, 1); cent = zeros(nc, size(X, 2));
for c = 1:nc, m = lab == c; sz(c) = nnz(m); cent(c, :) = mean(X(m, :), 1); end
[sz, ord] = sort(sz, 'descend');
remap = zeros(nc, 1); remap(ord) = 1:nc;
C.label = remap(lab);  C.size = sz;  C.centroid = cent(ord, :);
C.ncomp = nc;  C.kthDist = Ds(:, k+1);
end

% ---------------------------------------------------------------------------
function M = safeSquareform(v)
if isempty(v), M = 0; else, M = squareform(v); end
end