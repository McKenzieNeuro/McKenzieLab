function STRESS = stressRandomSlope(varargin)
%STRESSRANDOMSLOPE  Is the random-slope LRT a calibrated test at 7-11 animals?
%
%   STRESS = stressRandomSlope('Widths',[10 60], 'NPerm',200, 'NSim',100, ...)
%
%   Compares four tests of the same contrast, each run through the pipeline's
%   own code (fitHazardGLMM):
%     plain  intercept-only LRT          rs   random-slope LRT (what is reported)
%     CR     Poisson GLM, animal-cluster-robust Wald, t(G-1)
%     anim   exact animal-level sign-flip on per-animal log rate ratios
%              (DiD: paired, animals with both arms)
%
%   PART 1 -- RANDOMIZATION NULL ON THE REAL TABLES (acute arms and DiD).
%     Each animal's pre/post labels are swapped (all its sessions together)
%     with probability 1/2; for the DiD the unit is animal x drug. Exposures,
%     censoring, sparsity and each animal's real MAGNITUDE of change are kept;
%     the sign is randomized, so the mean change is zero by symmetry -- the
%     null the sign-flip test assumes. Rejection rate at 0.05 over null worlds
%     is each method's true type I error on THIS design. "p cal" ranks the
%     observed rs p among the null-world rs p's: a calibrated p for the same
%     statistic. If the real data hold an effect it adds heterogeneity to
%     the null worlds, so this is conservative; the SAL arm is a real null.
%     Persistence is not randomized (Phase is collinear with days in study,
%     and flipping labels would break that); Part 2 covers its design.
%
%   PART 2 -- PARAMETRIC SWEEP ON THE REAL DESIGN (acute ASV arm, DiD, persistence).
%     Outcomes regenerated from a known model on the real animals, sessions,
%     exposures and phase labels: log hazard = baseline + animal + session
%     intercepts + (beta + tau*z_animal)*post [+ drift*days, persistence].
%     Each interval's observed length is its censoring limit and its time to
%     event is exponential, so event counts run slightly below the real
%     data's. Baseline and the animal/session SDs are taken from the fitted
%     model. Sweeps tau (heterogeneity of the effect), beta (effect), and
%     drift (persistence only: slow change with NO phase effect).
%       beta = 0 rows: rejection rate = type I error.  beta ~= 0: power.
%       tau-hat < 1e-3 = slope variance collapsed to the boundary (rs == plain).
%
%   Fit times are unknown until the first fit; the function times one and
%   prints an estimate. Total fits ~ 4*(NPerm*cells1 + NSim*conds*cells2).
%   Uses parfor; works serially without a pool.
%
%   Needs timescaleSweep_N<MatchN>.mat. If it was never written, from the
%   workspace that holds SWEEP:
%     save(fullfile(CACHE_ROOT,'timescaleSweep_N10.mat'),'SWEEP','SCALES')

p = inputParser;
p.addParameter('CacheRoot', 'C:\Users\samckenzie\syncCache');
p.addParameter('MatchN', 10);
p.addParameter('Widths', [10 60]);
p.addParameter('Bout', 'null95');
p.addParameter('Part', [1 2]);
p.addParameter('NPerm', 200);
p.addParameter('NSim', 100);
p.addParameter('Conds', [0 0 0; 0.5 0 0; 1 0 0; 0.5 log(1.5) 0; 0.5 0 log(1.04)]);  % [tau beta drift/day]
p.addParameter('Part1', {'ASV','entry'; 'ASV','exit'; 'SAL','entry'; 'SAL','exit'; 'DiD','entry'});
p.addParameter('Part2', {'ASV','entry'; 'DiD','entry'; 'persist','entry'; 'persist','exit'});
p.addParameter('Seed', 1);
p.parse(varargin{:}); o = p.Results;
rng(o.Seed);

f = fullfile(o.CacheRoot, sprintf('timescaleSweep_N%d.mat', o.MatchN));
if ~isfile(f)
    error('stressRandomSlope:noSweep', ['%s not found. From the workspace holding SWEEP: ' ...
        'save(''%s'',''SWEEP'',''SCALES'')'], f, f);
end
L = load(f, 'SWEEP', 'SCALES');
[found, wi] = ismember(o.Widths, L.SCALES(:,1));
if ~all(found), error('stressRandomSlope:width', 'Widths must be among %s', mat2str(L.SCALES(:,1)')); end

STRESS = struct('part1', {{}}, 'part2', {{}}, 'opts', o);
tag = matlab.lang.makeValidName(o.Bout);

%% ---- timing --------------------------------------------------------------
[t0, s0] = local_getCell(L.SWEEP{wi(1)}, tag, 'ASV', 'entry');
tic; local_pvals(t0, s0); tfit = toc;
n1 = numel(o.Widths) * size(o.Part1,1) * o.NPerm;
n2 = numel(o.Widths) * size(o.Part2,1) * size(o.Conds,1) * o.NSim;
fprintf('one pass of all four tests: %.1f s. Part 1: %d passes, Part 2: %d passes -> ~%.1f h serial (divide by workers).\n', ...
    tfit, n1, n2, tfit * (n1*any(o.Part==1) + n2*any(o.Part==2)) / 3600);

%% ---- Part 1 ----------------------------------------------------------------
if any(o.Part == 1)
    fprintf('\n=== PART 1: animal-level label-flip null on the real tables (%d null worlds per cell) ===\n', o.NPerm);
    fprintf('%-5s %-6s %-5s | %-6s | real p: %-7s %-7s %-7s %-7s | %-6s | null-world rejection at 0.05: %-5s %-5s %-5s %-5s | collapsed\n', ...
        'width','ctx','haz','G','plain','rs','CR','anim','p cal','plain','rs','CR','anim');
    for w = 1:numel(o.Widths)
        for c = 1:size(o.Part1, 1)
            [t, spec, ~] = local_getCell(L.SWEEP{wi(w)}, tag, o.Part1{c,1}, o.Part1{c,2});
            Rr = local_part1(t, spec, o);
            Rr.width = o.Widths(w); Rr.ctx = o.Part1{c,1}; Rr.hz = o.Part1{c,2};
            STRESS.part1{end+1} = Rr; %#ok<AGROW>
            fprintf('%-5g %-6s %-5s | %-6d | %15.4f %7.4f %7.4f %7.4f | %6.3f | %29.0f%% %4.0f%% %4.0f%% %4.0f%% | %.0f%%\n', ...
                Rr.width, Rr.ctx, Rr.hz, Rr.G, Rr.real.pPlain, Rr.real.pRS, Rr.real.pCR, Rr.real.pAnim, ...
                Rr.pCal, 100*Rr.rej(1), 100*Rr.rej(2), 100*Rr.rej(3), 100*Rr.rej(4), 100*Rr.collapsed);
        end
    end
    fprintf(['A test is calibrated if its null-world rejection is near 5%% (Monte Carlo SE %.1f%% at NPerm = %d).\n' ...
             'Above ~8%% is liberal; the sign-flip is exact by construction and should sit near 5%%.\n'], ...
        100*sqrt(0.05*0.95/o.NPerm), o.NPerm);
end

%% ---- Part 2 ----------------------------------------------------------------
if any(o.Part == 2)
    fprintf('\n=== PART 2: parametric sweep on the real design (%d simulations per condition) ===\n', o.NSim);
    for w = 1:numel(o.Widths)
        for c = 1:size(o.Part2, 1)
            [t, spec, F0] = local_getCell(L.SWEEP{wi(w)}, tag, o.Part2{c,1}, o.Part2{c,2});
            fprintf('\n-- width %g s, %s %s: G = %d, observed tau-hat in the real rs fit: %.2f --\n', ...
                o.Widths(w), o.Part2{c,1}, o.Part2{c,2}, local_G(t, spec.grp), local_slopeSD(F0));
            fprintf('   %-4s %-6s %-9s | rejection at 0.05: %-5s %-5s %-5s %-5s | rs fit ok  tau-hat<1e-3  mean tau-hat\n', ...
                'tau', 'beta', 'drift/d', 'plain', 'rs', 'CR', 'anim');
            Rr = local_part2(t, spec, F0, o);
            Rr.width = o.Widths(w); Rr.ctx = o.Part2{c,1}; Rr.hz = o.Part2{c,2};
            STRESS.part2{end+1} = Rr; %#ok<AGROW>
            for q = 1:size(o.Conds, 1)
                if ~strcmp(spec.type, 'persist') && o.Conds(q,3) ~= 0, continue; end
                fprintf('   %-4.2f %+6.2f %+9.3f | %17.0f%% %4.0f%% %4.0f%% %4.0f%% | %9.0f%% %11.0f%% %12.2f\n', ...
                    o.Conds(q,1), o.Conds(q,2), o.Conds(q,3), 100*Rr.rej(q,1), 100*Rr.rej(q,2), ...
                    100*Rr.rej(q,3), 100*Rr.rej(q,4), 100*Rr.rsOK(q), 100*Rr.collapsed(q), Rr.meanTau(q));
            end
        end
    end
    fprintf(['\nbeta = 0 rows are type I error (Monte Carlo SE %.1f%% at NSim = %d); beta ~= 0 rows are power.\n' ...
             'A drift row with beta = 0 that rejects is slow change being mistaken for a phase effect.\n'], ...
        100*sqrt(0.05*0.95/o.NSim), o.NSim);
end

save(fullfile(o.CacheRoot, 'stressRandomSlope.mat'), 'STRESS');
end


%% =========================================================================
function [t, spec, F0] = local_getCell(S, tag, ctx, hz)
if strcmp(ctx, 'persist')
    Q = S.persist.(tag); t = Q.([hz 'Tbl']); F0 = Q.(hz);
    spec = struct('type','persist', 'grp','Phase', 'fixed','Phase + Days', 'term','Phase', ...
        'coef','Phase_post', 'slope','Phase', 'group','Animal');
else
    Q = S.acute.(tag); t = Q.([hz 'Tbl']);
    if strcmp(ctx, 'DiD')
        F0 = Q.(hz);
        spec = struct('type','did', 'grp','Period', 'fixed','Drug*Period', 'term','Drug:Period', ...
            'coef','Drug_ASV:Period_post', 'slope','Period', 'group','AnimalDrug');
    else
        t = t(t.Drug == ctx, :); F0 = Q.arm.(ctx).(hz);
        spec = struct('type','arm', 'grp','Period', 'fixed','Period', 'term','Period', ...
            'coef','Period_post', 'slope','Period', 'group','Animal');
    end
end
end


function R = local_part1(t, spec, o)
real = local_pvals(t, spec);
if strcmp(spec.group, 'AnimalDrug'), key = strcat(cellstr(t.Animal), '_', cellstr(t.Drug));
else, key = cellstr(t.Animal); end
[~, ~, gU] = unique(key); nU = max(gU);
FL = rand(nU, o.NPerm) < 0.5;
out = nan(o.NPerm, 5);
parfor s = 1:o.NPerm
    r = local_pvals(local_flip(t, spec, gU, FL(:,s)), spec);
    out(s,:) = [r.pPlain r.pRS r.pCR r.pAnim r.tau];
end
R = struct('real', real, 'G', local_G(t, spec.grp), 'null', out);
R.rej = [mean(out(isfinite(out(:,1)),1) < 0.05), mean(out(isfinite(out(:,2)),2) < 0.05), ...
         mean(out(isfinite(out(:,3)),3) < 0.05), mean(out(isfinite(out(:,4)),4) < 0.05)];
fin = isfinite(out(:,2));
R.pCal = (1 + sum(out(fin,2) <= real.pRS)) / (1 + sum(fin));
R.collapsed = mean(~fin | out(:,5) < 1e-3);
end


function t2 = local_flip(t, spec, gU, f)
x = t.(spec.grp);
fl = f(gU); fl = fl(:);
newPost = ((x == 'post') & ~fl) | ((x == 'pre') & fl);
t2 = t;
t2.(spec.grp) = categorical(double(newPost), [0 1], {'pre','post'});
end


function R = local_part2(t, spec, F0, o)
n = height(t);
lam0 = sum(t.Event) / sum(t.Exposure);
[sdA, sdS] = local_reSD(F0);
gA = local_idx(cellstr(t.Animal)); gS = local_idx(cellstr(t.Session));
if strcmp(spec.group, 'AnimalDrug'), gB = local_idx(strcat(cellstr(t.Animal), '_', cellstr(t.Drug)));
else, gB = gA; end
% one set of random draws shared by every condition (common random numbers)
U = rand(n, o.NSim); ZA = randn(max(gA), o.NSim); ZS = randn(max(gS), o.NSim); ZB = randn(max(gB), o.NSim);
nq = size(o.Conds, 1);
R = struct('rej', nan(nq,4), 'rsOK', nan(nq,1), 'collapsed', nan(nq,1), 'meanTau', nan(nq,1), 'p', {cell(nq,1)});
for q = 1:nq
    tau = o.Conds(q,1); beta = o.Conds(q,2); drift = o.Conds(q,3);
    if ~strcmp(spec.type, 'persist') && drift ~= 0, continue; end
    out = nan(o.NSim, 5);
    parfor s = 1:o.NSim
        t2 = local_simulate(t, spec, lam0, sdA, sdS, tau, beta, drift, gA, gS, gB, ...
            U(:,s), ZA(:,s), ZS(:,s), ZB(:,s));
        r = local_pvals(t2, spec);
        out(s,:) = [r.pPlain r.pRS r.pCR r.pAnim r.tau];
    end
    for m = 1:4, R.rej(q,m) = mean(out(isfinite(out(:,m)), m) < 0.05); end
    R.rsOK(q) = mean(isfinite(out(:,2)));
    R.collapsed(q) = mean(isfinite(out(:,2)) & out(:,5) < 1e-3);
    R.meanTau(q) = mean(out(isfinite(out(:,5)), 5));
    R.p{q} = out;
end
end


function t2 = local_simulate(t, spec, lam0, sdA, sdS, tau, beta, drift, gA, gS, gB, u, za, zs, zb)
post = double(t.(spec.grp) == 'post');
lp = log(lam0) + sdA*za(gA) + sdS*zs(gS);
switch spec.type
    case 'did'
        lp = lp + (beta*double(t.Drug == 'ASV') + tau*zb(gB)) .* post;
    case 'persist'
        lp = lp + (beta + tau*zb(gB)) .* post + drift * (t.Days - mean(t.Days));
    otherwise
        lp = lp + (beta + tau*zb(gB)) .* post;
end
T = -log(u) ./ exp(lp);
E0 = t.Exposure;
t2 = t;
t2.Event = double(T <= E0);
t2.Exposure = max(min(T, E0), 1e-6);
end


function R = local_pvals(t, spec)
w = warning('off', 'all'); cleanup = onCleanup(@() warning(w)); %#ok<NASGU>
R = struct('pPlain',NaN, 'pRS',NaN, 'pCR',NaN, 'pAnim',NaN, 'tau',NaN, 'hr',NaN);
try
    F = fitHazardGLMM(t, spec.fixed, spec.term, spec.coef, spec.slope, spec.group);
    if F.ok, R.pPlain = F.pLRT; end
    if F.rs.ok, R.pRS = F.rs.pLRT; R.hr = F.rs.HR; R.tau = local_slopeSD(F); end
catch
end
R.pCR = local_crPoisson(t, spec);
R.pAnim = local_animalP(t, spec);
end


function p = local_crPoisson(t, spec)
% Poisson regression with offset, animal-cluster-robust (CR1) Wald, t(G-1)
p = NaN;
try
    post = double(t.(spec.grp) == 'post'); n = height(t); one = ones(n,1);
    switch spec.type
        case 'did',     asv = double(t.Drug == 'ASV'); X = [one, asv, post, asv.*post]; j = 4;
        case 'persist', X = [one, post, t.Days]; j = 2;
        otherwise,      X = [one, post]; j = 2;
    end
    y = t.Event; off = log(t.Exposure);
    b = glmfit(X, y, 'poisson', 'offset', off, 'constant', 'off');
    mu = exp(X*b + off);
    A = X' * (X .* mu);
    r = y - mu;
    [~, ~, gi] = unique(cellstr(t.Animal)); G = max(gi);
    if G < 3, return; end
    S = zeros(size(X,2));
    for i = 1:G
        s = X(gi==i,:)' * r(gi==i);
        S = S + s*s';
    end
    k = size(X,2);
    V = (A \ S / A) * (G/(G-1)) * ((n-1)/(n-k));
    p = 2 * tcdf(-abs(b(j)/sqrt(V(j,j))), G-1);
catch
end
end


function p = local_animalP(t, spec)
if strcmp(spec.type, 'did')
    [dA, iA] = local_animalLR(t(t.Drug == 'ASV', :), spec.grp);
    [dS, iS] = local_animalLR(t(t.Drug == 'SAL', :), spec.grp);
    [~, ia, is] = intersect(iA, iS);
    d = dA(ia) - dS(is);
else
    d = local_animalLR(t, spec.grp);
end
p = local_signflip(d);
end


function [d, ids] = local_animalLR(t, grp)
% per-animal log rate ratio post/pre, +0.5 correction (as runHazardAnalysis);
% animals with no events in either condition carry no information and are skipped
ua = unique(t.Animal); d = []; ids = {};
for i = 1:numel(ua)
    m = t.Animal == ua(i);
    ma = m & t.(grp) == 'pre'; mb = m & t.(grp) == 'post';
    if ~any(ma) || ~any(mb), continue; end
    if sum(t.Event(ma)) + sum(t.Event(mb)) == 0, continue; end
    ra = (sum(t.Event(ma)) + 0.5) / sum(t.Exposure(ma));
    rb = (sum(t.Event(mb)) + 0.5) / sum(t.Exposure(mb));
    d(end+1,1) = log(rb/ra); %#ok<AGROW>
    ids{end+1,1} = char(ua(i)); %#ok<AGROW>
end
end


function p = local_signflip(d)
d = d(isfinite(d)); n = numel(d);
if n < 2, p = NaN; return; end
obs = abs(mean(d));
if n <= 14
    S = 2*(dec2bin(0:2^n-1, n) == '1') - 1;
    p = mean(abs(S*d(:)/n) >= obs - 1e-12);
else
    nperm = 1e4; nul = zeros(nperm,1);
    for k = 1:nperm, nul(k) = mean((2*(rand(n,1)>0.5)-1).*d); end
    p = (1 + sum(abs(nul) >= obs)) / (1 + nperm);
end
end


function G = local_G(tb, grp)
ua = unique(tb.Animal); G = 0;
for i = 1:numel(ua)
    m = tb.Animal == ua(i);
    if any(tb.(grp)(m) == 'pre') && any(tb.(grp)(m) == 'post'), G = G + 1; end
end
end


function sd = local_slopeSD(F)
% SD of the animal-specific slope in a stored or fresh rs fit; NaN if unavailable
sd = NaN;
if ~isfield(F, 'rs') || ~F.rs.ok || ~isfield(F.rs, 'model') || isempty(F.rs.model), return; end
try
    psi = covarianceParameters(F.rs.model);
    for i = 1:numel(psi)
        if size(psi{i}, 1) == 2, sd = sqrt(max(psi{i}(2,2), 0)); return; end
    end
    sd = 0;      % no 2x2 block: slope term dropped by the fit
catch
end
end


function [sdA, sdS] = local_reSD(F)
% animal and session intercept SDs from the fitted intercept-only model;
% defaults if unavailable (printed by the caller via the tau line only)
sdA = 0.5; sdS = 0.3;
try
    psi = covarianceParameters(F.model);
    gv = F.model.Formula.GroupingVariableNames;
    for i = 1:numel(gv)
        v = sqrt(max(psi{i}(1,1), 0)); g = char(gv{i});
        if strcmp(g, 'Animal'), sdA = v; elseif strcmp(g, 'Session'), sdS = v; end
    end
catch
end
end


function g = local_idx(c)
[~, ~, g] = unique(c);
end
