%% reactivationPrecision.m
% Variance of the early ConA1 / ConB1 reactivation estimate across sessions, from the existing
% compareReactivation.mat (T.pairVec); no session is loaded and nothing is recomputed upstream.
%
% PER SESSION (pairs = [cT cPre cPost i j], i, j = units)
%   dEV   EV - REV, as in compareReactivation (difference of squared partial correlations)
%   sEV   signed alternative: atanh(r_T,POST|PRE) - atanh(r_T,PRE|POST) (Fisher z of the partial
%         correlations, sign kept, no squaring)
%   SE    unit-level jackknife (leave one unit and all its pairs out): pairs that share a unit are
%         not independent, so a pair bootstrap understates the error
% GROUP ESTIMATES (Control, SAL, ASV; ConB1, ConA1 and B1 - A1)
%   unweighted   mean over sessions (the current analysis)
%   RE           random-effects mean (DerSimonian-Laird): weights 1/(SE^2 + tau^2); tau^2 = the
%                between-session variance beyond each session's own error, so a session with
%                thousands of pairs does not take all the weight
%   minPairs     unweighted mean over sessions with >= MIN_PAIRS pairs
% BETWEEN GROUPS (SAL vs ASV)
%   unweighted: permutation over sessions;  RE: permutation of group labels with tau^2 re-estimated
%   each time;  weighted LME y ~ Group + (1|Animal), Weights = 1/(SE^2 + tau^2)

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
GROUPS = {'Control','SAL','ASV'};
MIN_PAIRS = 45;                 % minPairs estimate (~10 units)
NPERM = 1e4;

A = load(fullfile(CACHE_ROOT, 'compareReactivation.mat'), 'T'); T = A.T;
T = T(ismember(T.Group, GROUPS) & ismember(T.Context, {'ConA1','ConB1'}), :);
n = height(T);
[dEV, sEV, seD, seS, nU] = deal(nan(n, 1));
for r = 1:n
    PV = T.pairVec{r};
    if size(PV, 1) < 10, continue; end
    [dEV(r), sEV(r)] = local_stats(PV);
    uu = unique(PV(:, 4:5)); nU(r) = numel(uu);
    jd = nan(nU(r), 1); js = jd;
    for q = 1:nU(r)
        [jd(q), js(q)] = local_stats(PV(PV(:,4) ~= uu(q) & PV(:,5) ~= uu(q), :));
    end
    seD(r) = local_jkSE(jd); seS(r) = local_jkSE(js);
end
T.dEVc = dEV; T.sEV = sEV; T.seD = seD; T.seS = seS; T.nU = nU;

%% ---- per session ---------------------------------------------------------------------------
fprintf('per session (jackknife SE over units)\n');
fprintf('%-8s %-6s %-22s %-5s %6s %5s | %8s %7s | %8s %7s\n', 'group', 'animal', 'session', 'ctx', 'pairs', 'units', ...
    'dEV', 'SE', 'sEV', 'SE');
[~, o] = sortrows(T(:, {'Group','Context','Animal'}));
for r = o(:)'
    [~, sn] = fileparts(fileparts(T.Session{r}));
    fprintf('%-8s %-6s %-22s %-5s %6d %5d | %+8.3f %7.3f | %+8.3f %7.3f\n', T.Group{r}, T.Animal{r}, sn, T.Context{r}, ...
        size(T.pairVec{r}, 1), T.nU(r), T.dEVc(r), T.seD(r), T.sEV(r), T.seS(r));
end

%% ---- B1 - A1 within session -------------------------------------------------------------------
% SE of the difference: sqrt(SE_B^2 + SE_A^2) (the two templates' errors are treated as independent)
kb = find(strcmp(T.Context, 'ConB1'));
D = table();
for r = kb(:)'
    a = find(strcmp(T.Session, T.Session{r}) & strcmp(T.Context, 'ConA1'), 1);
    if isempty(a), continue; end
    D = [D; table(T.Group(r), T.Animal(r), T.Session(r), {'B1-A1'}, T.dEVc(r) - T.dEVc(a), T.sEV(r) - T.sEV(a), ...
        hypot(T.seD(r), T.seD(a)), hypot(T.seS(r), T.seS(a)), min(size(T.pairVec{r}, 1), size(T.pairVec{a}, 1)), ...
        'VariableNames', {'Group','Animal','Session','Context','dEVc','sEV','seD','seS','nPairs'})]; %#ok<AGROW>
end
S = T(:, {'Group','Animal','Session','Context','dEVc','sEV','seD','seS'});
S.nPairs = cellfun(@(x) size(x, 1), T.pairVec);
S = [S; D];

%% ---- group estimates and SAL vs ASV ---------------------------------------------------------
ms = {'dEVc','seD','dEV'; 'sEV','seS','sEV'};
for c = {'ConB1','ConA1','B1-A1'}
    for mm = 1:2
        y_ = ms{mm,1}; s_ = ms{mm,2};
        t = S(strcmp(S.Context, c{1}) & isfinite(S.(y_)) & isfinite(S.(s_)) & S.(s_) > 0, :);
        if isempty(t), continue; end
        fprintf('\n%s, %s\n', c{1}, ms{mm,3});
        fprintf('   %-8s | %-30s | %-38s | %-30s\n', 'group', 'unweighted mean +- SEM (n)', 'RE mean +- SE (tau), z p', ...
            sprintf('>= %d pairs: mean +- SEM (n)', MIN_PAIRS));
        for g = GROUPS
            k = strcmp(t.Group, g{1});
            if ~any(k), continue; end
            y = t.(y_)(k); se = t.(s_)(k);
            [mu, smu, tau] = local_re(y, se);
            kp = t.nPairs(k) >= MIN_PAIRS;
            fprintf('   %-8s | %+.3f +- %.3f (%2d)%10s | %+.3f +- %.3f (tau %.3f) p=%.3f%6s | %+.3f +- %.3f (%d)\n', g{1}, ...
                mean(y), std(y)/sqrt(numel(y)), numel(y), '', mu, smu, tau, 2*normcdf(-abs(mu/smu)), '', ...
                mean(y(kp)), std(y(kp))/sqrt(sum(kp)), sum(kp));
        end
        a = strcmp(t.Group, 'SAL'); b = strcmp(t.Group, 'ASV');
        if sum(a) < 2 || sum(b) < 2, continue; end
        ya = t.(y_)(a); yb = t.(y_)(b); sa = t.(s_)(a); sb = t.(s_)(b);
        kp = t.nPairs >= MIN_PAIRS;
        fprintf('   SAL vs ASV: unweighted %+.3f perm %.3f | RE %+.3f perm %.3f | weighted LME %.3f | >= %d pairs %+.3f perm %.3f\n', ...
            mean(yb) - mean(ya), local_perm(ya, yb, NPERM), ...
            local_reDiff(ya, sa, yb, sb), local_rePerm(ya, sa, yb, sb, NPERM), local_wlme(t(a | b, :), y_, s_), ...
            MIN_PAIRS, mean(t.(y_)(b & kp)) - mean(t.(y_)(a & kp)), local_perm(t.(y_)(a & kp), t.(y_)(b & kp), NPERM));
    end
end
fprintf(['\n(differences are ASV minus SAL; RE p within a group is a normal approximation; the RE permutation\n' ...
         ' keeps each session''s SE with its value; sessions from one animal are treated as independent in RE)\n']);
save(fullfile(CACHE_ROOT, 'reactivationPrecision.mat'), 'S');


%% =========================================================================
function [d, s] = local_stats(PV)
d = NaN; s = NaN;
if size(PV, 1) < 5, return; end
r = corr(PV(:, 1:3)); a = r(1,2); b = r(1,3); c = r(2,3);
pPost = (b - a*c) / sqrt((1 - a^2) * (1 - c^2));
pPre  = (a - b*c) / sqrt((1 - b^2) * (1 - c^2));
d = pPost^2 - pPre^2;
s = atanh(max(min(pPost, 0.999), -0.999)) - atanh(max(min(pPre, 0.999), -0.999));
end


function se = local_jkSE(th)
th = th(isfinite(th)); m = numel(th); se = NaN;
if m < 3, return; end
se = sqrt((m - 1) / m * sum((th - mean(th)).^2));
end


function [mu, se, tau] = local_re(y, s)
% DerSimonian-Laird random-effects mean
y = y(:); s = s(:); w = 1 ./ s.^2;
k = numel(y); mu = NaN; se = NaN; tau = NaN;
if k < 2, return; end
mfe = sum(w .* y) / sum(w);
Q = sum(w .* (y - mfe).^2);
C = sum(w) - sum(w.^2) / sum(w);
tau2 = max(0, (Q - (k - 1)) / C); tau = sqrt(tau2);
ws = 1 ./ (s.^2 + tau2);
mu = sum(ws .* y) / sum(ws); se = sqrt(1 / sum(ws));
end


function d = local_reDiff(ya, sa, yb, sb)
d = local_re(yb, sb) - local_re(ya, sa);
end


function p = local_rePerm(ya, sa, yb, sb, nperm)
y = [ya(:); yb(:)]; s = [sa(:); sb(:)]; na = numel(ya);
obs = abs(local_reDiff(ya, sa, yb, sb));
nul = zeros(nperm, 1);
for k = 1:nperm
    o = randperm(numel(y));
    nul(k) = abs(local_reDiff(y(o(1:na)), s(o(1:na)), y(o(na+1:end)), s(o(na+1:end))));
end
p = (1 + sum(nul >= obs - 1e-12)) / (1 + nperm);
end


function p = local_perm(x, y, nperm)
x = x(:); y = y(:); p = NaN;
if numel(x) < 2 || numel(y) < 2, return; end
z = [x; y]; na = numel(x); obs = abs(mean(y) - mean(x));
nul = arrayfun(@(k) local_pdiff(z(randperm(numel(z))), na), 1:nperm);
p = (1 + sum(nul >= obs - 1e-12)) / (1 + nperm);
end


function d = local_pdiff(z, na)
d = abs(mean(z(na+1:end)) - mean(z(1:na)));
end


function p = local_wlme(t, y_, s_)
% y ~ Group + (1|Animal), observation weights 1/(SE^2 + tau^2), tau^2 pooled over both groups
p = NaN;
try
    [~, ~, tau] = local_re(t.(y_), t.(s_));
    t.y = t.(y_); t.w = 1 ./ (t.(s_).^2 + tau^2);
    t.Group = categorical(t.Group); t.Animal = categorical(t.Animal);
    lme = fitlme(t, 'y ~ Group + (1|Animal)', 'Weights', t.w);
    [~, ~, C] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
    if isa(C, 'dataset'), C = dataset2table(C); end
    k = find(startsWith(string(C.Name), 'Group'), 1);
    if ~isempty(k), p = C.pValue(k); end
catch ME
    warning('prec:lme', '%s', ME.message);
end
end
