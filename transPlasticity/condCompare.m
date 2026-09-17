function out = condCompare(T, x, rootName, label, varargin)
%CONDCOMPARE  Compare conditions within animal (Baseline/ASV/ASV24/SAL/SAL24).
%
%   out = condCompare(T, x, rootName, label, 'Name', value, ...)
%
%   EXAMPLE (five conditions, Baseline as reference):
%     C = condCompare(a, log(dur1), 'TransPlasticity', 'log duration', ...
%           'Conditions', {'Baseline','ASV','ASV24','SAL','SAL24'});
%
%   T     : session struct array from syncBatch (e.g. `a`), .file must parse
%   x     : per-session metric, same length as T (NaNs allowed and dropped)
%   rootName : 'TransPlasticity'
%   label : name of the metric, for printing
%
%   NAME/VALUE
%     'Conditions' : cell array of condition names to KEEP, in the order you
%                    want them printed and contrasted. Everything else is
%                    dropped before any model is fit -- so a 'Baseline' folder
%                    (or any other stray directory that happens to sit at the
%                    condition level) never enters the omnibus test, the
%                    contrast count, or the Holm correction. Default: every
%                    condition found.
%     'Factorial'  : if true, ALSO fit Drug x Day (ASV/SAL x same-day/next-day).
%                    Default false.
%     'Contrasts'  : cell array of 2-element cells naming the ONLY pairwise
%                    comparisons to run, e.g.
%                      {{'ASV','SAL'}, {'ASV','ASV24'}}
%                    Holm then corrects across just those, which is a real
%                    power gain over correcting across all K(K-1)/2 -- but
%                    only legitimate if the contrasts were chosen BEFORE
%                    seeing the results. Default: all pairs.
%
%   DESIGN. Three levels, not two: animal -> condition -> repeated sessions.
%   Sessions in the same animal-condition cell (e.g. CP19\SAL24\260415 and
%   CP19\SAL24\260417) are more alike than sessions across cells, so the
%   random structure needs BOTH (1|Animal) and (1|Animal:Condition).
%   Omitting the nested term treats repeat sessions within a cell as
%   independent replicates of that condition, which inflates the condition
%   effect's precision -- the same pseudo-replication that made session-level
%   p-values untrustworthy in the between-dataset analysis, one level down.
%
%   MODEL:  Metric ~ Condition + (1|Animal) + (1|Animal:Condition)
%
%   If the nested variance collapses to zero (common when most cells hold a
%   single session), the function refits without it and says so -- a
%   degenerate term is worse than an absent one, because it looks like it
%   corrected for something it did not.
%
%   UNBALANCED CELLS. Animals missing a condition are NOT dropped. That is
%   the main reason to use a mixed model here rather than repeated-measures
%   ANOVA, which would require complete cases and would discard most of your
%   animals. But unbalance is not free: the omnibus test borrows strength
%   across animals, so read it alongside the printed design table rather
%   than on its own. A "significant" condition effect resting on two animals
%   that happen to have the rare condition is not a within-animal result.
%
%   TWO TESTS PER CONTRAST, same logic as nestingCheck/lmeCompare:
%     LME contrast   : uses all sessions, better powered, asymptotic.
%     Paired permutation : collapse each animal to its per-condition mean,
%                       take the within-animal difference for animals having
%                       BOTH conditions, sign-flip permute. Exact, makes no
%                       distributional assumption, uses only genuinely paired
%                       information, and reports how many animals it had.
%                       This is the number to trust when they disagree.
%
%   Pairwise p-values are Holm-corrected across the 6 contrasts.
%
%   See also ANIMALFROMPATH, CONDITIONFROMPATH, LMECOMPARE, NESTINGCHECK.

ip = inputParser;
ip.addParameter('Conditions', {}, @(c) iscell(c) || isempty(c));
ip.addParameter('Factorial', false, @(v) islogical(v) || isnumeric(v));
ip.addParameter('Contrasts', {}, @(c) iscell(c) || isempty(c));
ip.parse(varargin{:});
keepConds  = ip.Results.Conditions;
doFactorial = logical(ip.Results.Factorial);
wantPairs   = ip.Results.Contrasts;
if nargin < 4, label = 'metric'; end

animal = arrayfun(@(r) animalFromPath(r.file, rootName), T, 'UniformOutput', false);
cond   = arrayfun(@(r) conditionFromPath(r.file, rootName), T, 'UniformOutput', false);
animal = animal(:); cond = cond(:); x = x(:);

bad = cellfun(@isempty, animal) | cellfun(@isempty, cond) | isnan(x);
if any(bad)
    fprintf('%s: dropping %d/%d sessions (unparsed path or NaN metric)\n', ...
        label, sum(bad), numel(x));
end
animal = animal(~bad); cond = cond(~bad); x = x(~bad);

if ~isempty(keepConds)
    found = unique(cond);
    missing = setdiff(keepConds, found);
    if ~isempty(missing)
        warning('condCompare:missingCond', ...
            'requested condition(s) not present in the data: %s', strjoin(missing, ', '));
    end
    dropped = setdiff(found, keepConds);
    sel = ismember(cond, keepConds);
    if ~isempty(dropped)
        fprintf('%s: excluding %d sessions in non-requested condition(s): %s\n', ...
            label, sum(~sel), strjoin(dropped', ', '));
    end
    animal = animal(sel); cond = cond(sel); x = x(sel);
end

if isempty(x), error('condCompare:empty','nothing left after dropping.'); end

if isempty(keepConds)
    uc = unique(cond);
else
    uc = keepConds(ismember(keepConds, unique(cond)));   % keep requested order
    uc = uc(:)';
end
ua = unique(animal);

%% ---- design balance ----------------------------------------------------
fprintf('\n=== condition comparison: %s ===\n', label);
fprintf('%d sessions, %d animals, %d conditions (%s)\n', ...
    numel(x), numel(ua), numel(uc), strjoin(uc', ', '));
fprintf('\nsessions per animal x condition:\n');
w = max(10, max(cellfun(@numel, uc)) + 2);
fmtS = ['%' num2str(w) 's']; fmtD = ['%' num2str(w) 'd'];
fprintf('%-10s', 'animal');
for j = 1:numel(uc), fprintf(fmtS, uc{j}); end
fprintf(fmtS, 'total'); fprintf('\n');
nCell = zeros(numel(ua), numel(uc));
for i = 1:numel(ua)
    fprintf('%-10s', ua{i});
    for j = 1:numel(uc)
        nCell(i,j) = sum(strcmp(animal, ua{i}) & strcmp(cond, uc{j}));
        fprintf(fmtD, nCell(i,j));
    end
    fprintf(fmtD, sum(nCell(i,:))); fprintf('\n');
end
fprintf('%-10s', 'total');
for j = 1:numel(uc), fprintf(fmtD, sum(nCell(:,j))); end
fprintf(fmtD, sum(nCell(:))); fprintf('\n');

nComplete = sum(all(nCell > 0, 2));
fprintf('\n%d/%d animals have all %d conditions; %d have repeat sessions in >=1 cell\n', ...
    nComplete, numel(ua), numel(uc), sum(any(nCell > 1, 2)));
if nComplete < 3
    warning('condCompare:balance', ...
        ['Only %d animals have every condition. The omnibus test will lean on ', ...
         'between-animal comparison for the sparse conditions -- which is exactly ', ...
         'what a within-animal design is supposed to avoid. Read the pairwise ', ...
         'paired-permutation n values, not just the omnibus p.'], nComplete);
end

%% ---- LME ---------------------------------------------------------------
tbl = table(x, categorical(animal), categorical(cond, uc), ...
    'VariableNames', {'Metric','Animal','Condition'});

nested = true;
try
    lme = fitlme(tbl, 'Metric ~ Condition + (1|Animal) + (1|Animal:Condition)', ...
        'FitMethod','REML');
    psi = covarianceParameters(lme);
    if numel(psi) >= 2 && psi{2}(1,1) < 1e-10
        nested = false;
    end
catch
    nested = false;
end
if ~nested
    fprintf(['\nnested (1|Animal:Condition) variance is ~0 or failed to fit; ', ...
             'refitting without it.\n(Expected when most animal-condition cells ', ...
             'hold a single session.)\n']);
    lme = fitlme(tbl, 'Metric ~ Condition + (1|Animal)', 'FitMethod','REML');
end

st = anova(lme, 'DFMethod', 'Satterthwaite');
r = strcmp(st.Term, 'Condition');
fprintf('\nLME omnibus: F(%d,%.1f) = %.3f, p = %.4f  (nested term: %d)\n', ...
    st.DF1(r), st.DF2(r), st.FStat(r), st.pValue(r), nested);

out.lme = lme;
out.omnibusP = st.pValue(r);
out.nested = nested;

%% ---- pairwise: LME contrast + exact paired permutation ------------------
if isempty(wantPairs)
    % Every pairwise comparison. With K conditions this is K(K-1)/2
    % contrasts, and Holm across all of them is conservative when only a few
    % were ever of interest.
    pairs = nchoosek(1:numel(uc), 2);
else
    % Pre-specified contrasts only. This is the honest way to buy power:
    % correcting across 2 planned comparisons instead of 6 exploratory ones
    % is legitimate ONLY if the choice was made before seeing the results.
    % Naming the two that already looked smallest is not pre-specification.
    pairs = nan(numel(wantPairs), 2);
    for k = 1:numel(wantPairs)
        pk = wantPairs{k};
        if ~iscell(pk) || numel(pk) ~= 2
            error('condCompare:contrast', ...
                'Each entry of ''Contrasts'' must be a 2-element cell, e.g. {''ASV'',''SAL''}.');
        end
        i1 = find(strcmpi(uc, pk{1}), 1);
        i2 = find(strcmpi(uc, pk{2}), 1);
        if isempty(i1) || isempty(i2)
            error('condCompare:contrast', ...
                'Contrast {%s,%s}: condition not present after filtering (have: %s).', ...
                pk{1}, pk{2}, strjoin(uc', ', '));
        end
        pairs(k,:) = [i1 i2];
    end
    fprintf('\nusing %d PRE-SPECIFIED contrast(s); Holm corrects across these only.\n', ...
        size(pairs,1));
end
np = size(pairs,1);
pLME = nan(np,1); pPerm = nan(np,1); nPair = nan(np,1); dMean = nan(np,1);

cn = lme.CoefficientNames;
for k = 1:np
    A = uc{pairs(k,1)}; B = uc{pairs(k,2)};

    % --- LME contrast
    H = zeros(1, numel(cn));
    iA = find(strcmp(cn, ['Condition_' A]));
    iB = find(strcmp(cn, ['Condition_' B]));
    if ~isempty(iA), H(iA) = H(iA) + 1; end
    if ~isempty(iB), H(iB) = H(iB) - 1; end
    if any(H ~= 0)
        pLME(k) = coefTest(lme, H, 0, 'DFMethod', 'Satterthwaite');
    end

    % --- exact paired permutation on per-animal means
    d = [];
    for i = 1:numel(ua)
        vA = x(strcmp(animal, ua{i}) & strcmp(cond, A));
        vB = x(strcmp(animal, ua{i}) & strcmp(cond, B));
        if ~isempty(vA) && ~isempty(vB)
            d(end+1,1) = mean(vA) - mean(vB); %#ok<AGROW>
        end
    end
    nPair(k) = numel(d);
    if numel(d) >= 2
        dMean(k) = mean(d);
        pPerm(k) = local_signflip(d, 1e4);
    end
end

pLMEadj  = local_holm(pLME);
pPermadj = local_holm(pPerm);

fprintf('\npairwise (Holm-corrected across %d contrasts)\n', np);
fprintf('%-22s %9s %9s %9s %7s\n','contrast','meanDiff','p_LME','p_paired','n_anim');
for k = 1:np
    fprintf('%-22s %9.3f %9.4f %9.4f %7d\n', ...
        sprintf('%s-%s', uc{pairs(k,1)}, uc{pairs(k,2)}), ...
        dMean(k), pLMEadj(k), pPermadj(k), nPair(k));
end
fprintf(['\nn_anim is how many animals actually contributed a WITHIN-animal\n' ...
         'difference for that contrast. Where it is small, the paired test is\n' ...
         'underpowered and the LME contrast is leaning on between-animal\n' ...
         'information -- report the pair, not just the smaller p.\n']);

% A sign-flip test on n animals has 2^n possible sign assignments, so the
% smallest reachable two-sided p is 2/2^n -- BEFORE multiplicity correction.
% Below ~6 animals the exact test cannot reach 0.05 no matter how large the
% effect is, which is a property of the design, not evidence of absence.
for k = 1:np
    if nPair(k) >= 2
        floorP = 2 / 2^nPair(k);
        if floorP * np > 0.05
            fprintf(['  NOTE %s: n=%d animals -> smallest reachable paired p is ' ...
                     '%.3f (%.3f after Holm). This contrast CANNOT reach 0.05\n' ...
                     '       by the exact test; use the LME contrast and the effect size.\n'], ...
                sprintf('%s-%s', uc{pairs(k,1)}, uc{pairs(k,2)}), nPair(k), ...
                floorP, min(1, floorP*np));
        end
    end
end

out.pairs = uc(pairs);
out.pLME = pLMEadj; out.pPerm = pPermadj;
out.nPair = nPair; out.meanDiff = dMean;
out.nCell = nCell; out.animals = ua; out.conditions = uc;

%% ---- optional 2x2 factorial --------------------------------------------
if doFactorial
    % Baseline is not a drug x timepoint cell -- it has neither factor, so it
    % cannot enter a 2x2 and is dropped here (it stays in the omnibus and in
    % the pairwise contrasts above).
    keepF = ~strcmpi(cond, 'Baseline');
    if any(~keepF)
        fprintf('\n2x2: excluding %d Baseline sessions (no drug/timepoint level).\n', ...
            sum(~keepF));
    end
    condF = cond(keepF); xF = x(keepF); animalF = animal(keepF);
    drug = regexprep(condF, '24$', '');
    t24  = ~cellfun(@isempty, regexp(condF, '24$', 'once'));
    if numel(unique(drug)) ~= 2
        warning('condCompare:factorial', ...
            ['Drug parse gave %d levels (%s) -- expected 2. The 2x2 assumes ', ...
             'condition names are <drug> or <drug>24. Skipping.'], ...
            numel(unique(drug)), strjoin(unique(drug)', ','));
    else
        ftbl = table(xF, categorical(animalF), categorical(drug), ...
            categorical(double(t24)), 'VariableNames', ...
            {'Metric','Animal','Drug','T24'});

        % ASV and ASV24 (or SAL/SAL24) are the SAME drug administration
        % measured on consecutive days, so those two sessions share an
        % episode, not merely an animal. (1|Animal:Drug) captures that.
        % Without it, the day effect and the interaction are tested as if
        % each day were an independent draw from the animal, which
        % understates their standard errors -- the same pseudo-replication
        % issue as elsewhere in this pipeline, one level further down.
        fnested = true;
        try
            flme = fitlme(ftbl, 'Metric ~ Drug*T24 + (1|Animal) + (1|Animal:Drug)', ...
                'FitMethod','REML');
            fpsi = covarianceParameters(flme);
            if numel(fpsi) >= 2 && fpsi{2}(1,1) < 1e-10, fnested = false; end
        catch
            fnested = false;
        end
        if ~fnested
            fprintf(['  (1|Animal:Drug) variance ~0 or failed; refitting with\n' ...
                     '  (1|Animal) only. Expected when few animals have BOTH days\n' ...
                     '  of the same drug.\n']);
            flme = fitlme(ftbl, 'Metric ~ Drug*T24 + (1|Animal)', 'FitMethod','REML');
        end
        fst = anova(flme, 'DFMethod','Satterthwaite');
        fprintf('\n2x2 factorial: Drug (ASV/SAL) x Day (same-day / next-day)\n');
        for k = 1:height(fst)
            if ~strcmp(fst.Term{k},'(Intercept)')
                fprintf('  %-12s F(%d,%.1f) = %7.3f  p = %.4f\n', ...
                    fst.Term{k}, fst.DF1(k), fst.DF2(k), fst.FStat(k), fst.pValue(k));
            end
        end
        fprintf('  cell counts (sessions):\n');
        for dl = unique(drug)'
            for tl = [0 1]
                fprintf('    %-6s day%d : %3d sessions, %2d animals\n', dl{1}, tl, ...
                    sum(strcmp(drug,dl{1}) & t24==tl), ...
                    numel(unique(animalF(strcmp(drug,dl{1}) & t24==tl))));
            end
        end
        out.factorialLME = flme;
        out.factorialNested = fnested;

    end
end
end


function p = local_signflip(d, nperm)
d = d(:); d = d(~isnan(d));
n = numel(d); obs = mean(d);
null = zeros(nperm,1);
for k = 1:nperm
    null(k) = mean((2*(rand(n,1)>0.5)-1) .* d);
end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end


function padj = local_holm(p)
padj = nan(size(p));
ok = ~isnan(p);
pv = p(ok); m = numel(pv);
[ps, ord] = sort(pv);
adj = min(1, ps .* (m:-1:1)');
adj = cummax(adj);
tmp = nan(m,1); tmp(ord) = adj;
padj(ok) = tmp;
end