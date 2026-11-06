function [lme, tbl] = lmeCompare(a, b, x1, x2, rootA, rootB, label)
%LMECOMPARE  Mixed-effects test: dataset fixed effect, animal random intercept.
%
%   [lme, tbl] = lmeCompare(a, b, x1, x2, rootA, rootB, label)
%
%   The per-animal permutation test in nestingCheck.m is conservative by
%   construction: collapsing every animal to one mean throws away all
%   within-animal session-to-session variation, which is exactly the
%   information that tells you how PRECISELY each animal's mean is known.
%   An animal with 12 sessions and one with 1 session get equal weight in
%   that test even though the 12-session mean is far better estimated.
%
%   A mixed model uses that information instead of discarding it: dataset is
%   a fixed effect, animal is a random intercept NESTED WITHIN dataset (an
%   animal cannot belong to both groups here, but nesting is stated
%   explicitly via Group:Animal rather than assumed from naming, in case
%   animal IDs ever collide across datasets). This can recover power a
%   permutation-of-means test cannot, while still respecting -- not
%   pretending away -- the fact that sessions from the same animal are
%   correlated.
%
%   MODEL: Metric ~ Group + (1|Group:Animal)
%
%   DF METHOD: Satterthwaite, not the fitlme default (Residual). The default
%   uses N_sessions - (fixed-effect params) degrees of freedom, which treats
%   sessions as if they were as informative as independent animals -- the
%   same optimism the session-level permutation test has, just wearing a
%   parametric model instead of a nonparametric one. Satterthwaite discounts
%   the DF toward the number of CLUSTERS (animals), which is what a fixed
%   effect estimated from 10-12 animals can actually support. Even so, with
%   only 10-12 clusters per group, asymptotic p-values (of any flavor) are
%   approximate -- report this alongside, not instead of, the exact
%   permutation test in nestingCheck.m.
%
%   OUTPUT
%     lme : the fitted LinearMixedModel object (inspect with anova(lme),
%           lme.Coefficients, or plotResiduals(lme) for diagnostics)
%     tbl : the table the model was fit on, one row per session, with NaN
%           rows (sessions with no significant transients, e.g. for the
%           duration metric) already excluded
%
%   See also NESTINGCHECK, ANIMALFROMPATH.

if isempty(ver('stats'))
    error('lmeCompare:toolbox', ...
        'fitlme requires the Statistics and Machine Learning Toolbox.');
end
if nargin < 7, label = 'metric'; end

idA = arrayfun(@(r) animalFromPath(r.file, rootA), a, 'UniformOutput', false);
idB = arrayfun(@(r) animalFromPath(r.file, rootB), b, 'UniformOutput', false);

Group  = [repmat({rootA}, numel(a),1); repmat({rootB}, numel(b),1)];
Animal = [idA(:); idB(:)];
Metric = [x1(:); x2(:)];

tbl = table(Metric, categorical(Group), categorical(Animal), ...
    'VariableNames', {'Metric','Group','Animal'});

bad = isnan(tbl.Metric) | cellfun(@isempty, cellstr(tbl.Animal));
if any(bad)
    fprintf('%s: dropping %d/%d sessions (NaN metric or unparsed animal ID)\n', ...
        label, sum(bad), height(tbl));
    tbl(bad,:) = [];
end

lme = fitlme(tbl, 'Metric ~ Group + (1|Group:Animal)', ...
    'FitMethod', 'REML', 'DummyVarCoding', 'effects');

st = anova(lme, 'DFMethod', 'Satterthwaite');
pRow = strcmp(st.Term, 'Group');

fprintf('\n=== LME: %s ===\n', label);
fprintf('  Group fixed effect: F(%d,%.1f) = %.3f, p = %.4f\n', ...
    st.DF1(pRow), st.DF2(pRow), st.FStat(pRow), st.pValue(pRow));

[psi, mse, statsTbl] = covarianceParameters(lme); %#ok<ASGLU>
animalSD   = sqrt(psi{1}(1,1));
residualSD = sqrt(mse);
iccLME = animalSD^2 / (animalSD^2 + residualSD^2);
fprintf('  animal SD = %.4f, residual (within-animal) SD = %.4f, ICC = %.3f\n', ...
    animalSD, residualSD, iccLME);
fprintf('  (%d sessions, %d + %d animals)\n', height(tbl), numel(unique(idA)), numel(unique(idB)));

fprintf(['\n  compare this p to nestingCheck''s session-level and animal-level p:\n' ...
         '  session-level treats all sessions as independent (usually too liberal);\n' ...
         '  animal-level-mean treats every animal as equally informative regardless\n' ...
         '  of session count (usually too conservative); this LME sits between them,\n' ...
         '  weighting each animal by how much data it actually contributed.\n']);
end