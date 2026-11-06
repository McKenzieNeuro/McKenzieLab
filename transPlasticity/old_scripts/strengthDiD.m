function D = strengthDiD(St)
%STRENGTHDID  ASV-vs-SAL difference-in-differences of within-bout strength.
%
%   D = strengthDiD(St)
%
%   St : per-bout table from runHazardAnalysis (columns M, Animal, Session,
%        Period [pre/post], Drug [SAL/ASV], AnimalDrug).
%   D  : .b, .p    LME  M ~ Drug*Period + (1|Animal) + (1|Session)
%        .bRS,.pRS LME  M ~ Drug*Period + (1+Period|AnimalDrug) + (1|Session)
%   The interaction coefficient is located by name pattern (it contains both
%   'Drug' and 'Period'), so the order in which MATLAB writes the name does not
%   matter. p is the Satterthwaite F test for the Drug:Period term when found,
%   otherwise the coefficient t-test.
D = struct('b',NaN,'p',NaN,'bRS',NaN,'pRS',NaN);
if height(St) < 5, return; end
[D.b,   D.p]   = local_int(St, 'M ~ Drug*Period + (1|Animal) + (1|Session)');
[D.bRS, D.pRS] = local_int(St, 'M ~ Drug*Period + (1+Period|AnimalDrug) + (1|Session)');
end


function [b, p] = local_int(St, formula)
b = NaN; p = NaN;
try
    lme = fitlme(St, formula);
    cn = lme.CoefficientNames;
    i = find(contains(cn, 'Drug') & contains(cn, 'Period'), 1);
    if isempty(i), return; end
    b = lme.Coefficients.Estimate(i);
    p = lme.Coefficients.pValue(i);
    try
        av = anova(lme, 'DFMethod', 'Satterthwaite');
        j = find(contains(av.Term, 'Drug') & contains(av.Term, 'Period'), 1);
        if ~isempty(j), p = av.pValue(j); end
    catch
    end
catch ME
    warning('strengthDiD:fail', '%s: %s', formula, ME.message);
end
end
