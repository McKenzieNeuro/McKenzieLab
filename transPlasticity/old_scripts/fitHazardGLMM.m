% FITHAZARDGLMM  one piecewise-exponential Poisson GLMM contrast; moved out of
% runHazardAnalysis.m so the stress tests (stressRandomSlope.m) call exactly the
% code the pipeline uses. fitHazardGLMM() with no arguments returns a blank result.
function F = fitHazardGLMM(t, fixed, testTerm, coefName, slope, slopeGroup)
% Piecewise-exponential Poisson GLMM; log hazard ratios. LRT for the tested
% term via compare(), since Wald p on a GLMM with ~10 clusters is optimistic.
%
% RANDOM SLOPE. With intercept-only REs the pre->post effect is assumed
% identical in every animal, so its SE comes from the hundreds of bouts, not
% from the ~8 animals: p is anticonservative whenever animals differ in
% their response. The (1+slope|Animal) refit lets the effect vary by animal;
% its p is the one to quote, and it should sit near the animal-level
% sign-flip p.
%
% slopeGroup 'AnimalDrug' (DiD only): the slope varies by animal x drug,
% (1|Animal) + (1+slope|AnimalDrug). (1+Period|Animal) lets only the
% animal's SHARED pre->post change vary, so an animal-specific drug effect
% -- what the Drug:Period term is tested against -- was not in the error
% term, which is why the DiD model p sat far below the animal-level p.
%
% OFFSET. A table column OffsetAdd, if present, is added to log(exposure)
% (used to fix a drift slope estimated elsewhere).
if nargin == 0, F = local_blankFit(); return; end
if nargin < 5, slope = ''; end
if nargin < 6 || isempty(slopeGroup), slopeGroup = 'Animal'; end
F = local_blankFit();
if height(t) < 10 || sum(t.Event) < 5
    F.why = sprintf('%d obs, %d events', height(t), sum(t.Event)); return
end
nullFixed = local_dropTerm(fixed, testTerm);
off = log(t.Exposure);
if ismember('OffsetAdd', t.Properties.VariableNames), off = off + t.OffsetAdd; end
opts = {'Distribution','Poisson','Link','log','Offset',off,'FitMethod','Laplace'};
re = ' + (1|Animal) + (1|Session)';
try
    full = fitglme(t, ['Event ~ ' fixed re], opts{:});
    nul  = fitglme(t, ['Event ~ ' nullFixed re], opts{:});
catch ME1
    re = ' + (1|Animal)'; F.reSession = false;
    try
        full = fitglme(t, ['Event ~ ' fixed re], opts{:});
        nul  = fitglme(t, ['Event ~ ' nullFixed re], opts{:});
    catch ME2
        F.why = sprintf('fitglme failed: %s | without session RE: %s', ME1.message, ME2.message);
        return
    end
end
% match the tested coefficient by its parts, not an exact string: MATLAB can
% order interaction factors either way ('Drug_ASV:Period_post' vs
% 'Period_post:Drug_ASV')
parts = strsplit(coefName, ':');
cn = full.CoefficientNames;
k = find(cellfun(@(c) all(cellfun(@(q) contains(c,q), parts)) && ...
                      numel(strsplit(c,':')) == numel(parts), cn));
if numel(k) ~= 1
    F.why = sprintf('coefficient ''%s'' not found; model has: %s', coefName, strjoin(cn, ', '));
    return
end
ci = coefCI(full);
F.HR = exp(full.Coefficients.Estimate(k));
F.CI = exp(ci(k,:));
F.pWald = full.Coefficients.pValue(k);
c = compare(nul, full);
F.pLRT = c.pValue(2);
F.model = full; F.ok = true;

if isempty(slope), F.rs.why = 'not requested'; return; end
if strcmp(slopeGroup, 'Animal')
    reS = strrep(re, '(1|Animal)', sprintf('(1+%s|Animal)', slope));
else
    t.AnimalDrug = categorical(strcat(cellstr(t.Animal), '_', cellstr(t.Drug)));
    reS = strrep(re, '(1|Animal)', sprintf('(1|Animal) + (1+%s|AnimalDrug)', slope));
end
% session intercepts can be near-collinear with the extra slope grouping
% (one session per animal x drug); retry without them before giving up
reTry = {reS};
if contains(reS, '(1|Session)'), reTry{end+1} = strrep(reS, ' + (1|Session)', ''); end
msg = {};
for r = 1:numel(reTry)
    try
        fullS = fitglme(t, ['Event ~ ' fixed reTry{r}], opts{:});
        nulS  = fitglme(t, ['Event ~ ' nullFixed reTry{r}], opts{:});
        cnS = fullS.CoefficientNames;
        kS = find(cellfun(@(c) all(cellfun(@(q) contains(c,q), parts)) && ...
                               numel(strsplit(c,':')) == numel(parts), cnS));
        ciS = coefCI(fullS);
        F.rs.HR = exp(fullS.Coefficients.Estimate(kS));
        F.rs.CI = exp(ciS(kS,:));
        F.rs.pWald = fullS.Coefficients.pValue(kS);
        cS = compare(nulS, fullS);
        F.rs.pLRT = cS.pValue(2);
        F.rs.ok = true; F.rs.re = strtrim(reTry{r}); F.rs.model = fullS;
        return
    catch ME
        msg{end+1} = ME.message; %#ok<AGROW>
    end
end
F.rs.why = strjoin(msg, ' | ');
end


function F = local_blankFit()
F = struct('ok',false,'why','','HR',NaN,'CI',[NaN NaN],'pWald',NaN,'pLRT',NaN, ...
           'reSession',true,'model',[], ...
           'rs', struct('ok',false,'why','','HR',NaN,'CI',[NaN NaN],'pWald',NaN,'pLRT',NaN, ...
                        're','','model',[]));
end


function s = local_dropTerm(fixed, term)
% null-model fixed effects: the full formula minus the tested term
if strcmp(term, 'Drug:Period'), s = 'Drug + Period'; return; end
tok = strtrim(strsplit(fixed, '+'));
if ~any(strcmp(tok, term)), error('term %s not in %s', term, fixed); end
tok(strcmp(tok, term)) = [];
if isempty(tok), s = '1'; else, s = strjoin(tok, ' + '); end
end
