function out = prePostASV(T, x, rootName, label, varargin)
%PREPOSTASV  Compare everything before an animal's first ASV to everything after.
%
%   out = prePostASV(T, x, rootName, label, 'Name', value, ...)
%
%   Motivated by the observation that ASV can produce synchrony changes that
%   do not revert: if so, the right contrast is not drug-vs-vehicle on a
%   given day but the animal's whole recording history split at its first
%   ASV exposure. Within the five conditions of interest
%   (Baseline/ASV/ASV24/SAL/SAL24), a SAL or Baseline session recorded AFTER
%   the first ASV is post-phase, and is the most informative kind of session
%   in this design: same vehicle condition, no drug on board, so a difference
%   there is evidence of a persistent state change rather than acute
%   pharmacology.
%
%   NAME/VALUE
%     'Conditions'  : condition folders to KEEP. Everything else is dropped
%                     before the phase split. Default
%                       {'Baseline','ASV','ASV24','SAL','SAL24'}
%                     This matters more here than in condCompare: other
%                     folders (ANI, Thapsigargin, ...) are their own
%                     pharmacological manipulations, and letting them into
%                     either phase would attribute their effects to ASV
%                     exposure. TrackHabituation and donotanalyze are
%                     likewise excluded by not being on the list.
%     'ASVPattern'  : regexp identifying an ASV condition. Default '^ASV'
%                     (matches ASV and ASV24). Change if other folders count
%                     as exposure.
%
%   THE CONFOUND, STATED PLAINLY. Post-ASV is always later than pre-ASV, so
%   "irreversible ASV effect" and "electrode drift / habituation / implant
%   ageing" predict the SAME sign of difference. A significant Phase term on
%   its own therefore proves nothing. What distinguishes them is SHAPE:
%
%     persistent drug effect -> a STEP at first ASV, roughly flat either side
%     drift                  -> a gradual TREND with no special status for
%                               the ASV date
%
%   So three models are fitted and compared:
%     M1  Metric ~ Phase           + (1|Animal)   step only
%     M2  Metric ~ DaysInStudy     + (1|Animal)   drift only
%     M3  Metric ~ DaysInStudy + Phase + (1|Animal)  step ON TOP of drift
%
%   The test that matters is Phase in M3: does a step at first ASV explain
%   variance that a smooth time trend does not? If Phase is significant in
%   M1 but vanishes in M3, you have drift, not a persistent drug effect. If
%   Phase survives M3, that is genuine evidence for a step. AIC/BIC for all
%   three are reported so you can see which shape the data prefer overall.
%
%   Identification is only PARTIAL. If every animal's pre sessions are
%   tightly clustered before, and post sessions tightly after, DaysInStudy
%   and Phase are nearly collinear and M3 cannot separate them -- the
%   printed overlap diagnostic tells you whether there is enough temporal
%   spread on both sides for the comparison to mean anything.
%
%   See also SESSIONDATEFROMPATH, CONDCOMPARE, ANIMALFROMPATH.

ip = inputParser;
ip.addParameter('Conditions', {'Baseline','ASV','ASV24','SAL','SAL24'}, @iscell);
ip.addParameter('ASVPattern', '^ASV', @ischar);
ip.parse(varargin{:});
keepConds = ip.Results.Conditions;
asvPat    = ip.Results.ASVPattern;

animal = arrayfun(@(r) animalFromPath(r.file, rootName), T, 'UniformOutput', false);
cond   = arrayfun(@(r) conditionFromPath(r.file, rootName), T, 'UniformOutput', false);
dnum   = arrayfun(@(r) sessionDateFromPath(r.file), T);
animal = animal(:); cond = cond(:); dnum = dnum(:); x = x(:);

bad = cellfun(@isempty, animal) | cellfun(@isempty, cond) | isnan(x) | isnan(dnum);
if any(bad)
    fprintf('%s: dropping %d/%d sessions (unparsed path/date or NaN metric)\n', ...
        label, sum(bad), numel(x));
end
notWanted = ~ismember(lower(cond), lower(keepConds));
if any(notWanted & ~bad)
    fprintf('%s: excluding %d sessions in non-requested condition(s): %s\n', ...
        label, sum(notWanted & ~bad), strjoin(unique(cond(notWanted & ~bad))', ', '));
end
drop = bad | notWanted;
animal = animal(~drop); cond = cond(~drop); dnum = dnum(~drop); x = x(~drop);

% NOTE ON THE PHASE SPLIT AFTER FILTERING: the first-ASV date is computed
% from the RETAINED sessions only. That is correct as long as every ASV
% session is on the keep list -- if you ever narrow 'Conditions' so that
% some ASV exposure is filtered out, an animal's true first exposure could
% be earlier than the date used here, and post-phase sessions would be
% mislabelled pre. Keep ASV (and any other exposure folder) on the list.

ua = unique(animal);

%% ---- classify each session -------------------------------------------
phase = repmat({''}, numel(x), 1);
daysRel = nan(numel(x),1);
hasASV = false(numel(ua),1);

for i = 1:numel(ua)
    m = strcmp(animal, ua{i});
    isASV = m & ~cellfun(@isempty, regexp(cond, asvPat, 'once'));
    if ~any(isASV)
        phase(m) = {'pre'};              % never exposed: all pre
        daysRel(m) = dnum(m) - min(dnum(m));
        continue
    end
    hasASV(i) = true;
    t0 = min(dnum(isASV));               % first ASV date
    daysRel(m) = dnum(m) - t0;
    phase(m & dnum <  t0) = {'pre'};
    phase(m & dnum >= t0) = {'post'};
end

%% ---- design report -----------------------------------------------------
fprintf('\n=== pre/post first ASV: %s ===\n', label);
fprintf('%d sessions, %d animals (%d ever received ASV)\n', ...
    numel(x), numel(ua), sum(hasASV));
fprintf('\n%-10s %5s %5s   %-28s %s\n','animal','pre','post','pre conditions','span (d, pre|post)');
nBoth = 0;
for i = 1:numel(ua)
    m = strcmp(animal, ua{i});
    np = sum(m & strcmp(phase,'pre')); nq = sum(m & strcmp(phase,'post'));
    if np>0 && nq>0, nBoth = nBoth + 1; end
    pc = unique(cond(m & strcmp(phase,'pre')));
    sp = daysRel(m & strcmp(phase,'pre')); sq = daysRel(m & strcmp(phase,'post'));
    fprintf('%-10s %5d %5d   %-28s %s | %s\n', ua{i}, np, nq, ...
        strjoin(pc', ','), local_span(sp), local_span(sq));
end
fprintf('\n%d/%d animals have BOTH phases -- only these give a within-animal contrast.\n', ...
    nBoth, numel(ua));

% the key sub-population: vehicle/baseline sessions on each side
isVeh = ~cellfun(@isempty, regexp(cond, '^(SAL|Baseline)', 'once'));
fprintf('non-ASV-condition sessions: %d pre, %d post', ...
    sum(isVeh & strcmp(phase,'pre')), sum(isVeh & strcmp(phase,'post')));
fprintf('  <- post-phase vehicle sessions are the cleanest evidence of persistence\n');

if nBoth < 3
    warning('prePostASV:balance', ...
        'Only %d animals have both phases; the Phase term is mostly between-animal.', nBoth);
end

%% ---- models -----------------------------------------------------------
keep = ~cellfun(@isempty, phase);
tbl = table(x(keep), categorical(animal(keep)), ...
    categorical(phase(keep), {'pre','post'}), daysRel(keep), ...
    'VariableNames', {'Metric','Animal','Phase','DaysInStudy'});

M1 = fitlme(tbl, 'Metric ~ Phase + (1|Animal)', 'FitMethod','ML');
M2 = fitlme(tbl, 'Metric ~ DaysInStudy + (1|Animal)', 'FitMethod','ML');
M3 = fitlme(tbl, 'Metric ~ DaysInStudy + Phase + (1|Animal)', 'FitMethod','ML');
% ML (not REML) because AIC/BIC across models with DIFFERENT fixed effects
% is only comparable under ML.

a1 = anova(M1,'DFMethod','Satterthwaite');
a3 = anova(M3,'DFMethod','Satterthwaite');
p1 = a1.pValue(strcmp(a1.Term,'Phase'));
p3 = a3.pValue(strcmp(a3.Term,'Phase'));
pD = a3.pValue(strcmp(a3.Term,'DaysInStudy'));

fprintf('\nM1 step only    : Phase p = %.4f          AIC %.1f\n', p1, M1.ModelCriterion.AIC);
fprintf('M2 drift only   : DaysInStudy p = %.4f   AIC %.1f\n', ...
    a3.pValue(strcmp(a3.Term,'DaysInStudy')), M2.ModelCriterion.AIC);
fprintf('M3 step + drift : Phase p = %.4f, Days p = %.4f, AIC %.1f\n', p3, pD, M3.ModelCriterion.AIC);

ci = coefCI(M3);
ip3 = find(strcmp(M3.CoefficientNames,'Phase_post'));
if ~isempty(ip3)
    fprintf('\nM3 Phase_post estimate = %+.3f  [%.3f, %.3f]\n', ...
        M3.Coefficients.Estimate(ip3), ci(ip3,1), ci(ip3,2));
end

fprintf('\nverdict: ');
if p1 < 0.05 && p3 < 0.05
    fprintf(['step survives adjustment for a smooth time trend -- consistent with\n' ...
             '         a persistent change at first ASV rather than drift.\n']);
elseif p1 < 0.05 && p3 >= 0.05
    fprintf(['step is significant alone but NOT after adjusting for time --\n' ...
             '         this looks like drift, not a persistent ASV effect.\n']);
else
    fprintf('no phase effect detected either way.\n');
end

%% ---- vehicle-only sensitivity -----------------------------------------
% Repeat using ONLY non-ASV-condition sessions. Removes any acute drug-on-
% board contribution, so a surviving effect is persistence, not pharmacology.
vk = keep & isVeh;
if sum(vk & strcmp(phase,'pre')) >= 3 && sum(vk & strcmp(phase,'post')) >= 3
    vtbl = table(x(vk), categorical(animal(vk)), ...
        categorical(phase(vk), {'pre','post'}), daysRel(vk), ...
        'VariableNames', {'Metric','Animal','Phase','DaysInStudy'});
    V3 = fitlme(vtbl, 'Metric ~ DaysInStudy + Phase + (1|Animal)', 'FitMethod','ML');
    av = anova(V3,'DFMethod','Satterthwaite');
    fprintf('\nvehicle/baseline sessions only (n=%d): Phase p = %.4f (time-adjusted)\n', ...
        sum(vk), av.pValue(strcmp(av.Term,'Phase')));
    out.vehicleLME = V3;
else
    fprintf('\nvehicle-only model skipped: too few sessions on one side.\n');
end

out.M1 = M1; out.M2 = M2; out.M3 = M3;
out.pPhaseAlone = p1; out.pPhaseAdjusted = p3; out.pDays = pD;
out.table = tbl; out.phase = phase; out.daysRel = daysRel;
out.animal = animal; out.cond = cond;
end


function s = local_span(d)
if isempty(d), s = '-'; else, s = sprintf('%.0f..%.0f', min(d), max(d)); end
end
