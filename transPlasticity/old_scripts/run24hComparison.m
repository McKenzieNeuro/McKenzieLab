%% run24hComparison.m
% Drug-naive drug-free epochs vs drug-free epochs 24 h after an ASV exposure.
%
% Assumes runPersistenceAnalysis.m has been run and its results are saved:
%   load(fullfile(CACHE_ROOT,'persistenceResults.mat'))
% This script RECLASSIFIES those same epochs and runs the focused comparison.
%
% NOTE: runPersistenceAnalysis.m now applies the same corrected rule itself
% and saves eHasInf, so its own pre/post split and recency stratification are
% already correct. This script remains the focused 24 h contrast.
%
% THE CLASSIFICATION RULE (fixed in both scripts)
% A drug-free epoch from an infusion session is taken BEFORE that session's
% infusion. The earlier recency code counted exposures with
% expDates <= sessionDate, so a session's OWN same-day dose was treated as
% already delivered. Consequences:
%   - epochs at "0 days since last exposure" were mostly PRE-dose baselines
%   - the epoch from an animal's FIRST ASV session was labelled post-exposure
%     when it is in fact the cleanest drug-naive data that animal has
% Correct rule, applied below:
%   infusion session  -> prior exposures are those STRICTLY BEFORE this date
%   no-infusion session (ASV24, SAL24, ...) -> exposures on or before count
%
% This moves a substantial number of epochs from the post group into the
% naive group, which is where they belong and which is what makes the 24 h
% contrast worth running: the naive group stops being starved and the post
% group stops being diluted with pre-dose data.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
load(fullfile(CACHE_ROOT,'persistenceResults.mat'), ...
     'ea','ec','ed','eHasInf','rE','oE','dE','mE','daysInStudy');
EXPOSURE_PAT = '^ASV$';

% Whether each epoch precedes that session's infusion, saved by
% runPersistenceAnalysis. Earlier this was inferred from the condition name
% (no trailing "24", not Baseline), which would misclassify any ASV/SAL
% folder session that happened to lack an infusion marker.
hasInf = eHasInf(:);

%% ---- corrected exposure timing ------------------------------------------
uan = unique(ea);
dLast2 = nan(numel(ea),1); everExposedBefore = false(numel(ea),1);
for i = 1:numel(uan)
    m = find(strcmp(ea, uan{i}));
    expDates = unique(ed(m(~cellfun(@isempty, regexp(ec(m), EXPOSURE_PAT,'once')))));
    if isempty(expDates), continue; end
    for j = m'
        if hasInf(j)
            prior = expDates(expDates <  ed(j));   % strict: own dose not yet given
        else
            prior = expDates(expDates <= ed(j));
        end
        if ~isempty(prior)
            dLast2(j) = ed(j) - max(prior);
            everExposedBefore(j) = true;
        end
    end
end

isNaive = ~everExposedBefore;              % no exposure before this epoch
is24h   = everExposedBefore & dLast2 == 1;

fprintf('after correction:\n');
fprintf('  drug-naive epochs : %d (%d animals)\n', sum(isNaive), numel(unique(ea(isNaive))));
fprintf('  24h-post epochs   : %d (%d animals)\n', sum(is24h),   numel(unique(ea(is24h))));
both = intersect(unique(ea(isNaive)), unique(ea(is24h)));
fprintf('  animals with BOTH : %d (%s)\n', numel(both), strjoin(both', ', '));
fprintf('\n(previously "0 days since exposure" held %d epochs that were really\n', ...
        sum(hasInf & ~everExposedBefore));
fprintf(' pre-dose baselines.)\n');

%% ---- comparison ----------------------------------------------------------
MET = {'duty cycle', oE; 'transient rate', rE; 'log duration', log(dE); 'log magnitude', log(mE)};

fprintf('\n%-15s %9s %9s %9s %9s %8s\n','metric','naive','24h','diff','p_LME','p_paired');
for k = 1:size(MET,1)
    y = MET{k,2};
    g = (isNaive | is24h) & isfinite(y);
    if numel(unique(is24h(g))) < 2, continue; end

    grp = repmat({'naive'}, numel(y),1); grp(is24h) = {'x24h'};
    tb = table(y(g), categorical(ea(g)), categorical(grp(g),{'naive','x24h'}), daysInStudy(g), ...
        'VariableNames',{'M','Animal','Grp','Days'});
    mm = fitlme(tb,'M ~ Days + Grp + (1|Animal)','FitMethod','ML');
    aa = anova(mm,'DFMethod','Satterthwaite');
    pL = aa.pValue(strcmp(aa.Term,'Grp'));

    % exact within-animal paired test on animals having both
    d = [];
    for u = both'
        vN = y(isNaive & strcmp(ea,u{1}) & isfinite(y));
        v24 = y(is24h  & strcmp(ea,u{1}) & isfinite(y));
        if ~isempty(vN) && ~isempty(v24), d(end+1,1) = mean(v24) - mean(vN); end %#ok<AGROW>
    end
    if numel(d) >= 2
        nperm = 1e4; null = zeros(nperm,1);
        for q = 1:nperm, null(q) = mean((2*(rand(numel(d),1)>0.5)-1).*d); end
        pP = (1 + sum(abs(null) >= abs(mean(d)))) / (1 + nperm);
    else
        pP = NaN;
    end

    fprintf('%-15s %9.3f %9.3f %+9.3f %9.4f %8.4f  (paired n=%d)\n', MET{k,1}, ...
        mean(y(isNaive & isfinite(y))), mean(y(is24h & isfinite(y))), ...
        mean(y(is24h & isfinite(y))) - mean(y(isNaive & isfinite(y))), pL, pP, numel(d));
end

fprintf(['\np_LME uses all epochs with a time-in-study covariate; p_paired uses\n' ...
         'only animals contributing both, and is exact. With few paired animals\n' ...
         'the paired test has a hard p floor (2/2^n) -- check n before reading it\n' ...
         'as absence of an effect.\n']);