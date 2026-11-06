%% bootstrapPersistenceLME.m
% Is the MatchN=20 persistence p-value trustworthy, or an artifact of the
% Satterthwaite df approximation at low cluster count?
%
% The concern: with ~6-10 animals (clusters), the Satterthwaite-corrected
% p-value from anova(lme) can understate uncertainty. A parametric bootstrap
% likelihood-ratio test does not rely on that approximation -- it simulates
% the reference distribution of the LRT statistic under the null model
% directly, so it is valid at small cluster counts where the asymptotic
% chi-square (and Satterthwaite) are not.
%
% METHOD
%   full model : M ~ Days + Phase + (1|Animal)
%   null model : M ~ Days        + (1|Animal)   (Phase removed)
%   observed LRT = 2*(logLik_full - logLik_null)
%   Under H0 (no Phase effect), simulate NBOOT datasets FROM THE NULL MODEL,
%   refit both models to each, collect the LRT distribution, and compare the
%   observed statistic to it. p_boot = fraction of simulated LRTs >= observed.
%
% Both models are fit by ML (not REML): likelihoods are only comparable
% across models with different fixed effects under ML.
%
% Run after runPersistenceAnalysis.m has produced persistenceResults.mat at
% the MatchN of interest. Set MATCHN_LABEL for the printout only.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
NBOOT = 2000;
MATCHN_LABEL = 20;   % informational only; reflects whatever run produced the .mat

P = load(fullfile(CACHE_ROOT,'persistenceResults.mat'));
ea = P.ea(:); phase = P.phase(:); dLast = P.dLast(:);

% same day-in-study covariate the step analysis uses
uan = unique(ea);
daysInStudy = nan(numel(ea),1);
for i = 1:numel(uan)
    m = strcmp(ea, uan{i});
    daysInStudy(m) = P.ed(m) - min(P.ed(m));
end

MET = {'duty cycle', P.oE(:); 'transient rate', P.rE(:); ...
       'log duration', log(P.dE(:)); 'log magnitude', log(P.mE(:))};

sel0 = ismember(phase,{'pre','post'}) & (strcmp(phase,'pre') | dLast >= 1);

rng(7);  % reproducible bootstrap
fprintf('\nParametric bootstrap LRT for the Phase term (MatchN=%d data, %d sims)\n', ...
    MATCHN_LABEL, NBOOT);
fprintf('%-15s %8s %10s %10s %8s %8s\n', ...
    'metric','LRT','p_satt','p_boot','nAnim','nEpoch');

for k = 1:size(MET,1)
    y = MET{k,2};
    g = sel0 & isfinite(y);
    if numel(unique(phase(g))) < 2, continue; end

    tb = table(y(g), categorical(ea(g)), categorical(phase(g),{'pre','post'}), daysInStudy(g), ...
        'VariableNames', {'M','Animal','Phase','Days'});

    full = fitlme(tb, 'M ~ Days + Phase + (1|Animal)', 'FitMethod','ML');
    null = fitlme(tb, 'M ~ Days + (1|Animal)',         'FitMethod','ML');

    obsLRT = 2*(full.LogLikelihood - null.LogLikelihood);

    % Satterthwaite p for comparison
    aa = anova(full, 'DFMethod','Satterthwaite');
    pSatt = aa.pValue(strcmp(aa.Term,'Phase'));

    % --- parametric bootstrap from the NULL model ---
    bootLRT = nan(NBOOT,1);
    ok = 0;
    for b = 1:NBOOT
        ysim = random(null);                 % simulate under H0 (no Phase)
        tb.Msim = ysim;
        try
            f2 = fitlme(tb, 'Msim ~ Days + Phase + (1|Animal)', 'FitMethod','ML');
            n2 = fitlme(tb, 'Msim ~ Days + (1|Animal)',         'FitMethod','ML');
            bootLRT(b) = 2*(f2.LogLikelihood - n2.LogLikelihood);
            ok = ok + 1;
        catch
            % occasional singular fit on a simulated set; skip
        end
    end
    bootLRT = bootLRT(isfinite(bootLRT));
    pBoot = (1 + sum(bootLRT >= obsLRT)) / (1 + numel(bootLRT));

    fprintf('%-15s %8.3f %10.4f %10.4f %8d %8d\n', MET{k,1}, obsLRT, pSatt, pBoot, ...
        numel(unique(ea(g))), sum(g));
end

fprintf(['\np_satt is the Satterthwaite p from anova(lme); p_boot is the\n' ...
         'bootstrap LRT p, valid at low cluster count. If they agree, the\n' ...
         'Satterthwaite p is trustworthy. If p_boot is substantially larger,\n' ...
         'the small Satterthwaite p was an artifact of the df approximation.\n' ...
         'Bootstrap resamples that failed to fit were dropped (count reflected\n' ...
         'in the effective NBOOT).\n']);
