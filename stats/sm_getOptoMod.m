function [Y,beta_period] = sm_getOptoMod(cond1,cond2,varargin)

% FUNCTION: sm_getOptoMod
%
% PURPOSE:
%   This function performs a Poisson regression to compare two paired
%   conditions (e.g., spike counts under two different experimental
%   periods, such as control vs. optogenetic stimulation).
%
%   It tests whether the mean counts differ significantly between the two
%   periods and returns both the p-value and the estimated regression
%   coefficient for the period effect.
%
% INPUTS:
%   cond1 - [n x 1] numeric array
%           Counts for condition 1 (e.g., baseline or control)
%
%   cond2 - [n x 1] numeric array
%           Counts for condition 2 (e.g., stimulation or observation)
%
%           NOTE: cond1 and cond2 must have the same number of observations
%           (rows), since they are paired measurements.
%
% OUTPUTS:
%   Y            - p-value associated with the 'Period' coefficient in the
%                  Poisson regression (tests whether Period has a
%                  significant effect on count)
%
%   beta_period  - estimated regression coefficient for the 'Period'
%                  variable, representing the log change in counts between
%                  condition 1 and condition 2
%
% EXAMPLE:
%   >> [pval, beta] = sm_getOptoMod(spikes_pre, spikes_stim);
%   % pval gives significance of stimulation effect
%   % beta gives effect size on log scale
%
%==========================================================================

p = inputParser;


addParameter(p,'paired',true,@islogical) %method to correct for decaying signal



parse(p,varargin{:})


paired = p.Results.paired;

if size(cond1,1) ~= size(cond2,1) && paired
    
    error('counts must be paired')
    
end


nObs1 = size(cond1,1);
nObs2 = size(cond2,1);
period = [ones(nObs1,1);2*ones(nObs2,1)];

counts = [cond1;cond2];



% Fit Poisson regression with indicator for period (baseline vs. observation)

if paired
    pairID = [(1:nObs1)'; (1:nObs1)']; % Each pair appears twice (cond1 & cond2);
    
    % Create table
    tbl = table( period, counts, pairID , ...
        'VariableNames', {'Period', 'Count', 'Pair'});
    
    % Convert Period to categorical if preferred (not required here)
    tbl.Period = categorical(tbl.Period);
    
    mdl = fitglme(tbl, 'Count ~ Period + (1|Pair)', 'Distribution', 'Poisson', 'Link', 'log');
    
else
    % Create table
    tbl = table( period, counts, ...
        'VariableNames', {'Period', 'Count'});
    
    % Convert Period to categorical if preferred (not required here)
    tbl.Period = categorical(tbl.Period);
    
    mdl = fitglm(tbl, 'Count ~ Period', 'Distribution', 'poisson', 'Link', 'log');
end
Y = mdl.Coefficients.pValue(2);
beta_period =  mdl.Coefficients.Estimate(2);



end