function [Y,beta_period] = sm_getOptoMod(cond1,cond2)



if size(cond1,1) ~= size(cond2,1)
    
    error('counts must be paired')
    
end


nObs = size(cond1,1);

period = [ones(nObs,1);2*ones(nObs,1)];

counts = [cond1;cond2];

% Create table
tbl = table( period, counts, ...
    'VariableNames', {'Period', 'Count'});

% Convert Period to categorical if preferred (not required here)
 tbl.Period = categorical(tbl.Period);

% Fit Poisson regression with indicator for period (baseline vs. observation)
mdl = fitglm(tbl, 'Count ~ Period', 'Distribution', 'poisson', 'Link', 'log');
Y = mdl.Coefficients.pValue(2);
beta_period =  mdl.Coefficients.Estimate(2);



end