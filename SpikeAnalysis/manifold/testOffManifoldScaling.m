function result = testOffManifoldScaling(v_par, v_perp, gain)
% TESTOFFMANIFOLDSCALING
% Test the hypothesis that off-manifold movement is purely a scaling
% of the on-manifold pattern.
%
% Inputs:
%   v_par  : cell array {1 x N_B} of tangent vectors (F x 1)
%   v_perp : cell array {1 x N_B} of normal vectors (F x 1)
%   gain   : N_B x 1 vector of radial gains (||v_par||)
%
% Output:
%   result : struct with fields
%       cosThetaMean        - mean absolute cosine between tangent & normal
%       cosThetaAll          - vector of |cos(theta)| per point
%       gainSlope            - slope beta from ||v_perp|| ~ alpha + beta*||v_par||
%       gainRSquared         - R^2 of above regression
%       residualPCAFrac      - fraction of variance in residual after scaling
%       residualPCAExplained - vector of explained variance per PC
%       residuals            - residuals after removing best scaling
%       betaPerPoint         - best-fit scaling coefficient per point
%       summary              - textual summary of evidence

N_B = numel(v_par);
F = length(v_par{1});

cosThetaAll = nan(N_B,1);
normVperp   = nan(N_B,1);
betaPerPoint = nan(N_B,1);
residuals    = zeros(N_B,F);

% ----- 1. Angular alignment & scaling coefficients -----
for i = 1:N_B
    vp = v_par{i};
    vn = v_perp{i};
    nvp = norm(vp);
    nvn = norm(vn);
    normVperp(i) = nvn;
    
    if nvp>0 && nvn>0
        cosThetaAll(i) = abs(dot(vp,vn)/(nvp*nvn));
        % Best-fit scaling along tangent
        betaPerPoint(i) = dot(vp,vn)/(nvp^2);
        % Residual after removing tangent scaling
        residuals(i,:) = vn - betaPerPoint(i)*vp;
    else
        cosThetaAll(i) = NaN;
        betaPerPoint(i) = 0;
        residuals(i,:) = vn;
    end
end

cosThetaMean = nanmean(cosThetaAll);

% ----- 2. Gain vs normal regression -----
tbl = table(gain, normVperp, 'VariableNames', {'gain','normVperp'});
mdl = fitlm(tbl,'normVperp ~ gain');
gainSlope = mdl.Coefficients.Estimate(2);
gainRSquared = mdl.Rsquared.Ordinary;

% ----- 3. Residual PCA -----
[~,S,~] = svd(residuals,'econ');
singvals = diag(S);
residualPCAExplained = singvals.^2 / sum(singvals.^2);
residualPCAFrac = sum(residualPCAExplained(1:1)); % fraction explained by first PC

% ----- 4. Summary -----
summaryLines = {};
summaryLines{end+1} = sprintf('Mean |cos(theta)| between tangent & normal: %.3f', cosThetaMean);
summaryLines{end+1} = sprintf('Gain vs normal regression: slope=%.3f, R^2=%.3f', gainSlope, gainRSquared);
summaryLines{end+1} = sprintf('Residual variance fraction explained by first PC: %.3f', residualPCAFrac);
if cosThetaMean>0.8 && gainRSquared>0.7 && residualPCAFrac<0.2
    summaryLines{end+1} = 'Evidence supports hypothesis: off-manifold is mostly tangent scaling.';
else
    summaryLines{end+1} = 'Evidence does NOT fully support hypothesis: off-manifold likely introduces new directions.';
end

% ----- 5. Pack results -----
result.cosThetaMean        = cosThetaMean;
result.cosThetaAll         = cosThetaAll;
result.gainSlope           = gainSlope;
result.gainRSquared        = gainRSquared;
result.residualPCAFrac     = residualPCAFrac;
result.residualPCAExplained= residualPCAExplained;
result.residuals           = residuals;
result.betaPerPoint        = betaPerPoint;
result.summary             = summaryLines;

end
