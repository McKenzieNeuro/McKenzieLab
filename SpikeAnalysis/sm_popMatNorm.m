function [X_pv, gain_event, gain_axis] = sm_popMatNorm(R, opts)
%--------------------------------------------------------------------------
% Removes shared multiplicative population drive:
%   R(i,e) = b_i + w_i * g_e + noise
% Ensures mean(X_pv,1) is uncorrelated with mean(R,1)
%
% INPUT
%   R : [neurons x events] spike counts or rates
%   opts.use_spont_only = true/false
%   opts.spont_idx
%   opts.eps = ridge term (default 1e-6)
%
% OUTPUT
%   X_pv       : gain-corrected, variance-equalized PVs
%   gain_event : g_e   (shared drive per event)
%   gain_axis  : neuron gain axis
%--------------------------------------------------------------------------

if nargin < 2
    opts = struct();
end
if ~isfield(opts,'eps')
    opts.eps = 1e-6;
end

[N,E] = size(R);

%% -------------------------------------------------------------
% Determine events used for gain estimation
% -------------------------------------------------------------
if isfield(opts,'use_spont_only') && opts.use_spont_only
    idx = opts.spont_idx;
else
    idx = true(1,E);
end

%% -------------------------------------------------------------
% Remove neuron baselines
% -------------------------------------------------------------
b = mean(R(:,idx), 2);      % N x 1
R_center = R - b;           % centered data

%% -------------------------------------------------------------
% Compute first SVD vector as gain axis
% -------------------------------------------------------------
[U,~,~] = svd(R_center(:,idx), 'econ');
gain_axis = U(:,1);         % N x 1
gain_axis = gain_axis / norm(gain_axis);

%% -------------------------------------------------------------
% Project full data onto gain axis
% -------------------------------------------------------------
gain_event = gain_axis' * R_center;  % 1 x E

%% -------------------------------------------------------------
% Reconstruct gain component and residual
% -------------------------------------------------------------
R_gain = gain_axis * gain_event;  % N x E
R_resid = R_center - R_gain;

%% -------------------------------------------------------------
% Whiten neuron contributions
% -------------------------------------------------------------
sigma = sqrt(std(R_resid(:,idx), 0, 2).^2 + opts.eps);
X_pv = R_resid ./ sigma;

%% -------------------------------------------------------------
% Event-space orthogonalization: remove correlation with mean(R,1)
% AFTER whitening
% -------------------------------------------------------------
mR = mean(R,1);             % 1 x E
if std(mR) > 0
    X_mean = mean(X_pv,1);               % 1 x E
    proj = (X_mean * mR') / (mR * mR');  % scalar projection
    X_pv = X_pv - repmat(proj .* mR, N, 1); % subtract same component from all neurons
end

end
