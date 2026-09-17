function [posterior, info] = bayesianDecode(countMat, rateMaps, tau, prior)
% bayesianDecode  Memoryless (one-step) Poisson Bayesian decoder.
%   Zhang et al. 1998 J Neurophysiol; Davidson, Kloosterman & Wilson 2009 Neuron.
%
%   [posterior, info] = bayesianDecode(countMat, rateMaps, tau, prior)
%
% INPUTS
%   countMat : [nCells x nTimeBins] spike counts per cell per decoding bin.
%   rateMaps : [nCells x nPosBins]  firing rate (Hz), strictly positive
%              (use computeRateMaps1D with a rateFloor so log() is finite).
%   tau      : decoding bin width (s).
%   prior    : [1 x nPosBins] spatial prior (default uniform).
%
% OUTPUT
%   posterior : [nPosBins x nTimeBins], each column sums to 1.
%               Bins with zero spikes are set to the (flat) prior.
%   info      : struct with .mapBin (MAP position bin per time bin),
%               .emptyBins (logical), .nActiveCells, .nSpikeBins.
%
% Model:  P(x|n) ~ prior(x) * prod_i rate_i(x)^{n_i} * exp(-tau * sum_i rate_i(x))
% Computed in log space, max-subtracted for numerical stability, then normalized.

nPos = size(rateMaps, 2);
nT   = size(countMat, 2);
if nargin < 4 || isempty(prior), prior = ones(1, nPos); end
prior = prior(:).';

logR  = log(rateMaps);                         % nCells x nPos
logL  = countMat.' * logR;                     % nT x nPos
logL  = logL - tau * sum(rateMaps, 1);         % subtract tau*sum_i rate_i(x)  (1 x nPos)
logL  = logL + log(prior);                     % add log prior (1 x nPos)
logL  = logL - max(logL, [], 2);               % stabilize per time bin
P     = exp(logL);
P     = P ./ sum(P, 2);                        % normalize over position

empty = sum(countMat, 1) == 0;                 % 1 x nT
P(empty, :) = repmat(prior / sum(prior), sum(empty), 1);

posterior = P.';                               % nPos x nT
[~, mapBin] = max(posterior, [], 1);

info = struct('mapBin', mapBin, 'emptyBins', empty, ...
              'nActiveCells', sum(any(countMat > 0, 2)), ...
              'nSpikeBins', sum(~empty));
end
