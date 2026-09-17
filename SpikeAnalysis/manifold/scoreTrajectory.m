function S = scoreTrajectory(posterior, binCenters, tau, opts)
% scoreTrajectory  Quantify sequential structure in a decoded posterior.
%
%   S = scoreTrajectory(posterior, binCenters, tau, opts)
%
% INPUTS
%   posterior  : [nPosBins x nTimeBins], columns sum to 1 (from bayesianDecode).
%   binCenters : [1 x nPosBins] position bin centers (cm).
%   tau        : decoding bin width (s).
%   opts       : struct (optional):
%       .computeLine : also compute radon/line score   (default true)
%       .lineBand    : half-band width for line score  (default 15 cm)
%       .nSlope      : number of slopes searched        (default 120)
%
% OUTPUT struct S:
%   .weightedCorr : weighted Pearson r between time & position (Grosmark &
%                   Buzsaki 2016; Silva et al. 2015). Sign = decoded direction.
%   .lineScore    : max mean posterior within a band around the best line
%                   (Davidson et al. 2009; radon transform of the posterior),
%                   in [0,1]. NaN if computeLine=false.
%   .lineSlope    : best-line slope (cm/s); sign = direction.  NaN if skipped.
%   .lineIntercept: best-line position at first time bin (cm). NaN if skipped.

if nargin < 4, opts = struct(); end
if ~isfield(opts,'computeLine') || isempty(opts.computeLine), opts.computeLine = true; end
if ~isfield(opts,'lineBand')    || isempty(opts.lineBand),    opts.lineBand = 15;    end
if ~isfield(opts,'nSlope')      || isempty(opts.nSlope),      opts.nSlope   = 120;   end

[nPos, nT] = size(posterior);
x = binCenters(:);                         % nPos x 1 (cm)
t = ((0:nT-1) + 0.5) * tau;                % 1 x nT (s)

% ---- weighted correlation ----
Tg = repmat(t, nPos, 1);
Xg = repmat(x, 1, nT);
W  = sum(posterior(:));
mt = sum(sum(posterior .* Tg)) / W;
mx = sum(sum(posterior .* Xg)) / W;
cov_tx = sum(sum(posterior .* (Tg - mt) .* (Xg - mx))) / W;
var_t  = sum(sum(posterior .* (Tg - mt).^2)) / W;
var_x  = sum(sum(posterior .* (Xg - mx).^2)) / W;
S.weightedCorr = cov_tx / sqrt(var_t * var_x + eps);

% ---- radon / line score ----
S.lineScore = NaN; S.lineSlope = NaN; S.lineIntercept = NaN;
if opts.computeLine
    binSize = median(diff(binCenters));
    dbin = max(1, round(opts.lineBand / binSize));
    j = 0:nT-1;
    aMax = (nPos - 1) / max(1, (nT - 1));          % full track coverage, both signs
    slopes = linspace(-aMax, aMax, opts.nSlope);   % position-bins per time-bin
    csum = [zeros(1, nT); cumsum(posterior, 1)];   % (nPos+1) x nT, for fast band sums
    best = -Inf; bestA = NaN; bestB = NaN;
    for a = slopes
        for b = 1:nPos
            k  = round(a * j + b);                 % 1 x nT line position (bins)
            lo = max(1, k - dbin);
            hi = min(nPos, k + dbin);
            valid = hi >= lo;
            bandSum = zeros(1, nT);
            cols = find(valid);
            for jj = cols
                bandSum(jj) = csum(hi(jj)+1, jj) - csum(lo(jj), jj);
            end
            score = sum(bandSum) / nT;
            if score > best
                best = score; bestA = a; bestB = b;
            end
        end
    end
    S.lineScore     = best;
    S.lineSlope     = bestA * binSize / tau;                 % cm/s
    S.lineIntercept = binCenters(min(max(round(bestB),1),nPos));
    S.lineParams    = [bestA, bestB];                        % [bins/bin, startBin]
end
end
