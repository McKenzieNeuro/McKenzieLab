function [dA, dB, nA, nB] = splitEpochDuty(r, width, step, blockSec)
%SPLITEPOCHDUTY  Two duty-cycle estimates from disjoint blocks of one epoch.
%
%   [dA, dB, nA, nB] = splitEpochDuty(r, width, step, blockSec)
%
%   r        : one session struct from syncBatch (needs .tc and .transients)
%   width,step : detection parameters used to produce r
%   blockSec : length of the interleaved blocks (default 600 s)
%
%   WHY. When the same baseline epoch is used both to SELECT animals (large
%   acute change = post-infusion minus this baseline) and as the REFERENCE
%   the persistence step is measured against, the selection conditions on the
%   baseline's noise and regression to the mean manufactures a persistence
%   effect. Splitting the epoch into alternating blocks gives two estimates
%   from non-overlapping minutes: use A for selection, B for the test, and
%   the shared-noise path is cut without discarding any animal.
%
%   HOW. Detection is NOT re-run -- the transients and the window grid are
%   taken as computed over the whole epoch, which is what you want: the
%   surrogate band and cluster-significance threshold stay based on the full
%   epoch rather than on a short fragment with its own weaker threshold.
%   Only the ACCOUNTING is split: each window is assigned to block
%   floor(t/blockSec), even blocks to A and odd to B, and duty is the
%   fraction of that block set's windows falling inside a significant
%   transient.
%
%   CHOOSING blockSec. Blocks must be long relative to a transient, or the
%   same transient lands in both A and B and the two estimates share exactly
%   the signal you were trying to separate. Transients here run several
%   hundred seconds, so the default 600 s is a floor, not a safe margin --
%   with a 3600 s epoch that gives only 3 blocks per side. This trades
%   independence against the number of blocks and there is no setting that
%   makes both comfortable at this epoch length. Check nA/nB: if either is
%   tiny the estimate is noise.
%
%   Returns NaN when a side has no windows.
%
%   See also TRANSIENTSTATS, SYNCBATCH.

if nargin < 4, blockSec = 600; end

dA = NaN; dB = NaN; nA = 0; nB = 0;
if isempty(r.tc), return; end

tc = r.tc(:);
t0 = min(tc);
blk = floor((tc - t0) / blockSec);
isA = mod(blk,2) == 0;

inTrans = false(size(tc));
if ~isempty(r.transients)
    sig = r.transients([r.transients.significant]);
    for k = 1:numel(sig)
        inTrans = inTrans | (tc >= sig(k).tStart & tc <= sig(k).tEnd);
    end
end

valid = ~isnan(r.Ct(:));            % windows that were actually evaluated
nA = sum(isA & valid); nB = sum(~isA & valid);
if nA > 0, dA = sum(inTrans & isA & valid) / nA; end
if nB > 0, dB = sum(inTrans & ~isA & valid) / nB; end
end
