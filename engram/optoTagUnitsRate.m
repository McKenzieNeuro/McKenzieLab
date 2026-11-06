function tag = optoTagUnitsRate(ts, pulseTimes, dt, varargin)
%OPTOTAGUNITSRATE  Tag units by paired baseline-vs-stim rate modulation.
%
%   tag = optoTagUnitsRate(ts, pulseTimes, dt, ...)
%
%   Replaces the short-latency/jitter-permutation criterion in
%   optoTagUnits.m with sm_getOptoMod: a paired Poisson GLMM
%   (Count ~ Period + (1|Pair)) comparing baseline counts to stim counts,
%   pair-matched by pulse so pulse-to-pulse rate drift is absorbed by the
%   random intercept rather than inflating the Period test.
%
%   ts         : {N x 1} cell of spike times (s)
%   pulseTimes : K x 2 [on off]
%   dt         : bin duration (s) for BOTH baseline and stim windows
%
%   PER-PULSE PAIRING. For pulse k:
%     baseline bin = [on(k) - PreOffset - dt,  on(k) - PreOffset]
%     stim bin     = [on(k) + PostOffset,      on(k) + PostOffset + dt]
%   binnedPop_base(unit,k) and binnedPop(unit,k) are the resulting spike
%   COUNTS. sm_getOptoMod fits counts directly (Poisson), so no round trip
%   through rate is needed here -- counts are taken straight from the bins.
%
%   NAME/VALUE
%     'PreOffset'  s before pulse onset where the baseline bin ENDS
%                  (default 0 -- baseline bin is [on-dt, on]). Increase if
%                  a bin flush against onset risks anticipatory/ramping
%                  contamination.
%     'PostOffset' s after pulse onset where the stim bin STARTS
%                  (default 0 -- stim bin is [on, on+dt]).
%     'Paired'     passed to sm_getOptoMod (default true)
%     'Alpha'      significance threshold on the p-value (default 0.01)
%     'MinBeta'    minimum |beta_period| to call tagged on top of
%                  significance (default 0)
%     'RequireIncrease' if true (default), only a positive beta (higher
%                  rate during stim) counts as tagged; a unit significantly
%                  SUPPRESSED by light is flagged (.isSuppressed) but not
%                  tagged. Set false to tag on significant modulation in
%                  either direction.
%
%   OUTPUT: struct array, one per unit:
%     .p (opto_mod), .beta (opto_sign), .isTagged, .isSuppressed,
%     .nPulsesUsed, .meanBaseline, .meanStim  (Hz, for reference)
%
%   See also SM_GETOPTOMOD, TAGGEDSYNCHRONY.

p = inputParser;
p.addParameter('PreOffset', 0);
p.addParameter('PostOffset', 0);
p.addParameter('Paired', true, @islogical);
p.addParameter('Alpha', 0.01);
p.addParameter('MinBeta', 0);
p.addParameter('RequireIncrease', true, @islogical);
p.parse(varargin{:});
o = p.Results;

N = numel(ts);
pOn = pulseTimes(:,1);
K = numel(pOn);

tag = repmat(struct('p',NaN,'beta',NaN,'isTagged',false,'isSuppressed',false, ...
    'nPulsesUsed',0,'meanBaseline',NaN,'meanStim',NaN), N, 1);

if K < 5
    warning('optoTagUnitsRate:fewPulses', ...
        'only %d pulses -- paired GLMM will be unreliable with this few pairs.', K);
end

baseOn  = pOn - o.PreOffset - dt;
baseOff = pOn - o.PreOffset;
stimOn  = pOn + o.PostOffset;
stimOff = pOn + o.PostOffset + dt;

for n = 1:N
    x = sort(ts{n}(:));
    if isempty(x), continue; end

    cBase = local_binCounts(x, baseOn, baseOff);
    cStim = local_binCounts(x, stimOn, stimOff);

    good = isfinite(cBase) & isfinite(cStim);
    tag(n).nPulsesUsed = sum(good);
    if sum(good) < 5, continue; end

    tag(n).meanBaseline = mean(cBase(good)) / dt;
    tag(n).meanStim     = mean(cStim(good)) / dt;

    try
        [pv, beta] = sm_getOptoMod(cBase(good), cStim(good), 'paired', o.Paired);
    catch ME
        warning('optoTagUnitsRate:fit','unit %d: %s', n, ME.message);
        continue
    end
    tag(n).p = pv; tag(n).beta = beta;

    sig = pv < o.Alpha;
    if sig && beta > 0 && beta >= o.MinBeta
        tag(n).isTagged = true;
    elseif sig && beta < 0
        tag(n).isSuppressed = true;
        if ~o.RequireIncrease && abs(beta) >= o.MinBeta
            tag(n).isTagged = true;
        end
    end
end
end


function c = local_binCounts(x, onT, offT)
K = numel(onT);
c = nan(K,1);
for k = 1:K
    c(k) = sum(x >= onT(k) & x < offT(k));
end
end
