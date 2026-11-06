function tag = optoTagUnits(ts, pulseTimes, varargin)
%OPTOTAGUNITS  Classify units as light-responsive from pulse-locked spiking.
%
%   tag = optoTagUnits(ts, pulseTimes, ...)
%
%   ts         : {N x 1} cell of spike times (s), one per unit
%   pulseTimes : K x 2 [on off], from pulseInfo.time
%
%   CRITERION (stated explicitly -- these are defaults, not a standard;
%   check against your lab's usual protocol before trusting the tag):
%     - response window [0, LatWin] s after pulse ONSET
%     - a unit is tagged if BOTH:
%         (a) response probability (fraction of pulses with >=1 spike in
%             the window) exceeds MinProb, AND
%         (b) that probability is higher than chance, assessed by comparing
%             to NPerm response probabilities computed after jittering each
%             pulse onset by a random offset drawn uniformly from the
%             inter-pulse interval (so the null preserves each unit's own
%             firing rate and any slow rate covariation with the pulse
%             train, and tests only fine-timing locking, not "did rate
%             happen to be higher near pulses").
%     - mean first-spike latency and its jitter (SD) among responding
%       pulses are reported but NOT used to gate the tag by default (set
%       'MaxLatency'/'MaxJitter' to add that as a third gate; useful for
%       excluding likely network-mediated / indirect responses, which tend
%       to have longer and more variable latency than direct opsin-driven
%       spiking).
%
%   NAME/VALUE
%     'LatWin'     response window length, s (default 0.010 -- 10 ms)
%     'MinProb'    minimum response probability to even consider (default 0.1)
%     'NPerm'      jitter permutations for the null (default 1000)
%     'Alpha'      significance level (default 0.01; opto-tagging typically
%                  uses a stricter threshold than 0.05 given many units are
%                  tested per session -- multiplicity is on you to correct
%                  further if desired, e.g. Bonferroni by numel(ts))
%     'MaxLatency' if set, an additional gate: mean latency must be <= this (s)
%     'MaxJitter'  if set, an additional gate: latency SD must be <= this (s)
%     'ExcludeArtifact' s to blank immediately after pulse onset before
%                  counting a spike as a response, for light-artifact
%                  contamination on some systems (default 0)
%
%   OUTPUT: struct array, one element per unit, fields:
%     .responseProb, .pValue, .isTagged, .meanLatency, .latencyJitter,
%     .nPulses, .nResponses
%
%   See also OPTOPULSEINTERVALS.

p = inputParser;
p.addParameter('LatWin', 0.010);
p.addParameter('MinProb', 0.1);
p.addParameter('NPerm', 1000);
p.addParameter('Alpha', 0.01);
p.addParameter('MaxLatency', []);
p.addParameter('MaxJitter', []);
p.addParameter('ExcludeArtifact', 0);
p.parse(varargin{:});
o = p.Results;

N = numel(ts);
pOn = pulseTimes(:,1);
K = numel(pOn);
isi = diff(pOn);
if isempty(isi), isi = 5; end   % fallback if only one pulse

tag = repmat(struct('responseProb',NaN,'pValue',NaN,'isTagged',false, ...
    'meanLatency',NaN,'latencyJitter',NaN,'nPulses',K,'nResponses',0), N, 1);

if K == 0
    warning('optoTagUnits:noPulses','no pulses given; returning all untagged.');
    return
end

for n = 1:N
    x = sort(ts{n}(:));
    if isempty(x), continue; end

    [prob, lat] = local_responseStats(x, pOn, o.LatWin, o.ExcludeArtifact);
    nResp = round(prob * K);

    tag(n).responseProb = prob;
    tag(n).nResponses = nResp;
    if nResp > 0
        tag(n).meanLatency = mean(lat(~isnan(lat)));
        tag(n).latencyJitter = std(lat(~isnan(lat)));
    end

    if prob < o.MinProb, tag(n).pValue = 1; continue; end

    % jitter null: shift each pulse onset by U(0, local ISI), preserving
    % rate structure, breaking fine-timing lock to the pulse
    nullProb = zeros(o.NPerm,1);
    for r = 1:o.NPerm
        jitter = [isi; isi(end)] .* rand(K,1);
        pOnJ = pOn + jitter;
        nullProb(r) = local_responseStats(x, pOnJ, o.LatWin, o.ExcludeArtifact);
    end
    tag(n).pValue = (1 + sum(nullProb >= prob)) / (1 + o.NPerm);

    gate = tag(n).pValue < o.Alpha;
    if ~isempty(o.MaxLatency), gate = gate && tag(n).meanLatency <= o.MaxLatency; end
    if ~isempty(o.MaxJitter),  gate = gate && tag(n).latencyJitter <= o.MaxJitter; end
    tag(n).isTagged = gate;
end
end


function [prob, lat] = local_responseStats(x, pOn, latWin, artifactExclude)
K = numel(pOn);
lat = nan(K,1);
for k = 1:K
    w0 = pOn(k) + artifactExclude;
    w1 = pOn(k) + latWin;
    j = find(x >= w0 & x < w1, 1, 'first');
    if ~isempty(j), lat(k) = x(j) - pOn(k); end
end
prob = mean(~isnan(lat));
end
