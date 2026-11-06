function R = rippleFeaturesStored(spkFile, varargin)
%RIPPLEFEATURESSTORED  Per-ripple recruitment and spike phase from a ripples
%   file that already stores each event's phase time series -- no LFP needed.
%
%   R = rippleFeaturesStored(spkFile, <options>)
%
%   For datasets (e.g. R:\HuszarR) whose <session>.ripples.events.mat has
%     ripples.timestamps   N x 2 [on off] (s)
%     ripples.phase_ts(k)  .ts (s) and .phase (rad) sampled across event k
%   Returns the same fields as RIPPLEFEATURES with wins = [] (whole recording),
%   so runRippleExportControl can write the same RIPPLES struct.
%
%   SPIKE PHASE: the event's stored phase, unwrapped and linearly interpolated
%   at the spike time, then wrapped to (-pi, pi]. Spikes outside the stored
%   samples of their event get NaN. The phase CONVENTION is whatever the
%   detector used (which angle is 0): resultant length / PPC do not depend on
%   it, but a preferred phase is only comparable across datasets if the
%   conventions match.
%
%   PER RIPPLE: Dur from timestamps; Freq = mean instantaneous frequency of
%   the stored phase (envelope weighting not possible here, so unweighted);
%   NActive / FracActive / NSpk as in rippleFeatures; PeakT from ripples.peaks
%   when stored. Amp / Pow are NaN (no LFP); AmpZ is taken from
%   ripples.peakNormedPower when stored, else NaN.
%
%   OPTIONS
%     'RateRange'  [0 Inf]  unit rate band over the recording ([0 Inf] = all units)
%     'Duration'   []       recording length (s); default: .lfp length from the
%                           .xml if both exist, else the last spike
%
%   See also RIPPLEFEATURES, RUNRIPPLEEXPORTCONTROL.

p = inputParser;
addParameter(p, 'RateRange', [0 Inf]);
addParameter(p, 'Duration', []);
parse(p, varargin{:});
o = p.Results;

R = struct('ok',false,'why','','file',spkFile,'lfpFile','','chan',NaN,'fs',NaN, ...
    'band',[NaN NaN],'nUnits',0,'unitIdx',[],'wins',[],'ev',table(), ...
    'epochT',NaN,'epochTnrem',NaN,'hasSleep',false, ...
    'envMu',NaN,'envSd',NaN,'nDetected',NaN,'nDroppedOverlap',NaN, ...
    'spkRip',[],'spkUnit',[],'spkT',[],'spkPh',[],'source','stored phase_ts');
folder = fileparts(spkFile);

rf = dir(fullfile(folder, '*.ripples.events.mat'));
if isempty(rf), R.why = 'no *.ripples.events.mat'; return; end
L = load(fullfile(rf(1).folder, rf(1).name), 'ripples');
rip = L.ripples;
if ~isfield(rip, 'phase_ts'), R.why = 'no ripples.phase_ts'; return; end
ts = double(rip.timestamps(:, 1:2));
if numel(rip.phase_ts) ~= size(ts, 1)
    R.why = sprintf('phase_ts has %d entries but timestamps %d rows', numel(rip.phase_ts), size(ts,1)); return
end
R.nDetected = size(ts, 1);

% ---- events: sort, drop malformed and overlapping (as rippleFeatures) ----------
src = (1:size(ts,1))';
okT = ts(:,2) > ts(:,1); ts = ts(okT, :); src = src(okT);
[ts, so] = sortrows(ts, 1); src = src(so);
keep = true(size(ts,1), 1); lastOff = -Inf;
for k = 1:size(ts,1)
    if ts(k,1) <= lastOff, keep(k) = false; else, lastOff = ts(k,2); end
end
ts = ts(keep, :); src = src(keep);
edges = reshape(ts', [], 1);
okE = [true; diff(edges) > 0];
if ~all(okE)                                   % touching events
    bad = unique(ceil(find(~okE)/2)); kr = true(size(ts,1),1); kr(bad) = false;
    ts = ts(kr, :); src = src(kr); edges = reshape(ts', [], 1);
end
R.nDroppedOverlap = R.nDetected - size(ts, 1);
on = ts(:,1); off = ts(:,2); nR = numel(on);

% ---- units -------------------------------------------------------------------
S = load(spkFile, 'spikes');
allTs = cellfun(@(x) x(:), S.spikes.times(:), 'UniformOutput', false);
dur = o.Duration;
if isempty(dur)
    [~, bn] = fileparts(rf(1).name); bn = erase(bn, '.ripples.events');
    lf = fullfile(folder, [bn '.lfp']); xf = fullfile(folder, [bn '.xml']);
    if exist(lf, 'file') && exist(xf, 'file')
        try
            xml = LoadXml(xf); fs = 1250;
            if isfield(xml, 'lfpSampleRate') && ~isempty(xml.lfpSampleRate), fs = xml.lfpSampleRate; end
            fi = dir(lf); dur = fi.bytes / (2 * xml.nChannels) / fs;
        catch
        end
    end
    if isempty(dur)
        last = cellfun(@(x) max([x; -Inf]), allTs); dur = max(last(isfinite(last)));
    end
end
R.wins = [0 dur]; R.epochT = dur;
rate = cellfun(@numel, allTs) / dur;
R.unitIdx = find(rate >= o.RateRange(1) & rate <= o.RateRange(2));
N = numel(R.unitIdx); R.nUnits = N;
if N == 0, R.why = 'no units in rate range'; return; end

% ---- spikes inside ripples ------------------------------------------------------
nAct = zeros(nR,1); nSpk = zeros(nR,1);
spkK = []; spkT = []; spkU = [];
for u = 1:N
    t = allTs{R.unitIdx(u)};
    b = discretize(t, edges);
    in = ~isnan(b) & mod(b, 2) == 1;
    k = (b(in) + 1) / 2;
    if isempty(k), continue; end
    nSpk = nSpk + accumarray(k, 1, [nR 1]);
    uk = unique(k); nAct(uk) = nAct(uk) + 1;
    spkK = [spkK; k]; spkT = [spkT; t(in)]; spkU = [spkU; repmat(u, numel(k), 1)]; %#ok<AGROW>
end

% ---- phase at spikes, frequency per event, from the stored phase series ----------
spkPh = nan(numel(spkK), 1); Freq = nan(nR, 1);
for k = 1:nR
    P = rip.phase_ts(src(k));
    tk = double(P.ts(:)); pk = unwrap(double(P.phase(:)));
    ok = isfinite(tk) & isfinite(pk); tk = tk(ok); pk = pk(ok);
    if numel(tk) < 2, continue; end
    [tk, ord] = unique(tk); pk = pk(ord);
    Freq(k) = mean(diff(pk) ./ diff(tk)) / (2*pi);
    sel = spkK == k;
    if any(sel)
        spkPh(sel) = angle(exp(1i * interp1(tk, pk, spkT(sel), 'linear', NaN)));
    end
end

% ---- optional stored fields --------------------------------------------------------
PeakT = nan(nR,1); AmpZ = nan(nR,1);
if isfield(rip, 'peaks') && numel(rip.peaks) == R.nDetected, PeakT = double(rip.peaks(src)); PeakT = PeakT(:); end
if isfield(rip, 'peakNormedPower') && numel(rip.peakNormedPower) == R.nDetected
    AmpZ = double(rip.peakNormedPower(src)); AmpZ = AmpZ(:);
end
nremFlag = nan(nR, 1);
ss = dir(fullfile(folder, '*.SleepState.states.mat'));
if ~isempty(ss)
    try
        Q = load(fullfile(ss(1).folder, ss(1).name), 'SleepState');
        nrem = Q.SleepState.ints.NREMstate; mid = (on + off) / 2;
        nremFlag = double(any(mid >= nrem(:,1)' & mid <= nrem(:,2)', 2));
        ov = min(nrem(:,2), dur) - max(nrem(:,1), 0); R.epochTnrem = sum(ov(ov > 0));
        R.hasSleep = true;
    catch
    end
end

R.ev = table(ones(nR,1), on, off - on, nan(nR,1), AmpZ, nan(nR,1), Freq, ...
    nAct, nAct / N, nSpk, nremFlag, PeakT, src, ...
    'VariableNames', {'Epoch','On','Dur','Amp','AmpZ','Pow','Freq', ...
                      'NActive','FracActive','NSpk','NREM','PeakT','SrcIdx'});
R.spkRip = spkK(:); R.spkUnit = spkU(:); R.spkT = spkT(:); R.spkPh = spkPh(:);
R.ok = true;
end
