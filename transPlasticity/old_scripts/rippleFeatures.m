function R = rippleFeatures(spkFile, wins, varargin)
%RIPPLEFEATURES  Per-ripple LFP features, unit recruitment and spike phase for
%   one session, restricted to the epochs in WINS.
%
%   R = rippleFeatures(spkFile, wins, <options>)
%
%   spkFile : the session's *.spikes.cellinfo.mat (its folder must hold the
%             *.ripples.events.mat, the .lfp and the .xml)
%   wins    : K x 2 [start stop] epochs in s. A ripple belongs to epoch k if
%             it starts and ends inside it; ripples outside every epoch are
%             dropped. [] = the whole recording (length of the .lfp file).
%
%   OPTIONS
%     'Band'      ripple band (Hz) for filtering; default = the detector's
%                 detectionparms.passband if stored, else [130 200]
%     'RateRange' [0.1 10] Hz; units kept if their rate over 'UnitWin' is in range
%     'UnitWin'   M x 2 windows for the rate filter (default = wins), so pre and
%                 post epochs of one session use the SAME units
%     'ChunkSec'  600   LFP processed in chunks of this length
%     'Pad'       1     s of padding on each side of a chunk (filter edges)
%
%   LFP: detection channel (detectionparms.channel, 0-based) read with
%   LoadBinary at the .xml lfpSampleRate (1250 Hz if absent). Zero-phase
%   3rd-order Butterworth band-pass, Hilbert transform.
%
%   PER RIPPLE (R.ev, one row per ripple)
%     Epoch, On, Dur (s)
%     Amp        peak ripple-band envelope (raw LFP units)
%     AmpZ       Amp z-scored to the envelope over ALL samples of the session's
%                epochs (same normalisation for every epoch of the session)
%     Pow        mean squared envelope over the event
%     Freq       envelope^2-weighted mean instantaneous frequency (Hz) from the
%                unwrapped Hilbert phase
%     NActive    units with >= 1 spike in [On, On+Dur]
%     FracActive NActive / number of kept units
%     NSpk       spikes of kept units in the event
%     NREM       1/0 if a *.SleepState.states.mat exists (event midpoint in
%                NREMstate), NaN otherwise
%     PeakT      time of the ripple-band envelope peak (s)
%     SrcIdx     row of this event in ripples.timestamps
%   Events with offset <= onset, or overlapping an earlier event, are dropped
%   (R.nDetected, R.nDroppedOverlap).
%
%   PER SPIKE (every in-ripple spike of the kept units): R.spkRip (row of
%   R.ev), R.spkUnit (index into R.unitIdx), R.spkT (s), R.spkPh (rad).
%
%   SPIKE PHASE: phase of the ripple-band LFP (Hilbert angle, 0 = filtered
%   peak) at every spike of a kept unit inside a ripple. Stored as sufficient
%   statistics per epoch x unit: R.unitN, R.unitC (sum cos), R.unitS (sum sin).
%   Pairwise phase consistency follows as (C^2+S^2-n)/(n(n-1)), which is
%   unbiased by spike count (the resultant length is not).
%
%   PER EPOCH: R.epochT (s), R.epochTnrem (s, NaN without sleep scoring).

p = inputParser;
addParameter(p, 'Band', []);
addParameter(p, 'RateRange', [0.1 10]);
addParameter(p, 'UnitWin', []);
addParameter(p, 'ChunkSec', 600);
addParameter(p, 'Pad', 1);
addParameter(p, 'Channel', []);   % 1-based LFP channel; overrides the detector's
parse(p, varargin{:});
o = p.Results;

nW = size(wins, 1);
R = struct('ok',false,'why','','file',spkFile,'lfpFile','','chan',NaN,'fs',NaN, ...
    'band',[NaN NaN],'nUnits',0,'unitIdx',[],'wins',wins,'ev',table(), ...
    'epochT',diff(wins,1,2),'epochTnrem',nan(nW,1),'hasSleep',false, ...
    'unitN',[],'unitC',[],'unitS',[],'phaseEdges',[],'unitHist',[],'unitHistC',[], ...
    'nDetected',NaN,'nDroppedOverlap',NaN,'spkRip',[],'spkUnit',[],'spkT',[],'spkPh',[], ...
    'envMu',NaN,'envSd',NaN);
folder = fileparts(spkFile);

% ---- files -----------------------------------------------------------------
rf = dir(fullfile(folder, '*.ripples.events.mat'));
if isempty(rf), R.why = 'no *.ripples.events.mat'; return; end
pick = find(strcmp({rf.name}, 'amplifier_analogin_auxiliary_int16.ripples.events.mat'), 1);
if isempty(pick), pick = 1; end
L = load(fullfile(folder, rf(pick).name), 'ripples');
rip = L.ripples;
% detector parameters: buzcode-style ripples.detectorinfo.detectionparms when
% present; other datasets may store less, so every field is optional here
dp = struct();
if isfield(rip, 'detectorinfo') && isfield(rip.detectorinfo, 'detectionparms')
    dp = rip.detectorinfo.detectionparms;
end
if ~isempty(o.Channel)
    R.chan = o.Channel;
elseif isfield(dp, 'channel') && ~isempty(dp.channel)
    R.chan = dp.channel(1) + 1;                % stored 0-based
elseif isfield(dp, 'Channels') && ~isempty(dp.Channels)
    R.chan = dp.Channels(1) + 1;
elseif isfield(rip, 'detectorinfo') && isfield(rip.detectorinfo, 'detectionchannel') && ~isempty(rip.detectorinfo.detectionchannel)
    R.chan = rip.detectorinfo.detectionchannel(1) + 1;
else
    R.why = 'no detection channel stored: pass ''Channel'' (1-based)'; return
end
bnRaw = '';
if isfield(dp, 'basename'), bnRaw = dp.basename; end   % may be cell / string / empty
if iscell(bnRaw), if isempty(bnRaw), bnRaw = ''; else, bnRaw = bnRaw{1}; end, end
bnRaw = char(bnRaw);
if isempty(bnRaw)                              % fall back on the events file's own basename
    bnRaw = erase(rf(pick).name, '.ripples.events.mat');
end
[~, bn] = fileparts(bnRaw);
R.lfpFile = fullfile(folder, [bn '.lfp']);
if ~exist(R.lfpFile, 'file') && exist(fullfile(folder, [bn '.eeg']), 'file')
    R.lfpFile = fullfile(folder, [bn '.eeg']);   % older naming
end
xmlFile = fullfile(folder, [bn '.xml']);
if ~exist(R.lfpFile, 'file'), R.why = ['missing ' R.lfpFile]; return; end
if ~exist(xmlFile, 'file'),   R.why = ['missing ' xmlFile];   return; end
xml = LoadXml(xmlFile);
nCh = xml.nChannels;
R.fs = 1250;
if isfield(xml, 'lfpSampleRate') && ~isempty(xml.lfpSampleRate), R.fs = xml.lfpSampleRate; end
fs = R.fs;
if ~isempty(o.Band), R.band = o.Band;
elseif isfield(dp, 'passband') && numel(dp.passband) == 2, R.band = dp.passband(:)';
else, R.band = [130 200];
end
[bf, af] = butter(3, R.band/(fs/2), 'bandpass');

% wins = [] -> the whole recording, from the length of the .lfp file
if isempty(wins)
    fi = dir(R.lfpFile);
    wins = [0, fi.bytes / (2 * nCh) / fs];
    nW = 1; R.wins = wins; R.epochT = wins(2); R.epochTnrem = NaN;
end

% ---- events: sort, drop malformed and overlapping, assign epochs ------------
ts = double(rip.timestamps(:, 1:2));
src = (1:size(ts,1))';                         % row in ripples.timestamps, tracked through every filter
okT = ts(:,2) > ts(:,1);
ts = ts(okT, :); src = src(okT);
[ts, so] = sortrows(ts, 1); src = src(so);
keep = true(size(ts,1), 1); lastOff = -Inf;
for k = 1:size(ts,1)
    if ts(k,1) <= lastOff, keep(k) = false; else, lastOff = ts(k,2); end
end
ts = ts(keep, :); src = src(keep);
R.nDetected = size(rip.timestamps, 1); R.nDroppedOverlap = sum(~keep) + sum(~okT);
ep = zeros(size(ts,1), 1);
for w = 1:nW
    ep(ts(:,1) >= wins(w,1) & ts(:,2) <= wins(w,2)) = w;
end
ts = ts(ep > 0, :); src = src(ep > 0); ep = ep(ep > 0);
nR = size(ts, 1);
on = ts(:,1); off = ts(:,2);

% ---- units: fixed set per session ------------------------------------------
S = load(spkFile, 'spikes');
allTs = S.spikes.times(:);
uw = o.UnitWin; if isempty(uw), uw = wins; end
Tuw = sum(diff(uw, 1, 2));
rate = cellfun(@(t) sum(arrayfun(@(j) sum(t >= uw(j,1) & t < uw(j,2)), 1:size(uw,1))), allTs) / Tuw;
R.unitIdx = find(rate >= o.RateRange(1) & rate <= o.RateRange(2));
N = numel(R.unitIdx); R.nUnits = N;
if N == 0, R.why = 'no units in rate range'; return; end

% spikes inside ripples: discretize against interleaved [on off] edges
nAct = zeros(nR, 1); nSpk = zeros(nR, 1);
spkK = []; spkT = []; spkU = [];
if nR > 0
    edges = reshape([on off]', [], 1);
    okE = [true; diff(edges) > 0];             % ties between touching events
    if ~all(okE)
        bad = unique(ceil(find(~okE)/2));
        keepR = true(nR,1); keepR(bad) = false;
        on = on(keepR); off = off(keepR); ep = ep(keepR); src = src(keepR); nR = numel(on);
        R.nDroppedOverlap = R.nDroppedOverlap + sum(~keepR);
        nAct = zeros(nR,1); nSpk = zeros(nR,1);
        edges = reshape([on off]', [], 1);
    end
    for u = 1:N
        t = allTs{R.unitIdx(u)}(:);
        b = discretize(t, edges);
        in = ~isnan(b) & mod(b, 2) == 1;
        k = (b(in) + 1) / 2;
        if isempty(k), continue; end
        nSpk = nSpk + accumarray(k, 1, [nR 1]);
        nAct(unique(k)) = nAct(unique(k)) + 1;
        spkK = [spkK; k]; spkT = [spkT; t(in)]; spkU = [spkU; repmat(u, numel(k), 1)]; %#ok<AGROW>
    end
end
spkPh = nan(numel(spkK), 1);

% ---- sleep state (optional) -------------------------------------------------
nremFlag = nan(nR, 1);
ss = dir(fullfile(folder, '*.SleepState.states.mat'));
if ~isempty(ss)
    try
        Q = load(fullfile(folder, ss(1).name), 'SleepState');
        nrem = Q.SleepState.ints.NREMstate;
        mid = (on + off) / 2;
        nremFlag = double(any(mid >= nrem(:,1)' & mid <= nrem(:,2)', 2));
        for w = 1:nW
            ov = min(nrem(:,2), wins(w,2)) - max(nrem(:,1), wins(w,1));
            R.epochTnrem(w) = sum(ov(ov > 0));
        end
        R.hasSleep = true;
    catch ME
        warning('rippleFeatures:sleep', '%s: %s', ss(1).name, ME.message);
    end
end

% ---- LFP, chunk by chunk ----------------------------------------------------
Amp = nan(nR,1); Pow = nan(nR,1); Freq = nan(nR,1); PeakT = nan(nR,1);
sE = 0; sE2 = 0; nE = 0;
for w = 1:nW
    for a = wins(w,1):o.ChunkSec:wins(w,2)
        b = min(a + o.ChunkSec, wins(w,2));
        if b - a < 1, continue; end
        st = max(a - o.Pad, 0);
        x = LoadBinary(R.lfpFile, 'nChannels', nCh, 'channels', R.chan, ...
            'frequency', fs, 'start', st, 'duration', b + o.Pad - st);
        x = double(x(:));
        if numel(x) < 3*max(numel(bf), numel(af)), continue; end
        h = hilbert(filtfilt(bf, af, x));
        env = abs(h); ph = angle(h);
        i0 = max(1, round((a - st)*fs) + 1); i1 = min(numel(env), round((b - st)*fs));
        sE = sE + sum(env(i0:i1)); sE2 = sE2 + sum(env(i0:i1).^2); nE = nE + (i1 - i0 + 1);

        kc = find(ep == w & on >= a & on < b);
        for k = kc'
            j1 = max(1, floor((on(k) - st)*fs) + 1);
            j2 = min(numel(env), ceil((off(k) - st)*fs) + 1);
            if j2 - j1 < 2, continue; end
            e = env(j1:j2);
            [Amp(k), im] = max(e); Pow(k) = mean(e.^2);
            PeakT(k) = st + (j1 + im - 2) / fs;       % time of the envelope peak
            f = diff(unwrap(ph(j1:j2))) * fs / (2*pi);
            wt = e(1:end-1).^2;
            Freq(k) = sum(f .* wt) / sum(wt);
        end
        if ~isempty(kc) && ~isempty(spkK)
            sel = ismember(spkK, kc);
            idx = round((spkT(sel) - st)*fs) + 1;
            idx = min(max(idx, 1), numel(ph));
            spkPh(sel) = ph(idx);
        end
    end
end
R.envMu = sE / max(nE, 1);
R.envSd = sqrt(max(sE2 / max(nE, 1) - R.envMu^2, 0));

% ---- spike-phase sufficient statistics per epoch x unit ---------------------
R.unitN = zeros(nW, N); R.unitC = zeros(nW, N); R.unitS = zeros(nW, N);
if ~isempty(spkK)
    okP = isfinite(spkPh);
    eK = ep(spkK(okP)); uK = spkU(okP); thK = spkPh(okP);
    R.unitN = accumarray([eK uK], 1, [nW N]);
    R.unitC = accumarray([eK uK], cos(thK), [nW N]);
    R.unitS = accumarray([eK uK], sin(thK), [nW N]);
end

% ---- phase histograms per epoch x unit (36 bins of 10 deg) -----------------
% unitHist  : raw phase, 0 = ripple-band peak
% unitHistC : phase relative to that unit's own circular mean in that epoch,
%             so the WIDTH of each unit's distribution can be compared across
%             epochs independently of any shift in its preferred phase
nb = 36;
R.phaseEdges = linspace(-pi, pi, nb + 1);
R.unitHist = zeros(nW, N, nb); R.unitHistC = zeros(nW, N, nb);
if ~isempty(spkK) && any(okP)
    b = discretize(thK, R.phaseEdges);
    R.unitHist = accumarray([eK uK b], 1, [nW N nb]);
    mu = atan2(R.unitS, R.unitC);
    muK = reshape(mu(sub2ind([nW N], eK, uK)), [], 1);   % column: with nW = 1, mu is a row
    d = angle(exp(1i * (thK - muK)));                      % and mu(idx) would be a row too
    R.unitHistC = accumarray([eK uK discretize(d, R.phaseEdges)], 1, [nW N nb]);
end

R.ev = table(ep, on, off - on, Amp, (Amp - R.envMu) / R.envSd, Pow, Freq, ...
    nAct, nAct / N, nSpk, nremFlag, PeakT, src, ...
    'VariableNames', {'Epoch','On','Dur','Amp','AmpZ','Pow','Freq', ...
                      'NActive','FracActive','NSpk','NREM','PeakT','SrcIdx'});
% every in-ripple spike of the kept units: ripple row in R.ev, unit (index into
% R.unitIdx), time (s, recording clock), ripple-band phase (0 = filtered peak)
R.spkRip = spkK(:); R.spkUnit = spkU(:); R.spkT = spkT(:); R.spkPh = spkPh(:);
R.ok = true;
end
