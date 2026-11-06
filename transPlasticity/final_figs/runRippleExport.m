%% runRippleExport.m
% For every session: every detected ripple over the WHOLE recording with its
% features, unit recruitment and the ripple-band phase of every in-ripple
% spike of EVERY unit (no rate exclusion), saved next to the raw data as
% rippleEvents.mat (variable RIPPLES). Same time conventions as the synchrony
% exports:  ts1 = s from the start of the recording,  ts2 = s from the
% animal's first ASV infusion.
%
% RIPPLES
%   .ts1, .ts2           ripple onset          (nRip x 1)
%   .ts1Peak, .ts2Peak   ripple-band envelope peak
%   .ts1Off,  .ts2Off    ripple offset
%   .ts3, .ts3Peak, .ts3Off   the same, s from the animal's first SAL infusion
%   .ev     table, one row per ripple:
%       Dur (s), Amp (peak ripple-band envelope, raw LFP units),
%       AmpZ (Amp z-scored to the envelope over the whole recording),
%       Pow (mean squared envelope), Freq (envelope^2-weighted instantaneous
%       frequency, Hz), NActive (units with >= 1 spike), FracActive (NActive /
%       all units), NSpk, NREM (1/0, NaN without a SleepState file),
%       SrcIdx (row in ripples.timestamps)
%   .count  nRip x nUnits sparse spike counts per ripple and unit
%   .spikes table, one row per in-ripple spike:
%       Ripple (row of .ev), Unit (column of .count), ts1, ts2, Phase
%       (ripple-band Hilbert phase, rad, 0 = filtered peak)
%   .units  table, one row per unit: SpikesIdx (index into spikes.times),
%       UID (if stored), Rate (Hz, whole recording), nRipSpk, PrefPhase,
%       ResultantR, PPC (pairwise phase consistency over all its in-ripple
%       spikes; NaN below 2 spikes)
%   .meta   files, channel, band, detector parameters, drops, ts2 mode, ...
%
% Ripples whose offset is not after onset, or that overlap an earlier ripple,
% are dropped (meta.nDetected, meta.nDropped); .ev.SrcIdx maps back.

CACHE_ROOT = 'R:\TransPlasticity\syncCache';
ROOT = 'R:\TransPlasticity';
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
EXPOSURE_COND = 'ASV';
RF_OPTS = {};               % e.g. {'Band',[120 200]}; default = detector passband, else 130-200 Hz
OVERWRITE = true;
LIMIT = Inf;
OUT_NAME = 'rippleEvents.mat';

fils = getAllExtFiles(ROOT, 'mat', 1);
fils = fils(contains(fils, 'spikes.cellinfo.mat'));
cnd = cellfun(@(f) conditionFromPath(f, 'TransPlasticity'), fils, 'UniformOutput', false);
ani = cellfun(@(f) animalFromPath(f, 'TransPlasticity'),    fils, 'UniformOutput', false);
dnm = cellfun(@(f) sessionDateFromPath(f),                  fils);
keep = ismember(cnd, CONDS) & ~cellfun(@isempty, ani) & isfinite(dnm);
fils = fils(keep); cnd = cnd(keep); ani = ani(keep); dnm = dnm(keep);
inf_t = nan(numel(fils), 1);
for i = 1:numel(fils), inf_t(i) = infusionTime(fils{i}); end
recStart = cellfun(@recordingStartFromPath, fils);
ref = firstASVReference(fils, ani, cnd, dnm, inf_t, recStart, EXPOSURE_COND);
refSAL = firstASVReference(fils, ani, cnd, dnm, inf_t, recStart, 'SAL');   % ts3: first vehicle infusion
fprintf('%d sessions\n', numel(fils));

status = cell(numel(fils), 1); nDone = 0;
for i = 1:numel(fils)
    if nDone >= LIMIT, break; end
    f = fils{i};
    out = fullfile(fileparts(f), OUT_NAME);
    if ~OVERWRITE && exist(out, 'file')
        status{i} = 'exists'; fprintf('[%3d/%3d] exists  %s\n', i, numel(fils), f); continue
    end
    tic;
    try
        R = rippleFeatures(f, [], 'RateRange', [0 Inf], RF_OPTS{:});
    catch ME
        status{i} = ['error: ' ME.message];
        fprintf('[%3d/%3d] ERROR   %s: %s\n', i, numel(fils), f, ME.message); continue
    end
    if ~R.ok
        status{i} = ['skip: ' R.why]; fprintf('[%3d/%3d] skip    %s: %s\n', i, numel(fils), f, R.why); continue
    end
    tsF = @(t) ts2FromReference(t, ani{i}, dnm(i), recStart(i), ref);
    tsF3 = @(t) ts2FromReference(t, ani{i}, dnm(i), recStart(i), refSAL);
    ev = R.ev;
    nRip = height(ev); N = R.nUnits;

    RP = struct();
    RP.ts1 = ev.On;                  RP.ts2 = tsF(ev.On);
    RP.ts1Peak = ev.PeakT;           RP.ts2Peak = tsF(ev.PeakT);
    RP.ts1Off = ev.On + ev.Dur;      RP.ts2Off = tsF(ev.On + ev.Dur);
    RP.ts3 = tsF3(ev.On); RP.ts3Peak = tsF3(ev.PeakT); RP.ts3Off = tsF3(ev.On + ev.Dur);
    RP.ev = removevars(ev, {'Epoch','On','PeakT'});

    RP.count = sparse(R.spkRip, R.spkUnit, 1, nRip, N);
    RP.spikes = table(R.spkRip, R.spkUnit, R.spkT, tsF(R.spkT), tsF3(R.spkT), R.spkPh, ...
        'VariableNames', {'Ripple','Unit','ts1','ts2','ts3','Phase'});

    % per-unit summary over all in-ripple spikes
    S = load(f, 'spikes');
    uid = nan(N, 1);
    if isfield(S.spikes, 'UID') && numel(S.spikes.UID) >= max(R.unitIdx), uid = S.spikes.UID(R.unitIdx); uid = uid(:); end
    rate = cellfun(@numel, S.spikes.times(R.unitIdx)); rate = rate(:) / R.epochT(1);
    okP = isfinite(R.spkPh);
    n = accumarray(R.spkUnit(okP), 1, [N 1]);
    C = accumarray(R.spkUnit(okP), cos(R.spkPh(okP)), [N 1]);
    Sn = accumarray(R.spkUnit(okP), sin(R.spkPh(okP)), [N 1]);
    ppc = (C.^2 + Sn.^2 - n) ./ (n .* (n - 1)); ppc(n < 2) = NaN;
    Rres = sqrt(C.^2 + Sn.^2) ./ n; Rres(n == 0) = NaN;
    pref = atan2(Sn, C); pref(n == 0) = NaN;
    RP.units = table(R.unitIdx(:), uid, rate, n, pref, Rres, ppc, ...
        'VariableNames', {'SpikesIdx','UID','Rate','nRipSpk','PrefPhase','ResultantR','PPC'});

    [~, mode, dayOff, rf] = tsF(0);
    [~, mode3, day3, rf3] = tsF3(0);
    RP.meta = struct('file',f, 'folder',fileparts(f), 'animal',ani{i}, 'condition',cnd{i}, ...
        'sessionDate',dnm(i), 'sessionDateStr',datestr(dnm(i),'yyyy-mm-dd'), ...
        'recStart',recStart(i), 'infusionSec',inf_t(i), ...
        'firstASV',rf, 'dayFromFirstASV',dayOff, 'ts2Mode',mode, ...
        'firstSAL',rf3, 'dayFromFirstSAL',day3, 'ts3Mode',mode3, ...
        'ts2Note',['ts1 = s from recording start; ts2 = s from the animal''s first ASV infusion. ' ...
                   'ts3 = s from the animal''s first SAL (vehicle) infusion, same modes. ' ...
                   'ts2Mode clock = exact; date = assumes the same time-of-day start as the first-ASV recording (+-24 h); none = animal has no ASV session.'], ...
        'lfpFile',R.lfpFile, 'channel1based',R.chan, 'fs',R.fs, 'band',R.band, ...
        'filter','3rd-order Butterworth band-pass, zero-phase (filtfilt), Hilbert envelope and phase', ...
        'recordingSec',R.epochT(1), 'envMu',R.envMu, 'envSd',R.envSd, ...
        'nDetected',R.nDetected, 'nDropped',R.nDroppedOverlap, 'nUnits',N, ...
        'hasSleepState',R.hasSleep, 'rfOpts',{RF_OPTS}, ...
        'matlab',version, 'created',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    rr = dir(fullfile(fileparts(f), '*.ripples.events.mat'));
    if ~isempty(rr)
        Q = load(fullfile(rr(1).folder, rr(1).name), 'ripples');
        RP.meta.rippleFile = fullfile(rr(1).folder, rr(1).name);
        RP.meta.detectorinfo = Q.ripples.detectorinfo;
    end
    RP.meta.checksums = struct('rippleFeatures', local_md5('rippleFeatures'), 'runRippleExport', local_md5('runRippleExport'));

    RIPPLES = RP; %#ok<NASGU>
    README = ['RIPPLES: every ripple of the session with features (ev), per-unit spike counts (count), every in-ripple spike ' ...
              'with its ripple-band phase (spikes), per-unit phase summaries (units). ts1 = s from recording start, ' ...
              'ts2 = s from the animal''s first ASV infusion. See RIPPLES.meta. Written by runRippleExport.m.']; %#ok<NASGU>
    tmp = [out '.tmp'];
    save(tmp, 'RIPPLES', 'README', '-v7.3');
    movefile(tmp, out, 'f');
    nDone = nDone + 1; status{i} = 'ok';
    fprintf('[%3d/%3d] ok      %s  %d ripples (%d dropped), %d units, %d in-ripple spikes, ts2 %s, %.0f s\n', ...
        i, numel(fils), f, nRip, R.nDroppedOverlap, N, numel(R.spkT), mode, toc);
end
fprintf('\n%d written\n', nDone);
save(fullfile(CACHE_ROOT, 'runRippleExport_status.mat'), 'fils', 'status');


%% =========================================================================
function h = local_md5(name)
h = '';
p = which(name);
if isempty(p), return; end
try
    fid = fopen(p, 'r'); b = fread(fid, Inf, '*uint8'); fclose(fid);
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(b);
    h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
catch
end
end
