%% runRippleExportControl.m
% runRippleExport for CONTROL animals that never received an infusion: every
% ripple over the whole recording, the same features, per-unit counts and
% in-ripple spike phases (every unit, no rate exclusion), but no condition and
% no infusion reference. Output next to the raw data as
% rippleEvents_control.mat (variable RIPPLES), a different name so the
% infusion-aligned scripts never pick these up.
%
% RIPPLES
%   .ts1, .ts1Peak, .ts1Off   ripple onset / envelope peak / offset, s from
%                             recording start (the only time axis: no ts2/ts3)
%   .ev, .count, .units       as in runRippleExport
%   .spikes                   Ripple, Unit, ts1, Phase
%   .meta                     meta.condition = '', meta.control = true
%
% Sessions are every *.spikes.cellinfo.mat under CONTROL_ROOT that has a
% *.ripples.events.mat; any session with an Infusion_start marker is skipped
% and listed. Keep RF_OPTS identical to runRippleExport so the sets compare.

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
CONTROL_ROOT = 'R:\HuszarR';        % one folder per session: <root>\<session>\<session>.*
RF_OPTS = {};               % must match runRippleExport; add {'Channel', k} (1-based) if the
                            % events files store no detection channel
% per-session ripple channel (1-based), for events files that do not store it
% and sessions without session.channelTags.Ripple, e.g.
%   CHANNEL_MAP = containers.Map({'e13_26m1_211019','e13_26m1_211020'}, {42, 42});
CHANNEL_MAP = containers.Map('KeyType','char','ValueType','any');
OVERWRITE = false;
LIMIT = Inf;
OUT_NAME = 'rippleEvents_control.mat';

fils = getAllExtFiles(CONTROL_ROOT, 'mat', 1);
fils = fils(contains(fils, 'spikes.cellinfo.mat'));
[~, rootName] = fileparts(CONTROL_ROOT);
ani = cell(numel(fils), 1); dnm = nan(numel(fils), 1); inf_t = nan(numel(fils), 1);
for i = 1:numel(fils)
    ani{i} = local_animal(fils{i}, rootName);
    dnm(i) = sessionDateFromPath(fils{i});
    inf_t(i) = infusionTime(fils{i});
end
recStart = cellfun(@recordingStartFromPath, fils);
hasInf = isfinite(inf_t);
if any(hasInf)
    fprintf('skipping %d session(s) with an Infusion_start marker (not controls):\n', sum(hasInf));
    fprintf('   %s\n', fils{hasInf});
end
fils = fils(~hasInf); ani = ani(~hasInf); dnm = dnm(~hasInf); recStart = recStart(~hasInf);
fprintf('%d control sessions, %d animals\n', numel(fils), numel(unique(ani)));

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
        ro = RF_OPTS;
        [~, sess] = fileparts(fileparts(f));
        if isKey(CHANNEL_MAP, sess), ro = [ro, {'Channel', CHANNEL_MAP(sess)}]; end
        % events files that store each event's phase series (ripples.phase_ts):
        % use it, no LFP needed. Otherwise compute from the LFP as in runRippleExport.
        R = rippleFeaturesStored(f, 'RateRange', [0 Inf]);
        if ~R.ok && strcmp(R.why, 'no ripples.phase_ts')
            R = rippleFeatures(f, [], 'RateRange', [0 Inf], ro{:});
            R.source = 'LFP (rippleFeatures)';
        end
    catch ME
        status{i} = ['error: ' ME.message];
        fprintf('[%3d/%3d] ERROR   %s: %s\n', i, numel(fils), f, ME.message);
        for q = 1:min(3, numel(ME.stack))
            fprintf('            at %s (line %d)\n', ME.stack(q).name, ME.stack(q).line);
        end
        continue
    end
    if ~R.ok
        status{i} = ['skip: ' R.why]; fprintf('[%3d/%3d] skip    %s: %s\n', i, numel(fils), f, R.why); continue
    end
    ev = R.ev;
    nRip = height(ev); N = R.nUnits;

    RP = struct();
    RP.ts1 = ev.On; RP.ts1Peak = ev.PeakT; RP.ts1Off = ev.On + ev.Dur;
    RP.ev = removevars(ev, {'Epoch','On','PeakT'});
    RP.count = sparse(R.spkRip, R.spkUnit, 1, nRip, N);
    RP.spikes = table(R.spkRip, R.spkUnit, R.spkT, R.spkPh, ...
        'VariableNames', {'Ripple','Unit','ts1','Phase'});

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

    ds = ''; if isfinite(dnm(i)), ds = datestr(dnm(i), 'yyyy-mm-dd'); end
    RP.meta = struct('file',f, 'folder',fileparts(f), 'animal',ani{i}, 'condition','', 'control',true, ...
        'sessionDate',dnm(i), 'sessionDateStr',ds, 'recStart',recStart(i), ...
        'ts1Note','ts1 = s from recording start. Control recording: no infusion, so no ts2/ts3.', ...
        'lfpFile',R.lfpFile, 'channel1based',R.chan, 'fs',R.fs, 'band',R.band, ...
        'filter','3rd-order Butterworth band-pass, zero-phase (filtfilt), Hilbert envelope and phase', ...
        'recordingSec',R.epochT(1), 'envMu',R.envMu, 'envSd',R.envSd, ...
        'nDetected',R.nDetected, 'nDropped',R.nDroppedOverlap, 'nUnits',N, ...
        'hasSleepState',R.hasSleep, 'rfOpts',{RF_OPTS}, 'phaseSource',R.source, ...
        'matlab',version, 'created',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    rr = dir(fullfile(fileparts(f), '*.ripples.events.mat'));
    if ~isempty(rr)
        Q = load(fullfile(rr(1).folder, rr(1).name), 'ripples');
        RP.meta.rippleFile = fullfile(rr(1).folder, rr(1).name);
        RP.meta.detectorinfo = Q.ripples.detectorinfo;
    end
    if strcmp(R.source, 'stored phase_ts')
        RP.meta.filter = 'none: phase from ripples.phase_ts (detector''s own phase series and convention); Amp/Pow NaN, AmpZ = ripples.peakNormedPower if stored';
    end
    RP.meta.checksums = struct('rippleFeatures', local_md5('rippleFeatures'), 'rippleFeaturesStored', local_md5('rippleFeaturesStored'), ...
        'runRippleExportControl', local_md5('runRippleExportControl'));

    RIPPLES = RP; %#ok<NASGU>
    README = ['RIPPLES (control, no infusion): every ripple of the session with features (ev), per-unit spike counts (count), ' ...
              'every in-ripple spike with its ripple-band phase (spikes), per-unit phase summaries (units). ' ...
              'ts1 = s from recording start. See RIPPLES.meta. Written by runRippleExportControl.m.']; %#ok<NASGU>
    tmp = [out '.tmp'];
    save(tmp, 'RIPPLES', 'README', '-v7.3');
    movefile(tmp, out, 'f');
    nDone = nDone + 1; status{i} = 'ok';
    fprintf('[%3d/%3d] ok      %s  %d ripples (%d dropped), %d units, %d in-ripple spikes, %.0f s\n', ...
        i, numel(fils), f, nRip, R.nDroppedOverlap, N, numel(R.spkT), toc);
end
fprintf('\n%d written\n', nDone);
save(fullfile(CACHE_ROOT, 'runRippleExportControl_status.mat'), 'fils', 'status');


%% =========================================================================
function a = local_animal(f, ~)
% animal ID: the first CP<number> token in the path if there is one (TransPlasticity
% naming); otherwise the session folder name with a trailing _YYMMDD[_HHMMSS]
% removed, e.g. R:\HuszarR\e13_26m1_211019\...  ->  e13_26m1
tok = regexp(f, '(CP\d+)', 'tokens', 'once');
if ~isempty(tok), a = tok{1}; return; end
[d, ~] = fileparts(f); [~, a] = fileparts(d);
a = regexprep(a, '_\d{6}(_\d{6})?$', '');
end


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