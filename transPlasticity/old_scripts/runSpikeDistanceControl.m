%% runSpikeDistanceControl.m
% runSpikeDistance for CONTROL animals that never received an infusion: same
% SPIKE-distance profile, same interval-jitter nulls, same settings, but no
% condition and no infusion reference. Output next to the raw data as
% spikeDistanceProfile_control.mat (variable SDIST), a different name so the
% infusion-aligned scripts (plotExcessSyncAligned, addTs3) never pick these up.
%
% SDIST
%   .bin.ts1     bin centres, s from recording start (the only time axis:
%                there is no infusion, so no ts2 / ts3)
%   .bin.edges, .bin.S, .bin.Snull     as in runSpikeDistance
%   .exact.ts1, .exact.y1, .exact.y2   exact observed profile (when computed)
%   .D, .Dnull, .N, .meta              meta.condition = '', meta.control = true
%
% Sessions are every *.spikes.cellinfo.mat under CONTROL_ROOT. Any session
% that DOES carry an Infusion_start marker is skipped and listed, so a
% misplaced infusion recording cannot enter the control set.
%
% Settings must match runSpikeDistance for the two sets to be comparable.

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
CONTROL_ROOT = 'R:\HuszarR';   % <- set to where the control sessions live
JITTER = 1;              % s
NSURR = 100;
BINW = 0.1;              % s
MATCHN = [];             % [] = every unit in the rate band
RATE_RANGE = [0 inf];
KEEP_EXACT = 'auto';
MAX_PAIRS = 500;
OVERWRITE = false;
LIMIT = Inf;
OUT_NAME = 'spikeDistanceProfile_control.mat';

assert(spikeDistanceProfile(), 'spikeDistanceProfile self-test failed: do not run the batch');

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
        try
            L = load(out, 'SDIST'); q = L.SDIST.meta.params;
            if q.Jitter == JITTER && q.NSurr == NSURR && q.BinWidth == BINW && ...
                    isequal(q.MatchN, MATCHN) && isequal(q.RateRange, RATE_RANGE) && ...
                    isfield(q, 'MaxPairs') && isequal(q.MaxPairs, MAX_PAIRS)
                status{i} = 'exists'; fprintf('[%3d/%3d] exists  %s\n', i, numel(fils), f); continue
            end
        catch
        end
    end
    tic;
    try
        R = spikeDistanceSession(f, 'Jitter', JITTER, 'NSurr', NSURR, 'BinWidth', BINW, ...
            'MatchN', MATCHN, 'RateRange', RATE_RANGE, 'KeepExact', KEEP_EXACT, 'MaxPairs', MAX_PAIRS);
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
    ctr = 0.5 * (R.edges(1:end-1) + R.edges(2:end));

    SD = struct();
    SD.bin = struct('ts1',ctr, 'edges',R.edges, 'S',R.S, 'Snull',R.Snull);
    SD.exact = struct('ts1',R.E, 'y1',R.y1, 'y2',R.y2);
    SD.D = R.D; SD.Dnull = R.Dnull; SD.N = R.N;
    ds = ''; if isfinite(dnm(i)), ds = datestr(dnm(i), 'yyyy-mm-dd'); end
    m = struct('file',f, 'folder',fileparts(f), 'animal',ani{i}, 'condition','', 'control',true, ...
        'sessionDate',dnm(i), 'sessionDateStr',ds, 'recStart',recStart(i), ...
        'ts1Note','ts1 = s from recording start. Control recording: no infusion, so no ts2/ts3.', ...
        'dur',R.dur, 'keepIdx',R.keepIdx, 'detIdx',R.detIdx, 'uid',R.uid, 'pairs',R.pairs, ...
        'params',R.params, 'matlab',version, 'created',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    m.measure = ['Multivariate SPIKE-distance (Kreuz et al. 2011 J Neurosci Methods; rate-weighted form of Kreuz et al. 2013 J Neurophysiol; ' ...
        'PySpike 0.9.0 edge handling), mean over unit pairs (see params.nPairsUsed), Interval [0 dur]. Dissimilarity: 0 = identical timing.'];
    m.nullDefinition = ['Interval jitter, fixed grid anchored at t = 0 of the recording, bin width params.Jitter (s). ' ...
        'Each spike redrawn uniformly within its own bin; per-unit bin counts preserved. Surrogate s uses ' ...
        'RandStream(''mt19937ar'',''Seed'',params.Seed + s). Null profiles stored as exact bin means on bin.edges.'];
    m.checksums = struct('spikeDistanceProfile', local_md5('spikeDistanceProfile'), ...
        'spikeDistanceBinMean', local_md5('spikeDistanceBinMean'), 'spikeDistanceSession', local_md5('spikeDistanceSession'));
    SD.meta = m;
    SDIST = SD; %#ok<NASGU>
    README = ['SDIST (control, no infusion): SPIKE-distance profile over the whole session (lower = more synchronous). ' ...
              'bin.* on a fixed grid with interval-jitter nulls (bin.Snull); exact.* the observed piecewise-linear profile. ' ...
              'ts1 = s from recording start. See SDIST.meta. Written by runSpikeDistanceControl.m.']; %#ok<NASGU>
    tmp = [out '.tmp'];
    save(tmp, 'SDIST', 'README', '-v7.3');
    movefile(tmp, out, 'f');
    nDone = nDone + 1; status{i} = 'ok';
    fprintf('[%3d/%3d] ok      %s  N=%d, D=%.3f (null %.3f +- %.4f), %.0f s\n', i, numel(fils), f, ...
        R.N, R.D, mean(R.Dnull), std(R.Dnull), toc);
end
fprintf('\n%d written\n', nDone);
save(fullfile(CACHE_ROOT, 'runSpikeDistanceControl_status.mat'), 'fils', 'status', 'JITTER', 'NSURR', 'BINW', 'MATCHN', 'RATE_RANGE');


%% =========================================================================
function a = local_animal(f, rootName)
% animal ID: animalFromPath under the control root if it parses, else the
% first CP<number> token in the path, else the session folder name
a = '';
try, a = animalFromPath(f, rootName); catch, end %#ok<CTCH>
if ~isempty(a), return; end
tok = regexp(f, '(CP\d+)', 'tokens', 'once');
if ~isempty(tok), a = tok{1}; return; end
[d, ~] = fileparts(f); [~, a] = fileparts(d);
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
