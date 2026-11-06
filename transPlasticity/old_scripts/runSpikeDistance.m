%% runSpikeDistance.m
% For every session: the SPIKE-distance profile S(t) (Kreuz et al. 2011/2013;
% time-resolved, time-scale adaptive) over the WHOLE recording, plus
% interval-jitter nulls, saved next to the raw data as spikeDistanceProfile.mat
% (variable SDIST). Supersedes runSyncProfiles.m, which stored
% SPIKE-synchronization coincidence counts instead.
%
% SDIST
%   .bin.ts1     bin centres, s from recording start
%   .bin.ts2     bin centres, s from the animal's first ASV infusion
%   .bin.edges   bin edges (s from recording start); width .meta.params.BinWidth
%   .bin.S       observed S, exact mean over each bin
%   .bin.Snull   nBins x NSurr, the same for every jittered surrogate (single)
%   .exact.ts1   event times of the observed profile (s from recording start)
%   .exact.ts2   the same, s from the first ASV infusion
%   .bin.ts3, .exact.ts3   the same, s from the animal's first SAL infusion
%   .exact.y1/y2 value at the start / end of each interval [ts1(k) ts1(k+1));
%                S is linear in between, discontinuous at spikes
%   .D, .Dnull   time-averaged S, observed and per surrogate
%   .meta        everything needed to reproduce
%
% S is a DISSIMILARITY (0 = identical timing): lower = more synchronous.
% Excess synchrony at any resolution: average .bin.S and every column of
% .bin.Snull over the same bins, then compare (e.g. z or percentile).
%
% Before the batch: run spikeDistanceProfile() with no arguments. It checks
% the MATLAB port against a stored PySpike 0.9.0 result and should print PASS.
% Then test with LIMIT = 2 and NSURR = 5 for time and file size.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT = 'R:\TransPlasticity';
%ROOT = 'R:\TransPlasticity\Ephys\CP13\ASV24\CP13_250617\'

CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
EXPOSURE_COND = 'ASV';
JITTER = 1;              % s
NSURR = 1;
BINW = 0.1;              % s
MATCHN = 2;             % [] = every unit in the rate band
RATE_RANGE = [0.1 inf];
KEEP_EXACT = 'auto';     % true | false | 'auto': exact observed profile only if pairs x spikes <= 2e9
MAX_PAIRS = 1000;         % sessions with more unit pairs average over a seeded random subset of this many
                         % (unbiased for the all-pairs mean; same pairs for data and surrogates). Inf = all
OVERWRITE = true;
LIMIT = Inf;
OUT_NAME = 'spikeDistanceProfile.mat';

assert(spikeDistanceProfile(), 'spikeDistanceProfile self-test failed: do not run the batch');

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
        try
            L = load(out, 'SDIST'); q = L.SDIST.meta.params;
            if q.Jitter == JITTER && q.NSurr == NSURR && q.BinWidth == BINW && ...
                    isequal(q.MatchN, MATCHN) && isequal(q.RateRange, RATE_RANGE) && ...
                    (~isfield(q, 'MaxPairs') || isequal(q.MaxPairs, MAX_PAIRS))   % older files used all pairs: kept
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
        for q = 1:min(3, numel(ME.stack))                % where it happened
            fprintf('            at %s (line %d)\n', ME.stack(q).name, ME.stack(q).line);
        end
        continue
    end
    if ~R.ok
        status{i} = ['skip: ' R.why]; fprintf('[%3d/%3d] skip    %s: %s\n', i, numel(fils), f, R.why); continue
    end
    ctr = 0.5 * (R.edges(1:end-1) + R.edges(2:end));
    [ts2b, mode, dayOff, rf] = ts2FromReference(ctr, ani{i}, dnm(i), recStart(i), ref);
    ts2e = ts2FromReference(R.E, ani{i}, dnm(i), recStart(i), ref);
    [ts3b, mode3, day3, rf3] = ts2FromReference(ctr, ani{i}, dnm(i), recStart(i), refSAL);
    ts3e = ts2FromReference(R.E, ani{i}, dnm(i), recStart(i), refSAL);

    SD = struct();
    SD.bin = struct('ts1',ctr, 'ts2',ts2b, 'ts3',ts3b, 'edges',R.edges, 'S',R.S, 'Snull',R.Snull);
    SD.exact = struct('ts1',R.E, 'ts2',ts2e, 'ts3',ts3e, 'y1',R.y1, 'y2',R.y2);
    SD.D = R.D; SD.Dnull = R.Dnull; SD.N = R.N;
    m = struct('file',f, 'folder',fileparts(f), 'animal',ani{i}, 'condition',cnd{i}, ...
        'sessionDate',dnm(i), 'sessionDateStr',datestr(dnm(i),'yyyy-mm-dd'), ...
        'recStart',recStart(i), 'infusionSec',inf_t(i), ...
        'firstASV',rf, 'dayFromFirstASV',dayOff, 'ts2Mode',mode, ...
        'firstSAL',rf3, 'dayFromFirstSAL',day3, 'ts3Mode',mode3, ...
        'ts2Note',['ts1 = s from recording start; ts2 = s from the animal''s first ASV infusion. ' ...
                   'ts3 = s from the animal''s first SAL (vehicle) infusion, same modes. ' ...
                   'ts2Mode clock = exact; date = assumes the same time-of-day start as the first-ASV recording (+-24 h); none = animal has no ASV session.'], ...
        'dur',R.dur, 'keepIdx',R.keepIdx, 'detIdx',R.detIdx, 'uid',R.uid, 'pairs',R.pairs, ...
        'params',R.params, 'matlab',version, 'created',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    m.measure = ['Multivariate SPIKE-distance (Kreuz et al. 2011 J Neurosci Methods; rate-weighted form of Kreuz et al. 2013 J Neurophysiol; ' ...
        'PySpike 0.9.0 edge handling), mean over all unit pairs, Interval [0 dur]. Dissimilarity: 0 = identical timing.'];
    m.nullDefinition = ['Interval jitter, fixed grid anchored at t = 0 of the recording, bin width params.Jitter (s). ' ...
        'Each spike redrawn uniformly within its own bin; per-unit bin counts preserved. Surrogate s uses ' ...
        'RandStream(''mt19937ar'',''Seed'',params.Seed + s). Null profiles stored as exact bin means on bin.edges.'];
    m.checksums = struct('spikeDistanceProfile', local_md5('spikeDistanceProfile'), ...
        'spikeDistanceBinMean', local_md5('spikeDistanceBinMean'), 'spikeDistanceSession', local_md5('spikeDistanceSession'));
    SD.meta = m;
    SDIST = SD; %#ok<NASGU>
    README = ['SDIST: SPIKE-distance profile over the whole session (lower = more synchronous). bin.* on a fixed grid with ' ...
              'interval-jitter nulls (bin.Snull); exact.* the observed piecewise-linear profile. ts1 = s from recording start, ' ...
              'ts2 = s from the animal''s first ASV infusion. See SDIST.meta. Written by runSpikeDistance.m.']; %#ok<NASGU>
    tmp = [out '.tmp'];
    save(tmp, 'SDIST', 'README', '-v7.3');
    movefile(tmp, out, 'f');
    nDone = nDone + 1; status{i} = 'ok';
    fprintf('[%3d/%3d] ok      %s  N=%d, %d events, D=%.3f (null %.3f +- %.3f), ts2 %s, %.0f s\n', i, numel(fils), f, ...
        R.N, numel(R.E), R.D, mean(R.Dnull), std(R.Dnull), mode, toc);
end
fprintf('\n%d written\n', nDone);
save(fullfile(CACHE_ROOT, 'runSpikeDistance_status.mat'), 'fils', 'status', 'JITTER', 'NSURR', 'BINW', 'MATCHN', 'RATE_RANGE');


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