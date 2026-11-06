%% runSyncProfiles.m
% For every session: spike-resolved SPIKE-synchronization over the WHOLE
% recording plus interval-jitter nulls, saved next to the raw data as
% spikeSyncProfile.mat (variable PROFILE). One file per session folder.
%
% PROFILE
%   .ts1      spike time, s from the start of the recording
%   .ts2      spike time, s from the animal's first ASV infusion (see .meta.ts2Mode)
%   .unit     index (into .meta.detIdx) of the unit that fired each spike
%   .idx      rank of the spike within its unit's train
%   .count    observed coincidence count per spike (0..N-1);  C = count/(N-1)
%   .null     M x NSurr counts of the same spikes after interval jitter, on the
%             OBSERVED spike axis (row k = the spike that started as observed k)
%   .N        number of units
%   .meta     everything needed to reproduce (below)
%
% A reader, per spike:  C = P.count / P.nPartners;  Cnull = P.null / P.nPartners.
% Any readout (spike-count windows, kernels, time windows) should be applied to
% the observed column and to every null column identically.
%
% JITTER: interval jitter on a fixed grid anchored at t = 0 of the recording
% (see syncProfileSession). Default J = 1 s.
%
% ts2 needs each recording's start clock time (HHMMSS in the folder name);
% otherwise it is date-resolved and approximate -- see ts2FromReference. Both
% modes are recorded in meta.
%
% COST. One spikeSync pass over the whole session per surrogate. Test first:
% set LIMIT = 2 and NSURR = 5 and look at the time and file size.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';   % only for the session list log
ROOT = 'R:\TransPlasticity';
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
EXPOSURE_COND = 'ASV';
JITTER = 1;                 % s
NSURR = 100;
MATCHN = 10;                % [] = every unit in the rate band
RATE_RANGE = [0.1 10];
OVERWRITE = false;          % false: skip sessions whose file already has these settings
LIMIT = 1;                % process at most this many sessions (testing)
OUT_NAME = 'spikeSyncProfile.mat';

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

fprintf('%d sessions, %d animals with a first-ASV reference\n', numel(fils), ...
    sum(cellfun(@(s) ~isempty(s), {ref.item.session})));

status = cell(numel(fils), 1);
nDone = 0;
for i = 1:numel(fils)
    if nDone >= LIMIT, break; end
    f = fils{i};
    out = fullfile(fileparts(f), OUT_NAME);
    if ~OVERWRITE && exist(out, 'file')
        try
            L = load(out, 'PROFILE');
            q = L.PROFILE.meta.params;
            if q.Jitter == JITTER && q.NSurr == NSURR && isequal(q.MatchN, MATCHN) && isequal(q.RateRange, RATE_RANGE)
                status{i} = 'exists'; fprintf('[%3d/%3d] exists   %s\n', i, numel(fils), f); continue
            end
        catch
        end
    end
    tic;
    try
        R = syncProfileSession(f, 'Jitter', JITTER, 'NSurr', NSURR, 'MatchN', MATCHN, 'RateRange', RATE_RANGE);
    catch ME
        status{i} = ['error: ' ME.message]; fprintf('[%3d/%3d] ERROR    %s: %s\n', i, numel(fils), f, ME.message); continue
    end
    if ~R.ok
        status{i} = ['skip: ' R.why]; fprintf('[%3d/%3d] skip     %s: %s\n', i, numel(fils), f, R.why); continue
    end
    [ts2, mode, dayOff, rf] = ts2FromReference(R.t, ani{i}, dnm(i), recStart(i), ref);

    P = struct();
    P.ts1 = R.t; P.ts2 = ts2;
    P.unit = R.unit; P.idx = R.idx; P.count = R.count; P.null = R.null;
    P.N = R.N; P.nPartners = R.nPartners;
    P.C = R.C; P.CnullMean = R.CnullMean; P.CnullSd = R.CnullSd;
    m = struct('file',f, 'folder',fileparts(f), 'animal',ani{i}, 'condition',cnd{i}, ...
        'sessionDate',dnm(i), 'sessionDateStr',datestr(dnm(i),'yyyy-mm-dd'), ...
        'recStart',recStart(i), 'infusionSec',inf_t(i), ...
        'firstASV',rf, 'dayFromFirstASV',dayOff, 'ts2Mode',mode, ...
        'ts2Note',['ts1 = s from recording start; ts2 = s from the animal''s first ASV infusion. ' ...
                   'ts2Mode clock = exact; date = assumes the same time-of-day start as the first-ASV recording (+-24 h); none = animal has no ASV session.'], ...
        'dur',R.dur, 'keepIdx',R.keepIdx, 'detIdx',R.detIdx, 'uid',R.uid, ...
        'params',R.params, 'matlab',version, 'created',datestr(now,'yyyy-mm-dd HH:MM:SS'));
    m.nullDefinition = ['Interval jitter, fixed grid anchored at t = 0 of the recording, bin width params.Jitter (s). ' ...
        'Each spike redrawn uniformly within its own bin; per-unit bin counts preserved. Surrogate s uses ' ...
        'RandStream(''mt19937ar'',''Seed'',params.Seed + s). null(k,s) is on the observed spike axis (the i-th spike of a unit stays in its bin).'];
    m.coincidence = 'SPIKE-synchronization (Kreuz), spikeSync.m, Interval [0 dur]; count = C*(N-1), the number of the other N-1 units with a coincident spike';
    m.checksums = struct('spikeSync', local_md5('spikeSync'), 'syncProfileSession', local_md5('syncProfileSession'));
    P.meta = m;
    PROFILE = P; %#ok<NASGU>
    README = ['PROFILE: spike-resolved SPIKE-synchronization over the whole session. ts1 = s from recording start, ' ...
              'ts2 = s from the animal''s first ASV infusion, count/(N-1) = C(t_k), null = interval-jitter surrogates on the observed spike axis. ' ...
              'See PROFILE.meta. Written by runSyncProfiles.m.']; %#ok<NASGU>
    tmp = [out '.tmp'];
    save(tmp, 'PROFILE', 'README', '-v7.3');
    movefile(tmp, out, 'f');
    nDone = nDone + 1;
    status{i} = 'ok';
    fprintf('[%3d/%3d] ok       %s  M=%d spikes, N=%d, ts2 %s, %.0f s\n', i, numel(fils), f, numel(R.t), R.N, mode, toc);
end
fprintf('\n%d written, %d skipped/errors/existing\n', nDone, sum(~strcmp(status,'ok') & ~cellfun(@isempty,status)));
save(fullfile(CACHE_ROOT, 'runSyncProfiles_status.mat'), 'fils', 'status', 'JITTER', 'NSURR', 'MATCHN', 'RATE_RANGE');


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
