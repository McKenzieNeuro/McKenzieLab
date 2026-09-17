function T = syncBatch(fils, varargin)
%SYNCBATCH  Run SPIKE-Synchronization over a list of *.spikes.cellinfo.mat files.
%
%   T = syncBatch(fils, 'Duration', @myDurFcn, 'Exclude', @myExcludeFcn, ...)
%
%   Returns a struct array, one element per session, carrying everything the
%   downstream comparison needs -- including the quantities that confound C and
%   that a raw trace throws away (N, spike count, duration, rate distribution,
%   null level). Sessions that fail the unit criterion are RETURNED with
%   .ok = false rather than silently skipped, so indices stay aligned with fils.
%
%   OPTIONS
%     'RateRange'  [0.1 10]   unit inclusion band, Hz
%     'Duration'   []         [] = infer from last spike (biased low), or a
%                             scalar, or a function handle @(file) -> seconds.
%                             Pass the real recording length if you have it:
%                             it sets both the rate estimate and the coincidence
%                             window fallback for edge spikes.
%     'Width'      60         window width (s)
%     'Step'       5          window step (s). Step=1 with Width=60 is 98%
%                             overlap: 60x redundant, and the window count is
%                             not a sample size.
%     'MinSpikes'  10         windows below this return NaN
%     'Exclude'    []         @(file) -> K x 2 [start stop] intervals to drop.
%                             Any window OVERLAPPING an interval is NaN'd.
%     'NullMethod' 'shift'    'shift' | 'dither' | 'both' | 'none'
%     'Dither'     []         required if NullMethod uses dither
%     'NSurr'      50         surrogates per session (scalar null only)
%     'MatchN'     []         if set, also compute C on random unit subsets of
%                             this size, averaged over 'NSubsets' draws. This
%                             is the fix for 1/(N-1) dilution when sessions or
%                             datasets differ in unit yield. Without it, a
%                             dataset with more units reads as less synchronous
%                             for purely combinatorial reasons.
%     'NSubsets'   20
%     'Window'     []         restrict analysis to [t0 t1] seconds, or a
%                             function handle @(file) -> [t0 t1] for a
%                             per-session window (e.g. pre- vs post-infusion
%                             epochs, where the split time differs per
%                             recording). Spike times are re-zeroed to the
%                             window start; r.windowAbs keeps the absolute
%                             bounds. Sessions with no usable window are
%                             returned with .ok = false rather than erroring
%                             the batch.
%     'CacheDir'   []         per-session results cached here; reruns resume.
%                             Survives a crash between sessions but NOT a
%                             crash that kills the calling MATLAB session
%                             before syncBatch returns -- see CheckpointFile.
%     'CheckpointFile' []     path to a .mat file. If set, the accumulated T
%                             array is saved periodically DURING the run (see
%                             CheckpointEvery) via write-then-rename, so a
%                             mid-write crash never corrupts the checkpoint --
%                             worst case you lose back to the previous save,
%                             never the file itself. On the next call with
%                             this SAME CheckpointFile and an IDENTICAL option
%                             set + file list, already-completed sessions are
%                             loaded from the checkpoint and skipped entirely
%                             (no cache lookup, no recompute). A checkpoint
%                             from a different option set or file list is
%                             detected and ignored, not silently reused --
%                             same protection as the CacheDir key, at the top
%                             level instead of per-file.
%     'CheckpointEvery' 5     save the checkpoint every this many sessions
%                             (plus always on the last session). Lower values
%                             cost more disk I/O per session; higher values
%                             risk losing more progress to an interruption.
%     'Verbose'    true
%     'DetectTransients' false  find sync epochs above a shift-surrogate band.
%                             Adds r.transients: struct array with tStart,
%                             tEnd, nWin, peakDC (height above the surrogate
%                             median), significant (run length exceeds the
%                             (1-RunAlpha) percentile of the surrogate
%                             max-run-length distribution -- exact for the
%                             session's longest cluster, conservative for
%                             shorter ones in the same session). The exclusion
%                             mask from 'Exclude' is applied identically to
%                             every surrogate trace, so a transient adjacent to
%                             an excluded epoch is tested fairly.
%     'NSurrWindow' 50        surrogates for the windowed band (separate from
%                             'NSurr' -- windowing every surrogate is the
%                             expensive step; keep this lower for large batches
%                             and raise it for a session you're reporting on).
%     'RunAlpha'   0.05       significance level for the cluster run length.
%     'Parallel'   false      parfor across surrogate draws within each
%                             session (Parallel Computing Toolbox; falls back
%                             to serial with a warning if unavailable). Does
%                             NOT parallelize across sessions/files -- don't
%                             nest this inside an outer parfor over syncBatch
%                             calls.
%
%   OUTPUT fields per session
%     .ok .file .N .M .dur .rate .rateMed .rateRatio
%     .C .Cnull .Csd .dC .z          scalar synchrony, null, corrected
%     .Cmat .prof                     (dropped if 'Light', see below)
%     .Ct .tc .nSpk .nExcluded        windowed trace on its OWN time base
%     .CmatchN                        N-matched C (mean over subsets), if asked
%
%   See also SPIKESYNC, SPIKESYNCWINDOW, SPIKESYNCDIAGNOSE, SPIKESYNCCOMPARE.

p = inputParser;
p.addParameter('RateRange', [0.1 10]);
p.addParameter('Duration', []);
p.addParameter('MaxDuration', []);
p.addParameter('Window', []);
p.addParameter('Width', 60);
p.addParameter('Step', 5);
p.addParameter('MinSpikes', 10);
p.addParameter('Exclude', []);
p.addParameter('NullMethod', 'shift');
p.addParameter('Dither', []);
p.addParameter('NSurr', 50);
p.addParameter('MatchN', []);
p.addParameter('MinUnits', 2);
p.addParameter('NSubsets', 20);
p.addParameter('CacheDir', []);
p.addParameter('Light', false);
p.addParameter('Verbose', true);
p.addParameter('DetectTransients', false);
p.addParameter('NSurrWindow', 50);
p.addParameter('RunAlpha', 0.05);
p.addParameter('Parallel', false);
p.addParameter('CheckpointFile', []);
p.addParameter('CheckpointEvery', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.parse(varargin{:});
o = p.Results;
if o.Parallel && isempty(ver('parallel'))
    warning('syncBatch:parallel','Parallel Computing Toolbox not found; running serially.');
    o.Parallel = false;
end

if any(strcmpi(o.NullMethod, {'dither','both'})) && isempty(o.Dither)
    error('syncBatch:dither', ...
        ['NullMethod ''%s'' needs ''Dither''. Note that dither also destroys ', ...
         'within-train regularity, which inflates the apparent effect for ', ...
         'rhythmic or bursty units; ''shift'' preserves ISI structure exactly.'], ...
        o.NullMethod);
end
if ~isempty(o.CacheDir) && ~exist(o.CacheDir, 'dir'), mkdir(o.CacheDir); end

nF = numel(fils);
T = repmat(local_blank(), nF, 1);
done = false(nF,1);

% ---- resume from checkpoint, if one exists and matches this call's options
% Protects against losing the whole in-progress T array to something that
% kills the CALLING session before syncBatch returns -- an IdleTimeout on a
% parallel pool, a network drop, a closed laptop lid on an unattended run.
% The per-session CacheDir already survives that (each session is on disk
% individually), but nothing previously survived the top-level T array
% itself not yet having been returned when the crash happened.
%
% Same protection as the cache key: a checkpoint written under a different
% option set is NOT silently reused -- that would be the exact stale-result
% bug the cache-key fix addressed, just at the top level instead of per-file.
optHash = local_optHash(o, fils);
if ~isempty(o.CheckpointFile) && exist(o.CheckpointFile, 'file')
    try
        ckpt = load(o.CheckpointFile, 'T', 'optHash', 'nDone');
        if isfield(ckpt,'optHash') && strcmp(ckpt.optHash, optHash)
            nLoad = min(numel(ckpt.T), nF);
            T(1:nLoad) = ckpt.T(1:nLoad);
            % Only skip sessions that SUCCEEDED. Same reasoning as the cache
            % fix: a checkpointed failure may be a transient environment
            % problem (license contention, network drop), not a fact about
            % that session. Marking it 'done' here would mean resuming never
            % retries it, same bug as caching a failure, just at the
            % checkpoint level instead of the per-file cache level.
            done(1:nLoad) = arrayfun(@(r) r.ok, T(1:nLoad));
            if o.Verbose
                fprintf('resuming from checkpoint: %d/%d sessions already done\n', ...
                    sum(done), nF);
            end
        elseif o.Verbose
            fprintf(['checkpoint found but option set or file list differs from this call -- ', ...
                     'ignoring it and starting fresh (this is deliberate, not a bug: a stale ', ...
                     'checkpoint silently reused under new settings is worse than a slow start).\n']);
        end
    catch ME
        % Same failure mode as a corrupt per-session cache file: don't let an
        % unreadable checkpoint kill the whole call. Starting fresh costs
        % time; erroring out here costs the same time PLUS the debugging.
        warning('syncBatch:checkpointCorrupt', ...
            'checkpoint file unreadable (%s) -- starting fresh.', ME.message);
    end
end

for i = 1:nF
    f = fils{i};
    if done(i)
        continue
    end
    T(i).file = f;

    cf = local_cacheFile(o.CacheDir, f, o);
    if ~isempty(cf) && exist(cf, 'file')
        try
            S = load(cf, 'r'); T(i) = S.r;
            if o.Verbose, fprintf('[%3d/%3d] cached  %s\n', i, nF, local_short(f)); end
            local_maybeCheckpoint(o, T, optHash, i, nF, false);
            continue
        catch ME
            % A corrupt cache file used to be a hard error for the WHOLE
            % batch here (load() throwing mid-loop). Most likely cause: the
            % write that created it was interrupted -- previously a plain
            % save() with no atomicity, now write-then-rename below, but old
            % cache directories can still have files written the old way.
            % Delete the bad entry and fall through to recompute rather than
            % losing everything downstream of session i.
            warning('syncBatch:cacheCorrupt', ...
                '%s: cache file unreadable (%s) -- deleting and recomputing.', ...
                local_short(f), ME.message);
            try, delete(cf); catch, end %#ok<CTCH>
        end
    end

    try
        r = local_session(f, o);
    catch ME
        % Inside parfor, the exception that reaches here is often a
        % ParallelException whose .message is a generic wrapper -- the real
        % error is in .cause{1}. Unwrap it, and record identifier + where it
        % was thrown, so 'why' is actually diagnosable rather than "an error
        % occurred". Worker-side warning() output does not reliably reach the
        % client console, so this struct field is the only record you get.
        cause = ME;
        while strcmp(cause.identifier,'MATLAB:parallel_function:ParallelException') ...
                && ~isempty(cause.cause)
            cause = cause.cause{1};
        end
        where = '';
        if ~isempty(cause.stack)
            where = sprintf(' [%s:%d]', cause.stack(1).name, cause.stack(1).line);
        end
        r = local_blank(); r.file = f; r.ok = false;
        r.why = sprintf('%s: %s%s', cause.identifier, cause.message, where);
        r.errStack = cause.stack;
        warning('syncBatch:session','%s: %s', local_short(f), r.why);
    end
    r.file = f;
    T(i) = r;

    if ~isempty(cf) && r.ok
        % Only cache SUCCESSES. A caught exception (r.ok=false via the catch
        % above) is frequently a transient environment problem -- a floating
        % toolbox license someone else has checked out, a network drop to
        % R:\, an antivirus lock -- not a deterministic fact about the data.
        % Caching a failure used to mean a five-minute license blip became a
        % PERMANENT hole in the dataset: the next run would cache-hit on
        % "failed" and never retry. Deterministic exclusions (too few units,
        % bad duration) are cheap to redetect anyway -- the expensive part
        % never ran for those, so not caching them costs almost nothing.
        tmp = [cf '.tmp'];
        try
            save(tmp, 'r'); %#ok<USENS>
            movefile(tmp, cf, 'f');
        catch ME
            warning('syncBatch:cacheWrite','%s: cache write failed (%s)', ...
                local_short(f), ME.message);
        end
    end
    if o.Verbose
        if r.ok
            fprintf('[%3d/%3d] N=%3d/%3d  M=%7d  %6.0fs  C=%.3f  null=%.3f  dC=%+.3f  %s\n', ...
                i, nF, r.N, r.Nfull, r.M, r.dur, r.C, r.Cnull, r.dC, local_short(f));
        else
            fprintf('[%3d/%3d] skipped (%s)  %s\n', i, nF, r.why, local_short(f));
        end
    end

    local_maybeCheckpoint(o, T, optHash, i, nF, false);
end

local_maybeCheckpoint(o, T, optHash, nF, nF, true);  % final save, forced

if o.Verbose, local_summary(T, o); end

end


function local_maybeCheckpoint(o, T, optHash, i, nF, force)
if isempty(o.CheckpointFile), return; end
if ~force && mod(i, o.CheckpointEvery) ~= 0 && i ~= nF, return; end
nDone = sum(arrayfun(@(r) ~strcmp(r.why,'not run'), T));
% write-then-rename: a crash mid-write leaves the .tmp file corrupt but the
% previous good checkpoint (or none) intact, never a half-written .mat that
% 'load' fails on when you go to resume.
tmp = [o.CheckpointFile '.tmp'];
save(tmp, 'T', 'optHash', 'nDone', '-v7.3');
[status, msg] = movefile(tmp, o.CheckpointFile, 'f');
if ~status
    warning('syncBatch:checkpoint','checkpoint write failed: %s', msg);
end
end


function h = local_optHash(o, fils)
% Hash covers the option set AND the file list, so a checkpoint from a
% different (or reordered, or filtered) file list is also refused rather
% than silently misaligned against the current fils.
c = {o.RateRange, o.Width, o.Step, o.MinSpikes, o.NullMethod, o.Dither, ...
    o.NSurr, o.MatchN, o.NSubsets, o.DetectTransients, o.NSurrWindow, ...
    o.RunAlpha, o.MinUnits, o.MaxDuration, fils(:)};
s = char(jsonencode(c));
h = sprintf('%08x', mod(sum(double(s) .* (1:numel(s))), 2^32));
end


% =========================================================================
function r = local_session(f, o)
r = local_blank();

S = load(f, 'spikes');
if ~isfield(S, 'spikes') || ~isfield(S.spikes, 'times')
    error('no spikes.times in file');
end
ts = S.spikes.times(:);
ts = cellfun(@(x) sort(x(:)), ts, 'UniformOutput', false);

% ---- duration ------------------------------------------------------------
% max() of an empty unit returns [], which is why cellfun on raw times breaks.
lastSpk = cellfun(@(x) max([x; -Inf]), ts);
if isempty(o.Duration)
    dur = max(lastSpk(isfinite(lastSpk)));
    r.durInferred = true;
elseif isa(o.Duration, 'function_handle')
    dur = o.Duration(f); r.durInferred = false;
else
    dur = o.Duration; r.durInferred = false;
end
if isempty(dur) || ~isfinite(dur) || dur <= 0, error('bad duration'); end
r.dur = dur;

% ---- optional truncation for a duration-matched sensitivity check --------
% Caps BOTH the duration and the actual spike times to the first MaxDuration
% seconds. Not the same as passing a shorter 'Duration': that only moves the
% Interval boundary while spikes past it stay in ts, which makes spikeSync
% warn about out-of-interval spikes and biases the rate denominator without
% actually removing the extra spikes.
r.durTruncated = false;
if ~isempty(o.MaxDuration) && dur > o.MaxDuration
    dur = o.MaxDuration;
    ts = cellfun(@(x) x(x <= dur), ts, 'UniformOutput', false);
    r.dur = dur; r.durTruncated = true;
end

% ---- optional arbitrary analysis window ---------------------------------
% 'Window' restricts the analysis to [t0 t1] seconds -- e.g. the pre-infusion
% baseline epoch of an ASV session, or the post-infusion epoch. Supply either
% a fixed [t0 t1] or (usually) a function handle @(file) -> [t0 t1], since
% the infusion time differs in every recording.
%
% Spike times are RE-ZEROED to the window start, so downstream everything
% (Interval, transient tStart/tEnd, the windowed trace's time base) is
% expressed relative to window onset rather than to the recording. That
% keeps a 20-minute pre-infusion epoch and a 20-minute post-infusion epoch
% directly comparable instead of one of them sitting at a large time offset,
% and it means the edge-spike handling in spikeSync sees a full window
% rather than a sparse interval preceded by nothing.
%
% r.windowAbs keeps the original absolute [t0 t1] so results can be mapped
% back to the recording.
r.windowAbs = [];
if ~isempty(o.Window)
    if isa(o.Window, 'function_handle'), w = o.Window(f); else, w = o.Window; end
    if isempty(w) || numel(w) ~= 2 || any(~isfinite(w)) || w(2) <= w(1)
        error('no usable analysis window for this session');
    end
    w(1) = max(w(1), 0); w(2) = min(w(2), dur);
    if w(2) - w(1) <= 0, error('analysis window is empty after clipping to recording'); end
    ts = cellfun(@(x) x(x >= w(1) & x < w(2)) - w(1), ts, 'UniformOutput', false);
    dur = w(2) - w(1);
    r.dur = dur; r.windowAbs = w;
end

% ---- unit selection ------------------------------------------------------
M = cellfun(@numel, ts);
FR = M / dur;
keep = FR > o.RateRange(1) & FR < o.RateRange(2) & M > 0;
if sum(keep) < max(2, o.MinUnits)
    r.ok = false;
    r.why = sprintf('%d units pass rate band (need %d)', sum(keep), max(2, o.MinUnits));
    r.N = sum(keep); r.Nfull = sum(keep); return
end
ts = ts(keep);
r.Nfull = numel(ts); r.M = sum(M(keep));
r.rate = FR(keep); r.rateMed = median(r.rate);
r.rateRatio = max(r.rate) / max(min(r.rate), eps);
r.keepIdx = find(keep);

iv = [0 dur];

% ---- N-matching: everything downstream runs on a FIXED-SIZE subset -------
% Detection sensitivity scales with N: the surrogate band narrows as more
% units are averaged into C_i, so a high-yield session detects more
% transients from identical underlying physiology. Correcting the scalar C
% alone (the old MatchN) left the transient detector uncorrected, which is
% what showed up as rho(N, transient rate) ~ 0.43 across sessions.
%
% Fix: draw one subset of exactly MatchN units and compute the trace, the
% surrogate band and the transients from THAT. Sensitivity is then constant
% across sessions by construction rather than residualised after the fact.
% Seeded per file so the subset is reproducible across reruns.
%
% Side benefit: cost goes as N^2, so 10 units instead of 45 is ~20x less
% pairwise work per surrogate.
r.Cfull = spikeSync(ts, 'Interval', iv);   % full-N value, kept for reference

if ~isempty(o.MatchN)
    if r.Nfull < o.MatchN
        r.ok = false;
        r.why = sprintf('%d units < MatchN %d', r.Nfull, o.MatchN);
        r.N = r.Nfull; return
    end
    % mt19937ar, not threefry: threefry/philox are Parallel Computing
    % Toolbox generators. Using one here made every MatchN subset draw
    % silently depend on that toolbox being licensed AND available (a
    % floating seat can be checked out by someone else), which is a much
    % worse failure mode than it sounds -- see the cache-write note below.
    rs = RandStream('mt19937ar','Seed', local_seedFromName(f));
    sel = randperm(rs, r.Nfull, o.MatchN);
    r.detIdx = r.keepIdx(sel);
    ts = ts(sel);
end
r.N = numel(ts);

% ---- observed ------------------------------------------------------------
[C, Cmat, prof] = spikeSync(ts, 'Interval', iv);
r.C = C;
if ~o.Light, r.Cmat = Cmat; r.prof = prof; end

% ---- trace ---------------------------------------------------------------
[Ct, tc, nSpk] = spikeSyncWindow(prof, 'Width', o.Width, 'Step', o.Step, ...
    'Range', iv, 'MinSpikes', o.MinSpikes);

% exclusion intervals -> NaN every window that OVERLAPS one. The original
% histc approach tagged only the window whose edge preceded the event, which
% leaves ~Width/Step contaminated windows in the trace.
nEx = 0;
excludeMask = false(size(tc));
if ~isempty(o.Exclude)
    ex = o.Exclude(f);
    if ~isempty(ex)
        ex = reshape(ex, [], 2);
        bad = false(size(tc));
        for e = 1:size(ex,1)
            bad = bad | (tc + o.Width/2 > ex(e,1) & tc - o.Width/2 < ex(e,2));
        end
        Ct(bad) = NaN; nEx = sum(bad);
        excludeMask = bad;
    end
end
r.Ct = Ct; r.tc = tc; r.nSpk = nSpk; r.nExcluded = nEx;

% ---- null + transients: ONE shared surrogate loop ------------------------
% local_null (scalar) and local_transients (windowed) previously each drew
% their own independent set of shift-surrogates and ran a full pairwise
% spikeSync on every one -- redundant work when both use 'shift', which is
% the only method local_transients supports. For a 9-hour, ~20-unit session
% each full spikeSync is the expensive step (O(N^2) pairs x O(M log M) each),
% so duplicating the surrogate draws roughly doubled runtime for nothing.
% Now: one surrogate loop produces both the scalar C and, if requested, the
% windowed Ct, from the SAME surrogate spike trains.
needShiftNull = strcmpi(o.NullMethod,'shift') || strcmpi(o.NullMethod,'both');
nSurrShared = max(o.NSurr * needShiftNull, o.NSurrWindow * o.DetectTransients);

if nSurrShared > 0 && (needShiftNull || o.DetectTransients)
    [Cs, CtS] = local_sharedSurrogates(ts, iv, nSurrShared, o.Width, o.Step, ...
        o.MinSpikes, excludeMask, o.DetectTransients, o.Parallel);
    if needShiftNull
        r.Cnull = mean(Cs(1:o.NSurr)); r.Csd = std(Cs(1:o.NSurr));
    end
else
    Cs = []; CtS = [];
end

switch lower(o.NullMethod)
    case 'none'
        r.Cnull = NaN; r.Csd = NaN;
    case 'shift'
        % filled above
    case 'dither'
        [r.Cnull, r.Csd] = local_null(ts, iv, 'dither', o.NSurr, o.Dither);
    case 'both'
        [r.CnullDither, r.CsdDither] = local_null(ts, iv, 'dither', o.NSurr, o.Dither);
    otherwise
        error('unknown NullMethod %s', o.NullMethod);
end
r.dC = r.C - r.Cnull;
r.z  = r.dC / max(r.Csd, eps);

% ---- N-matched C ---------------------------------------------------------
% (superseded: with MatchN set, ts IS the matched subset and r.C is already
% the N-matched value. r.Cfull holds the all-units number.)
r.CmatchN = r.C;

r.ok = true; r.why = '';

% ---- transients ------------------------------------------------------
if o.DetectTransients
    r.transients = local_transientsFromSurrogates(Ct, tc, CtS, o.RunAlpha);
end
end


function [Cs, CtS] = local_sharedSurrogates(ts, iv, nsurr, width, step, minSpk, ...
    excludeMask, wantWindow, useParallel)
N = numel(ts); L = iv(2) - iv(1);
Cs = zeros(nsurr,1);
CtTmp = cell(nsurr,1);
if useParallel
    parfor s = 1:nsurr
        sur = cell(N,1);
        for n = 1:N
            x = ts{n};
            sur{n} = sort(iv(1) + mod(x - iv(1) + rand*L, L));
        end
        if wantWindow
            [Cs(s), ~, pS] = spikeSync(sur, 'Interval', iv);
            ct = spikeSyncWindow(pS, 'Width', width, 'Step', step, ...
                'Range', iv, 'MinSpikes', minSpk);
            ct(excludeMask) = NaN;
            CtTmp{s} = ct;
        else
            Cs(s) = spikeSync(sur, 'Interval', iv);
        end
    end
else
    for s = 1:nsurr
        sur = cell(N,1);
        for n = 1:N
            x = ts{n};
            sur{n} = sort(iv(1) + mod(x - iv(1) + rand*L, L));
        end
        if wantWindow
            [Cs(s), ~, pS] = spikeSync(sur, 'Interval', iv);
            ct = spikeSyncWindow(pS, 'Width', width, 'Step', step, ...
                'Range', iv, 'MinSpikes', minSpk);
            ct(excludeMask) = NaN;
            CtTmp{s} = ct;
        else
            Cs(s) = spikeSync(sur, 'Interval', iv);
        end
    end
end
CtS = [];
if wantWindow
    CtS = [CtTmp{:}];   % each ct is nW x 1; horzcat gives nW x nsurr
end
end


function [mu, sd] = local_null(ts, iv, method, nsurr, dither)
N = numel(ts); L = iv(2) - iv(1);
c = zeros(nsurr,1);
for s = 1:nsurr
    sur = cell(N,1);
    for n = 1:N
        x = ts{n};
        switch lower(method)
            case 'shift',  sur{n} = sort(iv(1) + mod(x - iv(1) + rand*L, L));
            case 'dither', sur{n} = sort(x + dither*(2*rand(numel(x),1) - 1));
            otherwise, error('unknown NullMethod %s', method);
        end
    end
    c(s) = spikeSync(sur, 'Interval', iv);
end
mu = mean(c); sd = std(c);
end


function out = local_transientsFromSurrogates(Ct, tc, CtS, alpha)
% Detect windows where observed C exceeds a shift-surrogate band, as
% contiguous clusters. Takes the OBSERVED trace and the SHARED surrogate
% matrix already computed by local_sharedSurrogates -- no recomputation.
%
% Significance: a cluster counts as a transient if its run length meets or
% exceeds the (1-alpha) percentile of the max-run-length distribution taken
% across surrogates. Exact for the session's longest cluster; conservative
% for shorter ones in the same session -- report the count as a lower bound.
if isempty(CtS)
    out = struct('tStart', {}, 'tEnd', {}, 'nWin', {}, 'peakDC', {}, 'significant', {});
    return
end
nW = numel(tc); nsurr = size(CtS,2);
band = nan(nW,3);
for w = 1:nW
    band(w,:) = local_pct(CtS(w,:), [2.5 50 97.5]);
end

supra = Ct > band(:,3) & ~isnan(Ct);
[starts, ends] = local_runs(supra);

surRun = zeros(nsurr,1);
for s = 1:nsurr
    supraS = CtS(:,s) > band(:,3) & ~isnan(CtS(:,s));
    surRun(s) = local_maxRunLen(supraS);
end
cutoff = max(local_pct(surRun, 100*(1-alpha)), 1);

out = struct('tStart', {}, 'tEnd', {}, 'nWin', {}, 'peakDC', {}, 'significant', {});
for k = 1:numel(starts)
    idxR = starts(k):ends(k);
    out(k).tStart = tc(starts(k));
    out(k).tEnd = tc(ends(k));
    out(k).nWin = ends(k) - starts(k) + 1;
    out(k).peakDC = max(Ct(idxR) - band(idxR,2));
    out(k).significant = out(k).nWin >= cutoff;
end
end


function [starts, ends] = local_runs(tf)
tf = [false; tf(:); false];
d = diff(tf);
starts = find(d == 1);
ends = find(d == -1) - 1;
end


function r = local_maxRunLen(tf)
tf = logical(tf(:));
r = 0; c = 0;
for i = 1:numel(tf)
    if tf(i), c = c + 1; if c > r, r = c; end, else, c = 0; end
end
end


function y = local_pct(x, q)
x = sort(x(~isnan(x(:))));
n = numel(x);
if n == 0, y = nan(size(q)); return; end
if n == 1, y = repmat(x, size(q)); return; end
pos = (0.5:n-0.5) / n * 100;
y = interp1(pos, x, q, 'linear');
y(q < pos(1)) = x(1); y(q > pos(end)) = x(end);
end


function r = local_blank()
r = struct('ok', false, 'why', 'not run', 'file', '', ...
    'N', NaN, 'Nfull', NaN, 'M', NaN, 'dur', NaN, 'durInferred', NaN, ...
    'keepIdx', [], 'detIdx', [], 'Cfull', NaN, 'durTruncated', NaN, ...
    'windowAbs', [], ...
    'errStack', [], ...
    'rate', [], 'rateMed', NaN, 'rateRatio', NaN, ...
    'C', NaN, 'Cnull', NaN, 'Csd', NaN, 'dC', NaN, 'z', NaN, ...
    'CnullDither', NaN, 'CsdDither', NaN, ...
    'Cmat', [], 'prof', [], ...
    'Ct', [], 'tc', [], 'nSpk', [], 'nExcluded', NaN, ...
    'transients', [], ...
    'CmatchN', NaN, 'CmatchNsd', NaN);
end


function cf = local_cacheFile(dir_, f, o)
if isempty(dir_), cf = []; return; end
% Cache key includes the options that change what gets computed. Keying on
% file path alone means a call with different Width/MatchN/DetectTransients
% silently returns a stale result computed under different settings.
relevant = {o.RateRange, o.Width, o.Step, o.MinSpikes, o.NullMethod, ...
    o.Dither, o.NSurr, o.MatchN, o.NSubsets, o.DetectTransients, ...
    o.NSurrWindow, o.RunAlpha, o.MinUnits, o.MaxDuration, ...
    func2str_safe(o.Window), ...
    func2str_safe(o.Duration), func2str_safe(o.Exclude)};
h = local_hash(relevant);
k = regexprep(f, '[^\w]', '_');
if numel(k) > 100, k = k(end-99:end); end
cf = fullfile(dir_, sprintf('%s_%s.mat', k, h));
end

function s = func2str_safe(x)
if isempty(x), s = 'none'; elseif isa(x,'function_handle'), s = func2str(x);
else, s = mat2str(x); end
end

function h = local_hash(c)
s = char(jsonencode(c)); %#ok<CHARTEN>
h = sprintf('%08x', mod(sum(double(s) .* (1:numel(s))), 2^32));
end


function v = local_seedFromName(f)
% deterministic per-file seed so subset draws reproduce across reruns
v = mod(sum(double(f) .* (1:numel(f))), 2^31 - 1);
end


function s = local_short(f)
[a, b] = fileparts(f); [~, a] = fileparts(a);
s = fullfile(a, b);
end


function local_summary(T, o)
ok = [T.ok];
fprintf('\n%d/%d sessions usable\n', sum(ok), numel(T));
if ~any(ok), return; end
N = [T(ok).N]; dur = [T(ok).dur]; dC = [T(ok).dC];
fprintf('  N units    : median %d  range %d-%d\n', round(median(N)), min(N), max(N));
fprintf('  duration   : median %.0f s  range %.0f-%.0f\n', median(dur), min(dur), max(dur));
fprintf('  dC         : mean %+.4f  sd %.4f\n', mean(dC), std(dC));
if range(N) / max(median(N),1) > 0.5 && isempty(o.MatchN)
    fprintf(['  >> unit count varies by >50%% across sessions and MatchN is off.\n' ...
             '     C dilutes as 1/(N-1); set ''MatchN'' before comparing groups.\n']);
end
if any([T(ok).durInferred])
    fprintf('  >> duration inferred from last spike for some sessions (biases rate high).\n');
end
end