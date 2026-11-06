%% exportSyncSeries.m
% Write the windowed synchrony time series of every analysed session into the
% session's own folder as  synchronyTimeSeries.mat  (variable SYNC).
%
% HOW. For each config (width, step, MatchN) and context (persist drug-free
% period, acute pre-infusion epoch, acute post-infusion epoch) this re-calls
% syncBatch with EXACTLY the options runHazardAnalysis used, so every session
% is a cache hit (no recomputation). A session that is not in the cache IS
% computed, so a config you never ran will take as long as the original run.
%
% SYNC(k) -- one element per (context, width, step, MatchN), merged into any
% existing file (an element with the same .id is replaced, others are kept):
%   .id        e.g. 'persist_W60S5_N20'
%   .context   'persist' | 'acutePre' | 'acutePost'
%   .series
%       ts1    window-centre time, s from the start of the recording
%       ts2    window-centre time, s from the animal's FIRST ASV infusion
%              (negative = before it). See .meta.ts2Mode.
%       Ct     SPIKE-synchronization C per window (NaN: excluded / < MinSpikes)
%       e      excess synchrony (Ct - p50)/(p97.5 - p50); e = 1 is the band edge
%       band   [p2.5 p50 p97.5] of the shift-surrogate C per window
%       nSpk   spikes per window
%       nullRunLen  [p50 p90 p95 p99] of surrogate supra-band run length (windows)
%   .summary   C, Cnull, Csd, z, dC, N units, rate stats of the whole epoch
%   .meta      session, animal, condition, dates, infusion times, first-ASV
%              reference, unit indices (keepIdx, detIdx), window bounds, every
%              syncBatch option, seed, file checksums, MATLAB version, and a
%              plain-text recipe to reproduce.
%
% ts2. Needs the clock time at which each recording started. That comes from a
% HHMMSS token in the session folder name (CP12_250604_075757). When both the
% session and the animal's first-ASV session have one, ts2 is exact
% (ts2Mode 'clock'). When either is missing only the DATE is known, and ts2 is
% approximated as  days*86400 + (ts1 - first-ASV infusion time)  -- i.e. it
% assumes the recording started at the same time of day as the first-ASV one
% (ts2Mode 'date'; error up to +-24 h). meta.dayFromFirstASV is the exact
% integer day offset either way. Animals with no ASV session get NaN.
%
% REPRODUCIBILITY CAVEAT. The unit subset (detIdx) is seeded from the file
% path, so it reproduces if the path is unchanged. The surrogate time shifts
% use the global, unseeded rand, so the band and nullRunLen are NOT bit-for-bit
% reproducible: a rerun gives statistically equivalent but not identical
% values. That is why the band is stored here instead of being recomputable.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT = 'R:\TransPlasticity';
CONFIGS = [10 2 10; 20 4 10; 30 5 10; 60 5 10; 60 5 20];   % [width step MatchN]
CONTEXTS = {'persist','acutePre','acutePost'};
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
EXPOSURE_COND = 'ASV';
RATE_RANGE = [0.1 10]; MINSPK = 10; MINUNITS = 10; NSURRW = 50; RUNALPHA = 0.05;
MIN_EPOCH = 1800; LAT = 4*3600; POSTLEN = 2*3600;
OUT_NAME = 'synchronyTimeSeries.mat';
USE_PARALLEL = ~isempty(ver('parallel'));

fils = getAllExtFiles(ROOT, 'mat', 1);
fils = fils(contains(fils, 'spikes.cellinfo.mat'));
cnd = cellfun(@(f) conditionFromPath(f, 'TransPlasticity'), fils, 'UniformOutput', false);
ani = cellfun(@(f) animalFromPath(f, 'TransPlasticity'),    fils, 'UniformOutput', false);
dnm = cellfun(@(f) sessionDateFromPath(f),                  fils);
keep = ismember(cnd, CONDS) & ~cellfun(@isempty, ani) & isfinite(dnm);
fils = fils(keep); cnd = cnd(keep); ani = ani(keep); dnm = dnm(keep);

inf_t = nan(numel(fils), 1);
for i = 1:numel(fils), inf_t(i) = infusionTime(fils{i}); end
hasInf = ~isnan(inf_t);
bad = ismember(cnd, {'ASV','SAL'}) & ~hasInf;          % drug state unknown: same drop as the analysis
fils = fils(~bad); cnd = cnd(~bad); ani = ani(~bad); dnm = dnm(~bad);
inf_t = inf_t(~bad); hasInf = hasInf(~bad);

% recording start clock (datenum) when the folder name carries HHMMSS
recStart = cellfun(@local_recStart, fils);
% first-ASV reference per animal
ref = local_firstASV(fils, ani, cnd, dnm, inf_t, recStart, EXPOSURE_COND);
refSAL = local_firstASV(fils, ani, cnd, dnm, inf_t, recStart, 'SAL');   % ts3: first vehicle infusion

% durations (screen run: same options and cache as runHazardAnalysis)
% duration + usability per session (same rule as the old syncBatch screen, without
% running spikeSync over every unit -- that ran out of memory on large sessions)
[durs, scOK] = sessionScreen(fils, RATE_RANGE, MINUNITS, fullfile(CACHE_ROOT, 'screenLite'));

% epochs per context
w = [zeros(numel(fils),1), durs]; w(hasInf,2) = min(inf_t(hasInf), durs(hasInf));
usePersist = scOK & (w(:,2) - w(:,1)) >= MIN_EPOCH;
pre  = [zeros(numel(fils),1), inf_t];
post = [inf_t + LAT, min(durs, inf_t + LAT + POSTLEN)];
useAcute = hasInf & ismember(cnd, {'ASV','SAL'}) & scOK & ...
    (pre(:,2)-pre(:,1)) >= MIN_EPOCH & (post(:,2)-post(:,1)) >= MIN_EPOCH;

sums = local_checksums({'syncBatch','spikeSync','spikeSyncWindow','spikeSyncSurrogate','boutHazard','exportSyncSeries'});

nWritten = 0;
for c = 1:size(CONFIGS, 1)
    WIDTH = CONFIGS(c,1); STEP = CONFIGS(c,2); MATCHN = CONFIGS(c,3);
    cfgB = { 'RateRange',RATE_RANGE,'Width',WIDTH,'Step',STEP,'MinSpikes',MINSPK, ...
             'MinUnits',MINUNITS,'MatchN',MATCHN,'NullMethod','none', ...
             'DetectTransients',true,'NSurrWindow',NSURRW,'RunAlpha',RUNALPHA, ...
             'Parallel', USE_PARALLEL, 'CheckpointEvery',5 };
    for x = 1:numel(CONTEXTS)
        ctx = CONTEXTS{x};
        switch ctx
            case 'persist'
                use = usePersist; W = w; UW = [];
                cdir = 'hzPersist';
            case 'acutePre'
                use = useAcute; W = pre; cdir = 'hzAcutePre';
            case 'acutePost'
                use = useAcute; W = post; cdir = 'hzAcutePost';
        end
        if ~any(use), continue; end
        uf = fils(use);
        mapW = containers.Map(uf, num2cell(W(use,:), 2));
        args = [cfgB, {'Window', @(f) mapW(f), 'CacheDir', fullfile(CACHE_ROOT, cdir)}];
        if ~strcmp(ctx, 'persist')
            mapU = containers.Map(uf, arrayfun(@(k) [pre(k,:); post(k,:)], find(use), 'UniformOutput', false));
            args = [args, {'UnitWindow', @(f) mapU(f)}];
        end
        fprintf('\n== %s, width %g s, step %g s, MatchN %d: %d sessions ==\n', ctx, WIDTH, STEP, MATCHN, numel(uf));
        T = syncBatch(uf, args{:});
        idxAll = find(use);
        for q = 1:numel(T)
            r = T(q);
            if ~r.ok, fprintf('  skip %s: %s\n', r.file, r.why); continue; end
            i = idxAll(q);
            el = local_element(r, ctx, WIDTH, STEP, MATCHN, cfgB, ref, ani{i}, cnd{i}, dnm(i), ...
                recStart(i), inf_t(i), sums, LAT, POSTLEN, MIN_EPOCH, CACHE_ROOT, cdir, refSAL);
            local_write(fullfile(fileparts(r.file), OUT_NAME), el);
            nWritten = nWritten + 1;
        end
    end
end
fprintf('\nwrote %d series into session folders (%s)\n', nWritten, OUT_NAME);


%% =========================================================================
function el = local_element(r, ctx, W, S, N, cfgB, ref, animal, cond, date_, rs, infSec, sums, LAT, POSTLEN, MIN_EPOCH, cacheRoot, cdir, refSAL)
t0 = 0; if ~isempty(r.windowAbs), t0 = r.windowAbs(1); end
ts1 = r.tc(:) + t0;

% ts2 against the animal's first ASV infusion
rf = struct('session','','date',NaN,'recStart',NaN,'infusionSec',NaN,'clock',false);
k = find(strcmp(ref.animal, animal), 1);
if ~isempty(k), rf = ref.item(k); end
ts2 = nan(size(ts1)); mode = 'none'; dayOff = NaN;
if ~isempty(rf.session) && isfinite(rf.infusionSec)
    dayOff = round(date_ - rf.date);
    if isfinite(rs) && rf.clock
        ts2 = ((rs + ts1/86400) - (rf.recStart + rf.infusionSec/86400)) * 86400;
        mode = 'clock';
    else
        ts2 = dayOff*86400 + (ts1 - rf.infusionSec);
        mode = 'date';
    end
end

band = nan(numel(ts1), 3); if isfield(r,'band') && ~isempty(r.band), band = r.band; end
den = band(:,3) - band(:,2);
e = (r.Ct(:) - band(:,2)) ./ den; e(~(den > 0)) = NaN;
nrl = []; if isfield(r,'nullRunLen'), nrl = r.nullRunLen; end

el = struct();
el.id = sprintf('%s_W%gS%g_N%d', ctx, W, S, N);
el.context = ctx;
[ts3, mode3, day3, rf3] = ts2FromReference(ts1, animal, date_, rs, refSAL);
el.series = struct('ts1',ts1, 'ts2',ts2, 'ts3',ts3, 'Ct',r.Ct(:), 'e',e, 'band',band, ...
    'nSpk',r.nSpk(:), 'nullRunLen',nrl);
el.summary = struct('C',r.C, 'Cnull',r.Cnull, 'Csd',r.Csd, 'dC',r.dC, 'z',r.z, ...
    'N',r.N, 'Nfull',r.Nfull, 'rateMed',r.rateMed, 'durEpoch',r.dur);

seedV = mod(sum(double(r.file) .* (1:numel(r.file))), 2^31 - 1);
p = struct('Width',W,'Step',S,'MatchN',N,'cfg',{cfgB}, 'windowAbs',r.windowAbs, ...
    'unitWindow', local_getfield(r,'unitWindow'), ...
    'epochRule', sprintf('LAT %g s, POSTLEN %g s, MIN_EPOCH %g s', LAT, POSTLEN, MIN_EPOCH), ...
    'keepIdx',r.keepIdx, 'detIdx',local_getfield(r,'detIdx'), ...
    'seed',seedV, 'seedRule','mod(sum(double(file).*(1:numel(file))), 2^31-1) on the full path string; RandStream mt19937ar; randperm(Nfull, MatchN)', ...
    'surrogateRNG','shift surrogates use the global unseeded rand: band and nullRunLen are not bit-reproducible');
m = struct('file',r.file, 'folder',fileparts(r.file), 'animal',animal, 'condition',cond, ...
    'sessionDate',date_, 'sessionDateStr',datestr(date_,'yyyy-mm-dd'), ...
    'recStart',rs, 'recStartStr',local_dstr(rs), ...
    'infusionSec',infSec, ...
    'firstASV',rf, 'dayFromFirstASV',dayOff, 'ts2Mode',mode, ...
    'firstSAL',rf3, 'dayFromFirstSAL',day3, 'ts3Mode',mode3, ...
    'ts2Note',['ts1 = s from recording start; ts2 = s from the animal''s first ASV infusion. ' ...
              'ts3 = s from the animal''s first SAL (vehicle) infusion, same modes. ' ...
              'ts2Mode clock = exact; date = assumes same time-of-day start as the first-ASV recording (+-24 h); none = animal has no ASV session.'], ...
    'params',p, 'cacheDir',fullfile(cacheRoot, cdir), 'checksums',sums, ...
    'matlab',version, 'created',datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
    'bandColumns','[p2.5 p50 p97.5] of the shift-surrogate C per window', ...
    'excessDef','e = (Ct - p50) / (p97.5 - p50)');
m.reproduce = sprintf(['syncBatch(file, cfg{:}, ''Window'', windowAbs [, ''UnitWindow'', unitWindow]); ' ...
    'cfg = params.cfg (Width %g, Step %g, MatchN %d, ...). Same path string => same unit subset. ' ...
    'e and bouts: boutHazard(r, step, ''null95'').'], W, S, N);
el.meta = m;
end


function ref = local_firstASV(fils, ani, cnd, dnm, inf_t, recStart, expCond)
ua = unique(ani(:));
ref = struct('animal',{ua}, 'item',repmat(struct('session','','date',NaN,'recStart',NaN,'infusionSec',NaN,'clock',false), numel(ua), 1));
for u = 1:numel(ua)
    m = find(strcmp(ani, ua{u}) & strcmp(cnd, expCond) & isfinite(inf_t));
    if isempty(m), continue; end
    key = dnm(m) + 1e-6*(1:numel(m))';
    key(isfinite(recStart(m))) = recStart(m(isfinite(recStart(m))));
    [~, o] = min(key); j = m(o);
    ref.item(u) = struct('session',fileparts(fils{j}), 'date',dnm(j), 'recStart',recStart(j), ...
        'infusionSec',inf_t(j), 'clock',isfinite(recStart(j)));
end
end


function d = local_recStart(f)
% datenum of recording start from CP12_250604_075757 (YYMMDD_HHMMSS), else NaN
d = NaN;
tok = regexp(f, '_(\d{6})_(\d{6})(?:\\|/|$)', 'tokens');
if isempty(tok), return; end
s = tok{end}{1}; h = tok{end}{2};
yy = str2double(s(1:2)); mm = str2double(s(3:4)); dd = str2double(s(5:6));
hh = str2double(h(1:2)); mi = str2double(h(3:4)); ss = str2double(h(5:6));
if mm < 1 || mm > 12 || dd < 1 || dd > 31 || hh > 23 || mi > 59 || ss > 59, return; end
d = datenum(2000 + yy, mm, dd, hh, mi, ss);
end


function s = local_dstr(d)
if isfinite(d), s = datestr(d, 'yyyy-mm-dd HH:MM:SS'); else, s = ''; end
end


function v = local_getfield(s, f)
if isfield(s, f), v = s.(f); else, v = []; end
end


function out = local_checksums(names)
out = struct();
for k = 1:numel(names)
    p = which(names{k});
    out.(names{k}) = struct('path',p,'md5','');
    if isempty(p), continue; end
    try
        fid = fopen(p, 'r'); b = fread(fid, Inf, '*uint8'); fclose(fid);
        md = java.security.MessageDigest.getInstance('MD5');
        md.update(b);
        out.(names{k}).md5 = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'))', 1, []));
    catch
    end
end
end


function local_write(fn, el)
% merge into an existing file: same id replaced, other ids kept
SYNC = el;
README = ['SYNC(k): one synchrony time-series set per (context, width, step, MatchN). ' ...
          'series.ts1 = s from recording start, series.ts2 = s from the animal''s first ASV infusion; ' ...
          'Ct, e, band, nSpk share these time axes. meta holds everything needed to reproduce. ' ...
          'Written by exportSyncSeries.m.'];
if exist(fn, 'file')
    try
        L = load(fn, 'SYNC');
        old = L.SYNC;
        old = old(~strcmp({old.id}, el.id));
        if ~isempty(old), SYNC = [local_conform(old, el), el]; end
    catch ME
        warning('exportSyncSeries:merge', '%s unreadable (%s): overwriting', fn, ME.message);
    end
end
[~, ix] = sort({SYNC.id}); SYNC = SYNC(ix); %#ok<NASGU>
tmp = [fn '.tmp'];
save(tmp, 'SYNC', 'README', '-v7');
movefile([tmp ''], fn, 'f');
end


function old = local_conform(old, el)
% older elements may lack newer fields; add them empty so the array concatenates
fo = fieldnames(el);
for k = 1:numel(fo)
    if ~isfield(old, fo{k}), [old.(fo{k})] = deal([]); end
end
old = orderfields(old, el);
end