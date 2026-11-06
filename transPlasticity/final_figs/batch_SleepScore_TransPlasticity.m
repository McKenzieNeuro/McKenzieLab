%% batch_SleepScore_TransPlasticity.m
% Batch buzcode SleepScoreMaster over TransPlasticity sessions.
%
%  Theta : pinned to the ripple-detection channel (0-based; buzcode is 0-based,
%          so no index shift).
%  SW    : auto-picked on each animal's first Baseline session, then pinned for
%          every other session of that animal (removes condition-correlated
%          channel-choice variance).
%  Naming: SleepScoreMaster takes basename = folder name and requires
%          basePath/basename.lfp + .xml. If the session folder is not named
%          BASENAME, a local junction named BASENAME -> session folder is made,
%          so all outputs land in the real session folder as
%          BASENAME.SleepState.states.mat etc. (next to the ripples file).
%
%  Every session is logged; requested vs. used channels are checked post hoc.
%  Requires MATLAB R2017b+ (isfile/isfolder), Windows (mklink).

%% ======================= config =======================
ROOT     = 'R:\TransPlasticity\Ephys';
CONDS    = {'Baseline','ASV','ASV24','SAL','SAL24'};
BASENAME = 'amplifier_analogin_auxiliary_int16';
RIPFILE  = [BASENAME '.ripples.events.mat'];
STAGE    = 'D:\data';          % local dir for basename-named links
BUZCODE  = 'R:\Analysis\buzcode';       % <-- edit

OPT.overwrite       = false;        % true: re-pick channels/recompute metrics AND discard TheStateEditor manual scoring
OPT.motionSource    = 'EMGfromLFP'; % or 'Accelerometer'
OPT.accelChans      = [];           % 0-based aux accel channels, needed only for 'Accelerometer'
OPT.pinSWtoBaseline = true;
OPT.dryRun          = false;        % true: resolve paths/channels and log only
OPT.ignoretimeFn    = @(basePath, cond) [];   % return [t0 t1; ...] (s) to exclude, e.g. infusion/handling
OPT.makeBlindQC     = true;         % stage condition-blind links for TheStateEditor review (needs links)
OPT.stageMode       = 'auto';       % 'link' | 'copy' | 'auto' (link, else local copy). R: network drive => copy

addpath(genpath(BUZCODE));
if ~isfolder(STAGE), mkdir(STAGE); end

%% ======================= session list =======================
fils = getAllExtFiles(ROOT, 'mat', 1);
fils = fils(contains(fils, 'spikes.cellinfo.mat'));
[~, rootName] = fileparts(ROOT);            % 'Ephys'
[ani, cnd, dnm] = parseTPPath(fils, rootName);
keep = ismember(cnd, CONDS) & ~cellfun(@isempty, ani) & isfinite(dnm);
if ~any(keep)
    error('No sessions survived the filter. Parsed conditions: {%s}', ...
          strjoin(unique(cnd(~cellfun(@isempty, cnd))), ', '));
end
fils = fils(keep); ani = ani(keep); cnd = cnd(keep); dnm = dnm(keep);
%%
basePaths = cellfun(@fileparts, fils, 'uni', 0);
[basePaths, iu] = unique(basePaths(:), 'stable');
ani = ani(iu); cnd = cnd(iu); dnm = dnm(iu);
nS  = numel(basePaths);

% Baselines first (by day) so the SW pin exists before drug sessions run
isBL = strcmp(cnd(:), 'Baseline');
iBL  = find(isBL);  [~, o] = sort(dnm(iBL));  iBL = iBL(o);
order = [iBL; find(~isBL)];

%% ======================= main loop =======================
rows     = repmat(emptyRow(), nS, 1);
swPinned = containers.Map('KeyType', 'char', 'ValueType', 'double');

for k = order(:)'
    r = emptyRow();
    r.animal = string(ani{k}); r.cond = string(cnd{k}); r.day = dnm(k);
    r.basePath = string(basePaths{k});
    fprintf('\n[%d/%d] %s | %s | day %g\n', find(order==k), nS, ani{k}, cnd{k}, dnm(k));
    if isfile(fullfile(basePaths{k}, RIPFILE))
         v    = load(fullfile(basePaths{k}, RIPFILE));
         BASENAME = 'amplifier_analogin_auxiliary_int16';
    else

        f = getAllExtFiles(basePaths{k},'mat',0);
        kp = contains(f,'ripples.events.mat');
        if sum(kp) ==1
        

        RIPFILE1  = f{find(kp)};
        end
         v    = load( RIPFILE1);
         BASENAME = bz_BasenameFromBasepath(basePaths{k});
    end

    if (OPT.overwrite || ~isfile(fullfile(basePaths{k}, [BASENAME '.SleepState.states.mat']))) ...
            && (isfile(fullfile(basePaths{k}, RIPFILE)) || isfile(RIPFILE1))
    %try
        % ---- 1. theta channel from ripple detector (0-based) ----
       
        thCh = v.ripples.detectorinfo.detectionparms.channel;
        assert(isscalar(thCh) && thCh == round(thCh), 'ripple channel is not a scalar integer');
        r.thReq = thCh;

        % ---- 2. basename-conformant path ----
        % A previous interrupted run may leave an onCleanup object alive.
        % Destroy it BEFORE recreating the same staging path.
        clear cleaner
        mode = OPT.stageMode;  if OPT.dryRun, mode = 'xml'; end   % dry run: don't move .lfp
        [ssmPath, cleanupFn, stg] = stageBasePath(basePaths{k}, BASENAME, STAGE, sprintf('s%04d', k), mode);
        cleaner = onCleanup(cleanupFn); %#ok<NASGU>
        r.stage = stg.mode;
        assert(isfolder(ssmPath), 'Staged directory disappeared: %s', ssmPath);
        assert(isfile(fullfile(ssmPath, [BASENAME '.xml'])), ...
               'XML is missing from staged directory: %s', ssmPath);
        if isfile(fullfile(basePaths{k}, [BASENAME '.sessionInfo.mat']))
            assert(isfile(fullfile(ssmPath, [BASENAME '.sessionInfo.mat'])), ...
                   'Source sessionInfo exists but staged copy/link is missing: %s', ssmPath);
        end

        % ---- 3. channel sets ----
        si     = bz_getSessionInfo(ssmPath, 'noPrompts', true);
        ampCh  = amplifierChannels(si);
        nonAmp = setdiff(si.channels, ampCh);          % analogin + aux in the combined file
        % Read Skip directly: cached sessionInfo may omit XML exclusions.
        xml = LoadXml(fullfile(basePaths{k}, [BASENAME '.xml']));
        bad = xmlSkippedChannels(xml);
        if isfield(si, 'badchannels')
            bad = unique([bad si.badchannels(:)']);
        end
        reject = unique([nonAmp(:)' bad]);
        r.nAmp = numel(ampCh); r.nTot = numel(si.channels);

        assert(ismember(thCh, ampCh), 'theta ch %d is not in an xml spike group', thCh);
        assert(~ismember(thCh, bad),  'theta ch %d is marked bad in the xml', thCh);
        if thCh == 0
            warning('theta ch is 0 (= SleepScoreMaster''s default); post-run THchanID check will confirm it was honored');
        end

        a = ani{k};
        if OPT.pinSWtoBaseline && isKey(swPinned, a)
            swCh = swPinned(a);              r.swSource = "pinned";
        else
            swCh = setdiff(ampCh, reject);   % candidate pool: good amplifier channels
            if strcmp(cnd{k}, 'Baseline'), r.swSource = "auto_baseline";
            else,                          r.swSource = "auto_noBaseline"; end
        end

        assert(~isempty(swCh), 'No good amplifier channels remain for SW scoring');
        assert(all(ismember(swCh, ampCh)) && ~any(ismember(swCh, reject)), ...
               'SW channel selection includes an excluded/non-amplifier channel: %s', mat2str(swCh));

        args = {'overwrite', OPT.overwrite, 'noPrompts', true, ...
                'ThetaChannels', thCh, 'SWChannels', swCh, ...
                'rejectChannels', reject, 'MotionSource', OPT.motionSource, ...
                'ignoretime', OPT.ignoretimeFn(basePaths{k}, cnd{k})};
        if strcmpi(OPT.motionSource, 'Accelerometer')
            assert(~isempty(OPT.accelChans), 'set OPT.accelChans for Accelerometer motion');
            args = [args, {'AccelerChans', OPT.accelChans}]; %#ok<AGROW>
        end

        if OPT.dryRun
            r.status = "dry"; rows(k) = r; clear cleaner; continue
        end

        % ---- 4. score ----
        % SleepScoreMaster derives basename from the final directory name.
        % Its default savedir is the parent, so outputs land in ssmPath.
        SleepScoreMaster(ssmPath, args{:});
        close all

        % ---- 5. verify channels actually used ----
        S = load(fullfile(ssmPath, [BASENAME '.SleepScoreLFP.LFP.mat']));
        r.thUsed = S.SleepScoreLFP.THchanID;
        r.swUsed = S.SleepScoreLFP.SWchanID;

        if r.swSource == "auto_baseline" && ~isKey(swPinned, a)
            swPinned(a) = r.swUsed;
        end

        ss = load(fullfile(ssmPath, [BASENAME '.SleepState.states.mat']));
        [r.fWAKE, r.fNREM, r.fREM, r.hrsScored] = stateFractions(ss.SleepState.ints);
        r.manual = isfield(ss.SleepState.detectorinfo, 'LastManualUpdate');

        if any(ismember([r.thUsed r.swUsed], reject))
            r.status = "BAD_CHANNEL_USED";
            r.msg = "Cached scoring used an excluded channel; recompute with overwrite=true (discards manual scoring).";
            warning('%s', char(r.msg));
        elseif r.thUsed ~= r.thReq
            r.status = "TH_MISMATCH";   % usually stale SleepScoreLFP with overwrite=false
        elseif r.swSource == "pinned" && r.swUsed ~= swCh
            r.status = "SW_MISMATCH";
        else
            r.status = "ok";
        end
        syncBack(stg, BASENAME);   % copy mode: outputs -> R: before local delete
  %  catch ME
  %      r.status = "error"; r.msg = string(ME.message);
  %      fprintf(2, '  ERROR: %s\n', ME.message);
  % end
    rows(k) = r;
    clear cleaner                          % remove link now, not at next iteration
    end
end

%% ======================= log + consistency checks =======================
T = struct2table(rows(order));
logf = fullfile(ROOT, ['sleepscore_log_' datestr(now, 'yyyymmdd_HHMM') '.csv']);
writetable(T, logf);
fprintf('\nLog: %s\n', logf);
disp(groupcounts(T, 'status'));

% ripple (=theta) channel should be constant within animal
ua = unique(T.animal);
for i = 1:numel(ua)
    th = unique(T.thReq(T.animal == ua(i) & isfinite(T.thReq)));
    if numel(th) > 1
        warning('%s: ripple channel differs across sessions: %s', ua(i), mat2str(th'));
    end
end

%% ======================= blinded QC staging =======================
% Junctions with neutral names so TheStateEditor review can't read condition
% off the path. Key saved separately; open with TheStateEditor(fullfile(link, BASENAME)).
if OPT.makeBlindQC && ~OPT.dryRun
    okIdx = find(T.status == "ok");
    perm  = okIdx(randperm(numel(okIdx)));
    qcDir = fullfile(STAGE, 'blindQC');  if ~isfolder(qcDir), mkdir(qcDir); end
    code  = strings(numel(perm), 1);
    try
        for j = 1:numel(perm)
            code(j) = sprintf('QC%03d', j);
            p = fullfile(qcDir, code(j));  if ~isfolder(p), mkdir(p); end
            stageBasePath(char(T.basePath(perm(j))), BASENAME, p, '', 'link');   % link persists
        end
    catch ME
        warning('Blind QC staging needs directory links; skipped (%s)', ME.message);
        code = code(1:j-1);  perm = perm(1:j-1);
    end
    key = table(code, T.basePath(perm), T.animal(perm), T.cond(perm), ...
                'VariableNames', {'code', 'basePath', 'animal', 'cond'});
    writetable(key, fullfile(ROOT, ['sleepscore_blindQC_KEY_' datestr(now, 'yyyymmdd_HHMM') '.csv']));
    fprintf('Blinded QC links in %s (key written to ROOT; keep it away from the scorer)\n', qcDir);
end

%% ======================= local functions =======================
function [ani, cnd, dnm] = parseTPPath(fils, rootName)
% Layout: <Ephys root>\<animal>\<condition>\<session>\<recording>\...
% session = <animal>_yymmdd; recording = <animal>_yymmdd_HHMMSS
%   R:\TransPlasticity\Ephys\CP19\ASV24\CP19_250604\... -> 'CP19', 'ASV24', datenum(2025,6,4)
% Unparseable paths return ''/'' /NaN and are dropped by the keep filter.
n   = numel(fils);
ani = repmat({''}, n, 1);  cnd = repmat({''}, n, 1);  dnm = nan(n, 1);
for i = 1:n
    parts = regexp(fils{i}, '[\\/]', 'split');
    parts = parts(~cellfun(@isempty, parts));
    k = find(strcmpi(parts, rootName), 1, 'first');
    if isempty(k) || k+3 > numel(parts), continue; end
    ani{i} = parts{k+1};
    cnd{i} = parts{k+2};
    tok = regexp(parts{k+3}, '(\d{6})$', 'tokens', 'once');
    if ~isempty(tok), dnm(i) = datenum(tok{1}, 'yymmdd'); end
end
end

function [p, cleanupFn, stg] = stageBasePath(src, BASENAME, parentDir, sub, mode)
% Return a path whose last folder is BASENAME and which resolves to (or mirrors) src.
%   mode 'link' : junction, then directory symlink; error if both fail
%        'copy' : local copy of .xml, .lfp/.eeg and existing SSM intermediates
%        'auto' : link, else copy
%        'xml'  : local copy of .xml/.sessionInfo only (dry run)
% Link removal uses 'rmdir' WITHOUT /S: deletes the link, never the target.
% Local copies are deleted only if they carry the marker file written here.
if nargin < 5, mode = 'auto'; end
cleanupFn = @() [];
stg = struct('mode', "direct", 'src', src, 'dst', src);
[~, folder] = fileparts(src);
assert(isfile(fullfile(src, [BASENAME '.xml'])), 'missing %s.xml in %s', BASENAME, src);
lfpExt = '.lfp';
if ~isfile(fullfile(src, [BASENAME lfpExt])), lfpExt = '.eeg'; end
assert(isfile(fullfile(src, [BASENAME lfpExt])), ...
       'missing %s.lfp/.eeg in %s (run bz_LFPfromDat on the staged path first)', BASENAME, src);
if strcmp(folder, BASENAME), p = src; return; end

parent = fullfile(parentDir, sub);
if ~isfolder(parent), mkdir(parent); end
p = fullfile(parent, BASENAME);
stg.dst = p;
system(sprintf('rmdir "%s" 2>nul', p));                 % stale link (no /S)
removeLocalCopy(p);                                     % stale copy (marker-guarded)

if any(strcmp(mode, {'link', 'auto'}))
    [st, m1] = system(sprintf('mklink /J "%s" "%s"', p, src));
    m2 = '';
    if st ~= 0, [st, m2] = system(sprintf('mklink /D "%s" "%s"', p, src)); end
    if st == 0
        stg.mode  = "link";
        cleanupFn = @() system(sprintf('rmdir "%s"', p));
        return
    end
    if strcmp(mode, 'link')
        error(['Could not link %s -> %s\n  junction: %s\n  symlink : %s\n' ...
               'Junctions need a local NTFS target; symlinks need admin or Developer Mode.'], ...
               p, src, strtrim(m1), strtrim(m2));
    end
end

% ---- copy staging ----
mkdir(p);
fclose(fopen(fullfile(p, stageMarker()), 'w'));
pats = {[BASENAME '.xml'], [BASENAME '.sessionInfo.mat']};
if ~strcmp(mode, 'xml')
    % intermediates carried over so overwrite=false / manual-scoring semantics match an in-place run
    pats = [pats, {[BASENAME lfpExt], [BASENAME '.SleepScoreLFP.LFP.mat'], ...
                   [BASENAME '.EMGFromLFP.LFP.mat'], [BASENAME '.SleepState*.mat']}];
end
t0 = tic;  nb = 0;
for i = 1:numel(pats)
    d = dir(fullfile(src, pats{i}));
    for j = 1:numel(d)
        copyfile(fullfile(src, d(j).name), p);  nb = nb + d(j).bytes;
    end
end
fprintf('  staged copy: %.2f GB in %.0f s\n', nb/2^30, toc(t0));
stg.mode  = string(mode);  if mode == "auto", stg.mode = "copy"; end
cleanupFn = @() removeLocalCopy(p);
end

function syncBack(stg, BASENAME)
% Copy SSM outputs from a local copy-stage back to the source folder.
if stg.mode ~= "copy", return; end              % link/direct: outputs already in place
pats = {[BASENAME '.SleepState*.mat'], [BASENAME '.SleepScoreLFP.LFP.mat'], ...
        [BASENAME '.EMGFromLFP.LFP.mat']};
for i = 1:numel(pats)
    d = dir(fullfile(stg.dst, pats{i}));
    for j = 1:numel(d), copyfile(fullfile(stg.dst, d(j).name), stg.src); end
end
% figure folder location not verified: check basePath and its parent
figDirs = {fullfile(stg.dst, 'StateScoreFigures'), fullfile(fileparts(stg.dst), 'StateScoreFigures')};
for i = 1:numel(figDirs)
    if isfolder(figDirs{i}), copyfile(figDirs{i}, fullfile(stg.src, 'StateScoreFigures')); end
end
end

function removeLocalCopy(p)
% Delete a copy-stage folder only if it carries our marker (never a link or source).
if isfolder(p) && isfile(fullfile(p, stageMarker()))
    rmdir(p, 's');
end
end

function m = stageMarker()
m = '_SSM_LOCAL_COPY_safe_to_delete';
end

function bad = xmlSkippedChannels(xml)
% LoadXml Channels are zero-based IDs; Skip entries correspond by position.
assert(isfield(xml, 'AnatGrps'), 'XML has no AnatGrps channel definitions');
bad = [];
for g = 1:numel(xml.AnatGrps)
    ch = xml.AnatGrps(g).Channels(:)';
    if ~isfield(xml.AnatGrps, 'Skip') || isempty(xml.AnatGrps(g).Skip)
        continue
    end
    skip = xml.AnatGrps(g).Skip(:)';
    assert(numel(skip) == numel(ch), ...
           'XML group %d: Channels and Skip lengths differ', g);
    assert(all(isfinite(skip)) && all(skip == 0 | skip == 1), ...
           'XML group %d: Skip must contain only 0 or 1', g);
    bad = [bad ch(skip == 1)]; %#ok<AGROW>
end
bad = unique(bad);
end

function ch = amplifierChannels(si)
% ASSUMPTION: xml spike groups contain amplifier channels only, so
% analogin/aux channels of the combined file are excluded.
if isfield(si, 'spikeGroups') && ~isempty(si.spikeGroups.groups)
    g = si.spikeGroups.groups;
else
    g = {si.AnatGrps.Channels};
end
ch = unique(cell2mat(cellfun(@(x) x(:)', g, 'uni', 0)));
end

function [fW, fN, fR, hrs] = stateFractions(ints)
nm  = {'WAKEstate', 'NREMstate', 'REMstate'};
dur = zeros(1, 3);
for i = 1:3
    if isfield(ints, nm{i}) && ~isempty(ints.(nm{i}))
        dur(i) = sum(diff(ints.(nm{i}), 1, 2));
    end
end
tot = sum(dur);
f = dur / tot;  fW = f(1); fN = f(2); fR = f(3);  hrs = tot / 3600;
end

function r = emptyRow()
r = struct('animal', "", 'cond', "", 'day', NaN, 'basePath', "", ...
           'thReq', NaN, 'thUsed', NaN, 'swSource', "", 'swUsed', NaN, ...
           'nAmp', NaN, 'nTot', NaN, 'fWAKE', NaN, 'fNREM', NaN, 'fREM', NaN, ...
           'hrsScored', NaN, 'manual', false, 'stage', "", 'status', "", 'msg', "");
end