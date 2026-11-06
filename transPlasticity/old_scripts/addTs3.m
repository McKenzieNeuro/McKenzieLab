%% addTs3.m
% Add ts3 = s from the animal's FIRST SAL (vehicle) infusion to the files
% already written into the session folders, without recomputing anything:
%   spikeDistanceProfile.mat  SDIST.bin.ts3, SDIST.exact.ts3
%   rippleEvents.mat          RIPPLES.ts3 / ts3Peak / ts3Off, RIPPLES.spikes.ts3
%   synchronyTimeSeries.mat   SYNC(k).series.ts3
% plus meta.firstSAL, meta.dayFromFirstSAL, meta.ts3Mode. Same rule as ts2
% (ts2FromReference): exact when both recordings carry a start clock time in
% the folder name, otherwise date-resolved (+-24 h); NaN for an animal with no
% SAL infusion. The runners now write ts3 themselves; this is for files that
% already exist.

ROOT = 'R:\TransPlasticity';
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
FILES = {'spikeDistanceProfile.mat','rippleEvents.mat','synchronyTimeSeries.mat'};

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
refSAL = firstASVReference(fils, ani, cnd, dnm, inf_t, recStart, 'SAL');
nRef = sum(cellfun(@(s) ~isempty(s), {refSAL.item.session}));
fprintf('%d sessions; %d of %d animals have a first-SAL reference\n', numel(fils), nRef, numel(refSAL.animal));

nUpd = zeros(1, numel(FILES));
for i = 1:numel(fils)
    d = fileparts(fils{i});
    F = @(t) ts2FromReference(t, ani{i}, dnm(i), recStart(i), refSAL);
    [~, mode3, day3, rf3] = F(0);
    for k = 1:numel(FILES)
        fn = fullfile(d, FILES{k});
        if ~exist(fn, 'file'), continue; end
        try
            switch k
                case 1
                    L = load(fn); X = L.SDIST;
                    X.bin.ts3 = F(X.bin.ts1); X.exact.ts3 = F(X.exact.ts1);
                    X.meta = local_meta(X.meta, rf3, day3, mode3);
                    SDIST = X; README = local_readme(L); %#ok<NASGU>
                    save([fn '.tmp'], 'SDIST', 'README', '-v7.3'); movefile([fn '.tmp'], fn, 'f');
                case 2
                    L = load(fn); X = L.RIPPLES;
                    X.ts3 = F(X.ts1); X.ts3Peak = F(X.ts1Peak); X.ts3Off = F(X.ts1Off);
                    if ~ismember('ts3', X.spikes.Properties.VariableNames)
                        X.spikes = addvars(X.spikes, F(X.spikes.ts1), 'After', 'ts2', 'NewVariableNames', 'ts3');
                    else
                        X.spikes.ts3 = F(X.spikes.ts1);
                    end
                    X.meta = local_meta(X.meta, rf3, day3, mode3);
                    RIPPLES = X; README = local_readme(L); %#ok<NASGU>
                    save([fn '.tmp'], 'RIPPLES', 'README', '-v7.3'); movefile([fn '.tmp'], fn, 'f');
                case 3
                    L = load(fn); X = L.SYNC;
                    for j = 1:numel(X)
                        X(j).series.ts3 = F(X(j).series.ts1);
                        X(j).meta = local_meta(X(j).meta, rf3, day3, mode3);
                    end
                    SYNC = X; README = local_readme(L); %#ok<NASGU>
                    save([fn '.tmp'], 'SYNC', 'README', '-v7'); movefile([fn '.tmp'], fn, 'f');
            end
            nUpd(k) = nUpd(k) + 1;
        catch ME
            fprintf('  FAILED %s: %s\n', fn, ME.message);
        end
    end
end
for k = 1:numel(FILES), fprintf('%-26s %d files updated\n', FILES{k}, nUpd(k)); end


%% =========================================================================
function m = local_meta(m, rf, day, mode)
m.firstSAL = rf; m.dayFromFirstSAL = day; m.ts3Mode = mode;
m.ts3Note = ['ts3 = s from the animal''s first SAL (vehicle) infusion. Same modes as ts2: clock = exact; ' ...
             'date = assumes the same time-of-day start as the first-SAL recording (+-24 h); none = no SAL session.'];
end

function r = local_readme(L)
add = ' ts3 = s from the animal''s first SAL infusion (added by addTs3.m).';
if ~isfield(L, 'README'), r = strtrim(add);
elseif contains(L.README, 'ts3 ='), r = L.README;
else, r = [L.README add]; end
end
