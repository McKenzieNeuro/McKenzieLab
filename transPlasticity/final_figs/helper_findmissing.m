% listASVWindows  Per-session table behind the ASV groups of compareRateWindowsByCellType.
% Prints, for every session of an animal that has a first-ASV reference: condition, ts2 mode, offset of the
% recording start from the first-ASV infusion, recording length, whether spikes + cell_metrics + SleepState
% exist, and the seconds it contributes to PreASV, ASVshort (2-8 h) and ASVlong (>=24 h).
% Sessions contributing >= 1800 s to ASVshort are starred; animals with two of them are flagged.
% Saves the table to listASVWindows.mat (variable R) next to the cache folder for diffing between runs.
ROOT = 'R:\TransPlasticity'; CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
if exist('SHARED', 'var'), ROOT = SHARED.root; CACHE_ROOT = SHARED.cacheRoot; end
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
SHORT_WIN = [2 8] * 3600; LONG_AFTER = 24 * 3600; MIN_EPOCH = 1800;

fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'spikes.cellinfo.mat'));
cnd = cellfun(@(f) conditionFromPath(f, 'TransPlasticity'), fa, 'UniformOutput', false);
ani = cellfun(@(f) animalFromPath(f, 'TransPlasticity'),    fa, 'UniformOutput', false);
dnm = cellfun(@sessionDateFromPath, fa);
k = ismember(cnd, CONDS) & ~cellfun(@isempty, ani) & isfinite(dnm);
fprintf('%d spikes files found, %d pass the condition / animal / date filter\n', numel(fa), sum(k));
fprintf('rejected by the filter:\n'); fprintf('   %s\n', fa{~k});
fa = fa(k); cnd = cnd(k); ani = ani(k); dnm = dnm(k);
inf_t = cellfun(@infusionTime, fa);
recStart = cellfun(@recordingStartFromPath, fa);
ref = firstASVReference(fa, ani, cnd, dnm, inf_t, recStart, 'ASV');

n = numel(fa);
mode = repmat({''}, n, 1); off = nan(n, 1); len = nan(n, 1); sPre = nan(n, 1); sShort = nan(n, 1); sLong = nan(n, 1);
hasCM = false(n, 1); hasSS = false(n, 1); first = false(n, 1);
for i = 1:n
    fld = fileparts(fa{i});
    hasCM(i) = ~isempty(dir(fullfile(fld, '*cell_metrics.cellinfo.mat')));
    hasSS(i) = ~isempty(dir(fullfile(fld, '*.SleepState.states.mat')));
    [o, md, ~, rf] = ts2FromReference(0, ani{i}, dnm(i), recStart(i), ref);
    mode{i} = md; off(i) = o / 3600; first(i) = ~isempty(rf.session) && strcmp(fld, rf.session);
    if isempty(rf.session), continue; end
    L = load(fa{i}, 'spikes'); ts = L.spikes.times(:);
    dur = max(cellfun(@(x) max([x(~isnan(x)); -Inf]), ts)); len(i) = dur / 3600;
    if ~isfinite(dur), continue; end
    g = (0.5:1:dur)'; t2 = ts2FromReference(g, ani{i}, dnm(i), recStart(i), ref);
    later = isfinite(inf_t(i)) && ~first(i); cut = Inf; if later, cut = inf_t(i); end
    sPre(i)   = sum(t2 < 0);
    sShort(i) = sum(t2 >= SHORT_WIN(1) & t2 < SHORT_WIN(2) & g < cut);
    sLong(i)  = sum(t2 >= LONG_AFTER & g < cut);
end
star = sShort >= MIN_EPOCH;
[~, o] = sortrows(table(ani, dnm));
R = table(ani, cnd, dnm, mode, first, off, len, inf_t / 3600, sPre / 3600, sShort / 3600, sLong / 3600, hasCM, hasSS, star, fa, ...
    'VariableNames', {'Animal','Cond','Date','Mode','FirstASV','OffsetH','LenH','InfusionH','PreH','ShortH','LongH','hasCM','hasSS','ShortOK','File'});
R = R(o, :);
R = R(~strcmp(R.Mode, 'none'), :);
fprintf('\n%-6s %-6s %-7s %-6s %5s %8s %6s %6s %6s %6s %6s %3s %3s %s\n', 'Animal','Cond','Mode','First','','OffsetH','LenH','InfH','PreH','ShortH','LongH','CM','SS','');
for r = 1:height(R)
    fprintf('%-6s %-6s %-7s %-6d %5s %8.2f %6.2f %6.2f %6.2f %6.2f %6.2f %3d %3d %s%s\n', R.Animal{r}, R.Cond{r}, R.Mode{r}, R.FirstASV(r), '', ...
        R.OffsetH(r), R.LenH(r), R.InfusionH(r), R.PreH(r), R.ShortH(r), R.LongH(r), R.hasCM(r), R.hasSS(r), ternary(R.ShortOK(r), '* ', '  '), R.File{r}(end-min(60,end-1):end));
end
[ua, ~, ia] = unique(R.Animal); nS = accumarray(ia, R.ShortOK);
fprintf('\nASVshort windows per animal (>= %d s): ', MIN_EPOCH); fprintf('%s=%d ', [ua'; num2cell(nS')]); fprintf('\n');
fprintf('total %d windows; animals with >1: %s\n', sum(R.ShortOK), strjoin(ua(nS > 1)', ', '));
fprintf('ts2 modes: %d clock, %d date\n', sum(strcmp(R.Mode, 'clock')), sum(strcmp(R.Mode, 'date')));
save(fullfile(CACHE_ROOT, ['listASVWindows_' datestr(now, 'yyyymmdd_HHMMSS') '.mat']), 'R');

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end