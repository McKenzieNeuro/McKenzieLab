%% compareRippleLevels.m
% Overall ripple properties in never-infused control animals (R:\HuszarR)
% against the TransPlasticity animals, with no alignment to any infusion in
% the controls (there is none).
%
% GROUPS (one value per session x group, then averaged within animal)
%   Control   whole recording                         rippleEvents_control.mat
%   Naive     TransPlasticity, drug-naive: whole recordings before the animal's
%             first ASV (pre-infusion part only on infusion days), plus the
%             pre-infusion part of the first-ASV session; animals never given
%             ASV contribute all their drug-free time
%   ASVacute  LAT..LAT+POSTLEN after the FIRST ASV infusion
%   Post      >= 1 day after the first ASV, drug-free time only (pre-infusion
%             part on any later infusion day)
%
% METRICS: rate (ripples/h), logDur, freq, fracActive, unitPPC (units with
% >= MINSPK in-ripple spikes), popPPC.
%
% TESTS (animals are the unit)
%   Control vs each TransPlasticity group: two-sample permutation on animal means
%   Naive -> ASVacute and Naive -> Post: paired sign-flip within TransPlasticity animals
%
% THESE ARE BETWEEN-DATASET LEVELS. Different lab, spike sorting, probes,
% ripple detector, and phase/frequency methods (Huszar: detector's stored phase,
% unweighted frequency; TransPlasticity: Hilbert phase of the band-passed LFP,
% envelope-weighted frequency). fracActive also depends on how many and which
% units each dataset has. A Control-vs-TransPlasticity difference is therefore
% not attributable to the drug; it says whether the infusion-study animals sit
% in the range of an independent un-manipulated dataset. Only the paired
% within-TransPlasticity rows are like-for-like.

CACHE_ROOT   = 'C:\Users\samckenzie\syncCache';
ROOT         = 'R:\TransPlasticity';
CONTROL_ROOT = 'R:\HuszarR';
if exist('SHARED', 'var')     % set by runAllAnalyses.m
    CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root; CONTROL_ROOT = SHARED.controlRoot;
end
LAT = 4*3600; POSTLEN = 2*3600; MIN_EPOCH = 1800;
MINSPK = 20; NPERM = 1e4;
CONTROL_UNIT = 'session';   % 'session': each control recording is one replicate (N = sessions)
                            % 'animal' : control sessions averaged within animal first

fa = getAllExtFiles(ROOT, 'mat', 1);         fa = fa(endsWith(fa, 'rippleEvents.mat'));
fc = getAllExtFiles(CONTROL_ROOT, 'mat', 1); fc = fc(endsWith(fc, 'rippleEvents_control.mat'));
fprintf('%d infusion-study exports, %d control exports\n', numel(fa), numel(fc));

rows = {};
% ---- controls: whole recording ---------------------------------------------------
for i = 1:numel(fc)
    L = load(fc{i}, 'RIPPLES'); RP = L.RIPPLES;
    unitID = RP.meta.animal;
    if strcmp(CONTROL_UNIT, 'session'), [~, unitID] = fileparts(RP.meta.folder); end   % one replicate per recording
    rows(end+1,:) = [{'Control', unitID, fc{i}}, local_epoch(RP, [0 RP.meta.recordingSec], MINSPK)]; %#ok<AGROW>
end
% ---- TransPlasticity: windows by exposure state -------------------------------------
for i = 1:numel(fa)
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    dur = m.recordingSec; infS = m.infusionSec; hasInf = isfinite(infS);
    drugFree = [0, dur]; if hasInf, drugFree = [0, infS]; end
    noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
    isFirst = ~noASV && strcmp(m.folder, m.firstASV.session);
    if noASV || (isfinite(m.dayFromFirstASV) && m.dayFromFirstASV < 0)
        W = {'Naive', drugFree};
    elseif isFirst
        W = {'Naive', [0, infS]; 'ASVacute', [infS + LAT, min(dur, infS + LAT + POSTLEN)]};
    elseif m.dayFromFirstASV >= 1
        W = {'Post', drugFree};
    else
        continue                                  % same day as the first ASV, not that session
    end
    for k = 1:size(W, 1)
        w = W{k,2};
        if diff(w) < MIN_EPOCH, continue; end
        rows(end+1,:) = [{W{k,1}, m.animal, fa{i}}, local_epoch(RP, w, MINSPK)]; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Group','Animal','Session','Event','Hours','logDur','freq', ...
    'fracActive','ampZ','unitPPC','popPPC','popPhase'});
G = {'Control','Naive','ASVacute','Post'};
for g = G
    m = strcmp(T.Group, g{1});
    fprintf('%-9s %3d sessions, %2d animals, %.0f h\n', g{1}, sum(m), numel(unique(T.Animal(m))), sum(T.Hours(m)));
end

%% ---- per-animal means and tests ----------------------------------------------------
metrics = {'rate','logDur','freq','fracActive','unitPPC','popPPC'};
fprintf('\nanimal means (mean +- SEM across animals)\n');
fprintf('%-10s | %-18s | %-18s | %-18s | %-18s || %-26s | %-22s\n', 'metric', G{:}, ...
    'Control vs Naive/ASVac/Post p', 'Naive->ASVac, ->Post p');
OUT = struct();
for mm = metrics
    nm = mm{1};
    M = struct();
    for g = G, [M.(g{1}).v, M.(g{1}).ids] = local_animalMean(T(strcmp(T.Group, g{1}), :), nm); end
    pc = nan(1, 3);
    for j = 2:4
        if numel(M.Control.v) >= 2 && numel(M.(G{j}).v) >= 2, pc(j-1) = permtest(M.Control.v, M.(G{j}).v, NPERM); end
    end
    pp = nan(1, 2);
    for j = 3:4
        [~, ia, ib] = intersect(M.Naive.ids, M.(G{j}).ids);
        if numel(ia) >= 2, pp(j-2) = local_signflip(M.(G{j}).v(ib) - M.Naive.v(ia)); end
    end
    OUT.(nm) = struct('M', M, 'pControl', pc, 'pPaired', pp);
    c = cellfun(@(g) local_ms(M.(g).v), G, 'UniformOutput', false);
    fprintf('%-10s | %-18s | %-18s | %-18s | %-18s || %.3f / %.3f / %.3f | %.3f / %.3f\n', nm, c{:}, pc, pp);
end
fprintf('control replicate = %s (TransPlasticity replicate = animal)\n', CONTROL_UNIT);
fprintf(['rate in ripples/h over the whole window (state-dependent: no sleep scoring).\n' ...
         'Control vs TransPlasticity differences mix dataset/method with biology (see header).\n']);
save(fullfile(CACHE_ROOT, 'compareRippleLevels.mat'), 'T', 'OUT', 'LAT', 'POSTLEN');

%% ---- figure: per-animal values by group ------------------------------------------------
fig = figure('Color','w','Position',[80 80 1400 320]);
tiledlayout(1, numel(metrics), 'TileSpacing','compact');
col = {[0.4 0.4 0.4], [0 0 0], [0.85 0.2 0.1], [0.55 0.1 0.6]};
for k = 1:numel(metrics)
    nexttile; hold on
    for j = 1:4
        v = OUT.(metrics{k}).M.(G{j}).v;
        x = j + 0.12 * (rand(size(v)) - 0.5);
        plot(x, v, 'o', 'Color', col{j}, 'MarkerFaceColor', col{j}, 'MarkerSize', 4);
        if ~isempty(v), plot(j + [-0.25 0.25], mean(v) * [1 1], '-', 'Color', col{j}, 'LineWidth', 2); end
    end
    % paired lines within TransPlasticity animals
    Mn = OUT.(metrics{k}).M;
    for j = 3:4
        [~, ia, ib] = intersect(Mn.Naive.ids, Mn.(G{j}).ids);
        for q = 1:numel(ia), plot([2 j], [Mn.Naive.v(ia(q)) Mn.(G{j}).v(ib(q))], '-', 'Color', [0.8 0.8 0.8]); end
    end
    xlim([0.5 4.5]); xticks(1:4); xticklabels(G); xtickangle(30);
    title(metrics{k}, 'Interpreter','none');
end
exportgraphics(fig, fullfile(CACHE_ROOT, 'compareRippleLevels.pdf'), 'ContentType','vector');


%% =========================================================================
function y = local_epoch(RP, w, minSpk)
% {nRip, hours, logDur, freq, fracActive, ampZ, unitPPC, popPPC, popPhase} for one window
in = RP.ts1 >= w(1) & RP.ts1Off <= w(2);
ev = RP.ev(in, :);
rk = find(in);
sp = RP.spikes(ismember(RP.spikes.Ripple, rk) & isfinite(RP.spikes.Phase), :);
nU = height(RP.units); th = sp.Phase; u = sp.Unit;
n = accumarray(u, 1, [nU 1]);
Cs = accumarray(u, cos(th), [nU 1]); Sn = accumarray(u, sin(th), [nU 1]);
ppc = (Cs.^2 + Sn.^2 - n) ./ (n .* (n - 1));
ok = n >= minSpk;
nP = sum(n); CP = sum(Cs); SP = sum(Sn);
popPPC = NaN; popPh = NaN;
if nP >= 2, popPPC = (CP^2 + SP^2 - nP) / (nP * (nP - 1)); popPh = atan2(SP, CP); end
y = {height(ev), diff(w) / 3600, mean(log(ev.Dur), 'omitnan'), mean(ev.Freq, 'omitnan'), ...
     mean(ev.FracActive, 'omitnan'), mean(ev.AmpZ, 'omitnan'), mean(ppc(ok)), popPPC, popPh};
end


function [v, ids] = local_animalMean(T, nm)
% per animal: rate pooled over its windows (events / hours); others mean over windows
ua = unique(T.Animal); v = nan(numel(ua), 1);
for i = 1:numel(ua)
    m = strcmp(T.Animal, ua{i});
    if strcmp(nm, 'rate'), v(i) = sum(T.Event(m)) / sum(T.Hours(m));
    else, v(i) = mean(T.(nm)(m), 'omitnan'); end
end
ok = isfinite(v); v = v(ok); ids = ua(ok);
end


function s = local_ms(v)
if isempty(v), s = 'n/a'; return; end
s = sprintf('%.3g +- %.2g (n=%d)', mean(v), std(v) / sqrt(numel(v)), numel(v));
end


function p = local_signflip(d)
d = d(isfinite(d)); n = numel(d);
if n < 2, p = NaN; return; end
obs = abs(mean(d));
if n <= 14
    S = 2*(dec2bin(0:2^n-1, n) == '1') - 1;
    p = mean(abs(S*d(:)/n) >= obs - 1e-12);
else
    nul = arrayfun(@(k) mean((2*(rand(n,1) > 0.5) - 1) .* d), 1:1e4);
    p = (1 + sum(abs(nul) >= obs)) / (1 + 1e4);
end
end
