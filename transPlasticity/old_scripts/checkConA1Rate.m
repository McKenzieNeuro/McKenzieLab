%% checkConA1Rate.m
% Pre-infusion firing rate during the ConA1 exposure, paired within animal, SAL vs ASV.
%
% WHY. The sampling table in driftInclusionASV.m implies that pyramidal cells fire faster
% during ConA1 in saline sessions than in ASV sessions - before anything is infused. The draft
% reports pre-infusion home-cage rates of 1.98 (SAL) and 1.80 Hz (ASV), which are matched, so
% either the difference is confined to the exposure or the group means are carried by the three
% ASV animals that have no saline partner. Only the within-animal pairing distinguishes these.
%
% Run it after driftInclusionASV.m, which leaves T in the workspace; otherwise it loads the
% saved table. Nothing is recomputed.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
if ~exist('T', 'var')
    L = load(fullfile(CACHE_ROOT, 'driftInclusionASV.mat'), 'T'); T = L.T; %#ok<NASGU>
end

for ctx = {'ConA1', 'ConB1'}
    k = strcmp(T.Context, ctx{1}) & ismember(T.Group, {'SALpre', 'ASV1'});
    S = T(k,:);
    [~, ia] = unique(S.Session); S = S(ia,:);          % one row per session
    % nSpkCtx is summed over every recorded pyramidal cell and nUnits is that count, so this is
    % the mean pyramidal rate during the exposure, before any cell subsampling.
    S.rate = S.nSpkCtx ./ (S.nUnits .* S.ctxMin * 60);

    fprintf('\n===== %s: mean pyramidal rate during the exposure =====\n', ctx{1});
    for g = {'SALpre', 'ASV1'}
        kk = strcmp(S.Group, g{1});
        if ~any(kk), continue; end
        fprintf('  %-7s %5.2f Hz   (n = %d sessions, %.2f to %.2f, units %.0f to %.0f)\n', ...
            g{1}, mean(S.rate(kk)), sum(kk), min(S.rate(kk)), max(S.rate(kk)), ...
            min(S.nUnits(kk)), max(S.nUnits(kk)));
    end

    an = unique(S.Animal); d = []; lab = {};
    fprintf('  %-6s %8s %8s %8s %12s\n', 'animal', 'SAL', 'ASV', 'ASV/SAL', 'units S/A');
    for q = 1:numel(an)
        ka = strcmp(S.Animal, an{q}) & strcmp(S.Group, 'ASV1');
        ks = strcmp(S.Animal, an{q}) & strcmp(S.Group, 'SALpre');
        if ~any(ka) || ~any(ks), continue; end
        rs = mean(S.rate(ks)); ra = mean(S.rate(ka));
        d(end+1) = log2(ra / rs); lab{end+1} = an{q}; %#ok<AGROW,SAGROW>
        fprintf('  %-6s %8.2f %8.2f %8.2f %6.0f /%-5.0f\n', an{q}, rs, ra, ra / rs, ...
            mean(S.nUnits(ks)), mean(S.nUnits(ka)));
    end
    if numel(d) >= 2
        fprintf('  paired log2(ASV/SAL) = %+.2f  (x%.2f), p = %.4f, n = %d animals\n', ...
            mean(d), 2^mean(d), local_signflip(d), numel(d));
        fprintf('  animals contributing only one arm are excluded here and are what the group\n');
        fprintf('  means above include, so a gap between the two is the unpaired animals.\n');
    else
        fprintf('  too few animals with both arms to pair\n');
    end
end

fprintf('\nIf the paired ratio is near 1, the group difference is the unpaired ASV animals and\n');
fprintf('nothing further is needed. If ASV sessions really are slower within animal, the next\n');
fprintf('test is behavioural rather than electrophysiological: running speed or path length\n');
fprintf('during ConA1 on each day. ConA1 is met again on every session, so by the ASV day it is\n');
fprintf('more familiar than it was on the saline day, and less exploration would lower firing\n');
fprintf('during the exposure while leaving home-cage rates matched, which is what the draft\n');
fprintf('reports. Matched speed with a persisting rate difference would instead be a baseline\n');
fprintf('imbalance in excitability, which is a harder problem.\n');


function p = local_signflip(x)
% exact two-sided sign-flip test on the mean, floored at 2/2^n
x = x(isfinite(x)); n = numel(x);
if n < 2, p = NaN; return; end
if n > 20
    s = sign(randn(n, 2000)); s(:,1) = 1;
else
    s = double(dec2bin(0:2^n - 1) - '0')' * 2 - 1;
end
m = abs(x(:)' * s / n);
p = max(mean(m >= abs(mean(x)) - 1e-12), 2 / 2^n);
end
