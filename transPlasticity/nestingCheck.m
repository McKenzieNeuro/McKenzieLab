function R = nestingCheck(a, b, x1, x2, rootA, rootB, label)
%NESTINGCHECK  Is the session-level test overstating n because sessions cluster within animals?
%
%   R = nestingCheck(a, b, x1, x2, rootA, rootB, label)
%
%   a, b     : the .ok session structs from syncBatch (must have .file)
%   x1, x2   : the per-session metric being compared (e.g. occ1, occ2)
%   rootA/B  : dataset root folder names for animalFromPath ('TransPlasticity','OptoRAM')
%   label    : name of the metric, for printing
%
%   WHY THIS MATTERS. Every p-value computed so far permuted SESSION labels,
%   which assumes sessions are independent draws. If one animal contributes
%   several sessions, they are not: sessions from the same animal share an
%   implant, a unit population, a surgery, a behavioural history. The test
%   then has fewer effective degrees of freedom than it thinks, and reported
%   p-values are optimistic -- sometimes by a lot.
%
%   Three things get reported:
%
%   1. STRUCTURE: how many animals, sessions per animal. If every animal has
%      one session there is no nesting and nothing to correct.
%
%   2. ICC and design effect. ICC is the fraction of total variance in the
%      metric that sits BETWEEN animals rather than within. Design effect
%      DE = 1 + (mbar-1)*ICC, and n_eff = n/DE is roughly how many
%      independent sessions you actually have. ICC near 0 means sessions
%      within an animal are as different as sessions across animals, and the
%      session-level p-value is approximately fine. ICC near 1 means each
%      animal effectively contributes one data point.
%
%   3. ANIMAL-LEVEL TEST: collapse to one value per animal (mean over that
%      animal's sessions), then permute animal labels. This is the
%      conservative, assumption-light version -- it cannot be accused of
%      pseudo-replication because the unit of randomisation is the unit of
%      biological independence. If the effect survives here, it is not a
%      nesting artifact. If it does not, the session-level result was
%      leaning on repeated sampling of a few animals.
%
%   The animal-level test is the number to report. The session-level test is
%   supporting detail, not the headline, whenever ICC is non-trivial.
%
%   See also ANIMALFROMPATH, PERMTEST.

if nargin < 7, label = 'metric'; end

idA = arrayfun(@(r) animalFromPath(r.file, rootA), a, 'UniformOutput', false);
idB = arrayfun(@(r) animalFromPath(r.file, rootB), b, 'UniformOutput', false);
idA = idA(:); idB = idB(:);

bad = sum(cellfun(@isempty, idA)) + sum(cellfun(@isempty, idB));
if bad > 0
    warning('nestingCheck:unparsed', ...
        '%d session paths did not yield an animal ID -- check rootA/rootB.', bad);
end

[uA, ~, gA] = unique(idA);
[uB, ~, gB] = unique(idB);

fprintf('\n=== nesting: %s ===\n', label);
fprintf('TransPlasticity: %d sessions from %d animals\n', numel(idA), numel(uA));
local_printCounts(uA, gA);
fprintf('OptoRAM        : %d sessions from %d animals\n', numel(idB), numel(uB));
local_printCounts(uB, gB);

% ---- ICC / design effect, computed within each dataset -------------------
[iccA, deA, neffA] = local_icc(x1(:), gA);
[iccB, deB, neffB] = local_icc(x2(:), gB);
fprintf('\nICC (between-animal variance fraction)\n');
fprintf('  TransPlasticity: ICC=%.3f  design effect=%.2f  n_eff=%.1f (of %d sessions)\n', ...
    iccA, deA, neffA, numel(x1));
fprintf('  OptoRAM        : ICC=%.3f  design effect=%.2f  n_eff=%.1f (of %d sessions)\n', ...
    iccB, deB, neffB, numel(x2));

% ---- per-animal collapse + animal-level permutation ----------------------
mA = accumarray(gA, x1(:), [], @(v) mean(v, 'omitnan'));
mB = accumarray(gB, x2(:), [], @(v) mean(v, 'omitnan'));

pSession = permtest(x1, x2);
pAnimal  = permtest(mA, mB);

fprintf('\n%s: session-level  %.3f vs %.3f   p = %.4f  (n = %d, %d sessions)\n', ...
    label, mean(x1,'omitnan'), mean(x2,'omitnan'), pSession, numel(x1), numel(x2));
fprintf('%s: ANIMAL-level    %.3f vs %.3f   p = %.4f  (n = %d, %d animals)\n', ...
    label, mean(mA,'omitnan'), mean(mB,'omitnan'), pAnimal, numel(mA), numel(mB));
if any(isnan(mA)) || any(isnan(mB))
    fprintf('  (%d/%d TransPlasticity and %d/%d OptoRAM animals have NO significant\n', ...
        sum(isnan(mA)), numel(mA), sum(isnan(mB)), numel(mB));
    fprintf('   transients in ANY session -- excluded from the animal-level mean and test.)\n');
end

if pAnimal > 0.05 && pSession < 0.05
    fprintf(['>> the effect does NOT survive collapsing to animals. The session-level\n' ...
             '   p-value was relying on repeated sessions from the same animals.\n' ...
             '   Report the animal-level result.\n']);
elseif pAnimal < 0.05
    fprintf(['>> effect survives at the animal level -- not a pseudo-replication\n' ...
             '   artifact. Report this p-value as primary.\n']);
end

R = struct('idA',{idA},'idB',{idB},'animalsA',{uA},'animalsB',{uB}, ...
    'perAnimalA',mA,'perAnimalB',mB, 'iccA',iccA,'iccB',iccB, ...
    'neffA',neffA,'neffB',neffB,'pSession',pSession,'pAnimal',pAnimal);
end


function local_printCounts(u, g)
c = accumarray(g, 1);
[c, ord] = sort(c, 'descend');
for i = 1:numel(c)
    fprintf('    %-14s %d session(s)\n', u{ord(i)}, c(i));
end
end


function [icc, de, neff] = local_icc(x, g)
% One-way random-effects ICC(1,1) from the standard MS decomposition.
% No Statistics Toolbox needed.
ok = ~isnan(x);
x = x(ok); g = g(ok);
[ug, ~, gi] = unique(g);
k = numel(ug); n = numel(x);
if k < 2 || n <= k
    icc = NaN; de = NaN; neff = n; return
end
ni = accumarray(gi, 1);
grand = mean(x);
gm = accumarray(gi, x, [], @mean);

SSB = sum(ni .* (gm - grand).^2);
SSW = sum((x - gm(gi)).^2);
MSB = SSB / (k - 1);
MSW = SSW / (n - k);

% m0: average cluster size corrected for unequal sizes
m0 = (n - sum(ni.^2)/n) / (k - 1);
icc = (MSB - MSW) / (MSB + (m0 - 1)*MSW);
icc = max(min(icc, 1), 0);          % negative ICC estimates -> 0

mbar = mean(ni);
de = 1 + (mbar - 1) * icc;
neff = n / de;
end