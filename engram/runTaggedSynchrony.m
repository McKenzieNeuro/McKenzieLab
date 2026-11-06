%% runTaggedSynchrony.m
% Are opto-tagged units more synchronous with each other than chance,
% light off, across OptoRAM sessions?
%
% Per session: taggedSynchrony gives an observed z-score (tagged-group C vs
% a null of size-matched random untagged draws, same session). This script
% collects those z-scores and tests whether they are systematically
% positive across sessions/animals -- the animal-level test, not a pooled
% p-value, since sessions share animals. CHECK whether tagged units are
% re-identified across days in this dataset: if the same tagged cells recur
% across a given animal's sessions, that animal's per-session z-scores are
% not just nested but could be near-duplicates of each other, which the
% animal-level averaging below only partly addresses.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
% Options passed straight to taggedSynchrony (empty = its defaults). Every
% result stores the parsed options in r.params, so the saved file records
% what was run. SAVE_TAG distinguishes runs, e.g. a Pad sweep:
%   TS_OPTS = {'Pad',0.5};  SAVE_TAG = '_pad500ms';
if ~exist('TS_OPTS','var'),  TS_OPTS  = {}; end
if ~exist('SAVE_TAG','var'), SAVE_TAG = ''; end


fils2 = [ ...
       {'R:\WSun\ePhys\OptoRAM\OptoRAM13\OptoRAM13_250105_122857_positive\amplifier_analogin_auxiliary_int16.spikes.cellinfo.mat'}; ...
{'R:\WSun\ePhys\OptoRAM\OptoRAM16\OptoRAM16_250508_095103_Pos\amplifier_analogin_auxiliary_int16.spikes.cellinfo.mat'}; ...
];

fprintf('%d candidate OptoRAM sessions\n', numel(fils2));

R = repmat(taggedSynchronyBlank(), numel(fils2), 1);

for i = 1:numel(fils2)
    try
        r = taggedSynchrony(fils2{i}, TS_OPTS{:});
    catch ME
        r = taggedSynchronyBlank();
        r.why = ME.message;
    end
    r.file = fils2{i};
    R(i) = r;
    if R(i).ok
        fprintf('[%3d/%3d] nTagged=%2d nUntaggedPool=%3d  C=%.3f  z=%+.2f  p=%.4f\n', ...
            i, numel(fils2), R(i).nTagged, R(i).nUntaggedPool, R(i).Cobs, R(i).z, R(i).p);
    else
        fprintf('[%3d/%3d] skipped (%s)\n', i, numel(fils2), R(i).why);
    end
end

ok = [R.ok];
fprintf('\n%d/%d sessions usable\n', sum(ok), numel(R));
U = R(ok);

anim = arrayfun(@(r) animalFromPath(r.file,'OptoRAM'), U, 'UniformOutput', false);
ua = unique(anim);
fprintf('%d animals\n', numel(ua));
fprintf('%-10s %8s %8s\n','animal','nSess','meanZ');
zAnim = nan(numel(ua),1);
for i = 1:numel(ua)
    m = strcmp(anim, ua{i});
    zAnim(i) = mean([U(m).z]);
    fprintf('%-10s %8d %8.2f\n', ua{i}, sum(m), zAnim(i));
end

p_sess = local_signflip([U.z]', 1e4);
fprintf('\nsession-level (one-sample sign-flip vs 0): mean z = %.2f, p = %.4f\n', ...
    mean([U.z]), p_sess);

p_anim = local_signflip(zAnim, 1e4);
fprintf('animal-level  (one-sample sign-flip vs 0): mean z = %.2f, %d/%d animals positive, p = %.4f\n', ...
    mean(zAnim), sum(zAnim>0), numel(zAnim), p_anim);
if numel(zAnim) < 6
    fprintf('  (n=%d animals -> exact floor %.3f; cannot reach 0.05)\n', ...
        numel(zAnim), 2/2^numel(zAnim));
end

fprintf(['\nReport the ANIMAL-level result. The session-level one treats\n' ...
         'sessions from the same animal as independent, which they are not\n' ...
         '(shared implant, and possibly the same tagged units recurring\n' ...
         'across days -- check whether units are tracked/re-identified\n' ...
         'across sessions in this dataset before trusting session count as n).\n']);

outFile = fullfile(CACHE_ROOT, ['taggedSynchronyResults' SAVE_TAG '.mat']);
save(outFile, 'R','U','anim','zAnim','TS_OPTS');
fprintf('saved %s\n', outFile);

%% =========================================================================
function p = local_signflip(d, nperm)
d = d(:); d = d(isfinite(d));
n = numel(d); obs = mean(d);
null = zeros(nperm,1);
for k = 1:nperm, null(k) = mean((2*(rand(n,1)>0.5)-1).*d); end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end
