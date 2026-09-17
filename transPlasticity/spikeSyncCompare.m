function R = spikeSyncCompare(A, B, varargin)
%SPIKESYNCCOMPARE  Compare C between two conditions without fooling yourself.
%
%   R = spikeSyncCompare(A, B, 'Dither', 0.02, ...)
%
%   INPUT
%     A, B : either {N x 1} cell of spike-time vectors (one session), or
%            {S x 1} cell whose elements are themselves {N x 1} cells (S
%            sessions). Sessions are PAIRED: A{s} and B{s} must be the same
%            units in the same order. Unit identity must be preserved across
%            conditions -- C dilutes as 1/(N-1), so a difference in sorted
%            unit count produces a C difference mechanically.
%
%   REQUIRED
%     'Dither' : surrogate dither half-width. The comparison is made on
%                dC = C - C_null with the null computed SEPARATELY WITHIN EACH
%                CONDITION. That is the whole point: the within-condition null
%                carries that condition's own rates, unit count and background
%                level, so dC isolates timing. Raw C differences do not.
%
%   OPTIONS
%     'IntervalA','IntervalB' : [t0 t1], or S x 2 per session.
%     'NSurr'  : surrogates per condition per session (default 200).
%     'NPerm'  : sign-flip permutations for the across-session test (default 1e4).
%     'LOO'    : leave-one-unit-out influence (default true when N <= 32).
%
%   OUTPUT
%     .session      S x 1 struct: C, Cnull, dC per condition; N, M, rate,
%                   matchedFrac (fraction of spikes matching above chance --
%                   a proxy for how much of the train is event-locked rather
%                   than background)
%     .dC           S x 2 [A B] surrogate-corrected synchrony
%     .delta        S x 1 dC(B) - dC(A)
%     .stat         mean delta, sign-flip permutation p, and the same for RAW
%                   C so you can see how much of the raw effect was confound
%     .confound.rateDelta      S x 1 median rate change B - A
%     .confound.rhoRate        rank corr across sessions between delta and
%                              rate change. If this is large, your synchrony
%                              effect is a rate effect.
%     .confound.matchedDelta   change in event-locked fraction
%     .confound.rhoMatched     rank corr of delta with it. Large => the effect
%                              is background suppression, not tighter timing.
%     .localize     single-session only: per-unit dCi difference, per-pair
%                   Cmat difference, and leave-one-out influence on delta.
%                   A real population effect is distributed; concentration on
%                   one or two units is an assembly or a sorting result.
%
%   The statistical unit is the SESSION, not the spike and not the pair. The
%   N(N-1)/2 entries of Cmat share trains and are not independent samples.
%   With S = 1 no across-session test is possible and none is reported.
%
%   See also SPIKESYNC, SPIKESYNCDIAGNOSE, SPIKESYNCSURROGATE.

p = inputParser;
p.addParameter('Dither', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('IntervalA', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('IntervalB', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('NSurr', 200, @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.addParameter('NPerm', 1e4, @(x) isnumeric(x) && isscalar(x) && x >= 100);
p.addParameter('LOO', [], @(x) isempty(x) || islogical(x));
p.parse(varargin{:});
o = p.Results;
if isempty(o.Dither)
    error('spikeSyncCompare:dither','''Dither'' is required.');
end

A = local_asSessions(A); B = local_asSessions(B);
S = numel(A);
if numel(B) ~= S
    error('spikeSyncCompare:pairing','A and B must have the same number of sessions.');
end

R.session = struct([]);
dC = zeros(S,2); rawC = zeros(S,2);
rateMed = zeros(S,2); matchedFrac = zeros(S,2);

for s = 1:S
    if numel(A{s}) ~= numel(B{s})
        error('spikeSyncCompare:units', ...
            'Session %d has %d units in A and %d in B. C dilutes as 1/(N-1); use the intersection.', ...
            s, numel(A{s}), numel(B{s}));
    end
    ivA = local_iv(o.IntervalA, s); ivB = local_iv(o.IntervalB, s);
    a = local_one(A{s}, ivA, o.Dither, o.NSurr);
    b = local_one(B{s}, ivB, o.Dither, o.NSurr);

    if abs(a.dur - b.dur) / max(a.dur, b.dur) > 0.2
        warning('spikeSyncCompare:duration', ...
            'Session %d: epoch durations differ by >20%% (%.1f vs %.1f s). C variance scales with spike count.', ...
            s, a.dur, b.dur);
    end

    rawC(s,:)  = [a.C b.C];
    dC(s,:)    = [a.C - a.Cnull, b.C - b.Cnull];
    rateMed(s,:)     = [median(a.rate) median(b.rate)];
    matchedFrac(s,:) = [a.matchedFrac b.matchedFrac];
    R.session(s).A = a; R.session(s).B = b; %#ok<AGROW>
end

R.dC = dC; R.rawC = rawC;
R.delta = dC(:,2) - dC(:,1);
R.confound.rateDelta    = rateMed(:,2) - rateMed(:,1);
R.confound.matchedDelta = matchedFrac(:,2) - matchedFrac(:,1);

% ---- across-session test ------------------------------------------------
if S >= 2
    [R.stat.p, R.stat.mean] = local_signflip(R.delta, o.NPerm);
    [R.stat.pRaw, R.stat.meanRaw] = local_signflip(rawC(:,2) - rawC(:,1), o.NPerm);
    R.confound.rhoRate    = local_spearman2(R.delta, R.confound.rateDelta);
    R.confound.rhoMatched = local_spearman2(R.delta, R.confound.matchedDelta);
else
    R.stat = struct('p', NaN, 'mean', R.delta, 'pRaw', NaN, ...
                    'meanRaw', rawC(2) - rawC(1));
    R.confound.rhoRate = NaN; R.confound.rhoMatched = NaN;
end

% ---- localization (single session) --------------------------------------
R.localize = [];
if S == 1
    N = numel(A{1});
    doLOO = o.LOO; if isempty(doLOO), doLOO = N <= 32; end
    R.localize.dCmat = R.session(1).B.Cmat - R.session(1).A.Cmat;
    R.localize.dUnit = R.session(1).B.unitC - R.session(1).A.unitC;
    if doLOO && N > 2
        argsA = local_syncArgs(local_iv(o.IntervalA,1));
        argsB = local_syncArgs(local_iv(o.IntervalB,1));
        infl = nan(N,1);
        for n = 1:N
            k = setdiff(1:N, n);
            ca = spikeSync(A{1}(k), argsA{:});
            cb = spikeSync(B{1}(k), argsB{:});
            infl(n) = (cb - ca) - (rawC(2) - rawC(1));
        end
        R.localize.looInfluence = infl;
    end
end

% ---- report -------------------------------------------------------------
if S > 1, plural = 's'; else, plural = ''; end
fprintf('\nspikeSyncCompare  (%d session%s)\n', S, plural);
fprintf('  %-8s %8s %8s %8s %8s %8s\n','sess','C_A','C_B','dC_A','dC_B','delta');
for s = 1:S
    fprintf('  %-8d %8.4f %8.4f %8.4f %8.4f %+8.4f\n', ...
        s, rawC(s,1), rawC(s,2), dC(s,1), dC(s,2), R.delta(s));
end
if S >= 2
    fprintf('\n  surrogate-corrected  delta = %+.4f   p = %.4f (sign-flip)\n', ...
        R.stat.mean, R.stat.p);
    fprintf('  RAW                  delta = %+.4f   p = %.4f\n', ...
        R.stat.meanRaw, R.stat.pRaw);
    if abs(R.stat.meanRaw) > 2*abs(R.stat.mean)
        fprintf('  >> most of the raw effect did not survive the within-condition null.\n');
    end
    fprintf('  rho(delta, rate change)          = %+.2f %s\n', R.confound.rhoRate, ...
        local_flag(abs(R.confound.rhoRate) > 0.6, '<-- RATE EFFECT'));
    fprintf('  rho(delta, locked-fraction change) = %+.2f %s\n', R.confound.rhoMatched, ...
        local_flag(abs(R.confound.rhoMatched) > 0.6, '<-- BACKGROUND EFFECT'));
else
    fprintf('\n  S = 1: no across-session test. The pairs in Cmat are not\n');
    fprintf('  independent samples and must not be used as n.\n');
end
fprintf('\n');

end


% =========================================================================
function out = local_one(spk, iv, dither, nsurr)
args = local_syncArgs(iv);
[C, Cmat, prof, Sm] = spikeSync(spk, args{:});
N = numel(spk);
M = cellfun(@numel, spk(:));
dur = Sm.interval(2) - Sm.interval(1);

Cs = zeros(nsurr,1); CiS = [];
for s = 1:nsurr
    sur = cell(N,1);
    for n = 1:N
        x = spk{n}(:);
        sur{n} = sort(x + dither*(2*rand(numel(x),1) - 1));
    end
    [Cs(s), ~, pS] = spikeSync(sur, args{:});
    if s <= min(20, nsurr), CiS = [CiS; pS.C]; end %#ok<AGROW>
end

cut = local_pct2(CiS, 95);
out.C = C; out.Cmat = Cmat; out.prof = prof;
out.Cnull = mean(Cs); out.Csd = std(Cs);
out.N = N; out.M = sum(M); out.dur = dur;
out.rate = M / dur;
out.unitC = cellfun(@mean, Sm.Ci(:));
out.matchedFrac = mean(prof.C > cut);
end


function s = local_asSessions(X)
if ~iscell(X), error('spikeSyncCompare:input','Inputs must be cell arrays.'); end
if iscell(X{1}), s = X(:); else, s = {X}; end
end


function iv = local_iv(I, s)
if isempty(I), iv = []; elseif size(I,1) == 1, iv = I; else, iv = I(s,:); end
end


function a = local_syncArgs(iv)
if isempty(iv), a = {}; else, a = {'Interval', iv}; end
end


function [p, m] = local_signflip(d, nperm)
d = d(:); d = d(~isnan(d));
m = mean(d);
n = numel(d);
null = zeros(nperm,1);
for k = 1:nperm
    sgn = 2*(rand(n,1) > 0.5) - 1;
    null(k) = mean(sgn .* d);
end
p = (1 + sum(abs(null) >= abs(m))) / (1 + nperm);
end


function s = local_flag(tf, msg)
if tf, s = msg; else, s = ''; end
end


function y = local_pct2(x, q)
x = sort(x(~isnan(x(:))));
n = numel(x);
if n == 0, y = nan(size(q)); return; end
if n == 1, y = repmat(x, size(q)); return; end
pos = (0.5:n-0.5) / n * 100;
y = interp1(pos, x, q, 'linear');
y(q < pos(1)) = x(1); y(q > pos(end)) = x(end);
end


function rho = local_spearman2(x, y)
x = x(:); y = y(:);
k = isfinite(x) & isfinite(y);
x = x(k); y = y(k);
if numel(x) < 3, rho = NaN; return; end
rx = tied(x); ry = tied(y);
rx = rx - mean(rx); ry = ry - mean(ry);
rho = (rx'*ry) / max(sqrt((rx'*rx)*(ry'*ry)), eps);
end

function r = tied(x)
[xs, ord] = sort(x(:));
r = zeros(numel(x),1); r(ord) = 1:numel(x);
i = 1;
while i <= numel(xs)
    j = i;
    while j < numel(xs) && xs(j+1) == xs(i), j = j + 1; end
    if j > i, r(ord(i:j)) = mean(i:j); end
    i = j + 1;
end
end
