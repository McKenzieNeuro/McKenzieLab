function [E, y1, y2, D] = spikeDistanceProfile(trains, varargin)
%SPIKEDISTANCEPROFILE  Multivariate SPIKE-distance profile S(t) (Kreuz et al.).
%
%   [E, y1, y2, D] = spikeDistanceProfile(trains, 'Interval', [t0 t1], ...)
%
%   The time-resolved, time-scale-adaptive SPIKE-distance: a DISSIMILARITY
%   defined at every instant (0 = identical spike timing), piecewise linear
%   between the events of the pooled spike train. Bivariate profile, for
%   trains n, m at time t (Kreuz et al. 2013, J Neurophysiol, the rate-weighted
%   refinement of Kreuz et al. 2011, J Neurosci Methods):
%       S_n(t) = [dt_P^n x_F^n(t) + dt_F^n x_P^n(t)] / x_ISI^n(t)
%       S(t)   = [S_n x_ISI^m + S_m x_ISI^n] / ( 0.5 (x_ISI^n + x_ISI^m)^2 )
%   dt_P^n / dt_F^n: distance of n's previous / following spike to the nearest
%   spike of m; x_P, x_F: time since / until them; x_ISI: n's current ISI.
%   Multivariate profile = mean over all N(N-1)/2 pairs. D = time average.
%
%   OPTIONS
%     'Interval' [t0 t1]  (default [min max] of all spikes)
%     'RI'       false    true: rate-independent form (S_n + S_m)/(x^n + x^m),
%                         the original 2011 normalisation
%     'MRTS'     0        minimum relevant time scale (adaptive SPIKE-distance);
%                         0 = the standard measure
%
%   EDGES follow PySpike 0.9.0 exactly (auxiliary spikes extrapolated by one
%   ISI, the first/last spike's dt carried to the edge -- Kreuz et al. 2015).
%
%   OUTPUT
%     E   (K+1) x 1 event times (interval edges and every spike, unique)
%     y1  K x 1 profile value at the START of interval [E(k), E(k+1))
%     y2  K x 1 profile value at the END of that interval (left limit at E(k+1))
%         S is linear in between and generally discontinuous at spikes.
%     D   time-averaged SPIKE-distance
%
%   VALIDATION. The vectorised algorithm here was checked against
%   pyspike.spike_profile / spike_profile_multi on 2000 random spike-train sets
%   (2-6 trains, ties across trains, spikes on the interval edges, single-spike
%   trains, RI and MRTS variants): max |difference| 4.4e-16. Run
%   spikeDistanceProfile() with no arguments to check THIS MATLAB port against
%   a stored PySpike result.
%
%   FAST BINNED MODE. [m, ~, ~, D] = spikeDistanceProfile(trains, ..., 'Bins', g)
%   returns the exact mean of S(t) over each bin [g(k) g(k+1)] (and D over
%   [g(1) g(end)]) without forming the profile on the pooled grid. Same values
%   as spikeDistanceBinMean(E, y1, y2, g); far less time and memory when the
%   pooled spike count is large (cost per pair ~ its own spikes + numel(g)).
%
%   PAIRS. 'Pairs', P x 2 restricts the multivariate average to those unit
%   pairs (default: all N(N-1)/2). A random subset gives an unbiased estimate
%   of the all-pairs profile; use the SAME subset for data and surrogates.
%
%   See also SPIKEDISTANCEBINMEAN, SPIKEDISTANCESESSION.

if nargin == 0, E = local_selftest(); return; end
p = inputParser;
p.addParameter('Interval', []);
p.addParameter('RI', false);
p.addParameter('MRTS', 0);
p.addParameter('Bins', []);
p.addParameter('Pairs', []);          % P x 2 unit pairs to average over ([] = all)
p.parse(varargin{:});
o = p.Results;

trains = cellfun(@(x) unique(reshape(x(~isnan(x)), [], 1)), trains(:), 'UniformOutput', false);
N = numel(trains);
if N < 2, error('spikeDistanceProfile:input', 'need at least two trains'); end
if any(cellfun(@isempty, trains)), error('spikeDistanceProfile:empty', 'empty spike train'); end
allT = vertcat(trains{:});
if isempty(o.Interval), iv = [min(allT) max(allT)]; else, iv = sort(o.Interval(:))'; end
t0 = iv(1); t1 = iv(2);
if min(allT) < t0 || max(allT) > t1
    error('spikeDistanceProfile:interval', 'spikes outside the Interval');
end

if isempty(o.Pairs), PL = nchoosek(1:N, 2); else, PL = o.Pairs; end

if ~isempty(o.Bins)
    [E, D] = local_binned(trains, PL, t0, t1, o.Bins(:), o.MRTS, o.RI);   % E = bin means here
    y1 = []; y2 = [];
    return
end

E = unique([t0; t1; allT]);
E0 = E(1:end-1); E1 = E(2:end);
K = numel(E0);
y1 = zeros(K,1); y2 = zeros(K,1);
for q = 1:size(PL, 1)
    i = PL(q,1); j = PL(q,2);
    [s1s, s1e, x1] = local_terms(trains{i}, trains{j}, E0, E1, t0, t1);
    [s2s, s2e, x2] = local_terms(trains{j}, trains{i}, E0, E1, t0, t1);
    y1 = y1 + local_dist(s1s, s2s, x1, x2, o.MRTS, o.RI);
    y2 = y2 + local_dist(s1e, s2e, x1, x2, o.MRTS, o.RI);
end
y1 = y1 / size(PL, 1); y2 = y2 / size(PL, 1);
D = sum(0.5*(y1 + y2) .* diff(E)) / (t1 - t0);
end


function [m, D] = local_binned(trains, PL, t0, t1, g, mrts, ri)
% Exact bin means of the multivariate profile WITHOUT building it on the
% pooled grid: each pair is evaluated only on its own events (M_n + M_m + 2,
% not sum(M)), integrated exactly to the bin edges g, and the integrals are
% averaged over pairs. Same numbers as spikeDistanceBinMean on the exact
% profile; cost per pair ~ (M_n + M_m) + numel(g) instead of sum(M), and
% memory never exceeds one pair's events.
Fs = zeros(numel(g), 1);
for q = 1:size(PL, 1)
    i = PL(q,1); j = PL(q,2);
    E = unique([t0; t1; trains{i}; trains{j}]);
    E0 = E(1:end-1); E1 = E(2:end);
    [s1s, s1e, x1] = local_terms(trains{i}, trains{j}, E0, E1, t0, t1);
    [s2s, s2e, x2] = local_terms(trains{j}, trains{i}, E0, E1, t0, t1);
    y1 = local_dist(s1s, s2s, x1, x2, mrts, ri);
    y2 = local_dist(s1e, s2e, x1, x2, mrts, ri);
    F = [0; cumsum(0.5 * (y1 + y2) .* diff(E))];
    k = discretize(g, E);
    dx = g - E(k); h = E(k+1) - E(k);
    Fs = Fs + F(k) + y1(k) .* dx + (y2(k) - y1(k)) .* dx.^2 ./ (2 * h);
end
Fs = Fs / size(PL, 1);
m = diff(Fs) ./ diff(g);
D = (Fs(end) - Fs(1)) / (g(end) - g(1));
end


function [Ss, Se, isi] = local_terms(a, b, E0, E1, t0, t1)
% S_n at the start (t = E0) and end (t -> E1) of every global interval, and
% train a's ISI there, for train a against train b.
na = numel(a); K = numel(E0);
[bL, bR] = local_aux(b, t0, t1);
Dn = local_mindist(a, b, bL, bR);
p = discretize(E0, [-Inf; a; Inf]) - 1;          % number of spikes of a <= E0
Ss = zeros(K,1); Se = zeros(K,1); isi = zeros(K,1);

in = p >= 1 & p < na;
q = p(in);
tP = a(q); tF = a(q+1); dP = Dn(q); dF = Dn(q+1); x = tF - tP;
Ss(in) = (dP .* (tF - E0(in)) + dF .* (E0(in) - tP)) ./ x;
Se(in) = (dP .* (tF - E1(in)) + dF .* (E1(in) - tP)) ./ x;
isi(in) = x;

bef = p == 0;                                     % before the first spike
Ss(bef) = Dn(1); Se(bef) = Dn(1);
if na > 1, isi(bef) = max(a(1) - t0, a(2) - a(1)); else, isi(bef) = a(1) - t0; end

aft = p == na;                                    % after the last spike
if na == 1 && a(1) == t0
    % PySpike: a single spike AT t_start is followed by t_end as a real
    % following spike, whose distance to the other train is used
    dF = local_mindist(t1, b, bL, bR); x = t1 - a(1);
    Ss(aft) = (Dn(1) * (t1 - E0(aft)) + dF * (E0(aft) - a(1))) / x;
    Se(aft) = (Dn(1) * (t1 - E1(aft)) + dF * (E1(aft) - a(1))) / x;
    isi(aft) = x;
else
    Ss(aft) = Dn(end); Se(aft) = Dn(end);
    if na > 1, isi(aft) = max(t1 - a(end), a(end) - a(end-1)); else, isi(aft) = t1 - a(end); end
end
end


function [L, R] = local_aux(a, t0, t1)
% auxiliary edge spikes: one ISI beyond the first/last spike, at least the edge
if numel(a) > 1
    L = min(t0, a(1) - (a(2) - a(1)));
    R = max(t1, a(end) + (a(end) - a(end-1)));
else
    L = t0; R = t1;
end
end


function d = local_mindist(q, b, bL, bR)
% distance of each q to the nearest of {bL, b, bR}
B = [unique([bL; b(:); bR]); Inf];
k = discretize(q(:), B);
d = min(q(:) - B(k), B(k+1) - q(:));
end


function s = local_dist(s1, s2, i1, i2, mrts, ri)
m = 0.5 * (i1 + i2);
lim = max(mrts, m);
if ri
    s = 0.5 * (s1 + s2) ./ lim;
else
    s = 0.5 * (s1 .* i2 + s2 .* i1) ./ (m .* lim);
end
end


function ok = local_selftest()
% Reference: pyspike 0.9.0 spike_profile_multi on these four trains, [0 5].
tr = {[0.0 0.8 1.9 3.1 4.4], [0.5 1.9 2.7 4.95], [1.2 2.2 3.0], [0.3 2.6 3.3 3.9 4.6 5.0]};
Eref = [0 0.3 0.5 0.8 1.2 1.9 2.2 2.6 2.7 3.0 3.1 3.3 3.9 4.4 4.6 4.95 5.0]';
y1ref = [0.28571276040838733 0.29946470612608844 0.30910278006338765 0.30373737095183856 ...
    0.3180666697200601 0.28101108970494004 0.2656906940942237 0.3407758797578742 ...
    0.2641658783689904 0.1997473740596178 0.1945750019202544 0.20727126296911047 ...
    0.34354602050967104 0.25858401640080464 0.2382218315202117 0.22137679094005647]';
y2ref = [0.29946470612608844 0.30910278006338765 0.32348226142657754 0.2983343921356135 ...
    0.2555355573572034 0.251434146525055 0.24291775503125282 0.33464355332776385 ...
    0.2463453791387123 0.19792466011251567 0.20058972814867204 0.3608390262306138 ...
    0.25858401640080464 0.21851676525629515 0.22137679094005647 0.21908188641661164]';
[E, y1, y2, D] = spikeDistanceProfile(tr, 'Interval', [0 5]);
[~, ~, ~, Dri] = spikeDistanceProfile(tr, 'Interval', [0 5], 'RI', true);
g = [0 0.37 1.1 1.9 2.45 3.3 4.0 4.71 5.0]';
mb = spikeDistanceProfile(tr, 'Interval', [0 5], 'Bins', g);     % fast pair-local path
err = max([max(abs(E - Eref)), max(abs(y1 - y1ref)), max(abs(y2 - y2ref)), ...
    abs(D - 0.2752298734380775), abs(Dri - 0.25465030910191955), ...
    max(abs(mb - spikeDistanceBinMean(E, y1, y2, g)))]);
ok = err < 1e-12;
if ok, s = 'PASS'; else, s = 'FAIL'; end
fprintf('spikeDistanceProfile self-test vs PySpike 0.9.0: max |diff| = %.2g -> %s\n', err, s);
end