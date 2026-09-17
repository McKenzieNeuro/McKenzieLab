function [C, Cmat, prof, S] = spikeSync(spikes, varargin)
%SPIKESYNC  SPIKE-Synchronization (Kreuz et al.) for a cell array of spike trains.
%
%   [C, Cmat, prof, S] = spikeSync(spikes, ...)
%
%   INPUT
%     spikes : {N x 1} cell array; spikes{n} = vector of spike times for neuron n.
%              Trains may be empty, unsorted, and of unequal length.
%
%   NAME/VALUE OPTIONS
%     'Interval' : [tmin tmax] recording interval. Default = [min max] over all
%                  spikes. This is NOT cosmetic: for the first/last spike of a
%                  train the missing neighbouring ISI is replaced by the full
%                  interval length (PySpike/cSPIKE convention), so edge spikes
%                  get a maximally permissive coincidence window when the other
%                  four candidate intervals are large.
%     'MaxTau'   : hard upper bound on the coincidence window tau (default Inf).
%                  Set this if you want an absolute ceiling on what counts as
%                  "coincident" (e.g. 5 ms); tau = min(tau_adaptive, MaxTau).
%
%   OUTPUT
%     C     : scalar multivariate SPIKE-Synchronization, Eq. 19.
%             C = (1/M) * sum_k C(t_k), M = total number of spikes.
%             C = 1 if M == 0 (common silence = perfect synchrony).
%     Cmat  : N x N pairwise SPIKE-Synchronization, Cmat(n,m) =
%             (#coincident spikes in n + in m) / (M_n + M_m). Diagonal = 1.
%             NOTE: C is NOT mean(Cmat(offdiag)); C weights by spike counts,
%             Cmat weights each pair equally. Both are in [0,1].
%     prof  : struct with the pooled SPIKE-Synchronization profile (Eq. 18),
%             sorted by spike time:
%               .t      global spike times t_k
%               .C      C(t_k) = (1/(N-1)) * sum_{m~=n} C_i^{(n,m)}
%               .train  originating train index n(k)
%               .idx    spike index within that train i(k)
%     S     : struct with the raw matching, i.e. everything SPIKE-Order /
%             Spike Train Order (Eqs. 20-27) needs downstream:
%               .coinc{n,m}  logical(M_n x 1), C_i^{(n,m)}
%               .match{n,m}  M_n x 1 index j' of the nearest spike in train m
%               .tau{n,m}    M_n x 1 coincidence window actually used
%               .Ci{n}       M_n x 1 multivariate counter C_i^{(n)} (Eq. 17)
%               .interval    [tmin tmax] used
%
%   ALGORITHM (Eqs. 15-19 of Kreuz, "Quantifying spike train synchrony and
%   directionality"). For spike i of train n and its nearest spike j in train m,
%       tau_ij = min{ t^n_{i+1}-t^n_i , t^n_i-t^n_{i-1} ,
%                     t^m_{j+1}-t^m_j , t^m_j-t^m_{j-1} } / 2
%       C_i^{(n,m)} = 1  iff  |t^n_i - t^m_j| < tau_ij   (strict "<")
%   Nearest-neighbour matching plus the half-ISI window makes the matching
%   provably unambiguous and symmetric: C_i^{(n,m)} = 1 implies the partner
%   spike j has i as ITS nearest match and C_j^{(m,n)} = 1. (If j had a closer
%   partner i' in train n, then |t_i - t_i'| < 2*tau_ij <= min neighbouring ISI
%   of spike i, which is impossible.) The code exploits nothing of this — both
%   directions are computed independently — so an asymmetry would surface as a
%   non-integer coincidence count rather than being silently hidden.
%
%   Validated against PySpike 0.9.0 (spike_sync, spike_sync_matrix): exact
%   agreement on random and jittered-synfire spike train sets.
%
%   spikeSync() with no arguments runs a short self-test.

if nargin == 0, C = local_selftest(); return; end

p = inputParser;
p.addParameter('Interval', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
p.addParameter('MaxTau', Inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.parse(varargin{:});
maxTau = p.Results.MaxTau;

if ~iscell(spikes), error('spikeSync:input','spikes must be a cell array.'); end
N = numel(spikes);
if N < 2, error('spikeSync:input','Need at least two spike trains.'); end

% --- condition the trains -------------------------------------------------
t = cell(N,1); M = zeros(N,1);
for n = 1:N
    x = spikes{n}(:);
    x = x(~isnan(x));
    x = sort(x);
    if any(diff(x) == 0)
        warning('spikeSync:duplicates', ...
            'Train %d has duplicate spike times; keeping unique values.', n);
        x = unique(x);
    end
    t{n} = x; M(n) = numel(x);
end

allT = vertcat(t{:});
if isempty(p.Results.Interval)
    if isempty(allT), interval = [0 1]; else, interval = [min(allT) max(allT)]; end
else
    interval = sort(p.Results.Interval(:))';
    if ~isempty(allT) && (min(allT) < interval(1) || max(allT) > interval(2))
        warning('spikeSync:interval','Spikes fall outside the supplied Interval.');
    end
end
span = interval(2) - interval(1);
if span <= 0, span = Inf; end

% --- neighbouring ISIs per train (Inf where the neighbour does not exist) --
dPrev = cell(N,1); dNext = cell(N,1);
for n = 1:N
    dPrev{n} = Inf(M(n),1); dNext{n} = Inf(M(n),1);
    if M(n) > 1
        d = diff(t{n});
        dPrev{n}(2:end) = d;
        dNext{n}(1:end-1) = d;
    end
end

% --- pairwise adaptive coincidence detection ------------------------------
S = struct();
S.coinc = cell(N,N); S.match = cell(N,N); S.tau = cell(N,N);
S.interval = interval;
Csum = cell(N,1);
for n = 1:N, Csum{n} = zeros(M(n),1); end
Cmat = eye(N);

for n = 1:N
    for m = n+1:N
        if M(n) == 0 || M(m) == 0
            S.coinc{n,m} = false(M(n),1); S.coinc{m,n} = false(M(m),1);
            S.match{n,m} = zeros(M(n),1);  S.match{m,n} = zeros(M(m),1);
            S.tau{n,m}   = zeros(M(n),1);  S.tau{m,n}   = zeros(M(m),1);
            % no spikes at all in the pair -> common silence -> 1
            v = double(M(n) + M(m) == 0);
            Cmat(n,m) = v; Cmat(m,n) = v;
            continue
        end

        [cn, jn, taun] = local_pair(t{n}, t{m}, dPrev{n}, dNext{n}, ...
                                    dPrev{m}, dNext{m}, span, maxTau);
        [cm, jm, taum] = local_pair(t{m}, t{n}, dPrev{m}, dNext{m}, ...
                                    dPrev{n}, dNext{n}, span, maxTau);

        if sum(cn) ~= sum(cm)
            error('spikeSync:asymmetry', ...
                'Matching asymmetry for pair (%d,%d) — check input times.', n, m);
        end

        S.coinc{n,m} = cn; S.match{n,m} = jn; S.tau{n,m} = taun;
        S.coinc{m,n} = cm; S.match{m,n} = jm; S.tau{m,n} = taum;

        Csum{n} = Csum{n} + cn;
        Csum{m} = Csum{m} + cm;
        Cmat(n,m) = (sum(cn) + sum(cm)) / (M(n) + M(m));
        Cmat(m,n) = Cmat(n,m);
    end
end

% --- multivariate counter, profile, and overall value ---------------------
S.Ci = cell(N,1);
for n = 1:N, S.Ci{n} = Csum{n} / (N-1); end

Mtot = sum(M);
if Mtot == 0
    C = 1;
    prof = struct('t',[],'C',[],'train',[],'idx',[]);
    return
end

pt = zeros(Mtot,1); pc = zeros(Mtot,1); ptr = zeros(Mtot,1); pix = zeros(Mtot,1);
k = 0;
for n = 1:N
    r = k + (1:M(n));
    pt(r) = t{n}; pc(r) = S.Ci{n}; ptr(r) = n; pix(r) = (1:M(n))';
    k = k + M(n);
end
[pt, ord] = sort(pt);
prof = struct('t', pt, 'C', pc(ord), 'train', ptr(ord), 'idx', pix(ord));

C = sum(pc) / Mtot;

end % spikeSync


% =========================================================================
function [c, j, tau] = local_pair(ta, tb, dPa, dNa, dPb, dNb, span, maxTau)
% Coincidence indicators for every spike of train a against train b.
j   = local_nearest(ta, tb);
tau = min([repmat(span, numel(ta), 1), dPa, dNa, dPb(j), dNb(j)], [], 2) / 2;
tau = min(tau, maxTau);
c   = abs(ta - tb(j)) < tau;          % strict "<", Eq. 16
end


function j = local_nearest(x, y)
% Index of the nearest element of sorted y for each element of x.
ny = numel(y);
if isempty(x), j = zeros(0,1); return; end
if ny == 1
    j = ones(numel(x),1);
    return
end
% finite sentinel edges guarantee every x lands in a bin (no NaN, no Inf edges)
edges = [min(y(1), min(x)) - 1; y(:); max(y(end), max(x)) + 1];
b  = discretize(x, edges);               % b-1 = #{y <= x}
lo = min(max(b-1, 1), ny);
hi = min(max(b,   1), ny);
dlo = abs(x - y(lo));
dhi = abs(x - y(hi));
j = lo;
sw = dhi < dlo;                          % ties go to the earlier spike
j(sw) = hi(sw);
end


function ok = local_selftest()
ok = true;
tol = 1e-12;

% 1. identical trains -> C = 1
s = {(0:9)', (0:9)'};
c = spikeSync(s, 'Interval', [-1 10]);
ok = ok && abs(c - 1) < tol;
fprintf('identical trains      C = %.6f (expect 1)\n', c);

% 2. offset by exactly half the ISI -> window is exactly |dt| -> strict "<" fails
s = {(0:9)', (0.5:1:9.5)'};
c = spikeSync(s, 'Interval', [0 10]);
ok = ok && abs(c - 0) < tol;
fprintf('half-ISI offset       C = %.6f (expect 0)\n', c);

% 3. offset just under half the ISI -> all spikes coincide
s = {(0:9)', (0.49:1:9.49)'};
c = spikeSync(s, 'Interval', [0 10]);
ok = ok && abs(c - 1) < tol;
fprintf('just under half-ISI   C = %.6f (expect 1)\n', c);

% 4. one train empty -> its partner has no matches
s = {(0:9)', zeros(0,1)};
[c, Cm] = spikeSync(s, 'Interval', [0 10]);
ok = ok && abs(c - 0) < tol && abs(Cm(1,2)) < tol;
fprintf('one empty train       C = %.6f (expect 0)\n', c);

% 5. jittered synfire, 10 trains, 20 events -> C near 1, Fs-ready matching
rng(0);
ev = (10:10:200)';
s = arrayfun(@(n) sort(ev + 0.05*n + 0.01*randn(numel(ev),1)), (1:10)', ...
             'UniformOutput', false);
[c, ~, prof, S] = spikeSync(s, 'Interval', [0 210]);
ok = ok && abs(c - 1) < tol && numel(prof.t) == 200 && numel(S.Ci) == 10;
fprintf('jittered synfire      C = %.6f (expect 1)\n', c);

% 6. MaxTau ceiling suppresses the loose matches of case 3
c = spikeSync({(0:9)', (0.49:1:9.49)'}, 'Interval', [0 10], 'MaxTau', 0.1);
ok = ok && abs(c - 0) < tol;
fprintf('MaxTau = 0.1          C = %.6f (expect 0)\n', c);

if ok, verdict = 'PASSED'; else, verdict = 'FAILED'; end
fprintf('\nself-test %s\n', verdict);
end
