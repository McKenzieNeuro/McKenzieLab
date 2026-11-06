function D = spikeSyncDiagnose(spikes, varargin)
%SPIKESYNCDIAGNOSE  Work out what a SPIKE-Synchronization value actually means.
%
%   D = spikeSyncDiagnose(spikes, 'Dither', 0.02, ...)
%
%   C on its own is a disjunction: the same value is produced by weak diffuse
%   coupling, by a tight assembly diluted across a large population, by sparse
%   perfect synchrony sitting in a sea of background spikes, by propagation
%   delay, and by rate heterogeneity. This runs the checks that separate them
%   and returns the evidence rather than a verdict.
%
%   REQUIRED OPTION
%     'Dither'  : surrogate dither half-width, in the time units of the data.
%                 Sets the timescale the null is blind to. No default on
%                 purpose -- it is a scientific choice, not a nuisance knob.
%
%   OPTIONS
%     'Interval' : [t0 t1] recording window, passed to spikeSync.
%     'NSurr'    : surrogates (default 200).
%     'Width'    : if given, also compute a windowed trace (spikeSyncWindow).
%     'Step'     : window step (default Width/4).
%     'Plot'     : draw the four diagnostic panels (default true).
%
%   OUTPUT (fields, and what each one rules in or out)
%     .C, .Cmat, .prof     as returned by spikeSync
%     .null.C              surrogate mean C. Read .C against THIS, not against
%                          zero. Independent Poisson sits at 0.25.
%     .null.z              (C - mean_surr) / sd_surr
%     .null.p             (1 + #{surr >= obs}) / (1 + NSurr)
%
%     .spikes.aboveChance  fraction of observed spikes whose C_i exceeds the
%                          95th percentile of the pooled surrogate C_i. High
%                          C with a LOW value here means diffuse weak matching;
%                          low C with a HIGH value here means synchrony that is
%                          real but confined to a minority of spikes -- the
%                          background-dilution case, not asynchrony.
%     .spikes.deciles      deciles of prof.C. Unimodal near the null mean =
%                          diffuse. Bimodal = restricted but real.
%
%     .unit.C              per-unit mean C_i (N x 1)
%     .unit.rate           per-unit firing rate
%     .unit.rateRankCorr   Spearman rho between per-unit C and per-unit rate.
%                          Strongly negative => the rate-mismatch confound is
%                          driving the result, not timing.
%
%     .pair.rateRatio      N x N max/min rate ratio
%     .pair.rhoRateRatio   Spearman rho between off-diagonal Cmat and
%                          -log(rateRatio). Large positive => pairwise C is
%                          tracking rate similarity.
%
%     .assembly.loading    leading eigenvector of (Cmat - null Cmat), zero
%                          diagonal. Concentrated loading = a subpopulation
%                          carries the synchrony and the population C is
%                          diluted by 1/(N-1).
%     .assembly.eigRatio   lambda1 / lambda2. > ~2 indicates one dominant block.
%
%     .trace               [] unless 'Width' given: .Ct, .tc, .nSpk
%
%   No toolbox dependencies beyond base MATLAB.
%
%   See also SPIKESYNC, SPIKESYNCWINDOW, SPIKESYNCSURROGATE, SPIKESYNCCOMPARE.

p = inputParser;
p.addParameter('Dither', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Interval', [], @(x) isempty(x) || numel(x) == 2);
p.addParameter('NSurr', 200, @(x) isnumeric(x) && isscalar(x) && x >= 2);
p.addParameter('Width', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('Step', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('Plot', true, @(x) islogical(x) || isnumeric(x));
p.parse(varargin{:});
o = p.Results;
if isempty(o.Dither)
    error('spikeSyncDiagnose:dither', ...
        '''Dither'' is required: it sets the timescale the null destroys.');
end

N = numel(spikes);
syncArgs = {};
if ~isempty(o.Interval), syncArgs = {'Interval', o.Interval}; end

[C, Cmat, prof, S] = spikeSync(spikes, syncArgs{:});
D.C = C; D.Cmat = Cmat; D.prof = prof;

iv = S.interval; L = iv(2) - iv(1);
M = cellfun(@numel, spikes(:));
D.unit.rate = M / L;
D.unit.C = cellfun(@(c) mean(c), S.Ci(:));
D.unit.C(M == 0) = NaN;

% ---- surrogates ---------------------------------------------------------
wargs = {};
if ~isempty(o.Width)
    wargs = {'Width', o.Width};
    if ~isempty(o.Step), wargs = [wargs, {'Step', o.Step}]; end
end

Cs = zeros(o.NSurr,1);
CmatS = zeros(N,N,o.NSurr);
CiS = [];
CtS = [];
for s = 1:o.NSurr
    sur = cell(N,1);
    for n = 1:N
        x = spikes{n}(:);
        sur{n} = sort(x + o.Dither*(2*rand(numel(x),1) - 1));
    end
    [cS, cmS, pS] = spikeSync(sur, syncArgs{:});
    Cs(s) = cS; CmatS(:,:,s) = cmS;
    if s <= min(20, o.NSurr), CiS = [CiS; pS.C]; end %#ok<AGROW>
    if ~isempty(wargs)
        ct = spikeSyncWindow(pS, wargs{:});
        if isempty(CtS), CtS = nan(numel(ct), o.NSurr); end
        CtS(:,s) = ct; %#ok<AGROW>
    end
end
D.null.C  = mean(Cs);
D.null.sd = std(Cs);
D.null.z  = (C - mean(Cs)) / max(std(Cs), eps);
D.null.p  = (1 + sum(Cs >= C)) / (1 + o.NSurr);
D.null.Cmat = mean(CmatS, 3);

% ---- spike-level distribution ------------------------------------------
cut = local_pct(CiS, 95);
D.spikes.aboveChance = mean(prof.C > cut);
D.spikes.chanceCut   = cut;
D.spikes.deciles     = local_pct(prof.C, 10:10:90);

% ---- rate confounds -----------------------------------------------------
ok = M > 0;
D.unit.rateRankCorr = local_spearman(D.unit.C(ok), D.unit.rate(ok));

rr = max(D.unit.rate, D.unit.rate') ./ max(min(D.unit.rate, D.unit.rate'), eps);
D.pair.rateRatio = rr;
mask = triu(true(N), 1) & (M > 0) & (M' > 0);
D.pair.rhoRateRatio = local_spearman(Cmat(mask), -log(rr(mask)));

% ---- assembly structure -------------------------------------------------
A = Cmat - D.null.Cmat;
A(1:N+1:end) = 0;
A(~isfinite(A)) = 0;
A = (A + A')/2;
[V, E] = eig(A);
[ev, ord] = sort(abs(diag(E)), 'descend');
v = V(:, ord(1));
if sum(v) < 0, v = -v; end
D.assembly.loading  = v;
D.assembly.eigRatio = ev(1) / max(ev(2), eps);

% ---- optional trace -----------------------------------------------------
D.trace = [];
if ~isempty(wargs)
    [Ct, tc, nSpk] = spikeSyncWindow(prof, wargs{:});
    nW = numel(Ct);
    band = nan(nW,3);
    for w = 1:nW
        band(w,:) = local_pct(CtS(w,:), [2.5 50 97.5]);
    end
    pw = (1 + sum(CtS >= Ct, 2, 'omitnan')) ./ (1 + o.NSurr);
    pw(isnan(Ct)) = NaN;

    % cluster statistic: longest run of windows above the band, referred to
    % the same statistic computed on each surrogate. Pointwise p across
    % hundreds of overlapping windows is neither independent nor corrected.
    obsRun = local_maxRun(Ct > band(:,3));
    surRun = zeros(o.NSurr,1);
    for s = 1:o.NSurr
        surRun(s) = local_maxRun(CtS(:,s) > band(:,3));
    end
    D.trace = struct('Ct', Ct, 'tc', tc, 'nSpk', nSpk, 'null', band, ...
        'pPointwise', pw, 'clusterRun', obsRun, ...
        'clusterP', (1 + sum(surRun >= obsRun)) / (1 + o.NSurr));
    fprintf('  trace: longest supra-band run = %d windows, cluster p = %.4f\n', ...
        obsRun, D.trace.clusterP);
end

% ---- report -------------------------------------------------------------
fprintf('\nSPIKE-Synchronization diagnostics  (N = %d, M = %d, %.1f s)\n', N, sum(M), L);
fprintf('  C = %.4f   null = %.4f +/- %.4f   z = %+.2f   p = %.4f\n', ...
    C, D.null.C, D.null.sd, D.null.z, D.null.p);
fprintf('  spikes above chance C_i : %.3f  (chance cut %.3f)\n', ...
    D.spikes.aboveChance, cut);
fprintf('  rate:  median %.2f Hz, ratio max/min %.1f\n', ...
    median(D.unit.rate), max(D.unit.rate)/max(min(D.unit.rate), eps));
fprintf('  rho(unit C, unit rate)        = %+.2f  %s\n', D.unit.rateRankCorr, ...
    local_flag(abs(D.unit.rateRankCorr) > 0.5, 'RATE CONFOUND LIKELY'));
fprintf('  rho(pair C, rate similarity)  = %+.2f  %s\n', D.pair.rhoRateRatio, ...
    local_flag(D.pair.rhoRateRatio > 0.5, 'PAIRWISE C TRACKS RATE'));
fprintf('  assembly lambda1/lambda2      = %.2f  %s\n', D.assembly.eigRatio, ...
    local_flag(D.assembly.eigRatio > 2, 'SUBPOPULATION STRUCTURE'));
if D.C < D.null.C + 2*D.null.sd && D.spikes.aboveChance > 0.1
    fprintf('  NOTE: population C is at chance but %.0f%% of spikes match above\n', ...
        100*D.spikes.aboveChance);
    fprintf('        chance -- restricted synchrony diluted by N or by background.\n');
end
fprintf('\n');

if o.Plot, local_plot(D, N); end

end


% =========================================================================
function local_plot(D, N)
figure('Color','w','Name','spikeSyncDiagnose');

subplot(2,2,1);
edges = linspace(0, 1, 26);
histogram(D.prof.C, edges, 'Normalization','probability', ...
    'FaceColor',[.4 .4 .4], 'EdgeColor','none'); hold on
yl = ylim;
plot([D.spikes.chanceCut D.spikes.chanceCut], yl, 'r--', 'LineWidth', 1.2);
xlabel('C_i (per spike)'); ylabel('fraction of spikes');
title('spike-level distribution'); legend('observed','chance cut','Location','best');
box off

subplot(2,2,2);
imagesc(D.Cmat - D.null.Cmat); axis square; colorbar
xlabel('spike train'); ylabel('spike train');
title(sprintf('C_{obs} - C_{null}  (\\lambda_1/\\lambda_2 = %.1f)', D.assembly.eigRatio));

subplot(2,2,3);
plot(D.unit.rate, D.unit.C, 'ko', 'MarkerFaceColor',[.6 .6 .6]);
xlabel('firing rate (Hz)'); ylabel('mean C_i');
title(sprintf('rate confound  \\rho = %+.2f', D.unit.rateRankCorr));
set(gca,'XScale','log'); box off

subplot(2,2,4);
if ~isempty(D.trace)
    tc = D.trace.tc(:);
    lo = D.trace.null(:,1); hi = D.trace.null(:,3);
    k = ~isnan(lo) & ~isnan(hi);
    fill([tc(k); flipud(tc(k))], [lo(k); flipud(hi(k))], [.85 .85 .85], ...
        'EdgeColor','none'); hold on
    plot(tc, D.trace.null(:,2), 'r--');
    plot(tc, D.trace.Ct, 'k-', 'LineWidth', 1.1);
    xlabel('time'); ylabel('C (windowed)');
    title(sprintf('time-resolved C vs surrogate band (cluster p = %.3f)', ...
        D.trace.clusterP));
    box off
else
    stem(1:N, D.assembly.loading, 'k', 'filled', 'MarkerSize', 3);
    xlabel('spike train'); ylabel('leading eigenvector');
    title('assembly loading'); box off
end
end


function s = local_flag(tf, msg)
if tf, s = ['<-- ' msg]; else, s = ''; end
end


function y = local_pct(x, q)
% percentiles without the Statistics Toolbox (linear interpolation, MATLAB
% 'prctile' midpoint convention)
x = sort(x(~isnan(x(:))));
n = numel(x);
if n == 0, y = nan(size(q)); return; end
if n == 1, y = repmat(x, size(q)); return; end
pos = (0.5:n-0.5) / n * 100;
y = interp1(pos, x, q, 'linear');
y(q < pos(1))   = x(1);
y(q > pos(end)) = x(end);
end


function rho = local_spearman(x, y)
% rank correlation without the Statistics Toolbox
x = x(:); y = y(:);
k = ~isnan(x) & ~isnan(y) & isfinite(x) & isfinite(y);
x = x(k); y = y(k);
if numel(x) < 3, rho = NaN; return; end
rx = local_tiedrank(x); ry = local_tiedrank(y);
rx = rx - mean(rx); ry = ry - mean(ry);
rho = (rx' * ry) / max(sqrt((rx'*rx) * (ry'*ry)), eps);
end


function r = local_tiedrank(x)
[xs, ord] = sort(x(:));
r = zeros(numel(x),1);
r(ord) = 1:numel(x);
% average ranks within ties
i = 1;
while i <= numel(xs)
    j = i;
    while j < numel(xs) && xs(j+1) == xs(i), j = j + 1; end
    if j > i
        r(ord(i:j)) = mean(i:j);
    end
    i = j + 1;
end
end


function r = local_maxRun(tf)
% longest run of true values
tf = logical(tf(:)); tf(isnan(double(tf))) = false;
r = 0; c = 0;
for i = 1:numel(tf)
    if tf(i), c = c + 1; if c > r, r = c; end, else, c = 0; end
end
end
