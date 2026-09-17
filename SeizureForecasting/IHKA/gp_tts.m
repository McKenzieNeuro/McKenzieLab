function out = gp_tts(all_tim2seiz, all_pred, session, nTrain, kFold, dec)
%GP_TTS  GP regression: class scores -> log10 time-to-seizure.
%
%   out = gp_tts(all_tim2seiz, all_pred, session)
%
%   Session-held-out k-fold. Target is log10(seconds to seizure).
%   Training rows are subsampled evenly across half-decades of true time,
%   otherwise the fit collapses onto the grand mean (samples per decade
%   grow geometrically). Predictive SD comes straight from the GP.
%
%   nTrain : training rows per fold   (default 2500; exact GP is O(n^3))
%   kFold  : session folds            (default 5)
%   dec    : keep every dec-th test row for prediction (default 20)

if nargin < 4 || isempty(nTrain), nTrain = 2500; end
if nargin < 5 || isempty(kFold),  kFold  = 5;    end
if nargin < 6 || isempty(dec),    dec    = 20;   end

t = double(all_tim2seiz(:));
if any(isfinite(t)) && all(t(isfinite(t)) <= 0), t = -t; end   % negative-coded
X = double(all_pred(:, 2:7));
s = double(session(:));

% t < 1 s is below the 1-s window resolution: those values are the
% sub-second offset between seizure onset and the window grid, not
% resolvable time. They sit inside the final second before onset, are
% likely ictal-contaminated, and get wildly upweighted by the balanced
% subsampler below (a handful of rows filling a whole stratum).
keep = isfinite(t) & t >= 1 & all(isfinite(X),2) & isfinite(s);
t = t(keep); X = X(keep,:); s = s(keep);
y = log10(t);

% --- assign sessions to folds ---
us   = unique(s);
fmap = mod((1:numel(us))'-1, kFold) + 1;
foldOf = zeros(size(s));
for i = 1:numel(us), foldOf(s == us(i)) = fmap(i); end

% --- half-decade strata for balanced training subsample ---
edges = floor(min(y)*2)/2 : 0.5 : ceil(max(y)*2)/2;
strat = discretize(y, edges);

mu = nan(size(y)); sd = nan(size(y)); tested = false(size(y));

for f = 1:kFold
    te = find(foldOf == f);
    tr = find(foldOf ~= f);
    if isempty(te) || isempty(tr), continue; end

    % balanced subsample
    g  = strat(tr);
    ug = unique(g(~isnan(g)));
    nPer = max(20, floor(nTrain / numel(ug)));
    sel = [];
    for j = 1:numel(ug)
        idx = tr(g == ug(j));
        sel = [sel; idx(randperm(numel(idx), min(nPer, numel(idx))))]; %#ok<AGROW>
    end

    gp = fitrgp(X(sel,:), y(sel), ...
        'KernelFunction','ardsquaredexponential', ...
        'Standardize', true, ...
        'FitMethod','exact', 'PredictMethod','exact');

    teD = te(1:dec:end);                      % decimate test rows
    [m, ssd] = predict(gp, X(teD,:));
    mu(teD) = m; sd(teD) = ssd; tested(teD) = true;

    fprintf('fold %d/%d  n_train=%d  n_test=%d\n', f, kFold, numel(sel), numel(teD));
end

% --- output ---
ev = tested & isfinite(mu);
out.y_true = y(ev);
out.mu     = mu(ev);
out.sd     = sd(ev);
out.session = s(ev);
out.rho    = corr(out.y_true, out.mu, 'type','Spearman');
out.spreadRatio = std(out.mu) / std(out.y_true);
out.coverage95  = mean(abs(out.y_true - out.mu) <= 1.96*out.sd);

% --- one figure: decoded vs true, binned, with GP SD ---
figure('Color','w','Position',[100 100 760 320]);

subplot(1,2,1);
histogram2(out.y_true, out.mu, 50, 'DisplayStyle','tile','ShowEmptyBins','off');
hold on; L = [min(out.y_true) max(out.y_true)];
plot(L, L, 'w--','LineWidth',1.2);
xlabel('true log_{10} t (s)'); ylabel('GP mean log_{10} t (s)');
title(sprintf('\\rho = %.3f, spread %.2f', out.rho, out.spreadRatio));
axis square;

subplot(1,2,2);
be = linspace(min(out.y_true), max(out.y_true), 25);
b  = discretize(out.y_true, be);
c  = be(1:end-1) + diff(be)/2;
ok = ~isnan(b);
med = accumarray(b(ok), out.mu(ok),  [numel(c) 1], @median, NaN);
sdm = accumarray(b(ok), out.sd(ok),  [numel(c) 1], @median, NaN);
errorbar(c, med, 1.96*sdm, 'k-','LineWidth',1.2); hold on;
plot(c, c, 'r--');
xlabel('true log_{10} t (s)'); ylabel('GP mean \pm 1.96 SD');
title(sprintf('95%% coverage %.3f', out.coverage95));
axis square; grid on;
end