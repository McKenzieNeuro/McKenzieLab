function out = decode_continuous_tts(all_tim2seiz, all_pred, session, opts)
%DECODE_CONTINUOUS_TTS  Continuous time-to-seizure from RUSBoost class scores.
%
%   Maps the 6 class scores onto a continuous time-to-seizure estimate,
%   cross-validated leave-one-session-out, with empirical prediction
%   intervals. Two decoders are fit side by side:
%
%     ridge     : ridge regression from log class scores -> log10(t).
%                 Free to exploit the graded decay of the class-1 score,
%                 which is where the mid-range information actually lives.
%     centroid  : parameter-free posterior mean, sum_k w_k * mu_k, with
%                 mu_k the training-fold mean log10(t) of class k and
%                 w_k the scores renormalised over the TTE classes (1:4).
%                 This is the "convert the class posterior back to
%                 continuous time" decoder. It is the baseline the ridge
%                 fit has to beat to justify itself.
%
%   NOTE ON CV. The classifier scores are already out-of-session, but the
%   score -> time mapping is a second fit and must be cross-validated the
%   same way or the continuous decode is contaminated. That is what the
%   LOSO loop below is for. Do not fit this on pooled data.
%
%   INPUTS
%     all_tim2seiz : [N x 1] true seconds to next seizure (1-s bins)
%     all_pred     : [N x 7] col 1 = MAP class, cols 2:7 = class scores
%     session      : [N x 1] session / recording-day ID (REQUIRED)
%     opts         : struct, optional
%        .lambda    ridge penalty (per-sample scaled)        default 1e-3
%        .decimate  keep every k-th row for evaluation       default 60
%        .maxT      cap on t before taking log               default Inf
%        .exclude   [N x 1] logical, rows to drop            default none
%        .edges     TTE class edges, seconds            default [0 10 100 3600 Inf]
%        .plotOn    draw diagnostics                         default true
%
%   OUTPUT struct with per-row predictions, intervals and summary metrics.

if nargin < 4, opts = struct; end
def = struct('lambda',1e-3,'decimate',60,'maxT',Inf,'exclude',[], ...
             'edges',[0 10 100 3600 Inf],'plotOn',true, ...
             'classOrder',[4 3 2 1],'balance',true);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts,fn{i}) || isempty(opts.(fn{i}))
        opts.(fn{i}) = def.(fn{i});
    end
end

t    = double(all_tim2seiz(:));
S    = double(all_pred(:,2:7));
sess = double(session(:));
N    = numel(t);
assert(size(all_pred,1)==N && numel(sess)==N, ...
    'all_tim2seiz, all_pred and session must have the same number of rows.');

%% ---------- row selection -------------------------------------------
% t > 0 drops ictal rows. Postictal rows have a valid time to the NEXT
% seizure and are kept; pass opts.exclude to drop them if you want the
% decode to be pre-ictal only (recommended for the headline number, since
% postictal EEG is trivially decodable and carries renewal information,
% not prodromal information).
% Sign convention: some pipelines code pre-ictal time as NEGATIVE seconds
% (-3600 = one hour before onset). Detect and flip so t is a positive
% "seconds until the next seizure" throughout.
fin = isfinite(t);
if any(fin) && all(t(fin) <= 0)
    warning('decode_continuous_tts:signFlip', ...
        'All finite t <= 0; interpreting as negative-coded pre-ictal time and flipping sign.');
    t = -t;
end

keep = isfinite(t) & t > 0 & isfinite(sess) & all(isfinite(S),2);
if ~isempty(opts.exclude)
    keep = keep & ~logical(opts.exclude(:));
end

assert(sum(keep) > 100, ...
    ['Only %d rows survived selection. Check the sign convention of ' ...
     'all_tim2seiz (this function expects positive seconds-to-next-seizure, ' ...
     'and drops t <= 0 as ictal).'], sum(keep));
assert(numel(unique(sess(keep))) > 1, ...
    'Need >1 session for leave-one-session-out CV; found %d.', ...
    numel(unique(sess(keep))));

rowIdx = find(keep);
t = min(t(keep), opts.maxT);
S = S(keep,:);
sess = sess(keep);

%% ---------- target ---------------------------------------------------
% log10 because t spans four decades and squared error in linear time is
% entirely dominated by the far tail.
y = log10(t);

%% ---------- features -------------------------------------------------
tot = sum(S,2);
Sn  = S ./ max(tot,eps);                 % normalised scores
X   = [log(Sn + 1e-6), log(max(tot,eps))];
X   = X(:, std(X,0,1) > 0);              % drop constants (scores may sum to 1)
X   = [X, ones(size(X,1),1)];            % intercept LAST
p   = size(X,2);

%% ---------- true TTE class (derived from t) --------------------------
cls = discretize(t, opts.edges);         % 1 = 0-10s, 2 = 10-100s,
nCls = numel(opts.edges) - 1;            % 3 = 100-3600s, 4 = >3600s

% The classifier codes 1 = >1hr descending to 4 = 0-10s, but discretize()
% above returns 1 = 0-10s ascending. Reorder the score columns into
% ascending time so the centroid decoder pairs score k with centroid k.
Sasc = Sn(:, opts.classOrder);

%% ---------- leave-one-session-out ------------------------------------
us    = unique(sess);
yhatR = nan(size(y));   yhatC = nan(size(y));
lo    = nan(size(y));   hi    = nan(size(y));
ynull = nan(size(y));

P = eye(p); P(end,end) = 0;              % do not penalise the intercept

for k = 1:numel(us)
    te = sess == us(k);
    tr = ~te;
    if ~any(tr) || ~any(te), continue; end

    Xtr = X(tr,:);  ytr = y(tr);
    gtr = cls(tr);

    % --- inverse-density weights over log-decades -------------------
    % Without this, samples per decade grow geometrically and least
    % squares has no incentive to get short times right: the fit
    % collapses onto the grand mean. Same correction RUSBoost applies
    % to the classifier, in regression form.
    if opts.balance
        okg = ~isnan(gtr);
        cnt = accumarray(gtr(okg), 1, [nCls 1], @sum, 0);
        wv  = zeros(size(ytr));
        wv(okg) = 1 ./ max(cnt(gtr(okg)), 1);
        wv  = wv / mean(wv(okg));
    else
        wv = ones(size(ytr));
    end

    % --- ridge (weighted) ---
    XW = Xtr .* wv;
    w = (XW'*Xtr + opts.lambda*numel(ytr)*P) \ (XW'*ytr);
    fitTr     = Xtr*w;
    yhatR(te) = X(te,:)*w;

    % --- centroid / posterior-mean decoder ---
    ok = ~isnan(gtr);
    mu = accumarray(gtr(ok), ytr(ok), [nCls 1], @mean, NaN);
    w4 = Sasc ./ max(sum(Sasc,2), eps);
    good = ~isnan(mu);
    yhatC(te) = (w4(te,good) ./ max(sum(w4(te,good),2),eps)) * mu(good);

    % --- null: training-fold mean ---
    ynull(te) = mean(ytr);

    % --- empirical 90% interval from training residuals, conditioned on
    %     the fitted value (heteroscedastic by construction) ------------
    rtr = ytr - fitTr;
    qe  = unique(quantile(fitTr, linspace(0,1,11)));
    if numel(qe) > 2
        qe(1) = -inf; qe(end) = inf;
        btr = discretize(fitTr, qe);
        bte = discretize(yhatR(te), qe);
        ok2 = ~isnan(btr);
        nb  = numel(qe)-1;
        loB = accumarray(btr(ok2), rtr(ok2), [nb 1], @(v) quantile(v,0.05), NaN);
        hiB = accumarray(btr(ok2), rtr(ok2), [nb 1], @(v) quantile(v,0.95), NaN);
        bte(isnan(bte)) = 1;
        lo(te) = yhatR(te) + loB(bte);
        hi(te) = yhatR(te) + hiB(bte);
    end
end

%% ---------- evaluation -----------------------------------------------
% Adjacent 1-s rows are near-duplicates. Decimating does not fix the
% dependence but it stops n from being inflated by ~3 orders of magnitude.
% The honest unit of analysis is still the seizure, not the second.
ev = false(size(y));
ev(1:opts.decimate:end) = true;
ev = ev & isfinite(yhatR) & isfinite(ynull);

sse = @(yh) sum((y(ev)-yh(ev)).^2, 'omitnan');
out.R2_ridge     = 1 - sse(yhatR)/sse(ynull);
out.R2_centroid  = 1 - sse(yhatC)/sse(ynull);
out.MAE_log10    = median(abs(y(ev)-yhatR(ev)), 'omitnan');
out.foldError    = 10.^out.MAE_log10;      % multiplicative error
out.rho_spearman = corr(y(ev), yhatR(ev), 'type','Spearman','rows','complete');
out.coverage90   = mean(y(ev) >= lo(ev) & y(ev) <= hi(ev), 'omitnan');
out.nEval        = sum(ev);
out.nSessions    = numel(us);

% collapse diagnostic: if the decode has regressed to the mean, the
% predicted spread is a small fraction of the true spread and everything
% below is meaningless regardless of R^2
out.std_true  = std(y(ev), 'omitnan');
out.std_pred  = std(yhatR(ev), 'omitnan');
out.spreadRatio = out.std_pred / out.std_true;
out.range_pred  = [min(yhatR(ev)) max(yhatR(ev))];

% per-class bias and error: this is where you see whether the decode is
% just regressing everything to the grand mean
for i = 1:nCls
    m = ev & cls == i;
    out.byClass(i).edges = opts.edges(i:i+1);
    out.byClass(i).n     = sum(m);
    out.byClass(i).bias  = median(yhatR(m)-y(m), 'omitnan');
    out.byClass(i).mae   = median(abs(yhatR(m)-y(m)), 'omitnan');
    out.byClass(i).medPred = median(yhatR(m), 'omitnan');
end
% decade-balanced R^2: each decade contributes equally, so a mean-collapse
% scores ~0 here even when the unweighted R^2 looks respectable
out.R2_balanced = 1 - mean([out.byClass.mae].^2) / ...
                      mean(arrayfun(@(i) median(abs(y(ev & cls==i) - ...
                          mean(y(ev),'omitnan')),'omitnan').^2, 1:nCls), 'omitnan');

% per-row output, mapped back to original row indices
out.rowIdx    = rowIdx;
out.y_true    = y;
out.yhat      = yhatR;
out.yhat_cent = yhatC;
out.lo        = lo;
out.hi        = hi;
out.session   = sess;
out.opts      = opts;

%% ---------- diagnostics ----------------------------------------------
if opts.plotOn
    figure('Color','w','Position',[100 100 1100 340]);

    subplot(1,3,1);
    histogram2(y(ev), yhatR(ev), 60, 'DisplayStyle','tile', ...
        'ShowEmptyBins','off','Normalization','pdf');
    hold on; lim = [min(y(ev)) max(y(ev))];
    plot(lim, lim, 'w--', 'LineWidth', 1.2);
    xlabel('true log_{10} t (s)'); ylabel('decoded log_{10} t (s)');
    title(sprintf('R^2 = %.3f (centroid %.3f)', out.R2_ridge, out.R2_centroid));
    axis square; colorbar;

    subplot(1,3,2);
    edgesY = linspace(min(y(ev)), max(y(ev)), 25);
    b = discretize(y(ev), edgesY);
    yy = yhatR(ev); tt = y(ev); ok3 = ~isnan(b);
    ctr = edgesY(1:end-1) + diff(edgesY)/2;
    m50 = accumarray(b(ok3), yy(ok3), [numel(ctr) 1], @median, NaN);
    m05 = accumarray(b(ok3), yy(ok3), [numel(ctr) 1], @(v) quantile(v,0.05), NaN);
    m95 = accumarray(b(ok3), yy(ok3), [numel(ctr) 1], @(v) quantile(v,0.95), NaN);
    plot(ctr, m50, 'k-', 'LineWidth', 1.5); hold on;
    plot(ctr, m05, 'k:', ctr, m95, 'k:');
    plot(ctr, ctr, 'r--');
    xlabel('true log_{10} t (s)'); ylabel('decoded log_{10} t (s)');
    title('median + 5/95 pct'); axis square; grid on;

    subplot(1,3,3);
    cov = accumarray(b(ok3), double(tt(ok3) >= lo(subsref(find(ev),substruct('()',{ok3}))) & ...
                                    tt(ok3) <= hi(subsref(find(ev),substruct('()',{ok3})))), ...
                     [numel(ctr) 1], @mean, NaN);
    plot(ctr, cov, 'k-', 'LineWidth', 1.5); hold on;
    yline(0.90, 'r--');
    xlabel('true log_{10} t (s)'); ylabel('90% interval coverage');
    ylim([0 1]); title(sprintf('overall %.3f', out.coverage90));
    axis square; grid on;
end
end