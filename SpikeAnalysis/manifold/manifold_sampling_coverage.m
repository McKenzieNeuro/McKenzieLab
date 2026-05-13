function [out] = manifold_sampling_coverage(ses)
% MANIFOLD_SAMPLING_COVERAGE  Quantifies how broadly stimulated ripples
%   sample the replay manifold relative to spontaneous ripples.
%
% Two complementary statistics per block:
%
%   Effective cluster count (ECC):
%     ECC = exp(H(p))  where H = -sum(p_m * log(p_m))
%     Ranges from 1 (all ripples in one cluster) to K (uniform sampling)
%     Answers: how many clusters are effectively sampled?
%
%   Coverage fraction:
%     Proportion of clusters receiving >= N_min ripples (default N_min = 2)
%     Answers: does silencing leave any clusters completely unvisited?
%
% Comparisons (paired across blocks within session):
%   stim vs spont ECC        — does silencing change sampling breadth?
%   stim vs test ECC         — same comparison against held-out baseline
%   stim vs spont coverage   — does silencing exclude any clusters?
%
% ECC ratio = ECC_stim / ECC_spont:
%   ~1  : coverage preserved
%   <1  : stim concentrates in fewer clusters
%   >1  : stim spreads more broadly than spont
%
% Usage:
%   out = manifold_sampling_coverage(ses);

N_MIN = 1;   % minimum ripples to count a cluster as visited
             % (after size-matching, counts are low so threshold is 1)

nSes  = length(ses);
nFold = length(ses(1).tensor.block);

% ── Pooled storage (one value per block) ──────────────────────────────────
ECC_stim_all  = [];
ECC_spont_all = [];
ECC_test_all  = [];
cov_stim_all  = [];
cov_spont_all = [];
cov_test_all  = [];
ratio_all     = [];
ses_id        = [];

out.session = struct();

for i = 1:nSes

    ECC_stim_b  = nan(nFold, 1);
    ECC_spont_b = nan(nFold, 1);
    ECC_test_b  = nan(nFold, 1);
    cov_stim_b  = nan(nFold, 1);
    cov_spont_b = nan(nFold, 1);
    cov_test_b  = nan(nFold, 1);

    for j = 1:nFold

        blk      = ses(i).tensor.block(j);
        kmean_id = blk.kmean_id(:);
        train_ix = blk.training_ix(:);
        stim_ix  = blk.stim_ix(:);
        test_ix  = blk.test_ix(:);
        K        = blk.optimal_k;

        if isempty(stim_ix) || isempty(train_ix), continue; end

        % ── Size-matched subsampling ───────────────────────────────────────
        % ECC scales with N — subsample all conditions to the same count
        % so differences reflect genuine redistribution, not sample size
        N_match = min([length(stim_ix), length(test_ix)]);
        if N_match < 5, continue; end

        stim_ix_s  = stim_ix(randsample(length(stim_ix),   N_match));
        test_ix_s  = test_ix(randsample(length(test_ix),   N_match));
        train_ix_s = train_ix(randsample(length(train_ix), N_match));

        % ── Per-cluster counts ─────────────────────────────────────────────
        N_stim  = accumarray(kmean_id(stim_ix_s),  1, [K 1], @sum, 0);
        N_spont = accumarray(kmean_id(train_ix_s), 1, [K 1], @sum, 0);
        N_test  = accumarray(kmean_id(test_ix_s),  1, [K 1], @sum, 0);

        % ── Effective cluster count ────────────────────────────────────────
        ECC_stim_b(j)  = effective_count(N_stim);
        ECC_spont_b(j) = effective_count(N_spont);
        ECC_test_b(j)  = effective_count(N_test);

        % ── Coverage fraction ──────────────────────────────────────────────
        cov_stim_b(j)  = mean(N_stim  >= N_MIN);
        cov_spont_b(j) = mean(N_spont >= N_MIN);
        cov_test_b(j)  = mean(N_test  >= N_MIN);

    end

    out.session(i).ECC_stim   = ECC_stim_b;
    out.session(i).ECC_spont  = ECC_spont_b;
    out.session(i).ECC_test   = ECC_test_b;
    out.session(i).cov_stim   = cov_stim_b;
    out.session(i).cov_spont  = cov_spont_b;
    out.session(i).cov_test   = cov_test_b;
    out.session(i).ratio      = ECC_stim_b ./ (ECC_spont_b + 1e-12);

    vb = ~isnan(ECC_stim_b) & ~isnan(ECC_spont_b);
    ECC_stim_all  = [ECC_stim_all;  ECC_stim_b(vb) ];
    ECC_spont_all = [ECC_spont_all; ECC_spont_b(vb)];
    ECC_test_all  = [ECC_test_all;  ECC_test_b(vb) ];
    cov_stim_all  = [cov_stim_all;  cov_stim_b(vb) ];
    cov_spont_all = [cov_spont_all; cov_spont_b(vb)];
    cov_test_all  = [cov_test_all;  cov_test_b(vb) ];
    ratio_all     = [ratio_all;     ECC_stim_b(vb) ./ (ECC_spont_b(vb) + 1e-12)];
    ses_id        = [ses_id;        repmat(i, sum(vb), 1)];

end

% ── Session-level averages (average blocks within session first) ──────────
nSes_valid = 0;
ECC_stim_ses  = nan(nSes, 1);
ECC_spont_ses = nan(nSes, 1);
ECC_test_ses  = nan(nSes, 1);
cov_stim_ses  = nan(nSes, 1);
cov_spont_ses = nan(nSes, 1);
cov_test_ses  = nan(nSes, 1);
ratio_ses     = nan(nSes, 1);

for i = 1:nSes
    ECC_stim_ses(i)  = nanmedian(out.session(i).ECC_stim);
    ECC_spont_ses(i) = nanmedian(out.session(i).ECC_spont);
    ECC_test_ses(i)  = nanmedian(out.session(i).ECC_test);
    cov_stim_ses(i)  = nanmedian(out.session(i).cov_stim);
    cov_spont_ses(i) = nanmedian(out.session(i).cov_spont);
    cov_test_ses(i)  = nanmedian(out.session(i).cov_test);
    ratio_ses(i)     = ECC_stim_ses(i) / (ECC_spont_ses(i) + 1e-12);
end

% Keep only sessions with valid data
valid_ses = ~isnan(ECC_stim_ses) & ~isnan(ECC_spont_ses);
ECC_stim_ses  = ECC_stim_ses(valid_ses);
ECC_spont_ses = ECC_spont_ses(valid_ses);
ECC_test_ses  = ECC_test_ses(valid_ses);
cov_stim_ses  = cov_stim_ses(valid_ses);
cov_spont_ses = cov_spont_ses(valid_ses);
cov_test_ses  = cov_test_ses(valid_ses);
ratio_ses     = ratio_ses(valid_ses);
nSes_valid    = sum(valid_ses);

% ── Statistics (session-level, N = nSes_valid) ────────────────────────────
[~, p_stim_spont, ~, ts1] = ttest(ECC_stim_ses,  ECC_spont_ses);
[~, p_stim_test,  ~, ts2] = ttest(ECC_stim_ses,  ECC_test_ses);
[~, p_test_spont, ~, ts3] = ttest(ECC_test_ses,  ECC_spont_ses);
[~, p_cov,        ~, ts4] = ttest(cov_stim_ses,  cov_spont_ses);
[~, p_ratio,      ~, ts5] = ttest(ratio_ses, 1);

% ── Report ────────────────────────────────────────────────────────────────
fprintf('\n%s\n', repmat('=',1,68));
fprintf('Manifold sampling coverage — N sessions: %d  (%d blocks total)\n', ...
    nSes_valid, length(ECC_stim_all));
fprintf('%s\n', repmat('-',1,68));
fprintf('Effective cluster count (session median, mean ± SD across sessions):\n');
fprintf('  Spont (train):  %.2f ± %.2f\n', ...
    mean(ECC_spont_ses), std(ECC_spont_ses));
fprintf('  Spont (test):   %.2f ± %.2f  [baseline]\n', ...
    mean(ECC_test_ses),  std(ECC_test_ses));
fprintf('  Stim:           %.2f ± %.2f\n', ...
    mean(ECC_stim_ses),  std(ECC_stim_ses));
fprintf('\n');
fprintf('ECC ratio (stim/spont): %.3f ± %.3f\n', ...
    mean(ratio_ses), std(ratio_ses));
fprintf('  Ratio vs 1:     t(%d) = %.3f, p = %.4f\n', ...
    ts5.df, ts5.tstat, p_ratio);
fprintf('\n');
fprintf('Paired t-tests on ECC (N = %d sessions):\n', nSes_valid);
fprintf('  Stim vs spont:  t(%d) = %.3f, p = %.4f\n', ...
    ts1.df, ts1.tstat, p_stim_spont);
fprintf('  Stim vs test:   t(%d) = %.3f, p = %.4f\n', ...
    ts2.df, ts2.tstat, p_stim_test);
fprintf('  Test vs spont:  t(%d) = %.3f, p = %.4f  [should be ~0]\n', ...
    ts3.df, ts3.tstat, p_test_spont);
fprintf('%s\n', repmat('-',1,68));
fprintf('Coverage fraction (>= %d ripple, session median, mean ± SD):\n', N_MIN);
fprintf('  Spont:  %.3f ± %.3f\n', mean(cov_spont_ses), std(cov_spont_ses));
fprintf('  Test:   %.3f ± %.3f  [baseline]\n', ...
    mean(cov_test_ses), std(cov_test_ses));
fprintf('  Stim:   %.3f ± %.3f\n', mean(cov_stim_ses), std(cov_stim_ses));
fprintf('  Stim vs spont:  t(%d) = %.3f, p = %.4f\n', ...
    ts4.df, ts4.tstat, p_cov);
fprintf('%s\n', repmat('=',1,68));
fprintf('Interpretation:\n');
fprintf('  ECC ratio ~1, coverage ~spont => full manifold sampled\n');
fprintf('  ECC ratio <1, coverage <spont => silencing excludes some clusters\n');
fprintf('  ECC ratio >1, coverage ~spont => silencing spreads more broadly\n');

% ── Pack outputs ──────────────────────────────────────────────────────────
out.ECC_stim_all  = ECC_stim_all;
out.ECC_spont_all = ECC_spont_all;
out.ECC_test_all  = ECC_test_all;
out.cov_stim_all  = cov_stim_all;
out.cov_spont_all = cov_spont_all;
out.cov_test_all  = cov_test_all;
out.ratio_all     = ratio_all;
out.ses_id        = ses_id;

% Session-level summaries
out.ECC_stim_ses  = ECC_stim_ses;
out.ECC_spont_ses = ECC_spont_ses;
out.ECC_test_ses  = ECC_test_ses;
out.cov_stim_ses  = cov_stim_ses;
out.cov_spont_ses = cov_spont_ses;
out.cov_test_ses  = cov_test_ses;
out.ratio_ses     = ratio_ses;

out.stats.p_stim_spont = p_stim_spont;
out.stats.p_stim_test  = p_stim_test;
out.stats.p_test_spont = p_test_spont;
out.stats.p_cov        = p_cov;
out.stats.p_ratio      = p_ratio;
out.stats.ts_ratio     = ts5;

end

% ── Helper ────────────────────────────────────────────────────────────────
function ecc = effective_count(N)
% Effective cluster count from raw counts N (K x 1)
    N   = N(N > 0);
    if isempty(N), ecc = nan; return; end
    p   = N / sum(N);
    H   = -sum(p .* log(p + 1e-12));
    ecc = exp(H);
end