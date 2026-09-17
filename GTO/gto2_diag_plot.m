function gto2_diag_plot(R, covariate, covariate_name)
% GTO2_DIAG_PLOT  Does gate failure cluster, or is it uniform across probes?
%
%   gto2_diag_plot(R)                          % vs probe order (time) only
%   gto2_diag_plot(R, covariate, 'velocity')   % vs an external covariate too
%
% INPUTS
%   R               output struct from gto2_preflight
%   covariate       [n_prb x 1] optional, e.g. velocity, time-in-session,
%                   condition id -- anything matched 1:1 to the probes R
%                   was called on. Omit or pass [] to skip those panels.
%   covariate_name  string for axis labels (default 'covariate')
%
% Six panels:
%   1. q (leakage ratio) vs probe order, gate-1 pass/fail coloured
%   2. d_hat vs probe order
%   3. q vs covariate (if given)
%   4. C6 vs covariate, split by gate-1 pass/fail (checks whether gate-1
%      failures are concentrated in a C6 regime, i.e. confounded with
%      conditioning rather than independent of it)
%   5. histogram of d_hat, with the count at each integer value labelled
%   6. per-probe pass/fail across ALL evaluable gates simultaneously
%      (a probe failing q but nothing else looks different from one failing
%      everything -- this panel is the only place that distinction is visible)
%
% Reads directly from R.diag (table) and R.gate; does not recompute anything.

if nargin < 2, covariate = []; end
if nargin < 3, covariate_name = 'covariate'; end

n = height(R.diag);
x = (1:n)';
g1pass = R.gate(:,2) ~= 0;   % column 2 = gate 1 (dim)
q  = R.diag.q;
dh = R.diag.d_hat;
C6 = R.diag.C6;

figure('Name','gto2 gate diagnostics','Color','w');
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

% --- 1: q vs probe order -----------------------------------------------
nexttile;
scatter(x(g1pass), q(g1pass), 12, [0.2 0.6 0.2], 'filled'); hold on;
scatter(x(~g1pass), q(~g1pass), 12, [0.8 0.2 0.2], 'filled');
yline(1, '--', 'qmax default context', 'Color',[0.5 0.5 0.5]);
xlabel('probe order'); ylabel('q = \lambda_{d+1}/\lambda_d');
title('leakage ratio vs probe order'); legend({'gate1 pass','gate1 fail'},'Location','best');
box on;

% --- 2: d_hat vs probe order --------------------------------------------
nexttile;
scatter(x, dh, 12, C6, 'filled'); colorbar; colormap(gca, parula);
xlabel('probe order'); ylabel('d\_hat'); title('estimated local dimension (colour = C6)');
box on;

% --- 3: q vs covariate ---------------------------------------------------
nexttile;
if ~isempty(covariate)
    scatter(covariate(g1pass), q(g1pass), 12, [0.2 0.6 0.2], 'filled'); hold on;
    scatter(covariate(~g1pass), q(~g1pass), 12, [0.8 0.2 0.2], 'filled');
    xlabel(covariate_name); ylabel('q');
    title(sprintf('leakage ratio vs %s', covariate_name));
    box on;
else
    axis off; text(0.1,0.5,'no covariate supplied','Units','normalized');
end

% --- 4: C6 vs covariate, split by gate1 ----------------------------------
nexttile;
if ~isempty(covariate)
    scatter(covariate(g1pass), C6(g1pass), 12, [0.2 0.6 0.2], 'filled'); hold on;
    scatter(covariate(~g1pass), C6(~g1pass), 12, [0.8 0.2 0.2], 'filled');
    set(gca,'YScale','log');
    xlabel(covariate_name); ylabel('C6 (log scale)');
    title('conditioning vs covariate, by gate1 pass/fail');
    legend({'gate1 pass','gate1 fail'},'Location','best'); box on;
else
    axis off; text(0.1,0.5,'no covariate supplied','Units','normalized');
end

% --- 5: histogram of d_hat -----------------------------------------------
nexttile;
edges = min(dh)-0.5 : 1 : max(dh)+0.5;
histogram(dh, edges, 'FaceColor',[0.3 0.3 0.6]);
xlabel('d\_hat'); ylabel('count'); title('distribution of estimated dimension');
vals = unique(dh(~isnan(dh)));
for v = vals(:)'
    cnt = sum(dh==v);
    text(v, cnt, sprintf(' %d', cnt), 'VerticalAlignment','bottom');
end
box on;

% --- 6: full gate matrix --------------------------------------------------
nexttile;
gm = double(R.gate ~= 0);
gm(:,3) = NaN;                        % gate 2 not evaluable
imagesc(gm', [0 1]); colormap(gca, [0.8 0.2 0.2; 0.9 0.9 0.9; 0.2 0.6 0.2]);
set(gca,'YTick',1:7,'YTickLabel',{'0 chart','1 dim','2 est(N/A)','3 anchor','4 gain','5 book','6 C6'});
xlabel('probe order'); title('per-probe, per-gate pass (green) / fail (red)');

fprintf('--- diagnostic summary ---\n');
fprintf('  gate1 pass rate           : %.1f%%\n', 100*mean(g1pass));
if ~isempty(covariate) && any(~g1pass) && any(g1pass)
    [~, p] = ttest2(covariate(g1pass), covariate(~g1pass));
    fprintf('  %s differs pass vs fail : t-test p = %.4f (descriptive only, not a gate)\n', ...
        covariate_name, p);
end
qfin = q(isfinite(q));
fprintf('  q distribution            : median %.2f, IQR [%.2f, %.2f]\n', ...
    median(qfin), quantile(qfin,0.25), quantile(qfin,0.75));
end
