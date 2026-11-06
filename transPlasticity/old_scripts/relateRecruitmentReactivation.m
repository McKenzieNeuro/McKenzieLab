%% relateRecruitmentReactivation.m
% Does the change in deep-pyramidal ripple recruitment after an infusion track
% reactivation of the context coactivity in the same session?
% Run compareRippleWindowsByCellType and compareReactivation first (their caches).
%
% PER INFUSION SESSION (first SAL before ASV, first ASV)
%   dRec_<m> = metric m in the 2-8 h window after the infusion minus the same session's
%              pre-infusion part (SALshort - PreSAL, ASVshort - PreASV; same ripple file)
%              m = frac (recruitment), xfrac (beyond chance at the window's own rate),
%              gainW (log10), rRip (log10), all for PYRdeep
%   dEV        reactivation of ConB1 (novel, encoded on drug) and ConA1 (familiar, pre-drug):
%              early = first POST_WIN after the exposure, late = LATE_START-TC_MAX after ConB1
% TESTS (sessions; few animals, so read as descriptive)
%   Pearson and Spearman over sessions, p by permuting dRec (1e4)
%   within each group, and pooled with group as a covariate: LME dEV ~ dRec + Group + (1|Animal)

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
METRICS = {'frac','xfrac','gainW','rRip'};
STRAT = 'PYRdeep';
NPERM = 1e4;

R = load(fullfile(CACHE_ROOT, 'compareRippleWindowsByCellType.mat'), 'T'); W = R.T;
A = load(fullfile(CACHE_ROOT, 'compareReactivation.mat'), 'T', 'LB'); RT = A.T; LB = A.LB;

pairs = {'SAL', 'PreSAL', 'SALshort', 'SALpre';     % arm, pre window, post window, reactivation group
         'ASV', 'PreASV', 'ASVshort', 'ASV1'};
rows = {};
for a = 1:size(pairs, 1)
    sess = unique(W.Session(strcmp(W.Group, pairs{a,3})));
    for s = sess(:)'
        kPre = strcmp(W.Session, s{1}) & strcmp(W.Group, pairs{a,2});
        kPost = strcmp(W.Session, s{1}) & strcmp(W.Group, pairs{a,3});
        if ~any(kPre) || ~any(kPost), continue; end
        d = nan(1, numel(METRICS));
        for m = 1:numel(METRICS)
            nm = [METRICS{m} '_' STRAT];
            x0 = W.(nm)(find(kPre, 1)); x1 = W.(nm)(find(kPost, 1));
            if any(strcmp(METRICS{m}, {'gainW','rRip'})), x0 = log10(x0); x1 = log10(x1); end
            d(m) = x1 - x0;
        end
        % reactivation of the same session
        r = nan(1, 4);                              % B early, B late, A early, A late
        kB = strcmp(RT.Session, s{1}) & strcmp(RT.Context, 'ConB1') & strcmp(RT.Group, pairs{a,4});
        kA = strcmp(RT.Session, s{1}) & strcmp(RT.Context, 'ConA1') & strcmp(RT.Group, pairs{a,4});
        lB = strcmp(LB.Session, s{1}) & strcmp(LB.Template, 'ConB1');
        lA = strcmp(LB.Session, s{1}) & strcmp(LB.Template, 'ConA1');
        if any(kB), r(1) = RT.dEV(find(kB, 1)); end
        if any(lB), r(2) = LB.dEV(find(lB, 1)); end
        if any(kA), r(3) = RT.dEV(find(kA, 1)); end
        if any(lA), r(4) = LB.dEV(find(lA, 1)); end
        rows(end+1,:) = [{pairs{a,1}, W.Animal{find(kPost, 1)}, s{1}}, num2cell(d), num2cell(r)]; %#ok<AGROW>
    end
end
D = cell2table(rows, 'VariableNames', [{'Arm','Animal','Session'}, strcat('dRec_', METRICS), ...
    {'dEV_B_early','dEV_B_late','dEV_A_early','dEV_A_late'}]);
fprintf('%d infusion sessions with recruitment change: SAL %d, ASV %d; with ConB1 reactivation: %d\n', height(D), ...
    sum(strcmp(D.Arm, 'SAL')), sum(strcmp(D.Arm, 'ASV')), sum(isfinite(D.dEV_B_early)));

%% ---- correlations --------------------------------------------------------------------------
ys = {'dEV_B_early','dEV_B_late','dEV_A_early','dEV_A_late'};
fprintf('\n%-12s %-12s | %-34s | %-26s | %-26s | %s\n', 'recruitment', 'reactivation', 'pooled: r / rho (perm p), n', ...
    'SAL: r (p), n', 'ASV: r (p), n', 'LME slope (p), Group + (1|Animal)');
for m = 1:numel(METRICS)
    xn = ['dRec_' METRICS{m}];
    for j = 1:numel(ys)
        k = isfinite(D.(xn)) & isfinite(D.(ys{j}));
        if sum(k) < 4, continue; end
        [r, pr] = local_permCorr(D.(xn)(k), D.(ys{j})(k), 'Pearson', NPERM);
        [rho, ps] = local_permCorr(D.(xn)(k), D.(ys{j})(k), 'Spearman', NPERM);
        c = cell(1, 2);
        for a = 1:2
            kk = k & strcmp(D.Arm, pairs{a,1});
            c{a} = sprintf('n=%d', sum(kk));
            if sum(kk) >= 4
                [ra, pa] = local_permCorr(D.(xn)(kk), D.(ys{j})(kk), 'Pearson', NPERM);
                c{a} = sprintf('%+.2f (%.3f), n=%d', ra, pa, sum(kk));
            end
        end
        [bl, pl] = local_lme(D(k,:), xn, ys{j});
        fprintf('%-12s %-12s | %+.2f (%.3f) / %+.2f (%.3f), n=%-2d | %-26s | %-26s | %+.3f (%.3f)\n', METRICS{m}, ys{j}, ...
            r, pr, rho, ps, sum(k), c{:}, bl, pl);
    end
end
fprintf('recruitment change = 2-8 h after the infusion minus the same session before it (%s); gainW and rRip in log10.\n', STRAT);
save(fullfile(CACHE_ROOT, 'relateRecruitmentReactivation.mat'), 'D');

%% ---- figure ----------------------------------------------------------------------------------
col = containers.Map({'SAL','ASV'}, {[0.2 0.6 1], [0.85 0.2 0.1]});
fig = figure('Color','w','Position',[60 60 1400 640]);
tl = tiledlayout(2, numel(ys), 'TileSpacing','compact');
xm = {'frac','xfrac'};
for rr = 1:2
    xn = ['dRec_' xm{rr}];
    for j = 1:numel(ys)
        nexttile; hold on
        for a = {'SAL','ASV'}
            k = strcmp(D.Arm, a{1}) & isfinite(D.(xn)) & isfinite(D.(ys{j}));
            plot(D.(xn)(k), D.(ys{j})(k), 'o', 'Color', col(a{1}), 'MarkerFaceColor', col(a{1}), 'MarkerSize', 6);
        end
        k = isfinite(D.(xn)) & isfinite(D.(ys{j}));
        if sum(k) >= 3
            b = polyfit(D.(xn)(k), D.(ys{j})(k), 1); xx = [min(D.(xn)(k)) max(D.(xn)(k))];
            plot(xx, polyval(b, xx), 'k-', 'LineWidth', 1);
        end
        xline(0, 'k:'); yline(0, 'k:');
        xlabel(sprintf('\\Delta %s %s (2-8 h - pre)', xm{rr}, STRAT), 'Interpreter','tex');
        ylabel(strrep(ys{j}, '_', ' '));
        if rr == 1 && j == 1, legend({'SAL','ASV'}, 'Location','best', 'Box','off'); end
    end
end
title(tl, 'Deep pyramidal recruitment change after infusion vs reactivation in the same session');
exportgraphics(fig, fullfile(CACHE_ROOT, 'relateRecruitmentReactivation.pdf'), 'ContentType','vector');


%% =========================================================================
function [r, p] = local_permCorr(x, y, type, nperm)
r = corr(x(:), y(:), 'Type', type);
nul = arrayfun(@(k) corr(x(randperm(numel(x))), y(:), 'Type', type), 1:nperm);
p = (1 + sum(abs(nul) >= abs(r) - 1e-12)) / (1 + nperm);
end


function [b, p] = local_lme(D, xn, yn)
b = NaN; p = NaN;
if height(D) < 6 || numel(unique(D.Arm)) < 2, return; end
t = D; t.x = t.(xn); t.y = t.(yn); t.Arm = categorical(t.Arm); t.Animal = categorical(t.Animal);
try
    lme = fitlme(t, 'y ~ x + Arm + (1|Animal)');
    [~, ~, C] = fixedEffects(lme, 'DFMethod', 'satterthwaite');
    if isa(C, 'dataset'), C = dataset2table(C); end
    k = find(string(C.Name) == "x", 1);
    if ~isempty(k), b = C.Estimate(k); p = C.pValue(k); end
catch ME
    warning('rel:lme', '%s: %s', yn, ME.message);
end
end
