%% plotResultFigures.m
% Figures 1-4 of the results section. Nothing is recomputed: every panel reads a cache.
%   compareRateWindowsByCellType    -> compareRateWindowsByCellType.mat    (T)
%   compareRippleWindowsByCellType  -> compareRippleWindowsByCellType.mat  (T)
%   compareReactivation             -> compareReactivation.mat             (T, B, TC, LB, GN, gainS)
%
% FIGURE 1  firing rate by cell type      FIGURE 3  reactivation
% FIGURE 2  ripples and unit engagement   FIGURE 4  gain decomposition
%
% Conventions used in every panel
%   one dot = one session-window (rate, ripple) or one session (reactivation, gain)
%   a horizontal bar = the group mean (geometric for rates and gains, which are drawn on a log axis)
%   whiskers = SEM over those dots
%   paired panels: one line per animal, its mean before and after that arm's infusion
%   DiD panels: one dot per animal, (ASV change from PreASV) - (SAL change from PreSAL); p = exact
%               sign-flip over those animals, floor 2/2^n
% Set FIG to draw a subset, e.g. FIG = 3.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
if exist('SHARED', 'var'), CACHE_ROOT = SHARED.cacheRoot; end   % set by runAllAnalyses.m
REACT_MAT  = 'compareReactivation.mat';      % e.g. compareReactivation_spk50.mat for a re-thresholded run
FIG        = [1 2 3 4];
SAVE_PDF   = true;
GAIN_CELLS = 'PYR';                          % the stratum the cached gain analysis used

COL = struct('Control', [0.42 0.42 0.42], 'SAL', [0.16 0.44 0.78], 'ASV', [0.84 0.30 0.12]);
WIN = {'Control','PreASV','ASVshort','ASVlong','PreSAL','SALshort','SALlong'};
WCOL = [COL.Control; COL.ASV*0.45 + 0.55; COL.ASV; COL.ASV*0.7 + 0.3*[1 1 1]; ...
        COL.SAL*0.45 + 0.55; COL.SAL; COL.SAL*0.7 + 0.3*[1 1 1]];

R = load(fullfile(CACHE_ROOT, 'compareRateWindowsByCellType.mat'), 'T');      RT = R.T;
P = load(fullfile(CACHE_ROOT, 'compareRippleWindowsByCellType.mat'), 'T');    PT = P.T;
A = load(fullfile(CACHE_ROOT, REACT_MAT), 'T', 'B', 'TC', 'LB', 'GN', 'gainS');

%% ======================= FIGURE 1: firing rate by cell type =======================
if ismember(1, FIG)
    f1 = figure('Color', 'w', 'Position', [40 40 1250 700], 'Name', 'Fig 1 rate');
    t1 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % (a) baseline: control vs TransPlasticity (every session, before its own infusion), by stratum
    nexttile; hold on
    st = {'PYR','PYRdeep','PYRsup','INT'};
    for s = 1:numel(st)
        v1 = RT.(['rate_' st{s}])(strcmp(RT.Group, 'Control'));
        v2 = RT.(['rate_' st{s}])(strcmp(RT.Group, 'BaseAll'));   % cohort baseline: every session, pre-infusion
        local_strip(s - 0.17, v1(isfinite(v1)), COL.Control, true);
        local_strip(s + 0.17, v2(isfinite(v2)), [0.1 0.1 0.1], true);
    end
    set(gca, 'YScale', 'log'); xlim([0.5 numel(st) + 0.5]); xticks(1:numel(st)); xticklabels(st);
    ylabel('firing rate (Hz)'); title('a  baseline: control (grey) vs TransPlasticity (black)', 'FontWeight', 'normal');

    % (b, c) acute change within animals
    for s = {'PYR','PYRdeep'}
        nexttile; hold on
        local_pairedPanel(RT, ['rate_' s{1}], {'PreSAL','SALshort'; 'PreASV','ASVshort'}, COL, true);
        ylabel(sprintf('%s rate (Hz)', s{1}));
        title(sprintf('%s  2-8 h after infusion, %s', char('a' + find(strcmp(s{1}, {'PYR','PYRdeep'}))), s{1}), ...
            'FontWeight', 'normal');
    end

    % (d) difference in differences, acute
    nexttile; hold on
    local_didPanel(RT, strcat('rate_', {'PYR','PYRdeep','PYRsup','INT'}), {'PYR','PYRdeep','PYRsup','INT'}, ...
        'PreASV', 'ASVshort', 'PreSAL', 'SALshort', true);
    ylabel('\Deltalog_{10} rate: ASV - SAL'); title('d  acute drug specificity (one dot = one animal)', 'FontWeight', 'normal');

    % (e) persistent change within animals
    nexttile; hold on
    local_pairedPanel(RT, 'rate_PYR', {'PreSAL','SALlong'; 'PreASV','ASVlong'}, COL, true);
    ylabel('PYR rate (Hz)'); title('e  >= 24 h after infusion', 'FontWeight', 'normal');

    % (f) composition
    nexttile; hold on
    local_windowStrip(RT, 'fracINT', WIN, WCOL, false);
    ylabel('interneuron fraction of units'); title('f  recorded composition', 'FontWeight', 'normal');

    title(t1, 'Figure 1  Firing rate by cell type: saline raises it, ASV does not');
    if SAVE_PDF, exportgraphics(f1, fullfile(CACHE_ROOT, 'figure1_rate.pdf'), 'ContentType', 'vector'); end
end

%% ================= FIGURE 2: ripple properties and unit engagement =================
if ismember(2, FIG)
    f2 = figure('Color', 'w', 'Position', [40 40 1250 700], 'Name', 'Fig 2 ripples');
    t2 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % (a, b) the events themselves
    nexttile; hold on
    local_windowStrip(PT, 'rate', WIN, WCOL, true); ylabel('ripples / h');
    title('a  ripple occurrence', 'FontWeight', 'normal');
    nexttile; hold on
    local_windowStrip(PT, 'freq', WIN, WCOL, false); ylabel('peak frequency (Hz)');
    title('b  ripple frequency', 'FontWeight', 'normal');

    % (c) engagement at baseline, by stratum
    nexttile; hold on
    st = {'PYR','PYRdeep','PYRsup','INT'};
    for s = 1:numel(st)
        v1 = PT.(['gainW_' st{s}])(strcmp(PT.Group, 'Control'));
        v2 = PT.(['gainW_' st{s}])(strcmp(PT.Group, 'BaseAll'));
        local_strip(s - 0.17, v1(isfinite(v1)), COL.Control, true);
        local_strip(s + 0.17, v2(isfinite(v2)), [0.1 0.1 0.1], true);
    end
    set(gca, 'YScale', 'log'); yline(1, 'k:'); xlim([0.5 numel(st) + 0.5]); xticks(1:numel(st)); xticklabels(st);
    ylabel('gain (in-ripple / out-of-ripple)'); title('c  ripple modulation at baseline', 'FontWeight', 'normal');

    % (d) acute recruitment within animals
    nexttile; hold on
    local_pairedPanel(PT, 'frac_PYRdeep', {'PreSAL','SALshort'; 'PreASV','ASVshort'}, COL, false);
    ylabel('deep PYR recruitment'); title('d  2-8 h after infusion', 'FontWeight', 'normal');

    % (e) the key contrast: which engagement measures separate the arms
    nexttile; hold on
    local_didPanel(PT, {'frac_PYRdeep','xfrac_PYRdeep','rRip_PYRdeep','gainW_PYRdeep'}, ...
        {'recruit','excess','in-rip rate','gain_{window}'}, 'PreASV', 'ASVshort', 'PreSAL', 'SALshort', ...
        [false false true true]);
    ylabel('acute change: ASV - SAL'); title('e  recruitment moves, window-matched gain does not', 'FontWeight', 'normal');

    % (f) phase locking
    nexttile; hold on
    for s = 1:numel(st)
        for g = 1:3
            gn = {'Control','PreASV','ASVlong'};
            v = PT.(['ppc_' st{s}])(strcmp(PT.Group, gn{g}));
            cc = {COL.Control, [0.1 0.1 0.1], COL.ASV};
            local_strip(s + (g - 2)*0.24, v(isfinite(v)), cc{g}, false);
        end
    end
    xlim([0.5 numel(st) + 0.5]); xticks(1:numel(st)); xticklabels(st);
    ylabel('unit PPC'); title('f  phase locking: control, PreASV, ASVlong', 'FontWeight', 'normal');

    title(t2, 'Figure 2  Ripples are unchanged; recruitment follows firing rate');
    if SAVE_PDF, exportgraphics(f2, fullfile(CACHE_ROOT, 'figure2_ripples.pdf'), 'ContentType', 'vector'); end
end

%% ========================= FIGURE 3: reactivation =========================
if ismember(3, FIG)
    GR = {'Control','SAL','ASV'};
    f3 = figure('Color', 'w', 'Position', [40 40 1250 700], 'Name', 'Fig 3 reactivation');
    t3 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % (a) EV against REV, one dot per session
    nexttile; hold on
    for g = GR
        k = strcmp(A.T.Group, g{1}) & strcmp(A.T.Context, 'ConB1');
        plot(A.T.REV(k), A.T.EV(k), 'o', 'Color', COL.(g{1}), 'MarkerFaceColor', COL.(g{1}), 'MarkerSize', 6);
    end
    lim = [0 max([A.T.EV; A.T.REV], [], 'omitnan')*1.05];
    plot(lim, lim, 'k:'); xlim(lim); ylim(lim); axis square
    xlabel('REV (template vs PRE)'); ylabel('EV (template vs POST)');
    legend(GR, 'Location', 'northwest', 'Box', 'off');
    title('a  novel context, first 30 min', 'FontWeight', 'normal');

    % (b) EV - REV by context, early
    nexttile; hold on
    ctx = {'ConA1','ConB1'};
    for c = 1:2
        for g = 1:3
            k = strcmp(A.T.Group, GR{g}) & strcmp(A.T.Context, ctx{c});
            local_strip(c + (g - 2)*0.24, A.T.dEV(k), COL.(GR{g}), false);
        end
    end
    yline(0, 'k:'); xlim([0.5 2.5]); xticks(1:2); xticklabels({'familiar (A1)','novel (B1)'});
    ylabel('EV - REV'); title('b  early (first 30 min after each exposure)', 'FontWeight', 'normal');

    % (c) late block
    nexttile; hold on
    for c = 1:2
        for g = 1:3
            k = strcmp(A.LB.Group, GR{g}) & strcmp(A.LB.Template, ctx{c});
            local_strip(c + (g - 2)*0.24, A.LB.dEV(k), COL.(GR{g}), false);
        end
    end
    yline(0, 'k:'); xlim([0.5 2.5]); xticks(1:2); xticklabels({'familiar (A1)','novel (B1)'});
    ylabel('EV - REV'); title('c  late (2-4 h after the novel exposure)', 'FontWeight', 'normal');

    % (d, e) time course
    for c = 1:2
        nexttile; hold on
        nb = max(A.TC.Bin);
        for g = GR
            mu = nan(nb, 1); se = mu;
            for b = 1:nb
                k = strcmp(A.TC.Template, ctx{c}) & strcmp(A.TC.Group, g{1}) & A.TC.Bin == b;
                if any(k), mu(b) = mean(A.TC.dEV(k)); se(b) = std(A.TC.dEV(k))/sqrt(sum(k)); end
            end
            errorbar((1:nb)' - 0.5, mu, se, '-o', 'Color', COL.(g{1}), 'MarkerFaceColor', COL.(g{1}), ...
                'LineWidth', 1.5, 'CapSize', 0);
        end
        yline(0, 'k:'); xlabel('h after the novel exposure ends'); ylabel('EV - REV');
        title(sprintf('%s  %s template, time course', char('c' + c), ctx{c}), 'FontWeight', 'normal');
    end

    % (f) novel minus familiar, within session
    nexttile; hold on
    for g = 1:3
        k = strcmp(A.B.Group, GR{g});
        local_strip(g, A.B.d(k), COL.(GR{g}), false);
    end
    yline(0, 'k:'); xlim([0.5 3.5]); xticks(1:3); xticklabels(GR);
    ylabel('\Delta(EV - REV), B1 - A1'); title('f  novel minus familiar, same session', 'FontWeight', 'normal');

    title(t3, 'Figure 3  Reactivation of the novel context survives ASV');
    if SAVE_PDF, exportgraphics(f3, fullfile(CACHE_ROOT, 'figure3_reactivation.pdf'), 'ContentType', 'vector'); end
end

%% =================== FIGURE 4: gain decomposition ===================
if ismember(4, FIG)
    GR = {'SAL','ASV'};          % ASV vs saline only: same animals, contexts and recordings
    GRD = {'Control','SAL','ASV'};  % the drift check alone also shows the never-infused dataset
    f4 = figure('Color', 'w', 'Position', [40 40 1250 700], 'Name', 'Fig 4 gain');
    t4 = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    % (a) the effect itself: change in ripple firing against novel-context activity
    nexttile; hold on
    ed = [-Inf -1 -0.33 0.33 1 Inf];
    for g = GR
        kk = find(strcmp(A.gainS(:,1), g{1}));
        X = []; Y = []; W = [];
        for q = kk(:)'
            U = A.gainS{q,4};
            ok = U.(GAIN_CELLS) & U.preB_n >= 10 & U.early_n >= 10 & isfinite(U.rB);
            if sum(ok) < 8, continue; end
            z = zscore(log10(U.rB(ok) + 0.01));
            X = [X; z]; %#ok<AGROW>
            Y = [Y; U.early_lrip(ok) - U.preB_lrip(ok)]; %#ok<AGROW>
            W = [W; repmat(1/sum(ok), sum(ok), 1)]; %#ok<AGROW>
        end
        if isempty(X), continue; end
        [~, ~, bi] = histcounts(X, ed); v = bi > 0;
        mu = accumarray(bi(v), Y(v).*W(v), [5 1]) ./ accumarray(bi(v), W(v), [5 1]);
        plot(1:5, mu, '-o', 'Color', COL.(g{1}), 'MarkerFaceColor', COL.(g{1}), 'LineWidth', 1.5);
    end
    yline(0, 'k:'); xticks(1:5); xticklabels({'<-1','-1..-0.3','mid','0.3..1','>1'});
    xlabel('activity in the novel context (z, log rate)'); ylabel('\Deltalog_{10} in-ripple rate');
    legend(GR, 'Location', 'northwest', 'Box', 'off');
    title('a  cells active in the novel context gain the most', 'FontWeight', 'normal');

    % (b, c) the two components, per session
    for v = {'a','bB'}
        nexttile; hold on
        for y = 1:2
            yn = {'lrip','lout'};
            for g = 1:numel(GR)
                k = strcmp(A.GN.Group, GR{g}) & strcmp(A.GN.Base, 'preB') & strcmp(A.GN.Post, 'early') & ...
                    strcmp(A.GN.Y, yn{y});
                local_strip(y + (g - 1.5)*0.3, A.GN.(v{1})(k), COL.(GR{g}), false);
            end
        end
        yline(0, 'k:'); xlim([0.5 2.5]); xticks(1:2); xticklabels({'in ripples','out of ripples'});
        ylabel(sprintf('%s (log_{10} per SD)', v{1}));
        if strcmp(v{1}, 'a')
            title('b  shared component a: the saline rise ASV prevents', 'FontWeight', 'normal');
        else
            title('c  experience-tied component b_B: unchanged by ASV', 'FontWeight', 'normal');
        end
    end

    % (d) drift control: the same fit before the novel experience happened
    nexttile; hold on
    for g = 1:numel(GRD)
        tr = A.GN(strcmp(A.GN.Group, GRD{g}) & strcmp(A.GN.Base, 'preA') & strcmp(A.GN.Post, 'early') & ...
            strcmp(A.GN.Y, 'lrip'), {'Session','bB'});
        tp = A.GN(strcmp(A.GN.Group, GRD{g}) & strcmp(A.GN.Base, 'preA') & strcmp(A.GN.Post, 'preB') & ...
            strcmp(A.GN.Y, 'lrip'), {'Session','bB'});
        tr.Properties.VariableNames = {'Session','real'}; tp.Properties.VariableNames = {'Session','plc'};
        J = innerjoin(tr, tp, 'Keys', 'Session');
        if isempty(J), continue; end
        x = [g - 0.18, g + 0.18];
        plot(x, [J.plc, J.real]', '-', 'Color', COL.(GRD{g})*0.5 + 0.5);
        local_strip(x(1), J.plc, COL.(GRD{g}), false); local_strip(x(2), J.real, COL.(GRD{g}), false);
    end
    yline(0, 'k:'); xlim([0.5 numel(GRD) + 0.5]); xticks(1:numel(GRD)); xticklabels(GRD);
    ylabel('b_B (log_{10} per SD)');
    title('d  drift check: placebo (left) vs real (right); control = measure check only', 'FontWeight', 'normal');

    % (e) does gain-coded reactivation track pattern reactivation?
    nexttile; hold on
    tg = A.GN(strcmp(A.GN.Base, 'preB') & strcmp(A.GN.Post, 'early') & strcmp(A.GN.Y, 'lrip'), :);
    rb = A.T(strcmp(A.T.Context, 'ConB1'), {'Session','dEV'});
    J = innerjoin(tg, rb, 'Keys', 'Session');
    for g = GR
        k = strcmp(J.Group, g{1});
        plot(J.bB(k), J.dEV(k), 'o', 'Color', COL.(g{1}), 'MarkerFaceColor', COL.(g{1}), 'MarkerSize', 6);
    end
    xline(0, 'k:'); yline(0, 'k:');
    xlabel('b_B, in-ripple rate'); ylabel('EV - REV, novel context');
    title('e  the two measures of reactivation, per session', 'FontWeight', 'normal');

    % (f) participation beyond chance
    nexttile; hold on
    for g = 1:numel(GR)
        k = strcmp(A.GN.Group, GR{g}) & strcmp(A.GN.Base, 'preB') & strcmp(A.GN.Post, 'early') & ...
            strcmp(A.GN.Y, 'xp');
        local_strip(g, A.GN.bB(k), COL.(GR{g}), false);
    end
    yline(0, 'k:'); xlim([0.5 numel(GR) + 0.5]); xticks(1:numel(GR)); xticklabels(GR);
    ylabel('b_B, participation beyond chance'); title('f  the same gain in ripple participation', 'FontWeight', 'normal');

    title(t4, 'Figure 4  ASV removes the shared rate change and spares the experience-tied gain');
    if SAVE_PDF, exportgraphics(f4, fullfile(CACHE_ROOT, 'figure4_gain.pdf'), 'ContentType', 'vector'); end
end


%% =========================================================================
function local_strip(x, v, col, useLog)
% jittered dots + mean (geometric when useLog) and SEM
if nargin < 4, useLog = false; end
v = v(isfinite(v)); if isempty(v), return; end
pale = col*0.45 + 0.55;                              % no alpha: some renderers drop RGBA in vector PDFs
jit = 0.07 * (mod((1:numel(v))' * 0.618, 1) - 0.5);  % deterministic spread, no RNG side effect
plot(x + jit, v, 'o', 'Color', pale, 'MarkerFaceColor', pale, 'MarkerSize', 4);
if useLog
    l = log10(v); m = 10^mean(l); lo = 10^(mean(l) - std(l)/sqrt(numel(l))); hi = 10^(mean(l) + std(l)/sqrt(numel(l)));
else
    m = mean(v); se = std(v)/sqrt(numel(v)); lo = m - se; hi = m + se;
end
plot(x + [-0.14 0.14], [m m], '-', 'Color', col, 'LineWidth', 2.2);
plot([x x], [lo hi], '-', 'Color', col, 'LineWidth', 1.2);
end


function local_windowStrip(T, nm, wins, wcol, useLog)
% one column per window, in the order given
for k = 1:numel(wins)
    v = T.(nm)(strcmp(T.Group, wins{k}));
    local_strip(k, v, wcol(k,:), useLog);
end
xlim([0.5 numel(wins) + 0.5]); xticks(1:numel(wins)); xticklabels(wins); xtickangle(30);
if useLog, set(gca, 'YScale', 'log'); end
end


function local_pairedPanel(T, nm, arms, col, useLog)
% arms = {preWindow, postWindow; ...}; one line per animal, its mean in each window
for a = 1:size(arms, 1)
    [v0, id0] = local_animalMean(T, nm, arms{a,1});
    [v1, id1] = local_animalMean(T, nm, arms{a,2});
    [ids, i0, i1] = intersect(id0, id1);
    if isempty(ids), continue; end
    cc = col.SAL; if contains(arms{a,2}, 'ASV'), cc = col.ASV; end
    x = [2*a - 1, 2*a];
    plot(x, [v0(i0), v1(i1)]', '-', 'Color', cc*0.5 + 0.5, 'LineWidth', 1);
    local_strip(x(1), v0(i0), cc, useLog); local_strip(x(2), v1(i1), cc, useLog);
end
xlim([0.4 2*size(arms, 1) + 0.6]); xticks(1:2*size(arms, 1)); xticklabels(arms');
xtickangle(30); if useLog, set(gca, 'YScale', 'log'); end
end


function local_didPanel(T, vars, labels, preA, postA, preS, postS, useLog)
% one dot per animal: (ASV change from preA) - (SAL change from preS); p = exact sign-flip
if numel(useLog) == 1, useLog = repmat(useLog, 1, numel(vars)); end
for q = 1:numel(vars)
    d = local_did(T, vars{q}, preA, postA, preS, postS, useLog(q));
    if isempty(d), continue; end
    local_strip(q, d, [0.1 0.1 0.1], false);
    text(q, max(d, [], 'omitnan'), sprintf(' p=%.3f', local_signflip(d)), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 8);
end
yline(0, 'k:'); xlim([0.5 numel(vars) + 0.5]); xticks(1:numel(vars)); xticklabels(labels); xtickangle(20);
end


function d = local_did(T, nm, preA, postA, preS, postS, useLog)
f = @(w) local_animalMean(T, nm, w);
[a0, i0] = f(preA); [a1, i1] = f(postA); [s0, j0] = f(preS); [s1, j1] = f(postS);
ids = intersect(intersect(i0, i1), intersect(j0, j1));
d = nan(numel(ids), 1);
for k = 1:numel(ids)
    v = [a0(strcmp(i0, ids{k})), a1(strcmp(i1, ids{k})), s0(strcmp(j0, ids{k})), s1(strcmp(j1, ids{k}))];
    if useLog, v = log10(v); end
    d(k) = (v(2) - v(1)) - (v(4) - v(3));
end
d = d(isfinite(d));
end


function [v, ids] = local_animalMean(T, nm, win)
k = strcmp(T.Group, win) & isfinite(T.(nm));
ids = unique(T.Animal(k));
v = cellfun(@(a) mean(T.(nm)(k & strcmp(T.Animal, a))), ids);
end


function p = local_signflip(d)
d = d(isfinite(d)); n = numel(d); p = NaN;
if n < 2, return; end
obs = abs(mean(d));
if n <= 14
    Sg = 2*(dec2bin(0:2^n-1, n) == '1') - 1; p = mean(abs(Sg*d(:)/n) >= obs - 1e-12);
else
    nul = arrayfun(@(k) mean((2*(rand(n,1) > 0.5) - 1) .* d(:)), 1:1e4); p = (1 + sum(abs(nul) >= obs))/(1 + 1e4);
end
end
