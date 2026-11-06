%% rateDivergenceTiming.m
% WHEN does the ASV-versus-saline firing difference appear: at the infusion, or at the novel
% exposure? Both arms receive an infusion and both see ConB1, so the arm contrast is already
% drug-specific; what is in question is which event the divergence follows.
%
%   HC  ConA1  HC  Infusion  HC  ConB1  HC  ConA2  HC  ConB2  HC
%
% MEASURE. Out-of-ripple firing in NREM, so neither ripple occurrence nor time awake can carry
% the effect. Per neuron and per time bin, the rate is spikes in (bin AND NREM AND outside any
% ripple) over the duration of that intersection. Each neuron is expressed as
%   d = log10(rate + 0.01) - log10(baseRate + 0.01),
% its own pre-infusion NREM out-of-ripple rate as baseline, so the day-order imbalance (ASV days
% mostly later) cannot enter: every trajectory is a within-session, within-neuron change.
% A session's value in a bin is the median d over the neurons of a stratum (PYR, PYRdeep, INT).
%
% TWO CLOCKS. Each bin is placed on both, from the same spikes:
%   'inf'  time since Infusion_start
%   'B1'   time since ConB1 ends (the first NREM after the novel exposure is bin 1)
% On the B1 clock the negative bins are the post-infusion, pre-novelty window: that is where the
% two hypotheses differ. The exposures themselves are almost all wake, so bins covering them fail
% the coverage floor and drop out on their own.
%
% COVERAGE. A bin is scored only with MIN_COV s of NREM out-of-ripple time, and the coverage is
% reported per bin per arm, so a difference that is really a difference in how much the arms slept
% is visible rather than hidden. Two sensitivities: NREM bouts trimmed to their first TRIM_S s
% (equalises bout composition), and bouts shorter than MIN_BOUT s dropped.
%
% WHEN. The arm difference is scanned bin by bin outward from each anchor. The first bin that
% clears a calibrated threshold is the divergence. Calibration is single-step max-t: the sign of
% each animal's WHOLE trajectory is flipped (bins within an animal are correlated, so flipping
% bins independently would be wrong), the per-bin t is recomputed, and the calibrated p of a bin
% is the fraction of flips whose largest |t| over the scanned bins reaches the observed |t|.
% Bins before the infusion are a null check on the pipeline: a divergence there is baseline
% imbalance, not an effect.
%
% TWO ANCHORS AT ONCE. The infusion-to-ConB1 gap is not identical across sessions, and that jitter
% separates the clocks: a bin at a given time since infusion is pre-novelty in a long-gap session
% and post-novelty in a short-gap one. The model
%   d ~ Group * PostB + tSinceInfusion + (1|Animal),      PostB = bin falls after ConB1 ends
% asks whether the drug difference is larger after the novel exposure at matched time since the
% infusion. The Group:PostB interaction is the test the timing question reduces to; the gap's
% spread is printed first, because it is what makes the term identifiable.
%
% Output: rateDivergenceTiming.mat (S per session, Res) and rateDivergenceTiming.pdf.

ROOT       = 'R:\TransPlasticity';
CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
if exist('SHARED', 'var'), CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root; end
O = struct();
O.binH      = 0.5;      % h per bin
O.preH      = 2;        % h before the anchor
O.postInfH  = 8;        % h after the infusion
O.postB1H   = 6;        % h after ConB1 ends
O.minCov    = 120;      % s of NREM out-of-ripple time a bin needs
O.minBase   = 300;      % s of baseline NREM out-of-ripple time a session needs
O.minBaseSpk = 20;      % spikes in that baseline a neuron needs
O.minUnits  = 5;        % scored neurons a stratum needs in a bin
O.trimS     = 120;      % sensitivity: use only the first trimS s of each NREM bout
O.minBoutS  = 60;       % sensitivity: drop NREM bouts shorter than this
O.nFlipMax  = 4096;     % exact sign flips up to 12 animals, else this many random
STRAT = {'PYR','PYRdeep','PYRsup','INT'};
CLK   = {'inf','B1'};
VAR   = {'full','trim','longbout'};
SAVE_OUT = true;
rng(1);

eInf = (-O.preH:O.binH:O.postInfH);  eB1 = (-O.preH:O.binH:O.postB1H);
EDG  = struct('inf', eInf, 'B1', eB1);

S = struct('Group', {}, 'Animal', {}, 'Session', {}, 'gapH', {}, 'nPYR', {}, 'baseS', {}, ...
    'critS', {}, 'D', {}, 'C', {}, 'N', {}, 'P', {});
led = cell(0, 3);
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
for i = 1:numel(fa)
    fld = fileparts(fa{i}); [~, snL] = fileparts(fld);
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    if ~ismember(m.condition, {'SAL','ASV'})
        led(end+1,:) = {snL, '', 'not an infused session'}; continue %#ok<AGROW>
    end
    g = m.condition;
    ev = local_events(fld);
    if isempty(ev), led(end+1,:) = {snL, g, 'no event file'}; continue; end %#ok<AGROW>
    tOf = @(lab) ev.time(find(strcmp(ev.description, lab), 1));
    [tInf, infSrc, infNote] = local_infTime(ev, m);
    if isempty(tInf), tInf = tOf('ConA1_end'); infSrc = 'ConA1_end'; end
    if ~isempty(infNote), fprintf('   %-27s %-4s %s\n', snL, g, infNote); end
    tB0 = tOf('ConB1_start'); tB1 = tOf('ConB1_end');
    if isempty(tInf) || isempty(tB0) || isempty(tB1)
        led(end+1,:) = {snL, g, 'no infusion time or no ConB1'}; continue %#ok<AGROW>
    end
    if ~(tB0 > tInf)
        led(end+1,:) = {snL, g, sprintf('ConB1 (%.0f s) does not follow the infusion (%.0f s, from %s)', ...
            tB0, tInf, infSrc)}; continue %#ok<AGROW>
    end

    % NREM, minus ripples; spikes and labels
    nrem = local_nrem(fld);
    if isempty(nrem), led(end+1,:) = {snL, g, 'no SleepState NREM'}; continue; end %#ok<AGROW>
    rip = [RP.ts1(:), RP.ts1Off(:)];
    [ts, isP, isD, isS, isI] = local_units(fld);
    if isempty(ts), led(end+1,:) = {snL, g, 'no spikes or no cell_metrics match'}; continue; end %#ok<AGROW>
    strat = struct('PYR', isP, 'PYRdeep', isP & isD, 'PYRsup', isP & isS, 'INT', isI);
    if sum(isP) < O.minUnits, led(end+1,:) = {snL, g, sprintf('%d PYR units', sum(isP))}; continue; end %#ok<AGROW>

    rec = [0, max(cellfun(@(t) max([t; 0]), ts))];
    OK = struct();
    OK.full     = local_ivSub(local_ivInt(local_ivMerge(nrem), rec), local_ivMerge(rip));
    OK.trim     = local_ivSub(local_ivInt(local_ivTrim(nrem, O.trimS), rec), local_ivMerge(rip));
    OK.longbout = local_ivSub(local_ivInt(local_ivMerge(nrem(diff(nrem, 1, 2) >= O.minBoutS, :)), rec), ...
                              local_ivMerge(rip));

    % baseline: pre-infusion NREM out of ripples
    base = local_ivInt(OK.full, [rec(1), tInf]);
    baseS = local_ivDur(base);
    if baseS < O.minBase
        led(end+1,:) = {snL, g, sprintf('%.0f s pre-infusion NREM', baseS)}; continue %#ok<AGROW>
    end
    nb = cellfun(@(t) local_ivCount(t, base), ts);
    rb = nb / baseS;
    okU = nb >= O.minBaseSpk;
    lb = log10(rb + 0.01);

    D = struct(); C = struct(); N = struct(); P = struct();
    for v = VAR
        for c = CLK
            e = EDG.(c{1}); t0 = tInf; if strcmp(c{1}, 'B1'), t0 = tB1; end
            nBin = numel(e) - 1;
            cov = zeros(1, nBin); dm = nan(numel(STRAT), nBin); nn = zeros(numel(STRAT), nBin);
            fp = nan(1, nBin);                                   % fraction of the bin's NREM after ConB1
            for b = 1:nBin
                W = local_ivInt(OK.(v{1}), t0 + 3600 * [e(b), e(b+1)]);    % edges are in h, times in s
                cov(b) = local_ivDur(W);
                if cov(b) > 0                                   % share of this bin's NREM after ConB1
                    Wp = W(W(:,2) > tB1, :);
                    if ~isempty(Wp), Wp(:,1) = max(Wp(:,1), tB1); end
                    fp(b) = local_ivDur(Wp) / cov(b);
                end
                if cov(b) < O.minCov, continue; end
                r = cellfun(@(t) local_ivCount(t, W), ts) / cov(b);
                d = log10(r + 0.01) - lb;
                for s = 1:numel(STRAT)
                    k = okU & strat.(STRAT{s}) & isfinite(d);
                    nn(s, b) = sum(k);
                    if sum(k) >= O.minUnits, dm(s, b) = median(d(k)); end
                end
            end
            D.(v{1}).(c{1}) = dm; C.(v{1}).(c{1}) = cov; N.(v{1}).(c{1}) = nn; P.(v{1}).(c{1}) = fp;
        end
    end

    S(end+1) = struct('Group', g, 'Animal', char(string(m.animal)), 'Session', snL, ...
        'gapH', (tB1 - tInf) / 3600, 'nPYR', sum(isP), 'baseS', baseS, ...
        'critS', local_ivDur(local_ivInt(OK.full, [tInf, tB0])), 'D', D, 'C', C, 'N', N, 'P', P); %#ok<AGROW>
    fprintf('%-30s %-4s gap %.2f h | base %5.0f s | infusion-to-ConB1 NREM %5.0f s | %3d PYR | scored bins: inf %d/%d, B1 %d/%d\n', ...
        snL, g, S(end).gapH, baseS, S(end).critS, sum(isP), ...
        sum(isfinite(D.full.inf(1,:))), size(D.full.inf, 2), sum(isfinite(D.full.B1(1,:))), size(D.full.B1, 2));
end
if isempty(S), error('rdt:none', 'no session produced a trajectory'); end

Res = local_report(S, O, STRAT, CLK, VAR, EDG);
local_figure(S, O, STRAT, EDG, Res, CACHE_ROOT);
if ~isempty(led)
    fprintf('\nnot used:\n');
    for q = 1:size(led, 1), fprintf('   %-30s %-4s %s\n', led{q,:}); end
end
if SAVE_OUT
    save(fullfile(CACHE_ROOT, 'rateDivergenceTiming.mat'), 'S', 'Res', 'O', 'EDG', 'STRAT');
    fprintf('\nsaved rateDivergenceTiming.mat\n');
end


%% =============================================================================================
%  report
%% =============================================================================================
function Res = local_report(S, O, STRAT, CLK, VAR, EDG)
Res = struct();
grp = {S.Group}; ani = {S.Animal};
fprintf('\n================ WHEN DOES THE ARM DIFFERENCE APPEAR? ================\n');
fprintf('Out-of-ripple NREM firing, each neuron against its own pre-infusion rate; a session is\n');
fprintf('the median over neurons; an animal is the mean over its sessions in that arm.\n');

% ---- coverage audit: the gap jitter and the critical window ----
fprintf('\n---- coverage (this decides what the timing question can be asked of) ----\n');
gap = [S.gapH]; crit = [S.critS]; base = [S.baseS];
fprintf('  infusion to end of ConB1: median %.2f h, range %.2f-%.2f h, SD %.2f h over %d sessions\n', ...
    median(gap), min(gap), max(gap), std(gap), numel(gap));
fprintf('  -> the spread is what makes Group:PostB identifiable at matched time since infusion\n');
for g = {'SAL','ASV'}
    k = strcmp(grp, g{1});
    if ~any(k), continue; end
    fprintf('  %-4s %2d sessions, %2d animals | pre-infusion NREM %5.0f s | infusion-to-ConB1 NREM %5.0f s (median)\n', ...
        g{1}, sum(k), numel(unique(ani(k))), median(base(k)), median(crit(k)));
    fprintf('       sessions with >= %d s in that window: %d of %d\n', O.minCov, sum(crit(k) >= O.minCov), sum(k));
end
Res.gap = gap; Res.crit = crit;

for v = VAR
    fprintf('\n################ %s ################\n', local_varName(v{1}, O));
    for c = CLK
        e = EDG.(c{1}); ctr = e(1:end-1) + diff(e)/2;
        fprintf('\n======== clock: %s ========\n', local_clkName(c{1}));
        for s = 1:numel(STRAT)
            Y = local_stack(S, v{1}, c{1}, s);                       % sessions x bins
            [d, uaBoth] = local_armDiff(Y, grp, ani);                % animals x bins
            if isempty(d), continue; end
            fprintf('\n  %s (%d animals with both arms)\n', STRAT{s}, numel(uaBoth));
            for g = {'SAL','ASV'}
                x = local_animalMat(Y, grp, ani, g{1});
                fprintf('    %-4s per bin:', g{1});
                fprintf(' %s', local_fmtRow(mean(x, 1, 'omitnan')));
                fprintf('\n');
            end
            b0 = find(ctr > 0, 1);                                   % scan forward from the anchor
            [pcal, tObs, praw] = local_maxT(d(:, b0:end), O.nFlipMax);
            dm = mean(d, 1, 'omitnan');
            fprintf('    ASV - SAL  :'); fprintf(' %s', local_fmtRow(dm)); fprintf('\n');
            fprintf('    bins scanned from %.2f h; calibrated by max-t over whole-trajectory sign flips\n', ctr(b0));
            first = NaN;
            for q = 1:numel(pcal)
                b = b0 + q - 1;
                if ~isfinite(tObs(q)), continue; end
                fprintf('      %+5.2f to %+5.2f h  d = %+.3f  t = %+5.2f  p = %.3f  p_cal = %.3f  (n = %d)\n', ...
                    e(b), e(b+1), dm(b), tObs(q), praw(q), pcal(q), sum(isfinite(d(:, b))));
                if pcal(q) < 0.05 && ~isfinite(first), first = e(b); end
            end
            if isfinite(first)
                fprintf('    -> first calibrated divergence in the bin starting %+.2f h\n', first);
            else
                fprintf('    -> no bin reaches calibrated significance after the anchor\n');
            end
            % null check: the pre-anchor bins
            pre = 1:(b0 - 1);
            if ~isempty(pre)
                [pcP, tP] = local_maxT(d(:, pre), O.nFlipMax);
                [~, iw] = max(abs(tP));
                if isfinite(tP(iw))
                    fprintf('    null check, pre-anchor bins: strongest at %+.2f h, t = %+.2f, p_cal = %.3f\n', ...
                        e(pre(iw)), tP(iw), pcP(iw));
                end
            end
            Res.(v{1}).(c{1}).(STRAT{s}) = struct('d', d, 'animals', {uaBoth}, 'dm', dm, ...
                'pcal', pcal, 'tObs', tObs, 'praw', praw, 'scanFrom', b0, 'first', first, 'ctr', ctr);
        end
    end
end

% ---- both anchors in one model, using the gap jitter ----
fprintf('\n================ BOTH ANCHORS IN ONE MODEL ================\n');
fprintf('d ~ Group*PostB + tInf + (1|Animal) + (1|Session) on the infusion clock; PostB = the bin\n');
fprintf('falls after\n');
fprintf('ConB1 ends. Group:PostB > 0 in magnitude means the arm difference grows at the novel\n');
fprintf('exposure beyond what time since the infusion accounts for.\n');
for v = VAR
    for s = 1:numel(STRAT)
        T = local_binTable(S, v{1}, 'inf', s, EDG.inf);
        if isempty(T)
            fprintf('  %-10s %-8s no scored bins\n', local_varName(v{1}, O), STRAT{s}); continue
        end
        nPre = sum(T.PostB == 0); nPost = sum(T.PostB == 1);
        if numel(unique(T.Animal)) < 3 || nPre < 4 || nPost < 4
            fprintf('  %-10s %-8s %d animals, %d pre-ConB1 and %d post-ConB1 bin-rows: not enough to fit\n', ...
                local_varName(v{1}, O), STRAT{s}, numel(unique(T.Animal)), nPre, nPost);
            continue
        end
        txt = sprintf('  %-10s %-8s %4d rows (%d pre, %d post), %d sessions', local_varName(v{1}, O), ...
            STRAT{s}, height(T), nPre, nPost, numel(unique(T.Session)));
        if exist('fitlme', 'file') == 2
            try
                M = fitlme(T, 'd ~ Group*PostB + tInf + (1|Animal) + (1|Session)');
                [b, nm, st] = fixedEffects(M, 'DFMethod', 'satterthwaite');
                kI = contains(nm.Name, 'Group') & contains(nm.Name, 'PostB');
                kG = strcmp(nm.Name, 'Group_ASV'); kP = strcmp(nm.Name, 'PostB');
                txt = [txt sprintf(' | Group %+.3f p = %.3f | PostB %+.3f p = %.3f | Group:PostB %+.3f p = %.3f', ...
                    b(kG), st.pValue(kG), b(kP), st.pValue(kP), b(kI), st.pValue(kI))]; %#ok<AGROW>
                Res.anchor.(v{1}).(STRAT{s}) = struct('b', b, 'name', {nm.Name}, 'p', st.pValue);
            catch ME
                txt = [txt ' | LME failed: ' ME.message]; %#ok<AGROW>
            end
        end
        fprintf('%s\n', txt);
    end
end

% ---- the same contrast with the animal as the unit ----
fprintf('\n---- the same contrast, animal as the unit (this is the reportable version) ----\n');
fprintf('Per session, the mean d over its pre-ConB1 bins and over its post-ConB1 bins; then per\n');
fprintf('animal and arm, the mean over sessions. The LME above treats bins as independent given\n');
fprintf('animal, so its p-values track the bin count rather than the evidence; these do not.\n');
fprintf('  post-pre  within each arm: does that arm rise at ConB1 at all? The saline arm answers\n');
fprintf('            the timing question on its own, with no reference to the drug\n');
fprintf('  dPre   ASV - SAL before ConB1 (the null the timing claim needs)\n');
fprintf('  dPost  ASV - SAL after ConB1\n');
fprintf('  DiD    (post - pre | ASV) - (post - pre | SAL), exact sign-flip over animals\n');
for v = VAR
    for s = 1:numel(STRAT)
        [D_, ua] = local_didAnimal(S, v{1}, s);
        if isempty(ua) || numel(ua) < 2
            fprintf('  %-10s %-8s too few animals with both arms\n', local_varName(v{1}, O), STRAT{s}); continue
        end
        fprintf(['  %-10s %-8s n = %d | SAL post-pre %+.3f +- %.3f p = %.3f, ASV %+.3f +- %.3f p = %.3f | ' ...
                 'dPre %+.3f p = %.3f | dPost %+.3f p = %.3f | DiD %+.3f +- %.3f p = %.3f\n'], ...
            local_varName(v{1}, O), STRAT{s}, numel(ua), ...
            mean(D_.chS, 'omitnan'), local_sem(D_.chS), local_signflip(D_.chS), ...
            mean(D_.chA, 'omitnan'), local_sem(D_.chA), local_signflip(D_.chA), ...
            mean(D_.pre, 'omitnan'), local_signflip(D_.pre), ...
            mean(D_.post, 'omitnan'), local_signflip(D_.post), ...
            mean(D_.did, 'omitnan'), local_sem(D_.did), local_signflip(D_.did));
        Res.did.(v{1}).(STRAT{s}) = D_; Res.did.(v{1}).animals = {ua};
    end
end

fprintf('\n  Reading: a divergence already present in the bins between the infusion and ConB1, with a\n');
fprintf('  flat Group:PostB, favours an infusion-associated trigger. A divergence that starts in the\n');
fprintf('  first bins after ConB1, with the pre-ConB1 bins null and Group:PostB away from zero,\n');
fprintf('  favours novelty or an infusion-by-novelty interaction. Check the coverage rows before\n');
fprintf('  reading a null in the pre-ConB1 window as evidence: it may simply be unscored.\n');
fprintf('\n  CAVEAT on the first-crossing bin. The scan holds its nominal level (simulated family-wise\n');
fprintf('  false positive 0.050 at 8 animals and 12 bins) and almost never fires before the true\n');
fprintf('  change (0.005 of runs), but it is biased LATE: a step of 1 SD was detected in 0.65 of\n');
fprintf('  runs and first flagged about two bins after it began, and a step of 0.5 SD in only 0.18.\n');
fprintf('  So the first-crossing bin is an UPPER BOUND on when the divergence began, not an estimate\n');
fprintf('  of its onset. That bias runs toward the wrong conclusion here: a divergence that truly\n');
fprintf('  starts at the infusion can be first flagged a bin or two later, which at this gap is\n');
fprintf('  around ConB1. The claim that novelty is the proximal cause therefore rests on the\n');
fprintf('  pre-ConB1 bins being null WITH adequate coverage, and on Group:PostB, not on where the\n');
fprintf('  scan first crosses.\n');
end


function nm = local_varName(v, O)
switch v
    case 'full',     nm = 'all NREM';
    case 'trim',     nm = sprintf('first %ds', O.trimS);
    otherwise,       nm = sprintf('bouts>%ds', O.minBoutS);
end
end


function nm = local_clkName(c)
if strcmp(c, 'inf'), nm = 'time since the infusion'; else, nm = 'time since ConB1 ends'; end
end


function s = local_fmtRow(v)
s = strtrim(sprintf('%+.2f ', v));
end


function Y = local_stack(S, v, c, s)
% sessions x bins, one stratum
Y = cell2mat(arrayfun(@(x) x.D.(v).(c)(s, :), S(:), 'UniformOutput', false));
end


function X = local_animalMat(Y, grp, ani, g)
% animals x bins, the mean over that animal's sessions in one arm
k = strcmp(grp, g); ua = unique(ani(k));
X = nan(numel(ua), size(Y, 2));
for q = 1:numel(ua)
    X(q, :) = mean(Y(k & strcmp(ani, ua{q}), :), 1, 'omitnan');
end
end


function [d, uaBoth] = local_armDiff(Y, grp, ani)
% paired ASV - SAL per animal per bin, animals contributing to both arms
d = []; uaBoth = {};
uA = unique(ani(strcmp(grp, 'ASV'))); uS = unique(ani(strcmp(grp, 'SAL')));
uaBoth = intersect(uA, uS);
if isempty(uaBoth), return; end
d = nan(numel(uaBoth), size(Y, 2));
for q = 1:numel(uaBoth)
    a = mean(Y(strcmp(grp, 'ASV') & strcmp(ani, uaBoth{q}), :), 1, 'omitnan');
    s = mean(Y(strcmp(grp, 'SAL') & strcmp(ani, uaBoth{q}), :), 1, 'omitnan');
    d(q, :) = a - s;
end
end


function [D, ua] = local_didAnimal(S, v, s)
% Pre- and post-ConB1 means per session (so a session with more scored bins does not dominate),
% then per animal and arm, then the arm contrasts. The split uses each session's own ConB1, taken
% from the share of the bin's NREM that lies after it, so the gap jitter is handled per session.
grp = {S.Group}'; ani = {S.Animal}';
pre = nan(numel(S), 1); post = nan(numel(S), 1);
for q = 1:numel(S)
    y = S(q).D.(v).inf(s, :); p = S(q).P.(v).inf; k = isfinite(y);
    if any(k & p <= 0.5),  pre(q)  = mean(y(k & p <= 0.5)); end
    if any(k & p >  0.5),  post(q) = mean(y(k & p >  0.5)); end
end
ua = intersect(unique(ani(strcmp(grp, 'ASV'))), unique(ani(strcmp(grp, 'SAL'))));
n = numel(ua);
D = struct('pre', nan(n,1), 'post', nan(n,1), 'did', nan(n,1), 'chA', nan(n,1), 'chS', nan(n,1));
for q = 1:n
    kA = strcmp(grp, 'ASV') & strcmp(ani, ua{q});
    kS = strcmp(grp, 'SAL') & strcmp(ani, ua{q});
    pA = mean(pre(kA), 'omitnan'); qA = mean(post(kA), 'omitnan');
    pS = mean(pre(kS), 'omitnan'); qS = mean(post(kS), 'omitnan');
    D.pre(q) = pA - pS; D.post(q) = qA - qS;
    D.chA(q) = qA - pA; D.chS(q) = qS - pS;
    D.did(q) = D.chA(q) - D.chS(q);
end
end


function e = local_sem(v)
v = v(isfinite(v));
e = NaN; if numel(v) >= 2, e = std(v) / sqrt(numel(v)); end
end


function p = local_signflip(d)
d = d(isfinite(d)); n = numel(d); p = NaN;
if n < 2, return; end
obs = abs(mean(d));
if n <= 14
    Sg = 2*(dec2bin(0:2^n - 1, n) == '1') - 1; p = mean(abs(Sg*d(:)/n) >= obs - 1e-12);
else
    nul = arrayfun(@(k) mean((2*(rand(n,1) > 0.5) - 1) .* d(:)), 1:1e4); p = (1 + sum(abs(nul) >= obs))/(1 + 1e4);
end
end


function [pcal, tObs, praw] = local_maxT(d, nFlipMax)
% single-step max-t over the scanned bins. A flip reverses an animal's WHOLE trajectory, so the
% correlation between bins within an animal is preserved under the null.
[na, nb] = size(d);
pcal = nan(1, nb); tObs = nan(1, nb); praw = nan(1, nb);
if na < 2 || nb < 1, return; end
tObs = local_tRow(d);
if na <= 12
    Sg = (2*(dec2bin(0:2^na - 1, na) - '0') - 1)';
else
    Sg = sign(randn(na, nFlipMax)); Sg(:, 1) = 1;
end
nS = size(Sg, 2); mx = nan(nS, 1); ge = zeros(1, nb);
for q = 1:nS
    tb = local_tRow(d .* Sg(:, q));
    mx(q) = max(abs(tb), [], 'omitnan');
    ge = ge + (abs(tb) >= abs(tObs) - 1e-12);
end
for b = 1:nb
    if ~isfinite(tObs(b)), continue; end
    pcal(b) = max(mean(mx >= abs(tObs(b)) - 1e-12), 1/nS);
    praw(b) = max(ge(b)/nS, 1/nS);
end
end


function t = local_tRow(x)
t = nan(1, size(x, 2));
for b = 1:size(x, 2)
    v = x(:, b); v = v(isfinite(v));
    if numel(v) >= 2 && std(v) > 0, t(b) = mean(v) / (std(v)/sqrt(numel(v))); end
end
end


function T = local_binTable(S, v, c, s, e)
% one row per session x bin on the infusion clock
ctr = e(1:end-1) + diff(e)/2;
d = []; gp = {}; an = {}; se = {}; tI = []; pb = [];
for q = 1:numel(S)
    y = S(q).D.(v).(c)(s, :);
    k = find(isfinite(y));
    for b = k
        d(end+1, 1) = y(b); %#ok<AGROW>
        gp{end+1, 1} = S(q).Group; an{end+1, 1} = S(q).Animal; se{end+1, 1} = S(q).Session; %#ok<AGROW>
        % PostB from where the bin's NREM actually lies, not from the bin centre: a bin straddling
        % ConB1 holds NREM from both sides, and the exposure itself contributes almost none
        tI(end+1, 1) = ctr(b); pb(end+1, 1) = double(S(q).P.(v).(c)(b) > 0.5); %#ok<AGROW>
    end
end
T = table();
if isempty(d), return; end
T = table(d, categorical(gp, {'SAL','ASV'}), categorical(an), se, tI, pb, ...
    'VariableNames', {'d','Group','Animal','Session','tInf','PostB'});
end


%% =============================================================================================
%  figure
%% =============================================================================================
function local_figure(S, O, STRAT, EDG, Res, cacheRoot)
COL = struct('SAL', [0.16 0.44 0.78], 'ASV', [0.84 0.30 0.12]);
grp = {S.Group}; ani = {S.Animal};
sl = STRAT(1:min(2, numel(STRAT)));
f = figure('Color', 'w', 'Position', [40 40 1250 700], 'Name', 'rate divergence timing');
t = tiledlayout(numel(sl), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:numel(sl)
    is = find(strcmp(STRAT, sl{s}));
    for c = {'inf','B1'}
        e = EDG.(c{1}); ctr = e(1:end-1) + diff(e)/2;
        nexttile; hold on
        Y = local_stack(S, 'full', c{1}, is);
        for g = {'SAL','ASV'}
            X = local_animalMat(Y, grp, ani, g{1});
            mu = mean(X, 1, 'omitnan');
            n = sum(isfinite(X), 1); sd = std(X, 0, 1, 'omitnan');
            se = sd ./ sqrt(max(n, 1)); se(n < 2) = NaN;
            ok = isfinite(mu);
            errorbar(ctr(ok), mu(ok), se(ok), '-o', 'Color', COL.(g{1}), 'MarkerFaceColor', COL.(g{1}), ...
                'LineWidth', 1.4, 'CapSize', 0, 'MarkerSize', 4);
        end
        yline(0, 'k:'); xline(0, 'k-', 'LineWidth', 1.2);
        if strcmp(c{1}, 'inf')
            xline(median([S.gapH]), 'k--', 'LineWidth', 1);
            xlabel('h since the infusion'); lg = {'SAL','ASV'};
        else
            xlabel('h since ConB1 ends'); lg = {'SAL','ASV'};
        end
        ylabel(sprintf('\\Deltalog_{10} %s rate', sl{s}));
        legend(lg, 'Location', 'southwest', 'Box', 'off');
        title(sprintf('%s  %s, %s', char('a' + (s-1)*3 + find(strcmp(c, {'inf','B1'})) - 1), ...
            sl{s}, local_clkName(c{1})), 'FontWeight', 'normal');
    end
    % the difference on the ConB1 clock, with the calibrated divergence marked
    nexttile; hold on
    e = EDG.B1; ctr = e(1:end-1) + diff(e)/2;
    if isfield(Res, 'full') && isfield(Res.full, 'B1') && isfield(Res.full.B1, sl{s})
        Q = Res.full.B1.(sl{s});
        n = sum(isfinite(Q.d), 1); sd = std(Q.d, 0, 1, 'omitnan');
        se = sd ./ sqrt(max(n, 1)); se(n < 2) = NaN;
        ok = isfinite(Q.dm);
        errorbar(ctr(ok), Q.dm(ok), se(ok), '-o', 'Color', [0.1 0.1 0.1], 'MarkerFaceColor', [0.1 0.1 0.1], ...
            'LineWidth', 1.4, 'CapSize', 0, 'MarkerSize', 4);
        sig = false(1, numel(ctr));
        sig(Q.scanFrom:end) = Q.pcal < 0.05;
        if any(sig), plot(ctr(sig), Q.dm(sig), 'o', 'Color', [0.8 0 0], 'MarkerFaceColor', [0.8 0 0], 'MarkerSize', 7); end
        if isfinite(Q.first), xline(Q.first, '--', 'Color', [0.8 0 0], 'LineWidth', 1.2); end
    end
    yline(0, 'k:'); xline(0, 'k-', 'LineWidth', 1.2);
    xlabel('h since ConB1 ends'); ylabel('ASV - SAL');
    title(sprintf('%s  %s, arm difference (filled = p_{cal} < 0.05)', ...
        char('a' + (s-1)*3 + 2), sl{s}), 'FontWeight', 'normal');
end
title(t, sprintf(['Out-of-ripple NREM firing against each neuron''s pre-infusion rate. ' ...
    'Solid line = anchor, dashed on the left panels = median ConB1 (%.2f h after the infusion). ' ...
    'Bins need %d s of NREM.'], median([S.gapH]), O.minCov));
exportgraphics(f, fullfile(cacheRoot, 'rateDivergenceTiming.pdf'), 'ContentType', 'vector');

% coverage, as its own figure: a null in an unscored window is not a null result
f2 = figure('Color', 'w', 'Position', [60 60 1100 380], 'Name', 'NREM coverage');
t2 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = {'inf','B1'}
    e = EDG.(c{1}); ctr = e(1:end-1) + diff(e)/2;
    nexttile; hold on
    for g = {'SAL','ASV'}
        k = strcmp(grp, g{1});
        Cv = cell2mat(arrayfun(@(x) x.C.full.(c{1}), S(k)', 'UniformOutput', false));
        plot(ctr, median(Cv, 1), '-o', 'Color', COL.(g{1}), 'MarkerFaceColor', COL.(g{1}), ...
            'LineWidth', 1.4, 'MarkerSize', 4);
    end
    yline(O.minCov, 'k:'); xline(0, 'k-', 'LineWidth', 1.2);
    if strcmp(c{1}, 'inf'), xline(median([S.gapH]), 'k--', 'LineWidth', 1); end
    xlabel(sprintf('h since %s', local_clkName(c{1}))); ylabel('NREM out-of-ripple s per bin (median)');
    legend({'SAL','ASV'}, 'Location', 'northeast', 'Box', 'off');
    title(local_clkName(c{1}), 'FontWeight', 'normal');
end
title(t2, 'NREM coverage per bin; dotted line = the floor a bin must clear to be scored');
exportgraphics(f2, fullfile(cacheRoot, 'rateDivergenceTiming_coverage.pdf'), 'ContentType', 'vector');
end


%% =============================================================================================
%  data loading
%% =============================================================================================
function [tInf, src, note] = local_infTime(ev, m)
% The event file is the clock to trust, because ConB1 is read from it. Taking the infusion from
% RIPPLES.meta instead mixes two sources: if a meta value was left in ms while LoadEvents returns
% the .evt in s, the infusion lands far past the end of the recording and ConB1 appears to precede
% it, which silently drops the session. So the event file wins where it has the label, meta is the
% fallback, and a meta value on the wrong scale is rescaled and reported rather than used as-is.
tInf = []; src = ''; note = '';
span = max(ev.time);
tE = ev.time(find(strcmp(ev.description, 'Infusion_start'), 1));
tM = [];
if isfield(m, 'infusionSec'), tM = double(m.infusionSec); end
if ~isscalar(tM) || ~isfinite(tM)
    tM = [];
elseif tM > 1.5*span && tM/1000 <= 1.5*span
    note = sprintf('meta infusionSec %.0f is on a ms scale against events spanning %.0f s; divided by 1000', tM, span);
    tM = tM / 1000;
end
if ~isempty(tE)
    tInf = tE; src = 'event file';
    if ~isempty(tM) && abs(tM - tE) > 60
        note = strtrim(sprintf('%s | meta %.0f s and event file %.0f s differ by %.0f s; using the event file', ...
            note, tM, tE, abs(tM - tE)));
    end
elseif ~isempty(tM)
    tInf = tM; src = 'RIPPLES.meta';
end
end


function iv = local_nrem(fld)
iv = zeros(0, 2);
sf = dir(fullfile(fld, '*.SleepState.states.mat'));
if isempty(sf), return; end
try
    Ls = load(fullfile(sf(1).folder, sf(1).name), 'SleepState');
catch
    return
end
if ~isfield(Ls, 'SleepState') || ~isfield(Ls.SleepState, 'ints'), return; end
if ~isfield(Ls.SleepState.ints, 'NREMstate'), return; end
iv = Ls.SleepState.ints.NREMstate;
if size(iv, 2) ~= 2, iv = reshape(iv, [], 2); end
iv = iv(all(isfinite(iv), 2) & iv(:,2) > iv(:,1), :);
end


function [ts, isP, isD, isS, isI] = local_units(fld)
ts = {}; isP = []; isD = []; isS = []; isI = [];
sf = dir(fullfile(fld, '*spikes.cellinfo.mat')); cf = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(sf) || isempty(cf), return; end
Sp = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); sp = Sp.spikes;
Cm = load(fullfile(cf(1).folder, cf(1).name), 'cell_metrics'); cm = Cm.cell_metrics;
ts = cellfun(@(x) sort(reshape(x(~isnan(x)), [], 1)), sp.times(:), 'UniformOutput', false);
nU = numel(ts); nC = numel(cm.putativeCellType);
cmType = cellstr(string(cm.putativeCellType(:)));
cmLay = repmat({''}, nC, 1);
if isfield(cm, 'deepSuperficial') && numel(cm.deepSuperficial) == nC
    cmLay = cellstr(string(cm.deepSuperficial(:)));
end
clu = []; if isfield(sp, 'UID') && numel(sp.UID) == nU, clu = sp.UID(:); end
if isfield(sp, 'cluID') && numel(sp.cluID) == nU, clu = sp.cluID(:); end
loc = matchCellMetrics(cm, (1:nU)', clu, nU);
ct = repmat({''}, nU, 1); ly = ct; ok = loc > 0;
ct(ok) = cmType(loc(ok)); ly(ok) = cmLay(loc(ok));
isP = contains(ct, 'Pyramidal'); isI = contains(ct, 'Interneuron');
isD = strcmp(ly, 'Deep'); isS = strcmp(ly, 'Superficial');
if ~any(isP), ts = {}; end
end


%% =============================================================================================
%  interval algebra
%% =============================================================================================
function M = local_ivMerge(iv)
M = zeros(0, 2);
if isempty(iv), return; end
iv = sortrows(iv(all(isfinite(iv), 2) & iv(:,2) > iv(:,1), :));
if isempty(iv), return; end
M = iv(1, :);
for q = 2:size(iv, 1)
    if iv(q,1) <= M(end,2), M(end,2) = max(M(end,2), iv(q,2)); else, M(end+1,:) = iv(q,:); end %#ok<AGROW>
end
end


function M = local_ivInt(A, B)
% intersection of two interval sets
M = zeros(0, 2);
A = local_ivMerge(A); B = local_ivMerge(B);
i = 1; j = 1;
while i <= size(A, 1) && j <= size(B, 1)
    lo = max(A(i,1), B(j,1)); hi = min(A(i,2), B(j,2));
    if hi > lo, M(end+1,:) = [lo hi]; end %#ok<AGROW>
    if A(i,2) < B(j,2), i = i + 1; else, j = j + 1; end
end
end


function M = local_ivSub(A, B)
% A minus B
A = local_ivMerge(A); B = local_ivMerge(B);
if isempty(B), M = A; return; end
big = max([A(:); B(:)]) + 1; small = min([A(:); B(:)]) - 1;
Cp = [small, B(1,1)];
for q = 2:size(B, 1), Cp(end+1,:) = [B(q-1,2), B(q,1)]; end %#ok<AGROW>
Cp(end+1,:) = [B(end,2), big];
M = local_ivInt(A, Cp(Cp(:,2) > Cp(:,1), :));
end


function M = local_ivTrim(iv, trimS)
% keep only the first trimS s of every bout
M = zeros(0, 2);
iv = local_ivMerge(iv);
if isempty(iv), return; end
M = [iv(:,1), min(iv(:,2), iv(:,1) + trimS)];
end


function d = local_ivDur(M)
d = 0;
if ~isempty(M), d = sum(diff(M, 1, 2)); end
end


function n = local_ivCount(t, M)
% spikes inside any interval of a merged, non-overlapping set
n = 0;
if isempty(M) || isempty(t), return; end
e = reshape(M', [], 1);
if any(diff(e) <= 0), M = local_ivMerge(M); e = reshape(M', [], 1); end
if isempty(e) || any(diff(e) <= 0), return; end
c = histcounts(t, e);
n = sum(c(1:2:end));
end


%% =============================================================================================
%  events (same parser as the other stages)
%% =============================================================================================
function ev = local_events(fld)
ev = []; T = {}; X = []; src = {};
c = dir(fullfile(fld, '*.evt.sti'));
if ~isempty(c)
    try
        ws_ = warning('off', 'all');
        E = LoadEvents(fullfile(c(1).folder, c(1).name));
        warning(ws_);
        if isfield(E, 'time') && isfield(E, 'description')
            T = cellstr(E.description(:)); X = E.time(:); src{end+1} = c(1).name;
        end
    catch
        if exist('ws_', 'var'), warning(ws_); end
    end
end
[b, ~, ~] = cellfun(@local_parseLabel, T, 'UniformOutput', false);
if ~any(startsWith(b, 'Con'))
    [~, sn] = fileparts(fld);
    dirs = {fld, fileparts(fld), fileparts(fileparts(fld))}; found = false;
    for k = 1:numel(dirs)
        a = dir(fullfile(dirs{k}, '*transitions*.mat'));
        for q = 1:numel(a)
            try
                L = load(fullfile(a(q).folder, a(q).name), 'all_events'); A = L.all_events;
            catch
                continue
            end
            m = true(height(A), 1);
            if ismember('video', A.Properties.VariableNames)
                hit = contains(string(A.video), sn);
                if any(hit), m = hit; elseif k > 1, continue; end
            end
            T = [T(~startsWith(b, {'HC','Con'})); cellstr(string(A.event(m)))]; %#ok<AGROW>
            X = [X(~startsWith(b, {'HC','Con'})); double(A.Intan_ts(m))]; %#ok<AGROW>
            src{end+1} = fullfile(a(q).folder, a(q).name); %#ok<AGROW>
            found = true; break
        end
        if found, break; end
    end
end
if isempty(T), return; end
[X, o] = sort(X); T = T(o);
[b, num, edge] = cellfun(@local_parseLabel, T, 'UniformOutput', false);
keep = ~cellfun(@isempty, b);
unparsed = unique(T(~keep));
T = T(keep); X = X(keep); b = b(keep); num = num(keep); edge = edge(keep);
cnt = struct('ConA', 0, 'ConB', 0); lab = cell(size(T));
for k = 1:numel(T)
    if startsWith(b{k}, 'Con')
        if isempty(num{k})
            if strcmp(edge{k}, 'start'), cnt.(b{k}) = cnt.(b{k}) + 1; end
            num{k} = sprintf('%d', max(1, cnt.(b{k})));
        elseif strcmp(edge{k}, 'start')
            cnt.(b{k}) = str2double(num{k});
        end
        lab{k} = sprintf('%s%s_%s', b{k}, num{k}, edge{k});
    else
        lab{k} = sprintf('%s_%s', b{k}, edge{k});
    end
end
starts = find(endsWith(lab, '_start') & startsWith(lab, 'Con'));
for k = starts(:)'
    nm = erase(lab{k}, '_start');
    if ~any(strcmp(lab, [nm '_end']) & X > X(k))
        nx = find(X > X(k), 1);
        tEnd = X(k) + 1800; if ~isempty(nx), tEnd = X(nx); end
        lab{end+1, 1} = [nm '_end']; X(end+1, 1) = tEnd; %#ok<AGROW>
        unparsed{end+1, 1} = sprintf('(%s_end inferred)', nm); %#ok<AGROW>
    end
end
[X, o] = sort(X); lab = lab(o);
ev = struct('time', X, 'description', {lab}, 'source', strjoin(src, ' + '), 'unparsed', {unparsed});
end


function [b, num, edge] = local_parseLabel(s)
b = ''; num = ''; edge = '';
x = lower(strtrim(char(s)));
x = regexprep(x, '[\s_\-\.:]+', ' ');
t = regexp(x, ['^(con(?:text)? ?[ab]|familiar|novel|hc|home ?cage|infusion|inf) ?(\d*) ?' ...
               '(start|begin|onset|on|end|stop|offset|off)$'], 'tokens', 'once');
if isempty(t), return; end
w = t{1};
if startsWith(w, 'con'), b = ['Con' upper(w(end))];
elseif strcmp(w, 'familiar'), b = 'ConA';
elseif strcmp(w, 'novel'), b = 'ConB';
elseif startsWith(w, 'h'), b = 'HC';
else, b = 'Infusion';
end
num = t{2};
if any(strcmp(t{3}, {'start','begin','onset','on'})), edge = 'start'; else, edge = 'end'; end
end