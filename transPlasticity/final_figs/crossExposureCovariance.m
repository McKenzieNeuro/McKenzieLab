%% crossExposureCovariance.m
% Does the pairwise covariance structure expressed in ConA1 / ConB1 persist into ConA2 / ConB2,
% and is that persistence specific to the context?
%
%   HC  ConA1  HC  Infusion  HC  ConB1  HC  ConA2  HC  ConB2  HC
%   A = familiar context, B = novel context. ConA1 is encoded drug-free; ConB1, ConA2 and ConB2
%   all follow the infusion, so in ASV sessions the novel map is formed, and both maps are
%   re-expressed, with the drug on board.
%
% MEASURE
%   For each exposure X, c^X_jk = partial correlation of pyramidal pair (j,k) in 100-ms bins of the
%   running (else theta) part of X, with the summed count of the other units partialled out, so a
%   shared rate fluctuation cannot make every pair look alike. Units need MIN_SPK spikes in that
%   part of EVERY exposure, so all four exposures describe the same pairs.
%   Similarity r(X,Y) = Pearson correlation of c^X and c^Y across pairs (Fisher z for averaging).
%   Much of r is stationary (anatomy, cell identity), so only contrasts are interpreted.
%
% CONTRASTS (per session; each compares two templates on the SAME target, so the target's own
% sampling noise is common to both terms)
%   S_A2 = z(A1,A2) - z(B1,A2)     does the familiar map predict the familiar re-exposure best?
%   S_B2 = z(B1,B2) - z(A1,B2)     does the novel map predict the novel re-exposure best?
%   I    = S_A2 + S_B2             context specificity. Under a drift that is linear in elapsed
%                                  time the lag terms cancel exactly in I (B1 is nearer both
%                                  targets than A1, which helps S_B2 and hurts S_A2 equally).
%   lag model: over the six exposure pairs of a session,
%        z(X,Y) = b0 + b_A 1(A1,A2) + b_B 1(B1,B2) + b_lag |t_X - t_Y|
%      b_A, b_B = same-context persistence of the familiar and the novel map beyond elapsed time.
%   Drug: ASV vs SAL on S_B2, b_B and b_B - b_A (the novel map is the one formed under drug; the
%   familiar map, formed drug-free, is the within-session reference).
%
% RELIABILITY
%   Each exposure is split into alternating 1-s windows; Spearman-Brown split-half reliability of
%   c^X gives rel_X, and every similarity is also reported disattenuated, r / sqrt(rel_X rel_Y), so a
%   drug change in spike counts (and hence in how well c is estimated) cannot pose as a change in
%   persistence.
%
% YIELD AND PER-PAIR STABILITY (does the all-four criterion hide instability?)
%   The inclusion rule keeps only units active in every exposure, so a unit that carried ConB1
%   coordination and then fell silent is excluded by construction.
%   (1) DROPOUT. Units passing MIN_SPK in each exposure ALONE are counted as well, and
%       retB1 = (units in all four) / (units passing in ConB1) is compared between arms. Equal
%       retention leaves the criterion a stated restriction; lower retention after ASV is itself
%       a result. The counts are repeated on the first nBinMin bins of every exposure, nBinMin =
%       the shortest exposure, so a shorter exposure cannot lose units by duration alone.
%   (2) WHERE THE LOSS FALLS. The per-exposure pass fractions separate loss in ConA2 / ConB2
%       (expression under drug) from loss spread over all four (yield or sorting over the day).
%   (5) PER-PAIR STABILITY. Across pairs: the fraction changing sign and the median |change|
%       between ConB1 and ConB2, with ConA1-ConA2 as the within-session reference. The median
%       |change| expected from sampling noise alone is sqrt(s_X^2 + s_Y^2)/2, where s_X is the
%       median |difference| between the two interleaved halves of exposure X (each half holds
%       about half the data, so its difference is about sqrt(2) wider than two full estimates).
%       ratio = observed / expected: about 1 means the change is within noise, above 1 means the
%       pairs reorganised. An across-pair correlation can stay high while that happens, so a
%       similar S_B2 with a larger ratio after ASV would be similarity masking instability.
%
% FAMILIARITY OF A
%   Context A is the same arena on every recording day; B is a new arena each day (the same within
%   a day). ASV days mostly come after saline days, so A has been seen more often on ASV days.
%   nPrior = number of earlier recording days (any session type with a ConA1 marker) for that
%   animal, DayIdx = days since its first such day. Both enter the arm comparison as covariates,
%   and S_A2 / b_A are regressed on nPrior within each arm: if familiar-map persistence grows with
%   A's familiarity inside an arm, the arm difference is session order, not drug.
%
% ASSIMILATION (does the novel map get folded into the familiar one?)
%   Per session, across pairs, in raw correlation units:
%        c^A2 = a0 + alpha c^A1 + beta  c^B1      beta  = novel structure absorbed into the
%                                                  familiar map beyond what A1 shared with B1
%        c^B2 = a0' + alpha' c^B1 + beta' c^A1     beta' = the reverse (control: assimilation
%                                                  should be directional, and B2 is also drugged)
%   OLS on the full-exposure correlations is attenuated by sampling noise in the predictors
%   (errors in variables), and noise is larger on ASV days (lower reliability), which would push
%   ASV beta toward 0 by itself. The IV estimate fixes this: predictors from one interleaved half
%   of each exposure, instruments from the other half (independent noise, same true correlation),
%   both assignments averaged. IV needs >= 30 pairs.
%   U = z(B1,A2) - z(A1,B1)      model-free version: does A2 move toward B1 relative to A1?
%   Rate counterpart, per neuron (log running rate, +0.1 Hz):
%        logr_A2 - logr_A1[half b] = g0 + g_A logr_A1[half a] + g_B logr_B1
%      g_B > 0: neurons active in B1 gain firing in A2 (A1 split so regression to the mean cannot
%      enter); control g'_A from logr_B2 - logr_B1[b] on logr_B1[a], logr_A1.
%   PRIMARY (fixed before this run): beta_IV, ASV - SAL, LME with DayIdx and animal intercept.
%   Everything else is secondary. Running time per exposure is reported and log(run A2/run A1)
%   enters as a covariate; rate maps are not available, so occupancy cannot be matched.
%   If pairCoincidenceGLM.mat is in the cache, each session's post-B1 ripple reinstatement
%   (joint model, at the never-infused tau) is related to beta within each arm.
%
% Never-infused sessions have two exposures only and are not used.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT       = 'R:\TransPlasticity';
SESSION_SET = 'all';
if exist('SHARED', 'var')
    CACHE_ROOT = SHARED.cacheRoot; ROOT = SHARED.root;
    if isfield(SHARED, 'sessionSet'), SESSION_SET = SHARED.sessionSet; end
end
O = struct('bin', 0.1, 'speedMin', 3, 'thetaRatio', 2, 'minSpk', 20, 'minUnits', 5);
CTX = {'ConA1', 'ConB1', 'ConA2', 'ConB2'};
SAVE_OUT = true;

rows = {}; led = cell(0, 3); dayA = cell(0, 2);   % every recording day with a ConA1 marker
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'rippleEvents.mat'));
for i = 1:numel(fa)
    fld = fileparts(fa{i}); [~, snL] = fileparts(fld);
    ev = local_events(fld);
    if isempty(ev) || ~any(contains(ev.description, 'ConA1_start'))
        led(end+1,:) = {snL, '', 'no event file or no ConA1'}; continue %#ok<AGROW>
    end
    L = load(fa{i}, 'RIPPLES'); RP = L.RIPPLES; m = RP.meta;
    dayA(end+1, :) = {char(string(m.animal)), local_sessionDate(snL)}; %#ok<AGROW>
    g = local_group(RP, ev, SESSION_SET);
    if isempty(g), led(end+1,:) = {snL, '', 'not SAL or ASV in this session set'}; continue; end %#ok<AGROW>
    tOf = @(lab) ev.time(find(strcmp(ev.description, lab), 1));
    iv = nan(4, 2);
    for c = 1:4
        on = tOf([CTX{c} '_start']); off = tOf([CTX{c} '_end']);
        if ~isempty(on) && ~isempty(off), iv(c, :) = [on off]; end
    end
    if any(~isfinite(iv(:)))
        led(end+1,:) = {snL, g, ['missing ' strjoin(CTX(any(~isfinite(iv), 2)), ', ')]}; continue %#ok<AGROW>
    end
    if any(diff(iv(:, 1)) <= 0)
        led(end+1,:) = {snL, g, 'exposures not in the order A1 B1 A2 B2'}; continue %#ok<AGROW>
    end
    [~, ~, ~, nU, ts] = local_ripCounts(fld, RP);
    if nU < O.minUnits, led(end+1,:) = {snL, g, sprintf('%d PYR units', nU)}; continue; end %#ok<AGROW>
    % binned counts, full and split halves, for each exposure
    Z = cell(4, 3); src = cell(4, 1);
    for c = 1:4
        [mask, src{c}] = local_runMask(fld, iv(c,1), iv(c,2), O.speedMin);
        if isempty(src{c}), [mask, src{c}] = local_thetaMask(RP, iv(c,1), iv(c,2), O.thetaRatio); end
        mask = local_chunk(mask, 1);
        Z{c, 1} = local_binMask(ts, mask, O.bin);
        Z{c, 2} = local_binMask(ts, mask(1:2:end, :), O.bin);
        Z{c, 3} = local_binMask(ts, mask(2:2:end, :), O.bin);
    end
    if any(cellfun(@(z) size(z, 1), Z(:, 1)) < 100)
        led(end+1,:) = {snL, g, 'an exposure has < 10 s of running/theta'}; continue %#ok<AGROW>
    end
    okU = true(nU, 1); nPass = zeros(1, 4); nPassM = zeros(1, 4);
    nBinMin = min(cellfun(@(z) size(z, 1), Z(:, 1)));      % bins of the shortest exposure
    for c = 1:4
        okC = sum(Z{c, 1}, 1)' >= O.minSpk;
        nPass(c)  = sum(okC);                              % units passing in this exposure alone
        nPassM(c) = sum(sum(Z{c, 1}(1:nBinMin, :), 1)' >= O.minSpk);   % duration-matched
        okU = okU & okC;
    end
    U = find(okU);
    if numel(U) < O.minUnits, led(end+1,:) = {snL, g, sprintf('%d units with %d spikes in all four', numel(U), O.minSpk)}; continue; end %#ok<AGROW>
    [jj, kk] = find(triu(true(numel(U)), 1)); idx = sub2ind([numel(U) numel(U)], jj, kk);
    C = nan(numel(idx), 4); C1 = C; C2 = C; rel = nan(1, 4);
    for c = 1:4
        M = local_partialCorr(Z{c, 1}, U); C(:, c) = M(idx);
        h1 = local_partialCorr(Z{c, 2}, U); h2 = local_partialCorr(Z{c, 3}, U);
        a = h1(idx); b = h2(idx); C1(:, c) = a; C2(:, c) = b; ok = isfinite(a) & isfinite(b);
        if sum(ok) >= 5, r = corr(a(ok), b(ok)); rel(c) = 2 * r / (1 + r); end
    end
    ok = all(isfinite(C), 2);
    okH = ok & all(isfinite(C1), 2) & all(isfinite(C2), 2);
    as = local_assim(C(ok, :), C1(okH, :), C2(okH, :), C(okH, :));
    ins = local_instab(C, C1, C2);
    tRun = cellfun(@(z) size(z, 1), Z(:, 1))' * O.bin;
    lr = zeros(numel(U), 4, 3);                 % unit x exposure x {full, half a, half b}
    for c = 1:4
        for hh = 1:3
            lr(:, c, hh) = log(sum(Z{c, hh}(:, U), 1)' / max(size(Z{c, hh}, 1) * O.bin, eps) + 0.1);
        end
    end
    rt = local_rateAssim(lr);
    if sum(ok) < 10, led(end+1,:) = {snL, g, sprintf('%d pairs with finite correlations', sum(ok))}; continue; end %#ok<AGROW>
    R = corr(C(ok, :));                                    % 4 x 4 similarity
    Rd = R ./ sqrt(max(rel, eps)' * max(rel, eps)); Rd(:, rel <= 0) = NaN; Rd(rel <= 0, :) = NaN;
    tm = mean(iv, 2) / 3600; lagM = abs(tm - tm');
    rows(end+1, :) = {g, char(string(m.animal)), snL, numel(U), sum(ok), R, Rd, rel, lagM, strjoin(unique(src), '; '), as, rt, tRun, nPass, nPassM, nU, ins}; %#ok<AGROW>
    fprintf('%-30s %-4s %3d/%d units %5d pairs | pass %d %d %d %d | rel %.2f %.2f %.2f %.2f\n', ...
        snL, g, numel(U), nU, sum(ok), nPass, rel);
end
if isempty(rows), error('cec:none', 'no session had all four exposures with enough units'); end
% plain struct of columns: the per-session matrices stay cells, scalars become vectors
T = struct('Group', {rows(:,1)}, 'Animal', {rows(:,2)}, 'Session', {rows(:,3)}, ...
    'nUnits', cell2mat(rows(:,4)), 'nPairs', cell2mat(rows(:,5)), 'R', {rows(:,6)}, 'Rd', {rows(:,7)}, ...
    'rel', {rows(:,8)}, 'lag', {rows(:,9)}, 'src', {rows(:,10)}, 'as', {rows(:,11)}, 'rt', {rows(:,12)}, ...
    'tRun', {rows(:,13)}, 'nPass', {rows(:,14)}, 'nPassM', {rows(:,15)}, ...
    'nPYR', cell2mat(rows(:,16)), 'ins', {rows(:,17)});
% familiarity of A: earlier recording days with a ConA1 marker, per animal
T.nPrior = nan(numel(T.Session), 1); T.DayIdx = nan(numel(T.Session), 1);
for q = 1:numel(T.Session)
    d0 = local_sessionDate(T.Session{q});
    k = strcmp(dayA(:,1), T.Animal{q}); dd = cell2mat(dayA(k, 2)); dd = unique(dd(isfinite(dd)));
    if isfinite(d0) && ~isempty(dd), T.nPrior(q) = sum(dd < d0); T.DayIdx(q) = d0 - min(dd); end
end
fprintf('\nsession order (prior ConA1 days, days since first):\n');
for a = unique(T.Animal)'
    k = find(strcmp(T.Animal, a{1})); [~, o] = sort(T.DayIdx(k)); k = k(o);
    fprintf('   %-5s', a{1});
    for q = k', fprintf('  %s n%d d+%d', T.Group{q}, T.nPrior(q), T.DayIdx(q)); end
    fprintf('\n');
end
Res = local_report(T);
Res.yield = local_yieldReport(T, O);
Res.stab  = local_instabReport(T);
Res.assim = local_assimReport(T, CACHE_ROOT);
if ~isempty(led)
    fprintf('\nnot used:\n');
    for q = 1:size(led, 1), fprintf('   %-30s %-4s %s\n', led{q,:}); end
end
if SAVE_OUT
    save(fullfile(CACHE_ROOT, 'crossExposureCovariance.mat'), 'T', 'Res', 'O');
    fprintf('\nsaved crossExposureCovariance.mat\n');
end


%% =============================================================================================
function Res = local_report(T)
Res = struct();
pr = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; nmP = {'A1-B1','A1-A2','A1-B2','B1-A2','B1-B2','A2-B2'};
fprintf('\n================ PERSISTENCE OF EXPOSURE COVARIANCE INTO THE RE-EXPOSURES ================\n');
fprintf('similarity = correlation across pairs of 100-ms partial correlations (Fisher z); mean +- SEM\n');
fprintf('over animal means; p = exact sign-flip over animals. "disatt" divides each r by the square\n');
fprintf('root of the two exposures'' split-half reliabilities.\n');
for g = {'SAL', 'ASV'}
    k = strcmp(T.Group, g{1});
    if ~any(k), continue; end
    Tg = local_sub(T, k); [ua, ~, ai] = unique(Tg.Animal);
    fprintf('\n---- %s: %d sessions, %d animals | units %s, pairs %s ----\n', g{1}, numel(Tg.Session), numel(ua), ...
        local_med(Tg.nUnits), local_med(Tg.nPairs));
    lg = cat(3, Tg.lag{:});
    fprintf('  hours between exposure midpoints (median): A1-B1 %.2f, B1-A2 %.2f, A2-B2 %.2f, A1-B2 %.2f\n', ...
        median(lg(1,2,:)), median(lg(2,3,:)), median(lg(3,4,:)), median(lg(1,4,:)));
    rl = vertcat(Tg.rel{:});
    fprintf('  split-half reliability: A1 %.2f  B1 %.2f  A2 %.2f  B2 %.2f (session medians)\n', median(rl, 1, 'omitnan'));
    for v = {'raw', 'disatt'}
        if strcmp(v{1}, 'raw'), Rs = cat(3, Tg.R{:}); else, Rs = cat(3, Tg.Rd{:}); end
        Zs = atanh(max(min(Rs, 0.999), -0.999));
        fprintf('  [%s] similarity (z):', v{1});
        for q = 1:size(pr, 1)
            x = local_byAnimal(squeeze(Zs(pr(q,1), pr(q,2), :)), ai);
            fprintf('  %s %.3f', nmP{q}, mean(x, 'omitnan'));
        end
        fprintf('\n');
        zz = @(a, b) squeeze(Zs(a, b, :));
        S = struct();
        S.SA2 = zz(1,3) - zz(2,3);
        S.SB2 = zz(2,4) - zz(1,4);
        S.I = S.SA2 + S.SB2;
        [S.bA, S.bB, S.bLag] = local_lagModel(Zs, cat(3, Tg.lag{:}), pr);
        S.bBA = S.bB - S.bA;
        lab = {'SA2', 'S_A2 = z(A1,A2) - z(B1,A2)  familiar map, familiar target'; ...
               'SB2', 'S_B2 = z(B1,B2) - z(A1,B2)  novel map, novel target'; ...
               'I',   'I = S_A2 + S_B2             context specificity (linear drift cancels)'; ...
               'bA',  'b_A  lag model: A1-A2 beyond elapsed time'; ...
               'bB',  'b_B  lag model: B1-B2 beyond elapsed time'; ...
               'bBA', 'b_B - b_A'; ...
               'bLag','b_lag per h'};
        for q = 1:size(lab, 1)
            x = local_byAnimal(S.(lab{q,1}), ai);
            fprintf('     %-58s %+.4f +- %.4f  p = %.3f (n = %d)\n', lab{q,2}, mean(x, 'omitnan'), local_sem(x), ...
                local_signflip(x), sum(isfinite(x)));
        end
        Res.(g{1}).(v{1}) = S; Res.(g{1}).(v{1}).Animal = Tg.Animal; Res.(g{1}).(v{1}).Session = Tg.Session;
        Res.(g{1}).(v{1}).nPrior = Tg.nPrior; Res.(g{1}).(v{1}).DayIdx = Tg.DayIdx;
        % does persistence grow with A's familiarity inside this arm?
        for f = {'SA2', 'bA', 'SB2', 'bB'}
            y = S.(f{1}); xx = Tg.nPrior; ok = isfinite(y) & isfinite(xx);
            if sum(ok) >= 5 && numel(unique(xx(ok))) >= 3
                rho = corr(xx(ok), y(ok), 'type', 'Spearman');
                txt = sprintf('     %-4s vs prior A days: Spearman rho = %+.2f (%d sessions)', f{1}, rho, sum(ok));
                if exist('fitlme', 'file') == 2
                    try
                        Dd = table(y(ok), xx(ok), Tg.Animal(ok), 'VariableNames', {'y','n','Animal'});
                        M = fitlme(Dd, 'y ~ n + (1|Animal)');
                        [b, nm, st] = fixedEffects(M, 'DFMethod', 'satterthwaite'); kk = strcmp(nm.Name, 'n');
                        txt = [txt sprintf(' | LME slope %+.4f per day, p = %.3f', b(kk), st.pValue(kk))]; %#ok<AGROW>
                    catch
                    end
                end
                fprintf('%s\n', txt);
            end
        end
    end
end
% ---- ASV vs SAL ----
if all(isfield(Res, {'SAL', 'ASV'}))
    fprintf('\n---- ASV vs SAL (novel map formed under drug; familiar map drug-free) ----\n');
    for v = {'raw', 'disatt'}
        fprintf('  [%s]\n', v{1});
        for f = {'SB2', 'SA2', 'I', 'bB', 'bA', 'bBA'}
            a = Res.ASV.(v{1}); s = Res.SAL.(v{1});
            [uA, ~, iA] = unique(a.Animal); [uS, ~, iS] = unique(s.Animal);
            mA = local_byAnimal(a.(f{1}), iA); mS = local_byAnimal(s.(f{1}), iS);
            [an, ia, ib] = intersect(uA, uS);
            d = mA(ia) - mS(ib);
            lme = '';
            if exist('fitlme', 'file') == 2
                try
                    D = table([s.(f{1}); a.(f{1})], categorical([repmat({'SAL'}, numel(s.(f{1})), 1); repmat({'ASV'}, numel(a.(f{1})), 1)], {'SAL','ASV'}), ...
                        [s.Animal; a.Animal], 'VariableNames', {'y','Group','Animal'});
                    D = D(isfinite(D.y), :);
                    M = fitlme(D, 'y ~ Group + (1|Animal)');
                    [b, nm, st] = fixedEffects(M, 'DFMethod', 'satterthwaite');
                    kk = strcmp(nm.Name, 'Group_ASV');
                    lme = sprintf(' | LME %+.4f p = %.3f (%d sessions)', b(kk), st.pValue(kk), height(D));
                catch
                end
            end
            fprintf('     %-4s ASV %+.4f  SAL %+.4f | paired ASV - SAL %+.4f +- %.4f p = %.3f (%d animals)%s\n', f{1}, ...
                mean(mA, 'omitnan'), mean(mS, 'omitnan'), mean(d, 'omitnan'), local_sem(d), local_signflip(d), numel(an), lme);
            if exist('fitlme', 'file') == 2
                for cv = {'nPrior', 'DayIdx'}
                    try
                        D = table([s.(f{1}); a.(f{1})], categorical([repmat({'SAL'}, numel(s.(f{1})), 1); repmat({'ASV'}, numel(a.(f{1})), 1)], {'SAL','ASV'}), ...
                            [s.Animal; a.Animal], [s.(cv{1}); a.(cv{1})], 'VariableNames', {'y','Group','Animal','c'});
                        D = D(isfinite(D.y) & isfinite(D.c), :);
                        M = fitlme(D, 'y ~ Group + c + (1|Animal)');
                        [b, nm, st] = fixedEffects(M, 'DFMethod', 'satterthwaite');
                        k1 = strcmp(nm.Name, 'Group_ASV'); k2 = strcmp(nm.Name, 'c');
                        fprintf('          + %-6s ASV - SAL %+.4f p = %.3f | covariate %+.4f p = %.3f\n', cv{1}, b(k1), st.pValue(k1), b(k2), st.pValue(k2));
                    catch
                    end
                end
            end
        end
    end
end
fprintf('\n  Reading: I > 0 is context-specific persistence of exposure covariance into the re-exposure.\n');
fprintf('  A drug effect on stabilization of the new map predicts S_B2 and b_B lower after ASV with S_A2\n');
fprintf('  and b_A unchanged. A2 and B2 are both under drug in ASV sessions, so a drug effect on\n');
fprintf('  expression rather than stabilization would lower both.\n');
end


function as = local_assim(Cf, Ca, Cb, CfH)
% Cf: full-exposure pair correlations (pairs x 4: A1 B1 A2 B2); Ca/Cb: the two halves, CfH: full
% values on the pairs that have both halves. OLS on full predictors; IV with half predictors.
as = struct('alpha', NaN, 'beta', NaN, 'alphaP', NaN, 'betaP', NaN, ...
    'alphaIV', NaN, 'betaIV', NaN, 'alphaPIV', NaN, 'betaPIV', NaN, 'U', NaN, 'nIV', size(Ca, 1));
X = [ones(size(Cf, 1), 1), Cf(:, 1), Cf(:, 2)];
b = X \ Cf(:, 3); as.alpha = b(2); as.beta = b(3);
X = [ones(size(Cf, 1), 1), Cf(:, 2), Cf(:, 1)];
b = X \ Cf(:, 4); as.alphaP = b(2); as.betaP = b(3);
z = atanh(max(min(corr(Cf), 0.999), -0.999)); as.U = z(2, 3) - z(1, 2);
if size(Ca, 1) >= 30
    bA = (local_iv(CfH(:, 3), Ca(:, [1 2]), Cb(:, [1 2])) + local_iv(CfH(:, 3), Cb(:, [1 2]), Ca(:, [1 2]))) / 2;
    bB = (local_iv(CfH(:, 4), Ca(:, [2 1]), Cb(:, [2 1])) + local_iv(CfH(:, 4), Cb(:, [2 1]), Ca(:, [2 1]))) / 2;
    as.alphaIV = bA(2); as.betaIV = bA(3); as.alphaPIV = bB(2); as.betaPIV = bB(3);
end
end


function b = local_iv(y, Xp, Xi)
% just-identified IV: predictors Xp (one half), instruments Xi (other half), intercept in both
P = [ones(numel(y), 1), Xp]; Q = [ones(numel(y), 1), Xi];
b = (Q' * P) \ (Q' * y);
end


function rt = local_rateAssim(lr)
% lr: unit x {A1 B1 A2 B2} x {full, half a, half b}. Delta uses the other half of the baseline
% exposure than the predictor, so regression to the mean cannot produce a coefficient.
rt = struct('gB', NaN, 'gA', NaN, 'gAP', NaN, 'gBP', NaN);
if size(lr, 1) < 5, return; end
g = zeros(3, 2); h = zeros(3, 2);
for s = 1:2
    a = 1 + s; bb = 4 - s;            % half used as predictor / half used in Delta
    y = lr(:, 3, 1) - lr(:, 1, bb);   X = [ones(size(lr, 1), 1), lr(:, 1, a), lr(:, 2, 1)];
    g(:, s) = X \ y;                  % [g0 g_A g_B]
    y = lr(:, 4, 1) - lr(:, 2, bb);   X = [ones(size(lr, 1), 1), lr(:, 2, a), lr(:, 1, 1)];
    h(:, s) = X \ y;                  % [g0 g'_B g'_A]
end
g = mean(g, 2); h = mean(h, 2);
rt.gA = g(2); rt.gB = g(3); rt.gBP = h(2); rt.gAP = h(3);
end


function out = local_assimReport(T, cacheRoot)
out = struct();
n = numel(T.Session);
get = @(f, sub) arrayfun(@(q) T.(f){q}.(sub), (1:n)');
M = struct();
for f = {'alpha','beta','alphaIV','betaIV','alphaP','betaP','alphaPIV','betaPIV','U'}, M.(f{1}) = get('as', f{1}); end
for f = {'gB','gA','gAP','gBP'}, M.(f{1}) = get('rt', f{1}); end
tr = vertcat(T.tRun{:}); M.logRun = log(tr(:, 3) ./ tr(:, 1));
Rz = cellfun(@(r) atanh(max(min(r, 0.999), -0.999)), T.R, 'UniformOutput', false);
M.SA2 = cellfun(@(z) z(1,3) - z(2,3), Rz); zA1A2 = cellfun(@(z) z(1,3), Rz); zB1A2 = cellfun(@(z) z(2,3), Rz);
% ripple reinstatement after B1, from the pair GLM
M.Arip = nan(n, 1);
fp = fullfile(cacheRoot, 'pairCoincidenceGLM.mat');
if exist(fp, 'file')
    try
        P = load(fp, 'F', 'R', 'O');
        md = 'J'; tauC = P.R.Control.ConB1.J.tau; [~, iC] = min(abs(log(P.O.tauGrid) - log(tauC)));
        for q = 1:n
            k = find(strcmp({P.F.Session}, T.Session{q}) & strcmp({P.F.Context}, 'ConB1') & strcmp({P.F.CellSet}, 'all'), 1);
            if isempty(k), continue; end
            Fq = P.F(k); mm = md; if isempty(Fq.(mm)), mm = 'E'; end
            M.Arip(q) = Fq.(mm).th(iC, 2);
        end
    catch err
        fprintf('  (pairCoincidenceGLM.mat not used: %s)\n', err.message);
    end
end
fprintf('\n================ ASSIMILATION OF THE NOVEL MAP INTO THE FAMILIAR MAP ================\n');
fprintf('raw correlation units; mean +- SEM over animal means; p = exact sign-flip over animals\n');
lab = {'betaIV', 'beta  (IV)   c^A2 on c^B1 given c^A1        PRIMARY'; ...
       'beta',   'beta  (OLS)  attenuated by predictor noise'; ...
       'alphaIV','alpha (IV)   retention of A1 in A2'; ...
       'betaPIV','beta'' (IV)   c^B2 on c^A1 given c^B1 (control)'; ...
       'alphaPIV','alpha'' (IV)  retention of B1 in B2'; ...
       'U',      'U = z(B1,A2) - z(A1,B1)  (model-free)'; ...
       'gB',     'g_B   rate: B1-active units gain in A2'; ...
       'gAP',    'g''_A  rate: A1-active units gain in B2 (control)'; ...
       'logRun', 'log(run time A2 / A1)'};
for g = {'SAL', 'ASV'}
    k = strcmp(T.Group, g{1}); [ua, ~, ai] = unique(T.Animal(k));
    fprintf('\n---- %s: %d sessions, %d animals | IV pairs median %g ----\n', g{1}, sum(k), numel(ua), ...
        median(arrayfun(@(q) T.as{q}.nIV, find(k))));
    for q = 1:size(lab, 1)
        x = accumarray(ai, M.(lab{q,1})(k), [], @(z) mean(z, 'omitnan'));
        fprintf('  %-52s %+.4f +- %.4f  p = %.3f (n = %d)\n', lab{q,2}, mean(x, 'omitnan'), local_sem(x), local_signflip(x), sum(isfinite(x)));
    end
end
% ---- ASV vs SAL
fprintf('\n---- ASV vs SAL ----\n');
G = categorical(T.Group, {'SAL', 'ASV'});
for q = 1:size(lab, 1)
    y = M.(lab{q,1});
    kS = strcmp(T.Group, 'SAL'); kA = strcmp(T.Group, 'ASV');
    [uS, ~, iS] = unique(T.Animal(kS)); [uA, ~, iA] = unique(T.Animal(kA));
    mS = accumarray(iS, y(kS), [], @(z) mean(z, 'omitnan')); mA = accumarray(iA, y(kA), [], @(z) mean(z, 'omitnan'));
    [an, ia, ib] = intersect(uA, uS); d = mA(ia) - mS(ib);
    fprintf('  %-8s paired ASV - SAL %+.4f +- %.4f p = %.3f (%d animals)\n', lab{q,1}, mean(d, 'omitnan'), local_sem(d), local_signflip(d), numel(an));
    if exist('fitlme', 'file') ~= 2, continue; end
    for cv = {'', 'DayIdx', 'nPrior', 'logRun'}
        try
            D = table(y, G, T.Animal, T.DayIdx, T.nPrior, M.logRun, 'VariableNames', {'y','Group','Animal','DayIdx','nPrior','logRun'});
            if isempty(cv{1}), fm = 'y ~ Group + (1|Animal)'; else, fm = ['y ~ Group + ' cv{1} ' + (1|Animal)']; end
            D = D(isfinite(D.y), :);
            if ~isempty(cv{1}), D = D(isfinite(D.(cv{1})), :); end
            Mx = fitlme(D, fm);
            [b, nm, st] = fixedEffects(Mx, 'DFMethod', 'satterthwaite'); kk = strcmp(nm.Name, 'Group_ASV');
            tag = 'LME'; if ~isempty(cv{1}), tag = ['+ ' cv{1}]; end
            star = ''; if strcmp(lab{q,1}, 'betaIV') && strcmp(cv{1}, 'DayIdx'), star = '   <- PRIMARY'; end
            fprintf('           %-9s ASV - SAL %+.4f p = %.3f (%d sessions)%s\n', tag, b(kk), st.pValue(kk), height(D), star);
        catch
        end
    end
end
% ---- interaction: does A2 move toward B1 and away from A1 differently by arm?
if exist('fitlme', 'file') == 2
    try
        L = table([zA1A2; zB1A2], [G; G], categorical([repmat({'A1A2'}, n, 1); repmat({'B1A2'}, n, 1)]), ...
            [T.Animal; T.Animal], [T.Session; T.Session], [T.DayIdx; T.DayIdx], ...
            'VariableNames', {'z','Group','Type','Animal','Session','DayIdx'});
        Mx = fitlme(L, 'z ~ Group * Type + DayIdx + (1|Animal) + (1|Animal:Session)');
        [b, nm, st] = fixedEffects(Mx, 'DFMethod', 'satterthwaite'); kk = strcmp(nm.Name, 'Group_ASV:Type_B1A2');
        fprintf('\n  interaction z ~ Group x {A1-A2, B1-A2} + DayIdx: %+.4f p = %.3f\n', b(kk), st.pValue(kk));
        fprintf('  (negative = on ASV days A2 sits closer to A1 and farther from B1 than on saline days)\n');
    catch err
        fprintf('  interaction LME not fitted (%s)\n', err.message);
    end
end
% ---- per animal, with drug order
fprintf('\n---- per animal (arm means; order = arm of the animal''s first included session) ----\n');
fprintf('  %-5s %-6s %8s %8s %8s %8s %8s | %8s %8s %8s %8s %8s\n', 'anim', 'first', 'bIV SAL', 'bIV ASV', 'U SAL', 'U ASV', 'SA2 SAL', 'SA2 ASV', 'gB SAL', 'gB ASV', 'nS', 'nA');
for a = unique(T.Animal)'
    k = strcmp(T.Animal, a{1}); [~, o] = min(T.DayIdx + 1e6 * ~k);
    kS = k & strcmp(T.Group, 'SAL'); kA = k & strcmp(T.Group, 'ASV');
    f = @(v, kk) mean(v(kk), 'omitnan');
    fprintf('  %-5s %-6s %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f | %+8.3f %+8.3f %+8.3f %8d %8d\n', a{1}, T.Group{o}, ...
        f(M.betaIV, kS), f(M.betaIV, kA), f(M.U, kS), f(M.U, kA), f(M.SA2, kS), f(M.SA2, kA), f(M.gB, kS), f(M.gB, kA), sum(kS), sum(kA));
end
% ---- ripple link
if any(isfinite(M.Arip))
    fprintf('\n---- post-B1 ripple reinstatement (pair GLM, joint model) vs assimilation, within arm ----\n');
    for g = {'SAL', 'ASV'}
        k = strcmp(T.Group, g{1}) & isfinite(M.Arip) & isfinite(M.betaIV);
        if sum(k) < 5, continue; end
        rho = corr(M.Arip(k), M.betaIV(k), 'type', 'Spearman');
        fprintf('  %s: Spearman(A_ripple, beta_IV) = %+.2f (%d sessions)\n', g{1}, rho, sum(k));
    end
    if exist('fitlme', 'file') == 2
        try
            D = table(M.betaIV, G, T.Animal, M.Arip, 'VariableNames', {'y','Group','Animal','A'});
            D = D(isfinite(D.y) & isfinite(D.A), :);
            Mx = fitlme(D, 'y ~ A * Group + (1|Animal)');
            [b, nm, st] = fixedEffects(Mx, 'DFMethod', 'satterthwaite');
            for nmq = {'A', 'A:Group_ASV'}
                kk = strcmp(nm.Name, nmq{1});
                if any(kk), fprintf('  LME beta_IV ~ A x Group: %-12s %+.4f p = %.3f\n', nmq{1}, b(kk), st.pValue(kk)); end
            end
        catch
        end
    end
    fprintf('  H predicts a positive slope in SAL (stronger offline reinstatement, more assimilation) and a\n');
    fprintf('  flat one in ASV, where reinstatement is present but assimilation is not.\n');
end
out = M; out.Session = T.Session; out.Group = T.Group; out.Animal = T.Animal;
end


function out = local_yieldReport(T, O)
% CHECKS 1-2: is the all-four-exposure inclusion rule selective, and where does the loss fall?
out = struct();
fprintf('\n================ UNIT YIELD AND SELECTIVE DROPOUT (checks 1-2) ================\n');
fprintf('A unit needs >= %d spikes in the running/theta part of EVERY exposure.\n', O.minSpk);
fprintf('  nPYR     pyramidal units in the session\n');
fprintf('  pass_X   units reaching the threshold in exposure X alone\n');
fprintf('  all4     units reaching it in all four (the analysed set)\n');
fprintf('  retB1    all4 / pass_B1: the fraction of ConB1-active units the rule keeps (check 1)\n');
fprintf('  matched  pass_X recomputed on the first nBinMin bins of every exposure, so an\n');
fprintf('           exposure that was shorter cannot lose units by duration alone\n');
nmE = {'A1', 'B1', 'A2', 'B2'};
for g = {'SAL', 'ASV'}
    k = strcmp(T.Group, g{1});
    if ~any(k), continue; end
    Tg = local_sub(T, k); [ua, ~, ai] = unique(Tg.Animal);
    P = vertcat(Tg.nPass{:}); Pm = vertcat(Tg.nPassM{:});
    den = max(Tg.nPYR, eps);
    ret = Tg.nUnits ./ max(P(:,2), eps); retAll = Tg.nUnits ./ den;
    fr = P ./ den; frm = Pm ./ den;
    fprintf('\n---- %s: %d sessions, %d animals ----\n', g{1}, numel(Tg.Session), numel(ua));
    fprintf('   session means: nPYR %.1f | pass A1 %.1f B1 %.1f A2 %.1f B2 %.1f | all4 %.1f\n', ...
        mean(Tg.nPYR), mean(P, 1), mean(Tg.nUnits));
    fprintf('   pass fraction of nPYR:            A1 %.3f B1 %.3f A2 %.3f B2 %.3f\n', mean(fr, 1));
    fprintf('   pass fraction, duration-matched:  A1 %.3f B1 %.3f A2 %.3f B2 %.3f\n', mean(frm, 1));
    x = local_byAnimal(ret, ai); y = local_byAnimal(retAll, ai);
    fprintf('   retB1 %.3f +- %.3f | all4/nPYR %.3f +- %.3f (animal means)\n', ...
        mean(x, 'omitnan'), local_sem(x), mean(y, 'omitnan'), local_sem(y));
    out.(g{1}) = struct('Animal', {Tg.Animal}, 'Session', {Tg.Session}, 'nPYR', Tg.nPYR, ...
        'nAll4', Tg.nUnits, 'nPass', P, 'nPassM', Pm, 'retB1', ret, 'retAll', retAll, ...
        'frac', fr, 'fracM', frm);
end
if all(isfield(out, {'SAL', 'ASV'}))
    fprintf('\n---- ASV vs SAL (paired sign-flip over animals) ----\n');
    it = {'retB1', 'all4 / pass_B1'; 'retAll', 'all4 / nPYR'};
    for q = 1:size(it, 1)
        [d, an] = local_pairArm(out.ASV.(it{q,1}), out.ASV.Animal, out.SAL.(it{q,1}), out.SAL.Animal);
        fprintf('   %-18s ASV - SAL %+.3f +- %.3f  p = %.3f (%d animals)\n', it{q,2}, ...
            mean(d, 'omitnan'), local_sem(d), local_signflip(d), numel(an));
    end
    for c = 1:4
        [d, an] = local_pairArm(out.ASV.frac(:,c), out.ASV.Animal, out.SAL.frac(:,c), out.SAL.Animal);
        [dm, ~]  = local_pairArm(out.ASV.fracM(:,c), out.ASV.Animal, out.SAL.fracM(:,c), out.SAL.Animal);
        fprintf('   pass_%-13s ASV - SAL %+.3f +- %.3f  p = %.3f (%d animals) | matched %+.3f p = %.3f\n', ...
            nmE{c}, mean(d, 'omitnan'), local_sem(d), local_signflip(d), numel(an), ...
            mean(dm, 'omitnan'), local_signflip(dm));
    end
    fprintf('\n  Reading: equal retB1 between arms leaves the all-four rule a stated restriction.\n');
    fprintf('  Lower retention after ASV is itself a result, and a loss concentrated in A2 / B2\n');
    fprintf('  points to expression under drug rather than to yield over the day.\n');
end
end


function out = local_instabReport(T)
% CHECK 5: per-pair stability, against the sampling noise of the correlation estimates.
out = struct();
fprintf('\n================ PER-PAIR STABILITY (check 5) ================\n');
fprintf('Across the pairs of one session, between two exposures:\n');
fprintf('  flip    fraction of pairs whose partial correlation changes sign\n');
fprintf('  flipS   the same among pairs in the top quartile of |c| in the first exposure\n');
fprintf('  mad     median |c(second) - c(first)|\n');
fprintf('  madExp  mad expected from sampling noise alone, sqrt(s_X^2 + s_Y^2)/2, s = median\n');
fprintf('          |difference| between the two interleaved halves of that exposure\n');
fprintf('  ratio   mad / madExp; about 1 = within sampling noise, above 1 = pairs reorganised\n');
fprintf('  keep    of the top-quartile pairs, the fraction still above the second exposure''s\n');
fprintf('          median |c| (chance 0.5)\n');
fprintf('B1-B2 is the novel pair; A1-A2 is the within-session familiar reference.\n');
V = {'flip', 'flipS', 'mad', 'madExp', 'ratio', 'keep'};
for g = {'SAL', 'ASV'}
    k = strcmp(T.Group, g{1});
    if ~any(k), continue; end
    Tg = local_sub(T, k); [ua, ~, ai] = unique(Tg.Animal);
    fprintf('\n---- %s: %d sessions, %d animals ----\n', g{1}, numel(Tg.Session), numel(ua));
    S = struct('Animal', {Tg.Animal}, 'Session', {Tg.Session});
    S.nUnits = Tg.nUnits; S.nPairs = Tg.nPairs; S.tRun = cellfun(@(t) min(t), Tg.tRun);
    for f = {'B1B2', 'A1A2', 'A1B1'}
        for v = V, S.([f{1} '_' v{1}]) = cellfun(@(x) x.(f{1}).(v{1}), Tg.ins); end
    end
    fprintf('   yield: units %s, shortest exposure %.0f s (median)\n', local_med(Tg.nUnits), median(S.tRun));
    for f = {'B1B2', 'A1A2'}
        fprintf('   %-5s', f{1});
        for v = V
            x = local_byAnimal(S.([f{1} '_' v{1}]), ai);
            fprintf('  %s %.3f+-%.3f', v{1}, mean(x, 'omitnan'), local_sem(x));
        end
        r = local_byAnimal(S.([f{1} '_ratio']), ai) - 1;
        fprintf('\n          ratio vs 1: p = %.3f (%d animals)\n', local_signflip(r), sum(isfinite(r)));
    end
    out.(g{1}) = S;
end
if all(isfield(out, {'SAL', 'ASV'}))
    fprintf('\n---- ASV vs SAL (paired sign-flip over animals) ----\n');
    for f = {'B1B2', 'A1A2'}
        for v = V
            fn = [f{1} '_' v{1}];
            [d, an] = local_pairArm(out.ASV.(fn), out.ASV.Animal, out.SAL.(fn), out.SAL.Animal);
            fprintf('   %-14s ASV - SAL %+.4f +- %.4f  p = %.3f (%d animals)\n', fn, ...
                mean(d, 'omitnan'), local_sem(d), local_signflip(d), numel(an));
        end
    end
    % the ratio's SENSITIVITY grows with yield, so a matched ratio is weak evidence when the
    % arms differ in units or running time: both enter as covariates.
    if exist('fitlme', 'file') == 2
        for fn = {'B1B2_ratio', 'B1B2_mad', 'A1A2_ratio'}
            try
                a = out.ASV; s = out.SAL;
                D = table([s.(fn{1}); a.(fn{1})], ...
                    categorical([repmat({'SAL'}, numel(s.(fn{1})), 1); repmat({'ASV'}, numel(a.(fn{1})), 1)], {'SAL','ASV'}), ...
                    [s.Animal; a.Animal], [s.nUnits; a.nUnits], [s.tRun; a.tRun], ...
                    'VariableNames', {'y','Group','Animal','nU','tRun'});
                D = D(all(isfinite([D.y D.nU D.tRun]), 2), :);
                M = fitlme(D, 'y ~ Group + nU + tRun + (1|Animal)');
                [b, nm2, st] = fixedEffects(M, 'DFMethod', 'satterthwaite');
                k1 = strcmp(nm2.Name, 'Group_ASV'); k2 = strcmp(nm2.Name, 'nU');
                fprintf('   %-14s yield-adjusted ASV - SAL %+.4f p = %.3f | per unit %+.4f p = %.3f (%d sessions)\n', ...
                    fn{1}, b(k1), st.pValue(k1), b(k2), st.pValue(k2), height(D));
            catch
            end
        end
    end
    fprintf('\n  Reading: a similar S_B2 between arms together with a larger B1-B2 ratio after ASV\n');
    fprintf('  would mean the across-pair similarity holds while individual pairs reorganise.\n');
    fprintf('  Similar ratios leave the stability claim intact for the units the rule kept.\n');
    fprintf('  CAVEAT: ratio is unbiased (1.00 under no true change at every yield simulated) but\n');
    fprintf('  its SENSITIVITY grows with units and running time: a reorganisation that gives\n');
    fprintf('  ratio 2.4 at 20 units / 3000 bins gives only 1.3 at 12 units / 400 bins. Read a\n');
    fprintf('  null arm difference next to the yields printed above, not on its own.\n');
end
end


function ins = local_instab(C, C1, C2)
% per-pair stability between exposures, with the sampling-noise expectation taken from the two
% interleaved halves of each exposure. C: pairs x 4 (A1 B1 A2 B2); C1, C2: the same per half.
ins = struct();
nm = {'A1', 'B1', 'A2', 'B2'};
cc = {'A1A2', 1, 3; 'B1B2', 2, 4; 'A1B1', 1, 2};
for q = 1:size(cc, 1), ins.(cc{q,1}) = local_pairStab(C(:, cc{q,2}), C(:, cc{q,3})); end
for q = 1:4, ins.(['h' nm{q}]) = local_pairStab(C1(:, q), C2(:, q)); end
for q = 1:size(cc, 1)
    s1 = ins.(['h' nm{cc{q,2}}]).mad; s2 = ins.(['h' nm{cc{q,3}}]).mad;
    e = sqrt(s1^2 + s2^2) / 2;            % each half holds ~half the data: ~sqrt(2) too wide
    ins.(cc{q,1}).madExp = e;
    ins.(cc{q,1}).ratio = ins.(cc{q,1}).mad / max(e, eps);
end
end


function st = local_pairStab(a, b)
st = struct('n', 0, 'flip', NaN, 'flipS', NaN, 'mad', NaN, 'keep', NaN, 'madExp', NaN, 'ratio', NaN);
ok = isfinite(a) & isfinite(b);
if sum(ok) < 10, return; end
a = a(ok); b = b(ok); st.n = numel(a);
st.flip = mean(sign(a) ~= sign(b));
st.mad  = median(abs(b - a));
s = abs(a) >= quantile(abs(a), 0.75);
if sum(s) >= 5
    st.flipS = mean(sign(a(s)) ~= sign(b(s)));
    st.keep  = mean(abs(b(s)) >= median(abs(b)));
end
end


function [d, an] = local_pairArm(yA, anA, yS, anS)
% paired ASV - SAL over the animals contributing to both arms, on animal means
[uA, ~, iA] = unique(anA); [uS, ~, iS] = unique(anS);
mA = local_byAnimal(yA, iA); mS = local_byAnimal(yS, iS);
[an, ia, ib] = intersect(uA, uS);
d = mA(ia) - mS(ib);
end


function Tg = local_sub(T, k)
fn = fieldnames(T); Tg = struct();
for q = 1:numel(fn), Tg.(fn{q}) = T.(fn{q})(k); end
end


function [bA, bB, bLag] = local_lagModel(Zs, lag, pr)
n = size(Zs, 3); bA = nan(n, 1); bB = nan(n, 1); bLag = nan(n, 1);
for s = 1:n
    y = arrayfun(@(q) Zs(pr(q,1), pr(q,2), s), (1:size(pr, 1))');
    L = arrayfun(@(q) lag(pr(q,1), pr(q,2), s), (1:size(pr, 1))');
    X = [ones(6, 1), double(ismember(pr, [1 3], 'rows')), double(ismember(pr, [2 4], 'rows')), L];
    ok = isfinite(y);
    if sum(ok) < 5, continue; end
    b = X(ok, :) \ y(ok);
    bA(s) = b(2); bB(s) = b(3); bLag(s) = b(4);
end
end


function M = local_chunk(mask, w)
% split every [start end] row into consecutive pieces of at most w s, so alternating rows give
% two interleaved halves of the same exposure
M = zeros(0, 2);
for k = 1:size(mask, 1)
    e = mask(k,1):w:mask(k,2);
    if e(end) < mask(k,2), e(end+1) = mask(k,2); end %#ok<AGROW>
    M = [M; e(1:end-1)', e(2:end)']; %#ok<AGROW>
end
M = M(diff(M, 1, 2) > 0, :);
end


function x = local_byAnimal(v, ai)
x = accumarray(ai(:), v(:), [], @(z) mean(z, 'omitnan'));
end


function s = local_med(x)
s = sprintf('%g (%g-%g)', median(x), min(x), max(x));
end


function e = local_sem(x)
x = x(isfinite(x));
if numel(x) < 2, e = NaN; else, e = std(x) / sqrt(numel(x)); end
end


function p = local_signflip(x)
x = x(isfinite(x)); n = numel(x);
if n < 2, p = NaN; return; end
if n > 20
    s = sign(randn(n, 5000)); s(:,1) = 1;
else
    s = double(dec2bin(0:2^n - 1) - '0')' * 2 - 1;
end
m = abs(x(:)' * s / n);
p = max(mean(m >= abs(mean(x)) - 1e-12), 2 / 2^n);
end


%% =============================================================================================
%  shared with pairCoincidenceGLM.m (copied unchanged)
%% =============================================================================================
function Z = local_binMask(ts, mask, binw)
Z = zeros(0, numel(ts));
if isempty(mask), return; end
edges = [];
for k = 1:size(mask, 1), edges = [edges; (mask(k,1):binw:mask(k,2)-binw)']; end %#ok<AGROW>
if numel(edges) < 2, return; end
Z = zeros(numel(edges), numel(ts));
for u = 1:numel(ts)
    tt = ts{u}; tt = tt(tt >= edges(1) & tt < edges(end) + binw);
    if isempty(tt), continue; end
    [~, ~, b] = histcounts(tt, [edges; edges(end) + binw]);
    b = b(:); v = b > 0; v(v) = tt(v) < edges(b(v)) + binw;
    if any(v), Z(:, u) = accumarray(b(v), 1, [numel(edges) 1]); end
end
end


function C = local_partialCorr(Z, U)
% pair correlation of binned counts after removing, from each unit, its regression on the
% summed count of the OTHER units in Z (all recorded PYR units)
s = sum(Z, 2); n = numel(U);
R = zeros(size(Z, 1), n);
for q = 1:n
    z = Z(:, U(q)); x = s - z;
    A = [ones(numel(z), 1), x];
    R(:, q) = z - A * (A \ z);
end
sd = std(R, 0, 1); C = nan(n);
ok = sd > 0;
if sum(ok) >= 2, C(ok, ok) = corr(R(:, ok)); end
end


function [X, dur, tOn, nU, ts] = local_ripCounts(fld, RP)
X = []; dur = []; tOn = []; nU = 0; ts = {};
sf = dir(fullfile(fld, '*spikes.cellinfo.mat')); cf = dir(fullfile(fld, '*cell_metrics.cellinfo.mat'));
if isempty(sf) || isempty(cf), return; end
S = load(fullfile(sf(1).folder, sf(1).name), 'spikes'); sp = S.spikes;
C = load(fullfile(cf(1).folder, cf(1).name), 'cell_metrics'); cm = C.cell_metrics;
ts = cellfun(@(x) sort(reshape(x(~isnan(x)), [], 1)), sp.times(:), 'UniformOutput', false);
n = numel(ts);
clu = []; if isfield(sp, 'UID') && numel(sp.UID) == n, clu = sp.UID(:); end
if isfield(sp, 'cluID') && numel(sp.cluID) == n, clu = sp.cluID(:); end
loc = matchCellMetrics(cm, (1:n)', clu, n); okL = loc > 0;
typ = repmat({''}, n, 1); ct = cellstr(string(cm.putativeCellType(:))); typ(okL) = ct(loc(okL));
idx = find(contains(typ, 'Pyramidal'));
if isempty(idx), ts = {}; return; end
on = RP.ts1(:); off = RP.ts1Off(:);
[on, o] = sort(on); off = off(o);
keep = off > on & [true; on(2:end) > off(1:end-1)];
on = on(keep); off = off(keep);
if numel(on) < 2, ts = {}; return; end
edges = reshape([on'; off'], [], 1);
X = zeros(numel(on), numel(idx));
for u = 1:numel(idx)
    c = histcounts(ts{idx(u)}, edges);
    X(:,u) = c(1:2:end)';
end
dur = off - on; tOn = on; nU = numel(idx); ts = ts(idx);
end


function s = local_name(fld)
[~, s] = fileparts(fld);
end


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


function g = local_group(RP, ev, sessionSet)
% Group label of one TransPlasticity session: 'SAL', 'ASV' or '' (not used).
%   sessionSet 'all'   : every SAL session is 'SAL' and every ASV session is 'ASV' (first and repeat).
%   sessionSet 'first' : the earlier design -- SAL sessions before the animal's first ASV, and the
%                        first-ASV session; SAL after the first ASV and repeat ASV are ''.
% A session needs an infusion time (RIPPLES.meta, else Infusion_start, else ConA1_end) either way.
if nargin < 3 || isempty(sessionSet), sessionSet = 'all'; end
m = RP.meta; g = '';
inf_ = m.infusionSec;
if ~isfinite(inf_) && nargin > 1 && ~isempty(ev)
    k = find(strcmp(ev.description, 'Infusion_start'), 1);
    if isempty(k), k = find(strcmp(ev.description, 'ConA1_end'), 1); end
    if ~isempty(k), inf_ = ev.time(k); end
end
if ~isfinite(inf_), return; end
if strcmp(sessionSet, 'all')
    if ismember(m.condition, {'SAL', 'ASV'}), g = m.condition; end
    return
end
noASV = ~isfield(m, 'firstASV') || isempty(m.firstASV.session);
off = RP.ts2(:) - RP.ts1(:); off = median(off(isfinite(off)));
if noASV || isempty(off) || ~isfinite(off), after = false; else, after = inf_ + off > 0; end
isFirst = ~noASV && strcmp(m.folder, m.firstASV.session);
switch m.condition
    case 'SAL'
        if ~after, g = 'SAL'; end
    case 'ASV'
        if isFirst, g = 'ASV'; end
end
end


function [mask, src] = local_runMask(fld, a, b, vmin)
% 1 s windows inside [a b] whose mean speed > vmin; src = '' when there is no tracking file
mask = zeros(0, 2); src = '';
d = dir(fullfile(fld, '*.Behavior.mat'));
if isempty(d), return; end
[~, o_] = sort(endsWith({d.name}, 'Tracking.Behavior.mat')); d = d(o_);   % <session>.Behavior.mat first
try
    L = load(fullfile(d(1).folder, d(1).name), 'behavior'); bh = L.behavior;
    t = bh.timestamps(:); x = bh.position.x(:); y = bh.position.y(:);
    lin = []; if isfield(bh.position, 'lin'), lin = bh.position.lin(:); end
catch ME
    warning('react:track', '%s: %s', d(1).name, ME.message); return
end
inExp = t >= a & t <= b;
useLin = numel(lin) == numel(t) && sum(isfinite(lin(inExp))) > 10;   % linearised position (cm) when available
if useLin, ok = isfinite(t) & isfinite(lin); else, ok = isfinite(t) & isfinite(x) & isfinite(y); end
t = t(ok); x = x(ok); y = y(ok); if useLin, lin = lin(ok); end
if numel(t) < 3, return; end
[t, iu] = unique(t); x = x(iu); y = y(iu); if useLin, lin = lin(iu); end
if useLin
    v = [0; abs(diff(lin)) ./ diff(t)];
    v(v > 200) = NaN;                          % jumps where the linearisation wraps (track ends / arms)
else
    v = [0; hypot(diff(x), diff(y)) ./ diff(t)];
end
dt = median(diff(t)); v = movmean(v, max(1, round(0.25 / dt)), 'omitnan');   % 250 ms smoothing
w0 = (floor(a):1:ceil(b) - 1)'; w0 = w0(w0 >= a & w0 + 1 <= b);
if isempty(w0), src = sprintf('tracking speed > %g', vmin); return; end
[~, ~, bi] = histcounts(t, [w0; w0(end) + 1]);
okb = bi > 0;
ms = accumarray(bi(okb), v(okb), [numel(w0) 1], @mean, NaN);
% units: vmin is cm/s. Position stored in metres (or normalised track units) gives
% speeds < ~1; detect that from this exposure's own speeds and rescale.
scale = 1;
if ~useLin                                     % x/y: detect metres from this exposure's speeds
    v95 = prctile(ms(isfinite(ms)), 95);
    if isfinite(v95) && v95 < 1, scale = 100; end
end
k = find(ms * scale > vmin);
mask = [w0(k), w0(k) + 1];
if useLin, src = sprintf('tracking position.lin speed > %g cm/s', vmin);
else, src = sprintf('tracking x/y speed > %g cm/s', vmin); end
if scale ~= 1, src = [src ' (x/y x100: metres)']; end
end


function [mask, src] = local_thetaMask(RP, a, b, ratio)
% 1 s windows inside [a b] with theta/delta power ratio > ratio; whole [a b] when no LFP
mask = [a b]; src = 'whole exposure (no LFP)';
m = RP.meta;
if ~isfield(m, 'lfpFile') || ~exist(m.lfpFile, 'file') || ~isfield(m, 'channel1based'), return; end
[fd, bn] = fileparts(m.lfpFile); xf = fullfile(fd, [bn '.xml']);
if ~exist(xf, 'file'), return; end
try
    xml = LoadXml(xf); fs = m.fs;
    x = double(LoadBinary(m.lfpFile, 'nChannels', xml.nChannels, 'channels', m.channel1based, ...
        'frequency', fs, 'start', a, 'duration', b - a));
catch ME
    warning('react:lfp', '%s: %s', m.lfpFile, ME.message); return
end
fs = round(fs); nw = floor(numel(x) / fs);
if nw < 1, mask = zeros(0, 2); return; end
X = reshape(x(1:nw*fs), fs, nw);
F = abs(fft(X - mean(X, 1))).^2; f = (0:fs-1)';      % 1 s windows: 1 Hz resolution
th = sum(F(f >= 6 & f <= 10, :), 1); de = sum(F(f >= 1 & f <= 4, :), 1);
k = find(th ./ de > ratio);
mask = [a + (k(:) - 1), a + k(:)];
src = sprintf('LFP theta/delta > %g', ratio);
end


function d = local_sessionDate(sname)
% Recording date from the session folder name, as a datenum. Names look like CP12_250604_075757
% or CP5_230922: the first six-digit field after an underscore is YYMMDD.
d = NaN;
t = regexp(char(sname), '_(\d{6})(?:_|$)', 'tokens', 'once');
if isempty(t), return; end
x = t{1};
try
    d = datenum(2000 + str2double(x(1:2)), str2double(x(3:4)), str2double(x(5:6)));
catch
    d = NaN;
end
end