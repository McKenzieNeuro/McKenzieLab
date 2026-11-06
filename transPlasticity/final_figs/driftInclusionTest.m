%% driftInclusionTest.m
% Does including data from a PRESCRIBED period improve prediction of future ripples?
%
% THE CLAIM THIS IS BUILT TO SUPPORT
%   "Prediction of the fast and slow components of future ripple content was possible only when
%    data from the context explorations was included in the training set."
%
% WHY NOT THE NOVELTY SCAN
% driftNoveltyScan.m searches for epochs at which nothing in the past fits. That is the wrong
% instrument for this claim and it is strictly weaker, for two reasons. The exposure times are
% known by design, so rediscovering them spends all the power of the design on a search that was
% never necessary. And an unsupervised scan has to be calibrated against a null over all possible
% epochs, which is dominated by brain state: in the control sessions N(t) correlates +0.62 with
% ripple rate, so the largest excursions are sleep, not experience. Here the epoch is named in
% advance and the only question is what adding it buys.
%
% THE TEST
% Build a baseline library L0 of home-cage sub-windows that lie BEFORE the context, and holding
% that library fixed, add exactly one expert. Two versions of that one expert:
%       e_ctx   fitted on the context exploration itself
%       e_hc    fitted on a home-cage block of the same duration, taken
%               from immediately before the context, so it is at least as close to the target
%       G_ctx(t) = [ D(L0; t) - D(L0 + e_ctx; t) ] / [ D_t(split) - D_t(oracle) ]
%       G_hc(t)  = same with e_hc
%       inclusion effect = mean_t [ G_ctx(t) - G_hc(t) ]
% The library is nested either way, so both gains are non-negative and neither means anything on
% its own: adding any expert to any library helps. What the contrast removes is precisely that.
% Library size, nesting, duration, bin count and position relative to the target are all held
% equal, and the only thing that differs is whether the animal was exploring. That is the
% comparison the claim needs, and nothing weaker will carry it. Firing rate is NOT matched, since
% running and resting differ in rate by construction; the free gain removes the level, and the
% ConB1-against-ConA1 contrast, where both sides are explorations, is what matches the state.
%
% SLOW AND FAST, SEPARATELY
%   slow  the mean rate vector over the target window: which cells fire at what rate.
%         Scored by Poisson deviance with a free multiplicative gain, so a level difference
%         between the context and the home cage cannot masquerade as prediction.
%   fast  the ripple-to-ripple covariance within the target window, after the radial direction is
%         removed. Scored by the fraction of that variance captured by the added expert's own
%         basis, against the window's oracle basis and split-half floor.
% The two can dissociate: a context can set which cells are active without setting which of them
% co-fire, or the reverse. Reporting one of them alone would hide that.
%
% ConA1 (familiar, pre-exposed by design) and ConB1 (novel) are each tested against their own
% matched home-cage twin, so the familiar context doubles as an internal specificity control.
%
% TIME COURSE
% Everything is reported in bins of time since the context ended, because an inclusion effect
% that is present at 20 min and gone at 3 h is a different claim from one that persists, and the
% reactivation time course in the same sessions is already binned this way.

CONTROLS = { ...
    'R:\HuszarR\e13_26m1_211019', [15149.6 16749.2], [7605.7 9322]; ...
    'R:\HuszarR\e13_26m1_211116', [13608 13939],     [7469.8 8707.2]; ...
    'R:\HuszarR\e13_26m1_211119', [14911.9 16446],   [6181.3 7572]; ...
    'R:\HuszarR\e15_13f1_220118', [13137 15046],     [5636.7 7700]; ...
    'R:\HuszarR\e15_13f1_220119', [21219 22700],     [6836.6 8860]; ...
    'R:\HuszarR\e15_13f1_220120', [18017.9 19801],   [9610.5 11950]; ...
    'R:\HuszarR\e16_3m2_211211',  [9584.6 11203],    [5294.1 6511]; ...
    'R:\HuszarR\e16_3m2_211212',  [16416.3 17731],   [6109.3 7845]; ...
    'R:\HuszarR\e16_3m2_211213',  [13291 15323],     [7567.9 9359]};   % folder, NOVEL, FAMILIAR
% For the TransPlasticity sessions, build the same three columns from the .evt.sti markers
% (ConA1_start/_end and ConB1_start/_end, as compareReactivation.m reads them) and the rest of
% this script is unchanged. The arm comparison is then the difference of inclusion effects.

SUB_WIN     = 600;        % s, target window width
GUARD       = 600;        % s after the context ends before targets begin, against adjacency
TGT_MAX     = 4*3600;     % s of home cage after the context used as targets
BIN         = 3600;       % s, bin width for the time course
LIB_MIN     = 6;          % baseline library needs this many sub-windows
MIN_RIP_SUB = 25;         % ripples a sub-window needs
BEH_BIN     = 0.25;       % s, bin width for the behaviour-time expert (classical templates use
                          % 0.1 to 1 s; the covariance basis is what this width sets)
MIN_BIN_CTX = 200;        % bins a context needs to supply an expert (50 s at 0.25 s)
MIN_SPK_CTX = 200;        % spikes a context needs across all cells
MIN_UNITS   = 20;
GEO_RANK    = 4;          % covariance directions compared (see driftReverseControl.m section 10)
N_SHUF      = 20;         % cell-identity permutations for the GLM control
LAG_MIN     = 0.25;       % h, closest home-cage pseudo-template used for the decay curve
PSEUDO_STEP = 1800;       % s, spacing of the home-cage pseudo-templates across the session
SEED        = 13;
FREE_GAIN   = true;       % one multiplicative gain per predictor, fitted on the target
SAVE_OUT    = true;

rng(SEED, 'twister');
rows = {};
for i = 1:size(CONTROLS, 1)
    fld = CONTROLS{i,1}; nov = CONTROLS{i,2}; fam = CONTROLS{i,3};
    rf = fullfile(fld, 'rippleEvents_control.mat');
    if ~exist(rf, 'file'), warning('inc:file', 'missing %s', rf); continue; end
    L = load(rf, 'RIPPLES'); RP = L.RIPPLES;
    [X, dur, tOn, nU, ts] = local_ripCounts(fld, RP);
    if isempty(X) || nU < MIN_UNITS
        fprintf('%-34s skipped: %d PYR units\n', local_name(fld), nU); continue
    end
    CTX = [fam; nov];                                   % every context interval in this session
    inCtx = false(size(tOn));
    for c = 1:size(CTX, 1), inCtx = inCtx | (tOn >= CTX(c,1) & tOn < CTX(c,2)); end

    % ---- the drift decay curve, measured where no experience distinguishes anything ---------
    % The cross-template interaction cannot separate specificity from temporal proximity, and
    % the adversarial-subset shortcut fails in this design: ConB1 is the last context, so no
    % target after it has a nearer competitor, and the adversarial windows are all ConA1-own.
    % That subset cannot cancel the template main effect and simply reports it.
    % So measure the decay directly. Every home-cage sub-window in the session is a pseudo-
    % template: experience-neutral, same kind of object, spanning the whole lag range in both
    % directions. Regressing prediction error on |lag| over those gives the drift-only
    % expectation at any lag, and the two context templates are then scored as RESIDUALS from
    % that curve. Whatever is left is what proximity does not explain.
    % The pseudo-templates are built per context, inside the loop below, because they must match
    % that context's DURATION as well as its construction. Building them from ripple sub-windows
    % would reintroduce the expert-type confound the twin arm exists to remove: the context
    % expert is binned spikes, so its lag-matched comparison has to be binned spikes too.

    for which = 1:2
        if which == 1, iv = fam; nm = 'ConA1'; else, iv = nov; nm = 'ConB1'; end

        % ---- the prescribed expert, and its duration-matched home-cage twin ----------------
        % BOTH experts are built from CONTINUOUS BINNED SPIKES, not from ripples. A context
        % exploration is precisely the epoch in which sharp-wave ripples are suppressed: these
        % recordings yield 0 to 78 ripples per exploration, median about 27, which is far too
        % few to estimate a rate vector and hopeless for a covariance basis. That is physiology,
        % not a threshold to lower. The classical reactivation template was never built from
        % ripples during behaviour either: it comes from binned firing during the behaviour, and
        % only the TARGET is ripple-based. The same split is used here.
        % The expert therefore describes theta-state firing and the target describes ripple-time
        % firing, which differ in overall rate by a large factor; the free multiplicative gain
        % absorbs exactly that, so the test is whether the PATTERN carries across the states.
        durC = iv(2) - iv(1);
        [eC, ZC] = local_behVec(ts, iv, BEH_BIN);
        if isempty(eC) || size(ZC, 1) < MIN_BIN_CTX || sum(ZC(:)) < MIN_SPK_CTX
            fprintf('%-34s %s: %d bins, %d spikes, no expert\n', local_name(fld), nm, ...
                size(ZC, 1), round(sum(ZC(:)))); continue
        end
        % twin: the equal-DURATION home-cage block ending where the context starts, stepped back
        % if it would overlap an earlier context. Duration matching is what matters now, since
        % the expert is estimated from bins rather than from ripples.
        tw = [iv(1) - durC, iv(1)];
        for c = 1:size(CTX, 1)
            if tw(1) < CTX(c,2) && tw(2) > CTX(c,1), tw = [CTX(c,1) - durC, CTX(c,1)]; end
        end
        if tw(1) < min(tOn)
            fprintf('%-34s %s: no room for a %.0f min twin before it\n', local_name(fld), nm, durC/60); continue
        end
        [eT, ZT] = local_behVec(ts, tw, BEH_BIN);
        if isempty(eT) || size(ZT, 1) < MIN_BIN_CTX
            fprintf('%-34s %s: twin has %d bins, skipped\n', local_name(fld), nm, size(ZT, 1)); continue
        end
        BC = local_binBasis(ZC, GEO_RANK);               % context expert, fast
        BT = local_binBasis(ZT, GEO_RANK);

        % ---- the OTHER context's expert: the template-specificity control -------------------
        % The twin arm cannot separate two explanations. The twin is a home-cage block and so is
        % the library, so the twin is redundant with it; the context is a different brain state
        % and is therefore less redundant whatever it contains. "Carries this experience" and
        % "is merely less redundant with home cage" predict the same ctx-over-twin result.
        % The other context separates them. It is an exploration expert over the same cells, in
        % the same state, equally non-redundant with the library, and it differs from the real
        % template only in WHICH experience it encodes. Scoring both on the same targets gives a
        % within-session 2x2: each template on its own targets and on the other's. A template
        % that beats the other one on its own targets, in both directions, is experience
        % specificity, and every property of expert type, state, redundancy, cell count and
        % library size has cancelled. For the ConA1 targets the ConB1 template also lies in the
        % FUTURE, so it additionally tests that a template is not expressed before its experience.
        ivO = nov; if which == 2, ivO = fam; end
        [eO, ZO] = local_behVec(ts, ivO, BEH_BIN);
        BO = local_binBasis(ZO, GEO_RANK);
        hasO = ~isempty(eO) && size(ZO, 1) >= MIN_BIN_CTX;
        % PRECISION MATCHING, which the cross-template test needs and the twin test does not.
        % The two explorations differ in length (20-39 min against 6-34), so their experts are
        % estimated from different numbers of bins, and a better-estimated rate vector predicts
        % better whatever it encodes. That asymmetry would appear as apparent specificity in
        % whichever direction favours the longer context, so both templates are rebuilt from the
        % same number of bins, spread evenly across each exploration so that neither is reduced
        % to one portion of its own behaviour.
        eCm = eC; eOm = eO; BCm = BC; BOm = BO; nBm = size(ZC, 1);
        if hasO
            nBm = min(size(ZC, 1), size(ZO, 1));
            iCm = unique(round(linspace(1, size(ZC, 1), nBm)));
            iOm = unique(round(linspace(1, size(ZO, 1), nBm)));
            ZCm = ZC(iCm,:); ZOm = ZO(iOm,:);
            eCm = sum(ZCm, 1)' / max(size(ZCm, 1) * BEH_BIN, eps);
            eOm = sum(ZOm, 1)' / max(size(ZOm, 1) * BEH_BIN, eps);
            BCm = local_binBasis(ZCm, GEO_RANK);
            BOm = local_binBasis(ZOm, GEO_RANK);
            hasO = ~isempty(BCm) && ~isempty(BOm);
        end

        % the library is unchanged: home-cage ripple sub-windows before the context. Its type
        % does not enter the contrast, because the two added experts are of identical type and
        % are added to an identical L0.
        kPre = find(tOn < min(iv(1), tw(1)) & ~inCtx);
        if numel(kPre) < LIB_MIN * MIN_RIP_SUB
            fprintf('%-34s %s: only %d home-cage ripples before it\n', local_name(fld), nm, numel(kPre)); continue
        end
        kLib = kPre;
        nSpkC = sum(ZC(:)); nSpkT = sum(ZT(:));

        % ---- a FAMILY OF TWINS AT MANY LAGS -------------------------------------------------
        % The twin arm matched construction and duration but sat at one fixed position, so it
        % could not separate "exploration" from "nearer in time". The decay curve matched lag but
        % used ripple-based home-cage windows, a different kind of object from the context
        % expert. These do both: home-cage blocks of the SAME duration as this context, binned
        % the same way, placed across the whole session so they span the lag range on both
        % sides. Regressing their performance on lag gives what a home-cage block of this
        % construction achieves at any distance, and the context is then scored against that.
        % This is the lag-matched version of the twin test, and it is the one that decides
        % whether exploration contributes anything beyond being recent.
        PS = struct('lam', {}, 'B', {}, 'tc', {});
        for p0 = min(tOn):PSEUDO_STEP:(max(tOn) - durC)
            ivP = [p0, p0 + durC];
            bad = false;
            for c = 1:size(CTX, 1)
                if ivP(1) < CTX(c,2) && ivP(2) > CTX(c,1), bad = true; end
            end
            if bad, continue; end
            [lamP, ZP] = local_behVec(ts, ivP, BEH_BIN);
            if isempty(lamP) || size(ZP, 1) < MIN_BIN_CTX, continue; end
            BP = local_binBasis(ZP, GEO_RANK);
            if isempty(BP), continue; end
            PS(end+1) = struct('lam', lamP, 'B', BP, 'tc', mean(ivP)); %#ok<AGROW>
        end
        tPS = [PS.tc];

        % ---- baseline library: home-cage sub-windows before the context --------------------
        [ML, ~, ~] = local_subWinsIdx(X, dur, tOn, kLib, SUB_WIN, MIN_RIP_SUB);
        if size(ML, 2) < LIB_MIN
            fprintf('%-34s %s: library of %d sub-windows\n', local_name(fld), nm, size(ML, 2)); continue
        end

        % ---- targets: home-cage sub-windows after the context ------------------------------
        kPost = find(tOn >= iv(2) + GUARD & tOn < iv(2) + TGT_MAX & ~inCtx);
        [MT, tT, ripT] = local_subWinsIdx(X, dur, tOn, kPost, SUB_WIN, MIN_RIP_SUB);
        if isempty(ripT), continue; end

        for q = 1:numel(ripT)
            k = ripT{q};
            Xt = X(k,:); dt = dur(k);
            if size(Xt, 1) < 2 * MIN_RIP_SUB, continue; end

            % slow: nested gains, in units of this window's own floor
            [d0, iBest] = min(arrayfun(@(c) local_devG(Xt, dt, ML(:,c), FREE_GAIN), 1:size(ML, 2)));
            dC = min(d0, local_devG(Xt, dt, eC, FREE_GAIN));
            dT = min(d0, local_devG(Xt, dt, eT, FREE_GAIN));
            mt = sum(Xt, 1)' / max(sum(dt), eps);
            ora = local_poissDev(Xt, dt, mt);
            ha = 1:2:size(Xt,1); hb = 2:2:size(Xt,1);
            mA = sum(Xt(ha,:), 1)' / max(sum(dt(ha)), eps);
            mB = sum(Xt(hb,:), 1)' / max(sum(dt(hb)), eps);
            spl = (local_devG(Xt(hb,:), dt(hb), mA, FREE_GAIN) * numel(hb) + ...
                   local_devG(Xt(ha,:), dt(ha), mB, FREE_GAIN) * numel(ha)) / size(Xt, 1);
            flr = max(spl - ora, eps);
            gC = (d0 - dC) / flr; gT = (d0 - dT) / flr;

            % The nested min() above is degenerate for the slow component and is kept only as a
            % diagnostic: a theta-state rate vector essentially never beats a ripple-time rate
            % vector from the same animal, so min() returns the library every time and the gain
            % is identically zero. Two tests with actual resolution replace it.
            %
            % (1) HEAD TO HEAD. Each expert alone, gain-free, against the other. No library, so
            %     nothing can mask the comparison, and the units are the window's own floor.
            dCa = local_devG(Xt, dt, eC, FREE_GAIN);
            dTa = local_devG(Xt, dt, eT, FREE_GAIN);
            slowHH = (dTa - dCa) / flr;              % positive: the context predicts better
            %
            % (2) NESTED POISSON GLM. The expert enters as a REGRESSOR beside the library rather
            %     than having to beat it outright, so a partial contribution is detectable:
            %         log lambda_j = b0 + b1 log(lib_j) + b2 log(expert_j),  offset log(total dur)
            %     For the mean-rate component the per-unit totals are sufficient, so this is a
            %     regression over nU rows, and the deviance drop from adding b2 is chi2 on 1 df.
            cj = sum(Xt, 1)'; offs = log(max(sum(dt), eps));
            lb = log(max(ML(:,iBest), 1e-6));
            dv1 = local_poissGLM(cj, [ones(numel(cj),1), lb], offs);
            dv2C = local_poissGLM(cj, [ones(numel(cj),1), lb, log(max(eC, 1e-6))], offs);
            dv2T = local_poissGLM(cj, [ones(numel(cj),1), lb, log(max(eT, 1e-6))], offs);
            glmC = dv1 - dv2C; glmT = dv1 - dv2T;
            %
            % (3) CELL-IDENTITY CONTROL, which is what decides whether (2) means anything.
            %     With 212 units and counts this size the firing is heavily overdispersed, so
            %     the chi2(1) reference for a deviance drop is void: ANY regressor that differs
            %     in shape from the library will absorb misspecification in how library rates map
            %     onto ripple-time rates, and the context expert differs in shape because it was
            %     recorded in a different brain state. The twin cannot do this, being the same
            %     kind of object as the library, which is enough on its own to produce the whole
            %     ctx-over-twin difference with no experience-specific content anywhere.
            %     So: permute which cell carries which rate. That preserves the expert's rate
            %     DISTRIBUTION exactly, and with it every shape property that could absorb
            %     misspecification, while destroying the only thing that can carry experience,
            %     namely which cell is which. What survives the permutation is the real effect.
            shC = nan(1, N_SHUF); shT = nan(1, N_SHUF);
            on1 = ones(numel(cj), 1);
            for sh = 1:N_SHUF
                pc = randperm(numel(eC)); pt = randperm(numel(eT));
                shC(sh) = dv1 - local_poissGLM(cj, [on1, lb, log(max(eC(pc), 1e-6))], offs);
                shT(sh) = dv1 - local_poissGLM(cj, [on1, lb, log(max(eT(pt), 1e-6))], offs);
            end
            glmCs = glmC - mean(shC, 'omitnan');      % identity-specific part, context
            glmTs = glmT - mean(shT, 'omitnan');      % identity-specific part, twin

            % fast: captured variance of the target's per-ripple residuals
            [vC, vT, vO, vS, vR] = local_fastScore(Xt, dt, BC, BT, GEO_RANK);

            % cross-template: own against other, on these same targets, at a common dimension
            xSlow = NaN; xFast = NaN;
            xSlowAdj = NaN; xFastAdj = NaN; nLagFit = 0; resSlow = NaN; resFast = NaN;
            if hasO
                dOm = local_devG(Xt, dt, eOm, FREE_GAIN);
                dCm = local_devG(Xt, dt, eCm, FREE_GAIN);
                xSlow = (dOm - dCm) / flr;                  % positive: own template fits better
                [vOwn, vOth] = local_fastScore(Xt, dt, BCm, BOm, GEO_RANK);
                xFast = vOwn - vOth;                        % positive: own template captures more

                % lag-adjusted version: both templates scored as residuals from the home-cage
                % decay curve, so a template that is merely nearer gets no credit for it
                lg = abs(tT(q) - tPS) / 3600;
                useL = lg > LAG_MIN;
                if sum(useL) >= 8
                    iu = find(useL);
                    dH = nan(1, numel(iu)); vH = nan(1, numel(iu));
                    for z = 1:numel(iu)
                        dH(z) = local_devG(Xt, dt, PS(iu(z)).lam, FREE_GAIN);
                        vH(z) = local_captOne(Xt, dt, PS(iu(z)).B, GEO_RANK);
                    end
                    ok2 = isfinite(dH) & isfinite(vH);
                    if sum(ok2) >= 8
                        % QUADRATIC in lag, not linear. Prediction error against lag saturates,
                        % and a straight line through a saturating curve leaves a bias that runs
                        % AGAINST finding specificity: simulating an exponential decay with a
                        % 4 h constant, a linear fit returns -0.042 under a true null and only
                        % +0.009 when the true effect is +0.05, while a quadratic returns +0.003
                        % and +0.050. The linear version would have hidden a real effect.
                        Lq = lg(iu(ok2))';
                        GL = [ones(numel(Lq), 1), Lq, Lq.^2];
                        bD = GL \ dH(ok2)'; bV = GL \ vH(ok2)';
                        pD = @(LL) bD(1) + bD(2) * LL + bD(3) * LL.^2;
                        pV = @(LL) bV(1) + bV(2) * LL + bV(3) * LL.^2;
                        lOwn = (tT(q) - iv(2)) / 3600; lOth = abs(tT(q) - ivO(2)) / 3600;
                        xSlowAdj = ((dOm - pD(lOth)) - (dCm - pD(lOwn))) / flr;
                        vOwn1 = local_captOne(Xt, dt, BCm, GEO_RANK);
                        vOth1 = local_captOne(Xt, dt, BOm, GEO_RANK);
                        xFastAdj = (vOwn1 - pV(lOwn)) - (vOth1 - pV(lOth));
                        nLagFit = sum(ok2);
                        % the single-arm residuals: is THIS exploration better than a home-cage
                        % block of the same construction, duration and temporal distance?
                        resSlow = (pD(lOwn) - dCm) / flr;     % positive: context predicts better
                        resFast = vOwn1 - pV(lOwn);           % positive: context captures more
                    end
                end
            end

            rows(end+1,:) = {local_name(fld), nm, (tT(q) - iv(2))/3600, numel(k), nU, ...
                gC, gT, gC - gT, d0, flr, ora, spl, ...
                slowHH, glmC, glmT, glmC - glmT, ...
                glmCs, glmTs, glmCs - glmTs, mean(shC, 'omitnan'), ...
                (dTa - ora) / flr, (dCa - ora) / flr, ...
                vC, vT, vC - vT, vO, vS, vR, ...
                xSlow, xFast, nBm, (tT(q) - iv(2))/3600, abs(tT(q) - ivO(2))/3600, ...
                xSlowAdj, xFastAdj, nLagFit, resSlow, resFast, numel(PS), ...
                nSpkC, nSpkT, durC/60, size(ZC, 1)}; %#ok<AGROW>
        end
        fprintf('%-34s %s: expert %.0f min / %d bins / %d spikes, twin %d spikes, library %d, targets %d\n', ...
            local_name(fld), nm, durC/60, size(ZC, 1), round(nSpkC), round(nSpkT), ...
            size(ML, 2), numel(ripT));
    end
end
if isempty(rows), error('inc:none', 'no session produced a usable context expert and targets'); end
T = cell2table(rows, 'VariableNames', {'Session','Context','HoursAfter','nRip','nUnits', ...
    'slowCtx','slowTwin','slowInc','devLib','floorDev','devOracle','devSplit', ...
    'slowHH','glmCtx','glmTwin','glmInc', ...
    'glmCtxId','glmTwinId','glmIncId','glmShufCtx','slowTwinAbs','slowCtxAbs', ...
    'fastCtx','fastTwin','fastInc','fastOracle','fastSplit','fastChance', ...
    'xSlow','xFast','nBinMatch','lagOwn','lagOther','xSlowAdj','xFastAdj','nLagFit','resSlow','resFast','nPseudo', ...
    'nSpkCtx','nSpkTwin','ctxMin','nBinCtx'});

%% ---------------------------------------------------------------------------------------------
%  REPORT
%% ---------------------------------------------------------------------------------------------
fprintf('\n============ PRESCRIBED-PERIOD INCLUSION TEST ============\n');
fprintf('%d sessions, %d target windows, %g min windows, guard %g min, free gain %d\n', ...
    numel(unique(T.Session)), height(T), SUB_WIN/60, GUARD/60, FREE_GAIN);
fprintf('Every number below is what ONE added expert buys, over a home-cage block matched to it\n');
fprintf('in duration and taken from just before the same context. Both experts are built from\n');
fprintf('continuous binned spikes; only the targets are ripples.\n');

CX = unique(T.Context);
fprintf('\n--- SLOW (mean rate) ---\n');
fprintf('%-8s %-12s %10s %10s %8s\n', 'context', 'test', 'mean', 'p', 'n ses');
for c = 1:numel(CX)
    k = strcmp(T.Context, CX{c});
    x = local_sesMean(T.slowInc(k), T.Session(k));
    fprintf('%-8s %-12s %10.3f %10.4f %8d   (ctx %.3f twin %.3f)\n', CX{c}, 'nested min', ...
        mean(x, 'omitnan'), local_signflip(x), numel(x), ...
        mean(T.slowCtx(k), 'omitnan'), mean(T.slowTwin(k), 'omitnan'));
    x = local_sesMean(T.slowHH(k), T.Session(k));
    fprintf('%-8s %-12s %10.3f %10.4f %8d\n', CX{c}, 'head to head', ...
        mean(x, 'omitnan'), local_signflip(x), numel(x));
    x = local_sesMean(T.glmInc(k), T.Session(k));
    fprintf('%-8s %-12s %10.3f %10.4f %8d   (ctx %.1f twin %.1f)\n', CX{c}, 'nested GLM', ...
        mean(x, 'omitnan'), local_signflip(x), numel(x), ...
        mean(T.glmCtx(k), 'omitnan'), mean(T.glmTwin(k), 'omitnan'));
    x = local_sesMean(T.glmIncId(k), T.Session(k));
    fprintf('%-8s %-12s %10.3f %10.4f %8d   (ctx %.1f twin %.1f, shuffle %.1f)\n', CX{c}, 'GLM id-spec', ...
        mean(x, 'omitnan'), local_signflip(x), numel(x), ...
        mean(T.glmCtxId(k), 'omitnan'), mean(T.glmTwinId(k), 'omitnan'), ...
        mean(T.glmShufCtx(k), 'omitnan'));
    fprintf('%-8s %-12s %10.3f %10s %8s   (twin %.3f)\n', CX{c}, 'abs vs oracle', ...
        mean(T.slowCtxAbs(k), 'omitnan'), '-', '-', mean(T.slowTwinAbs(k), 'omitnan'));
end
fprintf('"nested min" is kept only as a diagnostic: it is structurally zero whenever the added\n');
fprintf('expert never beats the whole library, which a behaviour-time rate vector never does.\n');
fprintf('"head to head" is each expert alone in floor units; "abs vs oracle" is how far each sits\n');
fprintf('above the window''s own mean, so the two arms can be read separately rather than only\n');
fprintf('as a difference. "nested GLM" is the raw deviance drop from adding the expert beside\n');
fprintf('the library: do NOT read it against chi2(1). With %d units and counts this size the\n', ...
    round(mean(T.nUnits, 'omitnan')));
fprintf('firing is overdispersed, so a large drop measures misspecification absorbed as much as\n');
fprintf('information carried, and the context expert can absorb more purely by being a different\n');
fprintf('kind of object than the library. "GLM id-spec" subtracts the mean drop obtained after\n');
fprintf('permuting which cell carries which rate: that permutation preserves the rate\n');
fprintf('distribution and every shape property, and destroys only cell identity. THAT is the\n');
fprintf('number to report. If id-spec collapses to zero while the raw drop is large, the GLM\n');
fprintf('result was shape, not experience.\n');

fprintf('\n--- FAST (ripple-to-ripple covariance), captured variance ---\n');
fprintf('%-8s %8s %8s %8s %8s %8s %8s %8s\n', 'context', 'CHANCE', 'ctx', 'twin', 'incl', 'p', 'oracle', 'split');
for c = 1:numel(CX)
    k = strcmp(T.Context, CX{c});
    x = local_sesMean(T.fastInc(k), T.Session(k));
    fprintf('%-8s %8.3f %8.3f %8.3f %8.3f %8.4f %8.3f %8.3f\n', CX{c}, ...
        mean(T.fastChance(k), 'omitnan'), ...
        mean(T.fastCtx(k), 'omitnan'), mean(T.fastTwin(k), 'omitnan'), ...
        mean(x, 'omitnan'), local_signflip(x), ...
        mean(T.fastOracle(k), 'omitnan'), mean(T.fastSplit(k), 'omitnan'));
end
fprintf('\nCHANCE is k/(n-1) with n = %d cells on average: the variance a RANDOM basis of the same\n', ...
    round(mean(T.nUnits, 'omitnan')));
fprintf('dimension captures, exactly, by symmetry. Read every other column against it first.\n');
fprintf('A context basis sitting at chance carries no information about the target whatever its\n');
fprintf('inclusion effect, because the contrast would then be between two chance-level numbers.\n');
fprintf('Oracle is the window''s own basis and split an independent half of it.\n');
for c = 1:numel(CX)
    k = strcmp(T.Context, CX{c});
    xc = local_sesMean(T.fastCtx(k) - T.fastChance(k), T.Session(k));
    fprintf('  %-8s ctx above chance = %+.4f, p = %.4f\n', CX{c}, mean(xc, 'omitnan'), local_signflip(xc));
end

fprintf('\n--- TEMPLATE SPECIFICITY: own context vs the other one, on the same targets ---\n');
fprintf('Both templates are explorations over the same cells, so expert type, brain state,\n');
fprintf('redundancy with the library, cell count and library size all cancel. Positive means the\n');
fprintf('experience that actually preceded these ripples fits them better than the other one.\n');
fprintf('%-8s %-10s %10s %10s %8s\n', 'targets', 'component', 'own - other', 'p', 'n ses');
for c = 1:numel(CX)
    k = strcmp(T.Context, CX{c});
    for f = {'xSlow', 'xFast'}
        x = local_sesMean(T.(f{1})(k), T.Session(k)); x = x(isfinite(x));
        if numel(x) < 2, continue; end
        fprintf('%-8s %-10s %10.4f %10.4f %8d\n', CX{c}, f{1}, mean(x), local_signflip(x), numel(x));
    end
end
fprintf('\n  DECOMPOSITION. Writing D(template, targets), the two rows above are\n');
fprintf('      row A = D(B,A) - D(A,A)        row B = D(A,B) - D(B,B)\n');
fprintf('  Each template appears once in the cross terms and once in the own terms, and so does\n');
fprintf('  each target set, so the MEAN of the rows cancels both main effects exactly and is the\n');
fprintf('  interaction: the specificity test. HALF THE DIFFERENCE is the template main effect,\n');
fprintf('  i.e. one template predicting ripple rates better everywhere. A negative row therefore\n');
fprintf('  does NOT refute specificity; it means the main effect is larger than the interaction.\n');
fprintf('%-8s %12s %10s %12s %10s\n', 'component', 'interaction', 'p', 'main effect', 'p');
for f = {'xSlow', 'xFast'}
    s = unique(T.Session); inter = nan(numel(s), 1); mainE = nan(numel(s), 1);
    for q = 1:numel(s)
        k1 = strcmp(T.Session, s{q}) & strcmp(T.Context, 'ConA1');
        k2 = strcmp(T.Session, s{q}) & strcmp(T.Context, 'ConB1');
        if any(k1) && any(k2)
            a = mean(T.(f{1})(k1), 'omitnan'); b = mean(T.(f{1})(k2), 'omitnan');
            inter(q) = (a + b) / 2; mainE(q) = (a - b) / 2;
        end
    end
    fprintf('%-8s %12.4f %10.4f %12.4f %10.4f\n', f{1}, ...
        mean(inter, 'omitnan'), local_signflip(inter), ...
        mean(mainE, 'omitnan'), local_signflip(mainE));
end
fprintf('  Both templates are built from %d bins, so neither wins on estimation precision.\n', ...
    round(mean(T.nBinMatch, 'omitnan')));

fprintf('\n  THE CONFOUND THE INTERACTION DOES NOT CANCEL: temporal proximity. Each template is\n');
fprintf('  nearer in time to its own targets than to the other set, and with ongoing drift a\n');
fprintf('  nearer template predicts better whatever it encodes. That produces a positive\n');
fprintf('  interaction with no experience specificity anywhere. Observed lags:\n');
for c = 1:numel(CX)
    k = strcmp(T.Context, CX{c});
    fprintf('    %-6s targets: own template %.2f h away, other %.2f h away\n', CX{c}, ...
        mean(T.lagOwn(k), 'omitnan'), mean(T.lagOther(k), 'omitnan'));
end
fprintf('\n  DRIFT-ADVERSARIAL SUBSET. The lags above are not what the layout suggests: the ConA1\n');
fprintf('  targets run up to %g h past ConA1, so many fall AFTER ConB1 has happened, and there\n', TGT_MAX/3600);
fprintf('  the other template is the MORE RECENT experience. Those windows are a drift-controlled\n');
fprintf('  test on their own: drift predicts the closer template wins, so if the own template\n');
fprintf('  still wins where it is the FURTHER one, proximity cannot be the explanation. No decay\n');
fprintf('  curve is needed for this; the design already contains both signs of the confound.\n');
kAdv = T.lagOther < T.lagOwn;        % the other template is closer: drift favours the OTHER
kSup = ~kAdv;
fprintf('%-28s %10s %10s %10s %8s %8s\n', 'subset', 'xSlow', 'p', 'xFast', 'p', 'n win');
for s = 1:2
    if s == 1, kk = kAdv; lab = 'other CLOSER (adversarial)'; else, kk = kSup; lab = 'own closer (supportive)'; end
    if ~any(kk), continue; end
    xs = local_sesMean(T.xSlow(kk), T.Session(kk)); xs = xs(isfinite(xs));
    xf = local_sesMean(T.xFast(kk), T.Session(kk)); xf = xf(isfinite(xf));
    ps = NaN; pf = NaN;
    if numel(xs) > 1, ps = local_signflip(xs); end
    if numel(xf) > 1, pf = local_signflip(xf); end
    fprintf('%-28s %10.4f %10.4f %10.4f %8.4f %8d\n', lab, ...
        mean(xs, 'omitnan'), ps, mean(xf, 'omitnan'), pf, sum(kk));
end
fprintf('  mean lag advantage of the OTHER template in the adversarial subset: %.2f h\n', ...
    mean(T.lagOwn(kAdv) - T.lagOther(kAdv), 'omitnan'));
fprintf('  COMPOSITION, which decides whether this subset can be read at all:\n');
for c = 1:numel(CX)
    fprintf('    %-6s own: %d of %d adversarial windows\n', CX{c}, ...
        sum(kAdv & strcmp(T.Context, CX{c})), sum(kAdv));
end
fprintf('  If one context supplies nearly all of them, the subset cannot cancel the template\n');
fprintf('  main effect and will simply report it; compare its value against the main effect\n');
fprintf('  above before believing it. ConB1 is the last context, so no target after it has a\n');
fprintf('  nearer competitor and this imbalance is structural, not a sampling accident.\n');

fprintf('\n--- LAG-ADJUSTED SPECIFICITY (the decay curve, not a subset) ---\n');
fprintf('Every home-cage sub-window is an experience-neutral pseudo-template. Regressing its\n');
fprintf('prediction error on |lag| gives the drift-only expectation at any lag; both context\n');
fprintf('templates are then scored as residuals from that curve, so neither is credited for\n');
fprintf('merely being nearer. This needs no balanced subset and no assumed decay constant.\n');
fprintf('%-10s %12s %10s %10s\n', 'component', 'lag-adjusted', 'p', 'n ses');
for f = {'xSlowAdj', 'xFastAdj'}
    s = unique(T.Session); d = nan(numel(s), 1);
    for q = 1:numel(s)
        k1 = strcmp(T.Session, s{q}) & strcmp(T.Context, 'ConA1');
        k2 = strcmp(T.Session, s{q}) & strcmp(T.Context, 'ConB1');
        if any(k1) && any(k2)
            d(q) = (mean(T.(f{1})(k1), 'omitnan') + mean(T.(f{1})(k2), 'omitnan')) / 2;
        end
    end
    fprintf('%-10s %12.4f %10.4f %10d\n', f{1}, mean(d, 'omitnan'), local_signflip(d), ...
        sum(isfinite(d)));
end
fprintf('\n--- LAG-MATCHED TWIN TEST: does exploration beat home cage AT THE SAME DISTANCE? ---\n');
fprintf('This is the claim that survives if specificity does not. The pseudo-templates are\n');
fprintf('home-cage blocks of the same duration, binned the same way, differing from the context\n');
fprintf('only in that the animal was not exploring and in how far away they sit - and the lag is\n');
fprintf('then regressed out. Positive means the exploration contributes something a home-cage\n');
fprintf('block at the same remove does not. Chance is not the comparison here; home cage is.\n');
fprintf('%-8s %-10s %12s %10s %8s\n', 'context', 'component', 'ctx - HC(lag)', 'p', 'n ses');
for c = 1:numel(CX)
    k = strcmp(T.Context, CX{c});
    for f = {'resSlow', 'resFast'}
        x = local_sesMean(T.(f{1})(k), T.Session(k)); x = x(isfinite(x));
        if numel(x) < 2, continue; end
        fprintf('%-8s %-10s %12.4f %10.4f %8d\n', CX{c}, f{1}, mean(x), local_signflip(x), numel(x));
    end
end
fprintf('  %d pseudo-templates per context on average, spaced %g min apart.\n', ...
    round(mean(T.nPseudo, 'omitnan')), PSEUDO_STEP/60);
fprintf('  Note this remains conservative: the home-cage blocks CONTAIN ripples and the context\n');
fprintf('  does not, and ripple co-firing predicts later ripple co-firing well.\n');

fprintf('\n  decay curve fitted on %d home-cage windows per target on average, no closer than\n', ...
    round(mean(T.nLagFit, 'omitnan')));
fprintf('  %.2f h. This is the interaction again, so template and target main effects still\n', LAG_MIN);
fprintf('  cancel; what is new is that the lag advantage has been removed from both arms.\n');

fprintf('\n--- NOVEL vs FAMILIAR, paired within session ---\n');
for fld = {'slowHH', 'glmInc', 'glmIncId', 'glmCtxId', 'fastInc', 'fastCtx', 'xSlow', 'xFast'}
    [dd, ns] = local_pairCtx(T, fld{1}, 'ConB1', 'ConA1');
    fprintf('  %-8s ConB1 - ConA1 = %+.3f, p = %.4f (n = %d sessions with both)\n', ...
        fld{1}, mean(dd, 'omitnan'), local_signflip(dd), ns);
end
fprintf('  ConA1 is familiar by design and ConB1 novel, so an inclusion effect carried by the\n');
fprintf('  novel context alone is the specific result; one present in both is a statement about\n');
fprintf('  exploration rather than about novelty.\n');

fprintf('\n--- TIME COURSE: inclusion effect by hours after the context ended ---\n');
eb = 0:BIN/3600:TGT_MAX/3600;
fprintf('%-8s %10s %10s %10s %8s\n', 'context', 'bin (h)', 'slow', 'fast', 'n ses');
for c = 1:numel(CX)
    for b = 1:numel(eb) - 1
        k = strcmp(T.Context, CX{c}) & T.HoursAfter >= eb(b) & T.HoursAfter < eb(b+1);
        if ~any(k), continue; end
        xs = local_sesMean(T.slowInc(k), T.Session(k));
        xf = local_sesMean(T.fastInc(k), T.Session(k));
        fprintf('%-8s %4.1f-%4.1f %10.3f %10.3f %8d\n', CX{c}, eb(b), eb(b+1), ...
            mean(xs, 'omitnan'), mean(xf, 'omitnan'), numel(xs));
    end
end

fprintf('\n--- IS PREDICTION POSSIBLE AT ALL WITHOUT THE CONTEXT? ---\n');
fprintf('The contrast above is a difference of gains, so it supports "including the context\n');
fprintf('improved prediction". The stronger phrasing, that prediction was possible ONLY with it,\n');
fprintf('needs the baseline library to sit at or below the split-half floor on its own:\n');
fprintf('%-8s %12s %12s %12s %10s\n', 'context', 'D(library)', 'D(split)', 'D(oracle)', 'frac>floor');
for c = 1:numel(CX)
    k = strcmp(T.Context, CX{c});
    fprintf('%-8s %12.4f %12.4f %12.4f %10.2f\n', CX{c}, ...
        mean(T.devLib(k), 'omitnan'), mean(T.devSplit(k), 'omitnan'), ...
        mean(T.devOracle(k), 'omitnan'), mean(T.devLib(k) < T.devSplit(k)));
end
fprintf('  frac>floor is the share of target windows the library already predicts better than\n');
fprintf('  half the window predicts the other half. If that is near one, the library alone works\n');
fprintf('  and the honest claim is the weaker, incremental one.\n');

fprintf('\n--- match quality (the contrast is only as good as this) ---\n');
fprintf('  context %.0f min / %d bins / %d spikes; twin same duration, %d spikes (ratio %.2f)\n', ...
    mean(T.ctxMin, 'omitnan'), round(mean(T.nBinCtx, 'omitnan')), ...
    round(mean(T.nSpkCtx, 'omitnan')), round(mean(T.nSpkTwin, 'omitnan')), ...
    mean(T.nSpkTwin ./ max(T.nSpkCtx, eps), 'omitnan'));
fprintf('  duration and bin count are matched by construction; the spike ratio is not, because\n');
fprintf('  running and resting differ in firing rate. The free gain removes the level, but a\n');
fprintf('  large ratio still means the two experts are estimated with different precision, which\n');
fprintf('  favours whichever has more spikes. Read the ConB1 - ConA1 contrast as the primary\n');
fprintf('  endpoint for that reason: both are explorations, so state and precision are matched.\n');
fprintf('  Note also that the twin CONTAINS ripples while the context does not, so the twin is a\n');
fprintf('  strong baseline and this test is conservative against the context.\n');

%% ---------------------------------------------------------------------------------------------
%  FIGURES
%% ---------------------------------------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [60 60 1200 760]);
tl = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
cols = [0.20 0.45 0.75; 0.80 0.35 0.15];

for p = 1:2
    if p == 1, fc = 'slow'; else, fc = 'fast'; end
    nexttile(tl); hold on
    for c = 1:numel(CX)
        k = strcmp(T.Context, CX{c});
        xc = local_sesMean(T.([fc 'Ctx'])(k), T.Session(k));
        xt = local_sesMean(T.([fc 'Twin'])(k), T.Session(k));
        for s = 1:numel(xc)
            plot([1 2] + (c-1)*2.5, [xt(s) xc(s)], '-', 'Color', [0.7 0.7 0.7]);
        end
        plot((1:2) + (c-1)*2.5, [mean(xt,'omitnan') mean(xc,'omitnan')], '-o', ...
            'LineWidth', 2, 'Color', cols(min(c,2),:), 'MarkerFaceColor', cols(min(c,2),:));
    end
    set(gca, 'XTick', [1 2 3.5 4.5], 'XTickLabel', {'twin','ctx','twin','ctx'});
    for c = 1:numel(CX)
        text(1.5 + (c-1)*2.5, max(ylim), CX{c}, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom');
    end
    ylabel(sprintf('%s: gain from one added expert', fc));
    title(sprintf('%s component, paired within session', fc)); box off
end

for p = 1:2
    if p == 1, fc = 'slowInc'; else, fc = 'fastInc'; end
    nexttile(tl); hold on
    for c = 1:numel(CX)
        mu = nan(1, numel(eb)-1); se = nan(1, numel(eb)-1); ct = (eb(1:end-1) + eb(2:end))/2;
        for b = 1:numel(eb) - 1
            k = strcmp(T.Context, CX{c}) & T.HoursAfter >= eb(b) & T.HoursAfter < eb(b+1);
            if ~any(k), continue; end
            x = local_sesMean(T.(fc)(k), T.Session(k));
            mu(b) = mean(x, 'omitnan'); se(b) = std(x, 'omitnan') / max(sqrt(sum(isfinite(x))), 1);
        end
        errorbar(ct, mu, se, '-o', 'LineWidth', 1.5, 'Color', cols(min(c,2),:), ...
            'MarkerFaceColor', cols(min(c,2),:), 'DisplayName', CX{c});
    end
    yline(0, 'k:'); xlabel('h after the context ended');
    ylabel('inclusion effect (ctx - twin)');
    title(sprintf('%s, time course', fc)); legend('Location','best','Box','off'); box off
end
title(tl, 'Does including the context exploration improve prediction of later ripples?');

if SAVE_OUT
    save(fullfile(pwd, 'driftInclusionTest_out.mat'), 'T', 'SUB_WIN', 'GUARD', 'TGT_MAX', ...
        'BIN', 'GEO_RANK', 'FREE_GAIN', 'LIB_MIN');
    fprintf('\nsaved driftInclusionTest_out.mat\n');
end


%% ---------------------------------------------------------------------------------------------
%  helpers
%% ---------------------------------------------------------------------------------------------
function [M, tC, ripC] = local_subWinsIdx(X, dur, tOn, keep, subWin, minRip)
% sub-windows built only from the ripples in `keep`, binned on absolute time so that a window
% straddling an excluded stretch simply holds fewer ripples and is dropped if it holds too few
M = []; tC = []; ripC = {};
if isempty(keep), return; end
tk = tOn(keep);
ed = min(tk):subWin:max(tk) + subWin;
M = nan(size(X, 2), numel(ed) - 1); tC = nan(1, numel(ed) - 1); ripC = cell(1, numel(ed) - 1);
for b = 1:numel(ed) - 1
    k = keep(tk >= ed(b) & tk < ed(b+1));
    if numel(k) < minRip, continue; end
    M(:,b) = sum(X(k,:), 1)' / sum(dur(k));
    tC(b) = (ed(b) + ed(b+1)) / 2;
    ripC{b} = k;
end
ok = all(isfinite(M), 1) & ~cellfun(@isempty, ripC);
M = M(:, ok); tC = tC(ok); ripC = ripC(ok);
end


function [vC, vT, vO, vS, vR] = local_fastScore(Xt, dt, BC, BT, k)
% fraction of the target window's ripple-to-ripple variance captured by the context basis, the
% twin basis, the window's own basis and an independent half of it. The window's own radial
% direction is removed from the residuals and from every candidate basis, so none of the four
% can be inflated by a difference in overall rate.
vC = NaN; vT = NaN; vO = NaN; vS = NaN; vR = NaN;
mt = sum(Xt, 1)' / max(sum(dt), eps);
rho = mt / max(norm(mt), eps);
R = Xt ./ max(dt(:), eps) - mt';
R = R - (R * rho) * rho';
vTot = sum(R(:).^2);
if vTot <= eps || isempty(BC) || isempty(BT), return; end
ha = 1:2:size(R,1); hb = 2:2:size(R,1);
if numel(ha) < 3 || numel(hb) < 3, return; end
cand = {local_orthTo(BC, rho), local_orthTo(BT, rho), ...
        local_orthTo(local_topV(R, k), rho), local_orthTo(local_topV(R(ha,:), k), rho)};
kE = min(cellfun(@(b) size(b, 2), cand));
if kE < 1, return; end
vC = sum((R * cand{1}(:,1:kE)).^2, 'all') / vTot;
vT = sum((R * cand{2}(:,1:kE)).^2, 'all') / vTot;
vO = sum((R * cand{3}(:,1:kE)).^2, 'all') / vTot;
vS = sum((R(hb,:) * cand{4}(:,1:kE)).^2, 'all') / max(sum(R(hb,:).^2, 'all'), eps);
% CHANCE. A random kE-dimensional basis captures exactly kE/(n-1) of the residual variance in
% expectation, by symmetry, whatever the eigenvalue spectrum: the residuals live in n-1
% dimensions once the radial direction is removed. Without this on the page a captured variance
% of 0.04 reads as "4% explained" when with a hundred cells it is precisely chance.
vR = kE / max(size(R, 2) - 1, 1);
end


function [d, ns] = local_pairCtx(T, fld, c1, c2)
% within-session difference between two contexts, so animal and recording cancel
s = unique(T.Session); d = nan(numel(s), 1);
for q = 1:numel(s)
    k1 = strcmp(T.Session, s{q}) & strcmp(T.Context, c1);
    k2 = strcmp(T.Session, s{q}) & strcmp(T.Context, c2);
    if any(k1) && any(k2)
        d(q) = mean(T.(fld)(k1), 'omitnan') - mean(T.(fld)(k2), 'omitnan');
    end
end
ns = sum(isfinite(d));
end


function dev = local_poissGLM(y, G, offs)
% Poisson regression with a log link and a scalar offset, by IRLS. Returns the deviance.
% Used on per-unit total counts, which are sufficient for the mean-rate component.
y = y(:); b = zeros(size(G, 2), 1); b(1) = log(max(mean(y), eps)) - offs;
for it = 1:50
    eta = G * b + offs;
    eta = min(max(eta, -30), 30);
    mu = exp(eta);
    W = max(mu, 1e-10);
    z = eta - offs + (y - mu) ./ W;
    A = G' * (W .* G); rhs = G' * (W .* z);
    A = A + 1e-8 * eye(size(A)) * max(abs(diag(A)) + eps);   % ridge against collinear predictors
    bn = A \ rhs;
    if any(~isfinite(bn)), break; end
    if max(abs(bn - b)) < 1e-8, b = bn; break; end
    b = bn;
end
mu = exp(min(max(G * b + offs, -30), 30));
t = y .* log(max(y, 1e-12) ./ mu); t(y == 0) = 0;
dev = 2 * sum(t - (y - mu));
end


function v = local_captOne(Xt, dt, B, k)
% captured variance of one basis on the target window's per-ripple residuals, with the window's
% own radial direction removed from both. Used for the context templates AND for every home-cage
% pseudo-template, so the decay curve and the templates are measured on one scale.
v = NaN;
if isempty(B) || size(Xt, 1) < 4, return; end
mt = sum(Xt, 1)' / max(sum(dt), eps);
rho = mt / max(norm(mt), eps);
R = Xt ./ max(dt(:), eps) - mt';
R = R - (R * rho) * rho';
vTot = sum(R(:).^2);
if vTot <= eps, return; end
Bo = local_orthTo(B, rho);
if isempty(Bo), return; end
kE = min(k, size(Bo, 2));
v = sum((R * Bo(:,1:kE)).^2, 'all') / vTot;
end


function [lam, Z] = local_behVec(ts, iv, binw)
% CONTINUOUS binned spikes over an interval: the rate vector and the bin x unit count matrix.
% This is the behaviour-time template, built the way a reactivation template is built, and it
% does not depend on ripples occurring during the interval.
lam = []; Z = zeros(0, numel(ts));
if iv(2) <= iv(1), return; end
ed = iv(1):binw:iv(2);
if numel(ed) < 3, return; end
Z = zeros(numel(ed) - 1, numel(ts));
for u = 1:numel(ts)
    Z(:,u) = histcounts(ts{u}, ed)';
end
lam = sum(Z, 1)' / (iv(2) - iv(1));          % spikes per second per unit
end


function B = local_binBasis(Z, k)
% covariance directions of binned firing, with the radial direction projected out so the basis
% describes which cells co-vary rather than how strongly the population fired
B = [];
if size(Z, 1) < 5, return; end
m = mean(Z, 1)';
rho = m / max(norm(m), eps);
R = Z - m';
R = R - (R * rho) * rho';
B = local_topV(R, k);
end


function B = local_blockBasis(X, dur, k)
B = [];
if size(X, 1) < 3, return; end
m = sum(X, 1)' / max(sum(dur), eps);
rho = m / max(norm(m), eps);
R = X ./ max(dur(:), eps) - m';
R = R - (R * rho) * rho';
B = local_topV(R, k);
end


function Bo = local_orthTo(B, rho)
Bo = [];
if isempty(B), return; end
B = B - rho * (rho' * B);
[Q, Rq] = qr(B, 0);
keep = abs(diag(Rq)) > 1e-8 * max(abs(diag(Rq)));
Bo = Q(:, keep);
end


function d = local_devG(X, dur, lam, useGain)
% Poisson deviance with one free multiplicative gain fitted on the target in closed form, so a
% uniform rate difference between the training period and the target cannot be read as pattern
lam = max(lam(:), 0);
if useGain
    den = sum(dur(:)) * sum(lam);
    if den > eps, lam = max(sum(X(:)) / den, eps) * lam; end
end
d = local_poissDev(X, dur, lam);
end


function d = local_poissDev(X, dur, lam)
Lam = max(dur(:) * lam(:)', 1e-9);
T = X .* log(max(X, 1e-12) ./ Lam);
T(X == 0) = 0;
d = 2 * sum(T - (X - Lam), 'all') / size(X, 1);
end


function V = local_topV(R, k)
k = max(min([k, size(R, 2) - 1, size(R, 1) - 1]), 1);
[~, ~, V] = svd(R, 'econ');
V = V(:, 1:min(k, size(V, 2)));
end


function x = local_sesMean(y, grp)
g = unique(grp);
x = nan(numel(g), 1);
for q = 1:numel(g), x(q) = mean(y(strcmp(grp, g{q})), 'omitnan'); end
end


function p = local_signflip(x)
x = x(isfinite(x)); n = numel(x);
if n < 2, p = NaN; return; end
if n > 20
    s = sign(randn(n, 2000)); s(:,1) = 1;
else
    s = double(dec2bin(0:2^n - 1) - '0')' * 2 - 1;
end
m = abs(x(:)' * s / n);
p = max(mean(m >= abs(mean(x)) - 1e-12), 2 / 2^n);
end


function [X, dur, tOn, nU, ts] = local_ripCounts(fld, RP)
% per-ripple spike counts AND the pyramidal spike trains themselves, since the behaviour-time
% experts are built by binning the trains rather than by counting ripples
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
if isempty(idx), return; end
on = RP.ts1(:); off = RP.ts1Off(:);
[on, o] = sort(on); off = off(o);
keep = off > on & [true; on(2:end) > off(1:end-1)];
on = on(keep); off = off(keep);
if numel(on) < 2, return; end
edges = reshape([on'; off'], [], 1);
X = zeros(numel(on), numel(idx));
for u = 1:numel(idx)
    c = histcounts(ts{idx(u)}, edges);
    X(:,u) = c(1:2:end)';
end
dur = off - on; tOn = on; nU = numel(idx);
ts = ts(idx);                 % pyramidal trains only, in the same column order as X
end


function s = local_name(fld)
[~, s] = fileparts(fld);
end