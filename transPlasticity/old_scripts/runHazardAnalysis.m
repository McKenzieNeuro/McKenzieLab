%% runHazardAnalysis.m
% Synchrony-state dynamics as a two-state survival process: does ASV change
% how easily the population ENTERS a synchronous bout, how long it STAYS in
% one, or both? Fit for persistence (drug-naive vs >=24 h post-exposure,
% drug-free recordings) and for the acute effect (pre-infusion vs 4-6 h
% post-infusion, ASV against vehicle).
%
% WHY. Rate, duration, and duty cycle per epoch are ratios of a handful of
% bouts (often 1-5), so a single change in bout dynamics can surface in any
% of them depending on which animals are sampled -- which is what happened
% when two animals were added (duty significant -> duration significant).
% Here every bout and every inter-bout interval is an observation, right-
% censoring at epoch edges is handled, and epochs with no bouts at all
% contribute to the entry hazard instead of being dropped.
%
%   exit hazard   lower after exposure  -> bouts last longer
%   entry hazard  higher after exposure -> bouts start more often
%   duty cycle ~ entry / (entry + exit) under the two-state model, so a pure
%   exit-hazard decrease raises duration and duty while slightly LOWERING
%   bout rate per hour -- the pattern the ratio metrics kept disagreeing on.
%
% BOUT DEFINITION is independent of epoch length (see boutHazard.m), so the
% WHOLE drug-free period is used: all of the pre-infusion baseline on
% infusion days, the whole recording on non-infusion days. Several times
% more data per session than the 1 h LFIX epochs.
%
% COST. This is a full syncBatch pass over much longer windows than before,
% in fresh cache directories (the old caches have no stored surrogate band).
% Run PERSIST and ACUTE separately if needed; both checkpoint.
%
% STILL TRUE: post-exposure recordings are always later than drug-naive ones.
% Days-in-study stays in the persistence model. The hazard framing adds power
% and removes censoring bias; it does not remove that confound.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
WIDTH = 60; STEP = 5; MATCHN = 20; MINUNITS = 10; NSURRW = 50;
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
EXPOSURE_PAT = '^ASV$';
MIN_EPOCH = 1800;                 % s of usable window required
LAT = 4*3600; POSTLEN = 2*3600;   % acute treatment window: 4-6 h post-infusion
MINRUNS = {'null95', 'null99', 24};  % bout definitions; first is primary
RUN_PERSIST = true;
RUN_ACUTE   = true;
RUN_RATES   = true;   % firing rate in/out of bouts (primary bout definition); reloads spike files
SAVE_TAG    = '';

% Timescale sweep: runTimescaleSweep.m defines the struct SW before calling
% this script once per (width, step). Nothing here changes when SW is absent.
if exist('SW','var') && isstruct(SW)
    WIDTH = SW.width; STEP = SW.step; MATCHN = SW.matchN;
    RUN_PERSIST = SW.persist; RUN_ACUTE = true; RUN_RATES = false;
    SAVE_TAG = SW.tag;
    % fixed bout definitions in ELAPSED seconds, converted to a window count
    % at this scale: k windows span (k-1)*STEP + WIDTH seconds
    MINRUNS = {'null95','null99'};
    for q = SW.minRunSec(:)'
        MINRUNS{end+1} = max(1, ceil((q - WIDTH)/STEP) + 1); %#ok<SAGROW>
    end
    fprintf('\n[sweep] width %g s, step %g s, MatchN %d; fixed bouts = %s windows (= %s s elapsed)\n', ...
        WIDTH, STEP, MATCHN, mat2str(cell2mat(MINRUNS(3:end))), ...
        mat2str((cell2mat(MINRUNS(3:end))-1)*STEP + WIDTH));
end

fils = getAllExtFiles('R:\TransPlasticity','mat',1);
fils = fils(contains(fils,'spikes.cellinfo.mat'));
cnd = cellfun(@(f) conditionFromPath(f,'TransPlasticity'), fils, 'UniformOutput',false);
ani = cellfun(@(f) animalFromPath(f,'TransPlasticity'),    fils, 'UniformOutput',false);
dnm = cellfun(@(f) sessionDateFromPath(f),                 fils);
keep = ismember(cnd, CONDS) & ~cellfun(@isempty, ani) & isfinite(dnm);
fils = fils(keep); cnd = cnd(keep); ani = ani(keep); dnm = dnm(keep);

inf_t = nan(numel(fils),1);
for i = 1:numel(fils), inf_t(i) = infusionTime(fils{i}); end
hasInf = ~isnan(inf_t);

% infusion-day sessions with no marker: drug state of any window unknown
isInfCond = ismember(cnd, {'ASV','SAL'});
bad = isInfCond & ~hasInf;
if any(bad)
    fprintf('dropping %d %s-day session(s) with no Infusion_start marker (drug state unknown):\n', ...
        sum(bad), 'infusion');
    fprintf('   %s\n', fils{bad});
end
fils = fils(~bad); cnd = cnd(~bad); ani = ani(~bad); dnm = dnm(~bad);
inf_t = inf_t(~bad); hasInf = hasInf(~bad);

% durations (cheap screen; same options as runPersistenceAnalysis -> cache hit)
% duration + usability per session (same rule as the old syncBatch screen, without
% running spikeSync over every unit -- that ran out of memory on large sessions)
[durs, scOK] = sessionScreen(fils, [0.1 10], MINUNITS, fullfile(CACHE_ROOT, 'screenLite'));

cfgB = { 'RateRange',[0.1 10],'Width',WIDTH,'Step',STEP,'MinSpikes',10, ...
         'MinUnits',MINUNITS,'MatchN',MATCHN,'NullMethod','none', ...
         'DetectTransients',true,'NSurrWindow',NSURRW,'RunAlpha',0.05, ...
         'Parallel', ~isempty(ver('parallel')), 'CheckpointEvery',5 };

%% =========================================================================
%% PERSISTENCE
%% =========================================================================
if RUN_PERSIST
    fprintf('\n================ PERSISTENCE (drug-free, whole period) ================\n');
    w = [zeros(numel(fils),1), durs];
    w(hasInf,2) = min(inf_t(hasInf), durs(hasInf));
    use = scOK & (w(:,2) - w(:,1)) >= MIN_EPOCH;
    fprintf('%d/%d sessions supply >= %.0f min drug-free (median %.1f h)\n', ...
        sum(use), numel(use), MIN_EPOCH/60, median(w(use,2)-w(use,1))/3600);

    uf = fils(use);
    mapP = containers.Map(uf, num2cell(w(use,:), 2));
    TP = syncBatch(uf, cfgB{:}, 'Window', @(f) mapP(f), ...
        'CacheDir', fullfile(CACHE_ROOT,'hzPersist'), ...
        'CheckpointFile', fullfile(CACHE_ROOT,['ckpt_hzPersist' SAVE_TAG '.mat']));

    ok = [TP.ok]';
    E = TP(ok);
    ea = ani(use); ea = ea(ok);  ec = cnd(use); ec = ec(ok);
    ed = dnm(use); ed = ed(ok);  eInf = hasInf(use); eInf = eInf(ok);

    % phase, same rule as runPersistenceAnalysis (strict for infusion days)
    phase = repmat({'pre'}, numel(E), 1); dLast = nan(numel(E),1);
    days  = nan(numel(E),1);
    uan = unique(ea);
    for i = 1:numel(uan)
        m = find(strcmp(ea, uan{i}));
        days(m) = ed(m) - min(ed(m));
        expD = unique(ed(m(~cellfun(@isempty, regexp(ec(m), EXPOSURE_PAT,'once')))));
        for j = m'
            if eInf(j), prior = expD(expD < ed(j)); else, prior = expD(expD <= ed(j)); end
            if ~isempty(prior), phase{j} = 'post'; dLast(j) = ed(j) - max(prior); end
        end
    end
    tooSoon = strcmp(phase,'post') & dLast < 1;
    if any(tooSoon)
        fprintf('dropping %d post epoch(s) recorded on an exposure date\n', sum(tooSoon));
    end
    sel = ~tooSoon;
    E = E(sel); ea = ea(sel); phase = phase(sel); days = days(sel); dLast = dLast(sel);
    ec = ec(sel); ed = ed(sel); eInf = eInf(sel);
    fprintf('%d epochs: %d drug-naive, %d post (>=1 d); %d animals, %d with both\n', ...
        numel(E), sum(strcmp(phase,'pre')), sum(strcmp(phase,'post')), numel(unique(ea)), ...
        sum(arrayfun(@(k) any(strcmp(ea,uan{k}) & strcmp(phase,'pre')) && ...
                          any(strcmp(ea,uan{k}) & strcmp(phase,'post')), 1:numel(uan))));

    RES.persist = struct();
    for mr = 1:numel(MINRUNS)
        minRun = MINRUNS{mr};
        [Xt, Nt, epochS, St] = local_stack(E, STEP, minRun, ea, phase, days, [], []);
        tag = local_mrTag(minRun);
        fprintf('\n---- bout definition: minRun = %s (windows used: median %g) ----\n', ...
            tag, median(epochS.minRunUsed));
        fprintf('%d bouts (%d ended in view), %d out-intervals (%d ended by a bout)\n', ...
            height(Xt), sum(Xt.Event), height(Nt), sum(Nt.Event));

        fx = local_fitHazard(Xt, 'Phase + Days', 'Phase', 'Phase_post', 'Phase');
        fn = local_fitHazard(Nt, 'Phase + Days', 'Phase', 'Phase_post', 'Phase');
        local_report('EXIT  (leave bout)', fx, Xt, 'Phase');
        local_report('ENTRY (start bout)', fn, Nt, 'Phase');
        [ax, nx] = local_animalLogRatio(Xt, 'Phase', 'pre', 'post');
        [an, nn] = local_animalLogRatio(Nt, 'Phase', 'pre', 'post');
        fprintf('  animal-level log HR, exit : mean %+.3f, %d/%d negative, sign-flip p = %.4f\n', ...
            mean(ax), sum(ax<0), numel(ax), local_signflip(ax));
        fprintf('  animal-level log HR, entry: mean %+.3f, %d/%d positive, sign-flip p = %.4f\n', ...
            mean(an), sum(an>0), numel(an), local_signflip(an));
        if nx < 6, fprintf('  (exact floor at n=%d animals: %.3f)\n', nx, 2/2^nx); end

        % within-bout strength: mean e per bout, post vs pre adjusted for Days
        SS = local_strengthPersist(St);
        fprintf('  WITHIN-BOUT strength (mean e per bout, %d bouts): post-pre %+.3f, LME p=%.4f, random-slope p=%.4f; animal-level %d/%d up, p=%.4f\n', ...
            SS.n, SS.b, SS.p, SS.pRS, sum(SS.dAnimal>0), numel(SS.dAnimal), SS.pAnimal);

        RES.persist.(matlab.lang.makeValidName(tag)) = struct('exit',fx,'entry',fn, ...
            'strength',SS,'strengthTbl',St, ...
            'exitTbl',Xt,'entryTbl',Nt,'animalExit',ax,'animalEntry',an, ...
            'pExit',local_signflip(ax),'pEntry',local_signflip(an), ...
            'minRunWin',median(epochS.minRunUsed), ...
            'minElapsedSec',(median(epochS.minRunUsed)-1)*STEP + WIDTH);
    end

    % threshold-free: excess synchrony distribution, no bout definition at all
    H = arrayfun(@(r) boutHazard(r, STEP, MINRUNS{1}), E, 'UniformOutput', false);
    H = [H{:}];
    fprintf('\n---- threshold-free excess synchrony e = (Ct-p50)/(p97.5-p50) ----\n');
    TF = struct();
    for fld = {'eMeanPos','eQ90','fracAbove'}
        y = [H.(fld{1})]';
        tb = table(y, categorical(ea(:)), categorical(phase(:),{'pre','post'}), days(:), ...
            'VariableNames', {'M','Animal','Phase','Days'});
        tb = tb(isfinite(tb.M),:);
        lme = fitlme(tb, 'M ~ Phase + Days + (1|Animal)');
        av = anova(lme, 'DFMethod','Satterthwaite');
        b = lme.Coefficients.Estimate(strcmp(lme.CoefficientNames,'Phase_post'));
        pPh = av.pValue(strcmp(av.Term,'Phase'));
        dTF = local_animalDiff(tb.M, cellstr(tb.Animal), tb.Phase=='pre', tb.Phase=='post');
        fprintf('  %-10s pre %.3f  post %.3f  Phase(adj. Days) %+.3f  p = %.4f   animals %d/%d up, p=%.4f\n', fld{1}, ...
            mean(tb.M(tb.Phase=='pre')), mean(tb.M(tb.Phase=='post')), b, pPh, ...
            sum(dTF>0), numel(dTF), local_signflip(dTF));
        TF.(fld{1}) = struct('pre',mean(tb.M(tb.Phase=='pre')), 'post',mean(tb.M(tb.Phase=='post')), ...
            'b',b, 'p',pPh, 'dAnimal',dTF, 'pAnimal',local_signflip(dTF));
    end
    RES.persist.threshFree = H;
    RES.persist.threshFreeStats = TF;

    %% ---- days-in-study checks: drug-naive epochs only ---------------------
    % Post-exposure epochs are always later, so Phase and time in study are
    % collinear, and the Days term in the models above is estimated partly
    % from post-exposure data. The drug-naive epochs span days of recording
    % BEFORE any ASV, so slow drift can be measured where the drug cannot
    % have caused it:
    %  (a) DRIFT      hazard ~ Days within naive epochs only, random Days
    %                 slope by animal; per-animal Poisson slopes, sign-flip.
    %  (b) CORRECTED  the naive slope is fixed into the offset
    %                 (log exposure + b_naive*Days) and Phase refit with no
    %                 Days term; animal level = per-animal log HR minus
    %                 b_naive * that animal's post-minus-pre day gap. Assumes
    %                 naive drift continues linearly into the post period --
    %                 an extrapolation, so the day ranges are printed.
    %  (c) PLACEBO    vehicle always preceded ASV, so the first SAL infusion
    %                 is treated as a sham exposure among naive epochs, with
    %                 the same phase rule and Phase + Days model. A
    %                 "persistent effect" of saline is a time effect.
    isNaive = strcmp(phase, 'pre');
    pseudo = repmat({''}, numel(E), 1); dLastS = nan(numel(E), 1);
    for i = 1:numel(uan)
        m = find(strcmp(ea, uan{i}));
        salD = unique(ed(m(~cellfun(@isempty, regexp(ec(m), '^SAL$', 'once')))));
        for j = m'
            if ~isNaive(j), continue; end
            if eInf(j), prior = salD(salD < ed(j)); else, prior = salD(salD <= ed(j)); end
            if isempty(prior), pseudo{j} = 'pre';
            else, pseudo{j} = 'post'; dLastS(j) = ed(j) - max(prior); end
        end
    end
    selP = isNaive & ~cellfun(@isempty, pseudo) & ~(strcmp(pseudo,'post') & dLastS < 1);
    uaN = unique(ea(isNaive));
    nSpan = sum(arrayfun(@(k) numel(unique(days(isNaive & strcmp(ea,uaN{k})))) >= 2, 1:numel(uaN)));
    uaP = unique(ea(selP));
    nBothP = sum(arrayfun(@(k) any(selP & strcmp(ea,uaP{k}) & strcmp(pseudo,'pre')) && ...
                               any(selP & strcmp(ea,uaP{k}) & strcmp(pseudo,'post')), 1:numel(uaP)));

    fprintf('\n================ DAYS-IN-STUDY CHECKS (drug-naive epochs) ================\n');
    fprintf('days in study: naive median %.0f (range %.0f-%.0f), post median %.0f (range %.0f-%.0f)\n', ...
        median(days(isNaive)), min(days(isNaive)), max(days(isNaive)), ...
        median(days(~isNaive)), min(days(~isNaive)), max(days(~isNaive)));
    fprintf('%d naive epochs, %d animals with naive epochs on >=2 days (drift estimable)\n', ...
        sum(isNaive), nSpan);
    fprintf('placebo: %d pre-SAL, %d post-SAL (>=1 d) naive epochs; %d animals with both\n', ...
        sum(selP & strcmp(pseudo,'pre')), sum(selP & strcmp(pseudo,'post')), nBothP);

    RES.persist.drift = struct(); RES.persist.placebo = struct();
    for mr = 1:numel(MINRUNS)
        minRun = MINRUNS{mr}; tag = local_mrTag(minRun); vt = matlab.lang.makeValidName(tag);
        fprintf('\n---- %s ----\n', tag);
        [Xn, Nn] = local_stack(E(isNaive), STEP, minRun, ea(isNaive), phase(isNaive), days(isNaive), [], []);
        [Xp, Np] = local_stack(E(selP),    STEP, minRun, ea(selP),    pseudo(selP),  days(selP),    [], []);
        D = struct(); PL = struct();
        for hz = {'exit','entry'}
            if strcmp(hz{1},'exit')
                tn = Xn; tp = Xp; tr = RES.persist.(vt).exitTbl;
            else
                tn = Nn; tp = Np; tr = RES.persist.(vt).entryTbl;
            end
            % (a) naive drift
            Fd = local_fitHazard(tn, 'Days', 'Days', 'Days', 'Days');
            s = local_animalSlope(tn);
            fprintf('  (a) %-5s drift, naive only : HR/day %s;  animals %d/%d up, p = %.4f\n', ...
                hz{1}, local_hrStr(Fd), sum(s>0), numel(s), local_signflip(s));
            % (b) phase, corrected by the naive drift
            b = local_logHR(Fd);
            Fc = local_blankFit(); dAdj = []; gap = [];
            if isfinite(b)
                tc = tr; tc.OffsetAdd = b * tc.Days;
                Fc = local_fitHazard(tc, 'Phase', 'Phase', 'Phase_post', 'Phase');
                [dRaw, ~, ids] = local_animalLogRatio(tr, 'Phase', 'pre', 'post');
                gap = arrayfun(@(k) local_wDays(tr, ids{k}, 'post') - local_wDays(tr, ids{k}, 'pre'), ...
                    (1:numel(ids))');
                dAdj = dRaw - b * gap;
                fprintf('  (b) %-5s phase, drift-corrected: HR %s;  animals %d/%d up, p = %.4f  (median gap %.0f d)\n', ...
                    hz{1}, local_hrStr(Fc), sum(dAdj>0), numel(dAdj), local_signflip(dAdj), median(gap));
            else
                fprintf('  (b) %-5s phase, drift-corrected: n/a (no drift estimate)\n', hz{1});
            end
            % (c) placebo
            Fp = local_fitHazard(tp, 'Phase + Days', 'Phase', 'Phase_post', 'Phase');
            [dP, ~] = local_animalLogRatio(tp, 'Phase', 'pre', 'post');
            fprintf('  (c) %-5s PLACEBO (SAL sham)  : HR %s;  animals %d/%d up, p = %.4f\n', ...
                hz{1}, local_hrStr(Fp), sum(dP>0), numel(dP), local_signflip(dP));
            D.(hz{1}) = struct('fit',Fd, 'animalSlope',s, 'pAnimal',local_signflip(s), ...
                'corrected',struct('fit',Fc, 'dAnimal',dAdj, 'gap',gap, 'pAnimal',local_signflip(dAdj)));
            PL.(hz{1}) = struct('fit',Fp, 'dAnimal',dP, 'pAnimal',local_signflip(dP));
        end
        RES.persist.drift.(vt) = D; RES.persist.placebo.(vt) = PL;
    end

    fprintf('\n---- threshold-free: drift, drift-corrected phase, placebo ----\n');
    TFD = struct();
    for fld = {'eMeanPos','eQ90','fracAbove'}
        y = [H.(fld{1})]'; g = isfinite(y);
        % (a) naive drift
        tn = table(y(g & isNaive), categorical(ea(g & isNaive)), days(g & isNaive), ...
            'VariableNames', {'M','Animal','Days'});
        [bN, pN] = local_lmeTerm(tn, 'M ~ Days + (1+Days|Animal)', 'Days');
        sA = local_animalOLS(tn);
        % (b) corrected per-animal phase difference
        ph = phase(g); aa_ = ea(g);
        dRaw = local_animalDiff(y(g),    aa_, strcmp(ph,'pre'), strcmp(ph,'post'));
        gapE = local_animalDiff(days(g), aa_, strcmp(ph,'pre'), strcmp(ph,'post'));
        dAdj = dRaw - bN * gapE;
        % (c) placebo
        gP = g & selP;
        tq = table(y(gP), categorical(ea(gP)), categorical(pseudo(gP),{'pre','post'}), days(gP), ...
            'VariableNames', {'M','Animal','Phase','Days'});
        [bP, pP] = local_lmeTerm(tq, 'M ~ Phase + Days + (1|Animal)', 'Phase');
        dPl = local_animalDiff(tq.M, cellstr(tq.Animal), tq.Phase=='pre', tq.Phase=='post');
        fprintf(['  %-10s drift/day %+.4f (p=%.3f, %d/%d up p=%.3f) | corrected phase %+.3f, %d/%d up p=%.3f ', ...
                 '| placebo %+.3f (p=%.3f, %d/%d up p=%.3f)\n'], fld{1}, ...
            bN, pN, sum(sA>0), numel(sA), local_signflip(sA), ...
            mean(dAdj), sum(dAdj>0), numel(dAdj), local_signflip(dAdj), ...
            bP, pP, sum(dPl>0), numel(dPl), local_signflip(dPl));
        TFD.(fld{1}) = struct('driftB',bN, 'driftP',pN, 'driftAnimal',sA, 'driftPAnimal',local_signflip(sA), ...
            'corrD',dAdj, 'corrPAnimal',local_signflip(dAdj), ...
            'placeboB',bP, 'placeboP',pP, 'placeboD',dPl, 'placeboPAnimal',local_signflip(dPl));
    end
    RES.persist.threshFreeDrift = TFD;

    if RUN_RATES
        [ovE, inE, raE] = local_boutRates(E, WIDTH, STEP, MINRUNS{1});
        fprintf('\n---- firing rate, bouts = %s (matched units) ----\n', local_mrTag(MINRUNS{1}));
        fprintf('  %-17s %8s %8s %11s %8s %16s\n','quantity','naive','post','Phase(adj)','p_LME','animal-level');
        Q = {'overall (Hz)', ovE; 'in-bout (Hz)', inE; 'log in/overall', log(raE)};
        for k = 1:size(Q,1)
            y = Q{k,2}(:); g = isfinite(y);
            tb = table(y(g), categorical(ea(g)), categorical(phase(g),{'pre','post'}), days(g), ...
                'VariableNames', {'M','Animal','Phase','Days'});
            lme = fitlme(tb, 'M ~ Phase + Days + (1|Animal)');
            av = anova(lme, 'DFMethod','Satterthwaite');
            b = lme.Coefficients.Estimate(strcmp(lme.CoefficientNames,'Phase_post'));
            d = local_animalDiff(y, ea, strcmp(phase,'pre'), strcmp(phase,'post'));
            fprintf('  %-17s %8.3f %8.3f %+11.3f %8.4f   %d/%d up, p=%.4f\n', Q{k,1}, ...
                mean(y(g & strcmp(phase,'pre'))), mean(y(g & strcmp(phase,'post'))), b, ...
                av.pValue(strcmp(av.Term,'Phase')), sum(d>0), numel(d), local_signflip(d));
        end
        fprintf('  (in-bout and ratio exist only for epochs with >=1 bout: %d/%d)\n', ...
            sum(isfinite(inE)), numel(inE));
        RES.persist.rates = struct('overall',ovE,'inBout',inE,'ratio',raE);
    end
    RES.persist.meta = struct('animal',{ea},'phase',{phase},'days',days,'dLast',dLast);
end

%% =========================================================================
%% ACUTE (pre-infusion vs 4-6 h, ASV against vehicle)
%% =========================================================================
if RUN_ACUTE
    fprintf('\n================ ACUTE (pre vs %g-%g h post-infusion) ================\n', ...
        LAT/3600, (LAT+POSTLEN)/3600);
    isInf = hasInf & ismember(cnd, {'ASV','SAL'}) & scOK;
    pre  = [zeros(numel(fils),1), inf_t];
    post = [inf_t + LAT, min(durs, inf_t + LAT + POSTLEN)];
    use = isInf & (pre(:,2)-pre(:,1)) >= MIN_EPOCH & (post(:,2)-post(:,1)) >= MIN_EPOCH;
    fprintf('%d infusion sessions with >= %.0f min both before infusion and in the window\n', ...
        sum(use), MIN_EPOCH/60);
    fprintf('  ASV: %d sessions / %d animals;  SAL: %d sessions / %d animals\n', ...
        sum(use & strcmp(cnd,'ASV')), numel(unique(ani(use & strcmp(cnd,'ASV')))), ...
        sum(use & strcmp(cnd,'SAL')), numel(unique(ani(use & strcmp(cnd,'SAL')))));

    uf = fils(use);
    mapPre  = containers.Map(uf, num2cell(pre(use,:), 2));
    mapPost = containers.Map(uf, num2cell(post(use,:), 2));
    % Units chosen ONCE per session from pre+post spikes and used in both
    % epochs. Previously the rate band was applied inside each Window, so the
    % pre and post traces of a session could be built from different neurons
    % (e.g. CP7_241109: 15 units pass pre, 19 post) -- a drug-induced rate
    % change would then alter WHICH units are analysed, not just their sync.
    mapU = containers.Map(uf, arrayfun(@(k) [pre(k,:); post(k,:)], find(use), ...
        'UniformOutput', false));
    uw = {'UnitWindow', @(f) mapU(f)};
    TA = syncBatch(uf, cfgB{:}, uw{:}, 'Window', @(f) mapPre(f), ...
        'CacheDir', fullfile(CACHE_ROOT,'hzAcutePre'), ...
        'CheckpointFile', fullfile(CACHE_ROOT,['ckpt_hzAcutePre' SAVE_TAG '.mat']));
    TB = syncBatch(uf, cfgB{:}, uw{:}, 'Window', @(f) mapPost(f), ...
        'CacheDir', fullfile(CACHE_ROOT,'hzAcutePost'), ...
        'CheckpointFile', fullfile(CACHE_ROOT,['ckpt_hzAcutePost' SAVE_TAG '.mat']));

    ok = [TA.ok]' & [TB.ok]';
    fprintf('%d sessions usable in both epochs\n', sum(ok));
    same = arrayfun(@(k) isequal(TA(k).detIdx, TB(k).detIdx), find(ok));
    fprintf('identical unit subset pre/post: %d/%d sessions\n', sum(same), numel(same));
    assert(all(same), 'pre and post epochs use different units');
    aa = ani(use); aa = aa(ok); dd = cnd(use); dd = dd(ok);
    E2 = [TA(ok); TB(ok)];
    per = [repmat({'pre'}, sum(ok), 1); repmat({'post'}, sum(ok), 1)];
    an2 = [aa; aa]; dr2 = [dd; dd];

    RES.acute = struct();
    if RUN_RATES
        [ovA, inA, raA] = local_boutRates(TA(ok), WIDTH, STEP, MINRUNS{1});
        [ovB, inB, raB] = local_boutRates(TB(ok), WIDTH, STEP, MINRUNS{1});
        fprintf('\n---- firing rate, bouts = %s (same units pre and post) ----\n', local_mrTag(MINRUNS{1}));
        Q = {'overall (Hz)', ovA, ovB; 'in-bout (Hz)', inA, inB; 'log in/overall', log(raA), log(raB)};
        for k = 1:size(Q,1)
            fprintf('  %s\n', Q{k,1});
            dsess = Q{k,3}(:) - Q{k,2}(:);
            dArm = struct();
            for arm = {'ASV','SAL'}
                m = strcmp(dd, arm{1}) & isfinite(dsess);
                [dA, ids] = local_animalMean(dsess, aa, m);
                dArm.(arm{1}) = struct('d', dA, 'ids', {ids});
                fprintf('    %-4s pre %7.3f  post %7.3f  change %+7.3f   %d/%d animals up, p=%.4f\n', ...
                    arm{1}, mean(Q{k,2}(m)), mean(Q{k,3}(m)), mean(dsess(m)), ...
                    sum(dA>0), numel(dA), local_signflip(dA));
            end
            [both, ia, is] = intersect(dArm.ASV.ids, dArm.SAL.ids);
            if numel(both) >= 2
                dp = dArm.ASV.d(ia) - dArm.SAL.d(is);
                fprintf('    paired DiD (ASV-SAL): %+7.3f, %d/%d > 0, exact p=%.4f (floor %.3f)\n', ...
                    mean(dp), sum(dp>0), numel(dp), local_signflip(dp), 2/2^numel(dp));
            end
        end
        fprintf('  (in-bout/ratio need >=1 bout in BOTH epochs: %d/%d sessions)\n', ...
            sum(isfinite(inA) & isfinite(inB)), numel(inA));
        RES.acute.rates = struct('pre',[ovA inA raA],'post',[ovB inB raB], ...
            'animal',{aa},'drug',{dd});
    end
    for mr = 1:numel(MINRUNS)
        minRun = MINRUNS{mr};
        [Xt, Nt, eS, St] = local_stack(E2, STEP, minRun, an2, [], [], per, dr2);
        tag = local_mrTag(minRun);
        fprintf('\n---- bout definition: minRun = %s (median %g windows = %g s elapsed minimum) ----\n', ...
            tag, median(eS.minRunUsed), (median(eS.minRunUsed)-1)*STEP + WIDTH);
        fprintf('%d bouts, %d out-intervals\n', height(Xt), height(Nt));

        % random slope by animal x drug: each animal's ASV change and SAL
        % change vary independently -- the variation the DiD is tested against
        fx = local_fitHazard(Xt, 'Drug*Period', 'Drug:Period', 'Drug_ASV:Period_post', 'Period', 'AnimalDrug');
        fn = local_fitHazard(Nt, 'Drug*Period', 'Drug:Period', 'Drug_ASV:Period_post', 'Period', 'AnimalDrug');
        local_report('EXIT  (leave bout), ASV-vs-SAL change', fx, Xt, 'Drug:Period');
        local_report('ENTRY (start bout), ASV-vs-SAL change', fn, Nt, 'Drug:Period');

        % Within-arm pre->post. With few vehicle animals the DiD above has no
        % real control; the ASV-arm change is what the data can support, and
        % the SAL-arm change is reported so the time-course confound is
        % visible rather than assumed away.
        ARM = struct();
        for arm = {'ASV','SAL'}
            xa = Xt(Xt.Drug==arm{1},:); na = Nt(Nt.Drug==arm{1},:);
            nA = numel(unique(xa.Animal));
            fprintf('  -- %s arm only (%d animals) --\n', arm{1}, nA);
            fxa = local_fitHazard(xa, 'Period', 'Period', 'Period_post', 'Period');
            fna = local_fitHazard(na, 'Period', 'Period', 'Period_post', 'Period');
            local_report(sprintf('EXIT  %s pre->post', arm{1}), fxa, xa, 'Period');
            local_report(sprintf('ENTRY %s pre->post', arm{1}), fna, na, 'Period');
            [dx, nx] = local_animalLogRatio(xa, 'Period', 'pre', 'post');
            [dn, nn] = local_animalLogRatio(na, 'Period', 'pre', 'post');
            ARM.(arm{1}) = struct('exit',fxa,'entry',fna,'nAnimals',nA, ...
                'dExit',dx,'dEntry',dn,'pExit',local_signflip(dx),'pEntry',local_signflip(dn));
            fprintf('     animal-level log HR exit  %+.3f, %d/%d negative, p = %.4f\n', ...
                mean(dx), sum(dx<0), nx, local_signflip(dx));
            fprintf('     animal-level log HR entry %+.3f, %d/%d positive, p = %.4f\n', ...
                mean(dn), sum(dn>0), nn, local_signflip(dn));
            if nx < 6, fprintf('     (exact floor at n=%d: %.3f)\n', nx, 2/2^max(nx,1)); end
        end

        % animal-level DiD backstop: per animal, log(post/pre) rate within each
        % drug, two-sample permutation on the per-animal values. Stored so the
        % sweep can report it next to the (optimistic) model LRT.
        for hz = {'exit','entry'}
            if strcmp(hz{1},'exit'), tt = Xt; else, tt = Nt; end
            [dA, ~] = local_animalLogRatio(tt(tt.Drug=='ASV',:), 'Period', 'pre', 'post');
            [dS, ~] = local_animalLogRatio(tt(tt.Drug=='SAL',:), 'Period', 'pre', 'post');
            pD = NaN;
            if numel(dA) >= 2 && numel(dS) >= 2
                pD = permtest(dS, dA);
                fprintf('  animal-level %-5s DiD: ASV %+.3f (n=%d), SAL %+.3f (n=%d), p = %.4f\n', ...
                    hz{1}, mean(dA), numel(dA), mean(dS), numel(dS), pD);
            end
            ARM.(['did_' hz{1}]) = struct('dASV',dA,'dSAL',dS, ...
                'diff',mean(dA)-mean(dS),'p',pD);
        end

        % PAIRED DiD: animals with both a SAL and an ASV session. The two-
        % sample tests above treat the arms as independent animals, which
        % they mostly are not (SAL precedes ASV in the same mice).
        for hz = {'exit','entry'}
            if strcmp(hz{1},'exit'), tt = Xt; else, tt = Nt; end
            [dA, ~, iA] = local_animalLogRatio(tt(tt.Drug=='ASV',:), 'Period', 'pre', 'post');
            [dS, ~, iS] = local_animalLogRatio(tt(tt.Drug=='SAL',:), 'Period', 'pre', 'post');
            [both, ia, is] = intersect(iA, iS);
            dpair = dA(ia) - dS(is);
            if numel(dpair) >= 2
                fprintf('  paired %-5s DiD (ASV-SAL, same animal): mean %+.3f, %d/%d > 0, exact sign-flip p = %.4f  (floor %.3f)  [%s]\n', ...
                    hz{1}, mean(dpair), sum(dpair>0), numel(dpair), local_signflip(dpair), 2/2^numel(dpair), strjoin(both', ','));
                ARM.(['paired_' hz{1}]) = struct('d',dpair,'p',local_signflip(dpair),'ids',{both});
            end
        end
        % within-bout strength (mean e per bout): per arm, then DiD
        SS = local_strengthAcute(St);
        fprintf('  WITHIN-BOUT strength (mean e per bout, %d bouts)\n', SS.n);
        for arm = {'ASV','SAL'}
            A = SS.(arm{1});
            fprintf('     %-4s post-pre %+.3f, LME p=%.4f, random-slope p=%.4f; animal-level %d/%d up, p=%.4f\n', ...
                arm{1}, A.b, A.p, A.pRS, sum(A.dAnimal>0), numel(A.dAnimal), A.pAnimal);
        end
        fprintf('     DiD  %+.3f, LME p=%.4f, random-slope p=%.4f; animal-level p=%.4f; paired n=%d p=%.4f\n', ...
            SS.did.b, SS.did.p, SS.did.pRS, SS.did.pAnimal, numel(SS.did.dPair), SS.did.pPair);

        RES.acute.(matlab.lang.makeValidName(tag)) = struct('exit',fx,'entry',fn, ...
            'strength',SS,'strengthTbl',St, ...
            'exitTbl',Xt,'entryTbl',Nt,'arm',ARM, ...
            'minRunWin',median(eS.minRunUsed),'minElapsedSec',(median(eS.minRunUsed)-1)*STEP + WIDTH);
    end
end

save(fullfile(CACHE_ROOT,['hazardResults' SAVE_TAG '.mat']), 'RES', 'MINRUNS', 'MATCHN', 'WIDTH', 'STEP', 'LAT', 'POSTLEN');
fprintf('\nsaved -> %s\n', fullfile(CACHE_ROOT,['hazardResults' SAVE_TAG '.mat']));


%% =========================================================================
function [Xt, Nt, epochS, St] = local_stack(E, step, minRun, anim, phase, days, period, drug)
% Stack exit and entry observations from every epoch into two tables, plus St:
% one row per bout with its mean/peak excess synchrony e (within-bout strength).
ex = {}; en = {}; sb = {}; mru = nan(numel(E),1);
for i = 1:numel(E)
    Hh = boutHazard(E(i), step, minRun);
    mru(i) = Hh.minRunUsed;
    sess = sprintf('S%04d', i);
    if isempty(period), sessKey = sess; else, sessKey = E(i).file; end  % acute: pre & post share a session
    meta = {anim{i}, sessKey};
    if ~isempty(phase), meta = [meta, phase(i), {days(i)}]; end
    if ~isempty(period), meta = [meta, period(i), drug(i)]; end
    for k = 1:size(Hh.exit,1),  ex(end+1,:) = [{Hh.exit(k,1),  Hh.exit(k,2)},  meta]; end %#ok<AGROW>
    for k = 1:size(Hh.entry,1), en(end+1,:) = [{Hh.entry(k,1), Hh.entry(k,2)}, meta]; end %#ok<AGROW>
    for k = 1:size(Hh.boutE,1)
        nw = Hh.boutWin(k,2) - Hh.boutWin(k,1) + 1;
        sb(end+1,:) = [{Hh.boutE(k,1), Hh.boutE(k,2), log(nw)}, meta]; end %#ok<AGROW>
end
if ~isempty(phase)
    vn = {'Exposure','Event','Animal','Session','Phase','Days'};
    vs = {'M','Pk','LogLen','Animal','Session','Phase','Days'};
else
    vn = {'Exposure','Event','Animal','Session','Period','Drug'};
    vs = {'M','Pk','LogLen','Animal','Session','Period','Drug'};
end
Xt = local_tbl(ex, vn); Nt = local_tbl(en, vn); St = local_tbl(sb, vs);
if ~isempty(period) && height(St) > 0
    St.AnimalDrug = categorical(strcat(string(St.Animal), '_', string(St.Drug)));
end
epochS.minRunUsed = mru;
end


function t = local_tbl(c, vn)
if isempty(c), t = cell2table(cell(0,numel(vn)), 'VariableNames', vn); return; end
t = cell2table(c, 'VariableNames', vn);
t.Animal  = categorical(t.Animal);
t.Session = categorical(t.Session);
if ismember('Phase', vn),  t.Phase  = categorical(t.Phase,  {'pre','post'}); end
if ismember('Period', vn), t.Period = categorical(t.Period, {'pre','post'}); end
if ismember('Drug', vn),   t.Drug   = categorical(t.Drug,   {'SAL','ASV'});  end
end


function F = local_fitHazard(varargin)
F = fitHazardGLMM(varargin{:});
end


function F = local_blankFit()
F = fitHazardGLMM();
end


function b = local_logHR(F)
% log HR of the tested coefficient, random-slope fit when available
if F.rs.ok, b = log(F.rs.HR); elseif F.ok, b = log(F.HR); else, b = NaN; end
end


function s = local_hrStr(F)
if F.rs.ok,  s = sprintf('%.3f [%.3f %.3f] LRT p %.4f', F.rs.HR, F.rs.CI(1), F.rs.CI(2), F.rs.pLRT);
elseif F.ok, s = sprintf('%.3f [%.3f %.3f] LRT p %.4f (no random slope)', F.HR, F.CI(1), F.CI(2), F.pLRT);
else,        s = sprintf('n/a (%s)', F.why); end
end


function s = local_animalSlope(t)
% per-animal Poisson slope of the hazard on Days (log HR per day), animals
% with observations on >=2 distinct days and >=3 events
ua = unique(t.Animal); s = [];
for i = 1:numel(ua)
    m = t.Animal == ua(i);
    if numel(unique(t.Days(m))) < 2 || sum(t.Event(m)) < 3, continue; end
    try
        g = fitglm(t.Days(m), t.Event(m), 'Distribution','poisson', 'Offset',log(t.Exposure(m)));
        s(end+1,1) = g.Coefficients.Estimate(2); %#ok<AGROW>
    catch
    end
end
end


function d = local_wDays(t, animal, ph)
% exposure-weighted mean day of one animal's observations in one phase
m = t.Animal == animal & t.Phase == ph;
d = sum(t.Days(m) .* t.Exposure(m)) / sum(t.Exposure(m));
end


function [b, p] = local_lmeTerm(tb, formula, term, coef)
% LME coefficient and Satterthwaite p for one term ('Phase' -> Phase_post);
% coef overrides the coefficient name (e.g. 'Drug_ASV:Period_post')
b = NaN; p = NaN;
if height(tb) < 5, return; end
try
    lme = fitlme(tb, formula);
    if nargin >= 4, cn = coef;
    elseif strcmp(term, 'Phase'), cn = 'Phase_post'; else, cn = term; end
    b = lme.Coefficients.Estimate(strcmp(lme.CoefficientNames, cn));
    av = anova(lme, 'DFMethod', 'Satterthwaite');
    p = av.pValue(strcmp(av.Term, term));
catch ME
    warning('lme:fail', '%s: %s', formula, ME.message);
end
end


function s = local_animalOLS(tb)
% per-animal least-squares slope of M on Days, animals with >=2 distinct days
ua = unique(tb.Animal); s = [];
for i = 1:numel(ua)
    m = tb.Animal == ua(i);
    if numel(unique(tb.Days(m))) < 2, continue; end
    c = polyfit(tb.Days(m), tb.M(m), 1);
    s(end+1,1) = c(1); %#ok<AGROW>
end
end


function local_report(label, F, t, term)
if ~F.ok, fprintf('  %s: NOT FIT -- %s\n', label, F.why); return; end
if any(strcmp(term,{'Phase','Period'}))
    g = t.(term);  lv = {'pre','post'};
    r1 = sum(t.Event(g==lv{1}))/sum(t.Exposure(g==lv{1}))*3600;
    r2 = sum(t.Event(g==lv{2}))/sum(t.Exposure(g==lv{2}))*3600;
    fprintf('  %s: crude %.2f -> %.2f /h;  HR %.3f [%.3f %.3f]  Wald p %.4f  LRT p %.4f%s\n', ...
        label, r1, r2, F.HR, F.CI(1), F.CI(2), F.pWald, F.pLRT, local_re(F));
else
    fprintf('  %s: ratio of HRs %.3f [%.3f %.3f]  Wald p %.4f  LRT p %.4f%s\n', ...
        label, F.HR, F.CI(1), F.CI(2), F.pWald, F.pLRT, local_re(F));
end
if F.rs.ok
    fprintf('      + random slope [%s]: HR %.3f [%.3f %.3f]  Wald p %.4f  LRT p %.4f\n', ...
        F.rs.re, F.rs.HR, F.rs.CI(1), F.rs.CI(2), F.rs.pWald, F.rs.pLRT);
elseif ~strcmp(F.rs.why, 'not requested')
    fprintf('      + random slope by animal: NOT FIT -- %s\n', F.rs.why);
end
end


function s = local_re(F)
if F.reSession, s = ''; else, s = '  [no session RE: fit failed with it]'; end
end


function [d, n, ids] = local_animalLogRatio(t, grp, a, b)
% Per animal: log( events/exposure in b  /  events/exposure in a ). The +0.5
% continuity correction keeps animals with zero events in one arm finite.
ua = unique(t.Animal); d = []; ids = {};
for i = 1:numel(ua)
    m = t.Animal == ua(i);
    ma = m & t.(grp) == a; mb = m & t.(grp) == b;
    if ~any(ma) || ~any(mb), continue; end
    % no events in either condition: the +0.5 correction would make the
    % "ratio" log(exposure_a/exposure_b) -- recording length, not hazard
    if sum(t.Event(ma)) + sum(t.Event(mb)) == 0, continue; end
    ra = (sum(t.Event(ma)) + 0.5) / sum(t.Exposure(ma));
    rb = (sum(t.Event(mb)) + 0.5) / sum(t.Exposure(mb));
    d(end+1,1) = log(rb/ra); %#ok<AGROW>
    ids{end+1,1} = char(ua(i)); %#ok<AGROW>
end
n = numel(d);
end


function p = local_signflip(d)
% Two-sided sign-flip test on the mean. Exact enumeration for n <= 14
% (2^14 patterns), Monte Carlo above.
d = d(isfinite(d)); n = numel(d);
if n < 2, p = NaN; return; end
obs = abs(mean(d));
if n <= 14
    S = 2*(dec2bin(0:2^n-1, n) == '1') - 1;     % all sign patterns
    p = mean(abs(S*d(:)/n) >= obs - 1e-12);
else
    nperm = 1e4; nul = zeros(nperm,1);
    for k = 1:nperm, nul(k) = mean((2*(rand(n,1)>0.5)-1).*d); end
    p = (1 + sum(abs(nul) >= obs)) / (1 + nperm);
end
end


function [ov, in, ra] = local_boutRates(T, width, step, minRun)
% Firing rate of the matched detection units overall and inside bouts, with
% bouts defined exactly as in the hazard models (boutHazard), not by the old
% r.transients(.significant) criterion.
n = numel(T); [ov, in, ra] = deal(nan(n,1));
for i = 1:n
    try
        Hh = boutHazard(T(i), step, minRun);
        o  = boutRate(T(i), width, step, Hh.boutWin);
        ov(i) = o.rateOverall; in(i) = o.rateInBout; ra(i) = o.ratio;
    catch ME
        warning('boutRate:fail', '%s: %s', T(i).file, ME.message);
    end
end
end


function [d, ids] = local_animalMean(x, anim, m)
% per-animal mean of x over rows m
ids = unique(anim(m)); d = nan(numel(ids),1);
for i = 1:numel(ids), d(i) = mean(x(m & strcmp(anim, ids{i}))); end
end


function [d, ids] = local_animalDiff(y, anim, mA, mB)
% per animal: mean(y in B) - mean(y in A), animals with finite data in both
ua = unique(anim); d = []; ids = {};
for i = 1:numel(ua)
    a = y(mA & strcmp(anim,ua{i}) & isfinite(y));
    b = y(mB & strcmp(anim,ua{i}) & isfinite(y));
    if ~isempty(a) && ~isempty(b)
        d(end+1,1) = mean(b) - mean(a); ids{end+1,1} = ua{i}; %#ok<AGROW>
    end
end
end


function S = local_strengthPersist(St)
% Within-bout strength, persistence: mean e per bout, post vs drug-naive.
% LME adjusted for Days (intercept-only and random-slope-by-animal), and the
% animal-level sign-flip on per-animal differences of bout-mean e (NOT adjusted
% for Days). Conditional on bouts existing in both phases.
S = struct('n',height(St),'b',NaN,'p',NaN,'bRS',NaN,'pRS',NaN,'dAnimal',[],'pAnimal',NaN);
if height(St) < 5, return; end
[S.b, S.p]     = local_lmeTerm(St, 'M ~ Phase + Days + (1|Animal) + (1|Session)', 'Phase');
[S.bRS, S.pRS] = local_lmeTerm(St, 'M ~ Phase + Days + (1+Phase|Animal) + (1|Session)', 'Phase');
an = cellstr(St.Animal);
S.dAnimal = local_animalDiff(St.M, an, St.Phase == 'pre', St.Phase == 'post');
S.pAnimal = local_signflip(S.dAnimal);
end


function S = local_strengthAcute(St)
% Within-bout strength, acute: mean e per bout, post vs pre infusion, per arm
% (LME intercept-only / random slope, animal-level sign-flip) and ASV-vs-SAL DiD
% (LME with random slope by animal x drug, animal-level permutation, paired).
S = struct('n',height(St));
blank = struct('b',NaN,'p',NaN,'bRS',NaN,'pRS',NaN,'dAnimal',[],'ids',{{}},'pAnimal',NaN);
S.ASV = blank; S.SAL = blank;
S.did = struct('b',NaN,'p',NaN,'bRS',NaN,'pRS',NaN,'pAnimal',NaN,'dPair',[],'pPair',NaN);
if height(St) < 5, return; end
for arm = {'ASV','SAL'}
    T = St(St.Drug == arm{1}, :);
    A = blank;
    [A.b, A.p]     = local_lmeTerm(T, 'M ~ Period + (1|Animal) + (1|Session)', 'Period', 'Period_post');
    [A.bRS, A.pRS] = local_lmeTerm(T, 'M ~ Period + (1+Period|Animal) + (1|Session)', 'Period', 'Period_post');
    [A.dAnimal, A.ids] = local_animalDiff(T.M, cellstr(T.Animal), T.Period == 'pre', T.Period == 'post');
    A.pAnimal = local_signflip(A.dAnimal);
    S.(arm{1}) = A;
end
dl = strengthDiD(St);
S.did.b = dl.b; S.did.p = dl.p; S.did.bRS = dl.bRS; S.did.pRS = dl.pRS;
if numel(S.ASV.dAnimal) >= 2 && numel(S.SAL.dAnimal) >= 2
    S.did.pAnimal = permtest(S.SAL.dAnimal, S.ASV.dAnimal);
end
[~, ia, is] = intersect(S.ASV.ids, S.SAL.ids);
S.did.dPair = S.ASV.dAnimal(ia) - S.SAL.dAnimal(is);
S.did.pPair = local_signflip(S.did.dPair);
end


function s = local_mrTag(m)
if ischar(m), s = m; else, s = sprintf('fixed%d', m); end
end