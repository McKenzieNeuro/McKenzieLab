%% runRippleAnalysis.m
% Do sharp-wave ripples change after ASV? Same sessions and epochs as
% runHazardAnalysis:
%   ACUTE        pre-infusion vs 4-6 h post-infusion, ASV arm, SAL arm, DiD
%   PERSISTENCE  drug-naive vs >=1 d post-exposure, drug-free periods, adj. Days
%
% MEASURES (rippleFeatures.m does the per-session work)
%   rate        ripples per hour of epoch (Poisson GLMM, offset log time)
%   rateNREM    ripples during NREM per hour of NREM (only if every session
%               has a *.SleepState.states.mat; otherwise skipped). Ripple rate
%               per recording hour is confounded with time asleep.
%   amp, ampZ   peak ripple-band envelope, raw and z-scored within session
%   freq        envelope-weighted mean instantaneous frequency (Hz)
%   logDur      log duration (s)
%   fracActive  fraction of the session's units firing >= 1 spike in the event
%   unitPPC     mean pairwise phase consistency of ripple-band phase at spike
%               times, units with >= MINSPK in-ripple spikes in the epoch
%   popPPC      PPC of all in-ripple spikes pooled
%   popPhase    preferred phase of pooled in-ripple spikes (circular;
%               animal-level test on the wrapped pre->post difference only)
%
% UNITS. One set per session, chosen by rate over the epochs analysed (acute:
% pre+post together), so a drug-induced rate change cannot change which units
% are counted. Rate band as in the synchrony pipeline.
%
% STATS (per metric, as in the hazard analysis)
%   epoch-level LME: one row per session x epoch (per-ripple rows are
%   pseudo-replicates of the session), intercept-only and random slope;
%   animal-level sign-flip on per-animal changes; DiD with random slope by
%   animal x drug, animal-level permutation and paired sign-flip.
%   Persistence: Phase + Days in the model; the animal-level p is NOT adjusted
%   for days in study.
%
% Per-session results are cached; delete CACHE_ROOT\ripples to recompute.

CACHE_ROOT = 'R:\TransPlasticity\syncCache';
WIDTH = 60; STEP = 5; MINUNITS = 10;           % screen settings = runHazardAnalysis (cache hit)
CONDS = {'Baseline','ASV','ASV24','SAL','SAL24'};
EXPOSURE_PAT = '^ASV$';
MIN_EPOCH = 1800;
LAT = 4*3600; POSTLEN = 2*3600;
RATE_RANGE = [0.1 10];
MINSPK = 20;                                   % in-ripple spikes for a unit's PPC
RUN_PERSIST = true;
RUN_ACUTE   = true;
FIRST_EXPOSURE_ONLY = true;   % acute: first ASV per animal; SAL only if given before it
RF_OPTS = {};                                  % e.g. {'Band',[120 200]}
VERSION = 2;                                   % 2: phase histograms per unit (rippleFeatures)
cdir = fullfile(CACHE_ROOT, 'ripples');
if ~exist(cdir, 'dir'), mkdir(cdir); end

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
bad = ismember(cnd, {'ASV','SAL'}) & ~hasInf;
fils = fils(~bad); cnd = cnd(~bad); ani = ani(~bad); dnm = dnm(~bad);
inf_t = inf_t(~bad); hasInf = hasInf(~bad);

% duration + usability per session (same rule as the old syncBatch screen, without
% running spikeSync over every unit -- that ran out of memory on large sessions)
[durs, scOK] = sessionScreen(fils, RATE_RANGE, MINUNITS, fullfile(CACHE_ROOT, 'screenLite'));

METRICS = {'amp','ampZ','freq','logDur','fracActive','unitPPC','popPPC'};
RIP = struct();

%% =========================================================================
%% ACUTE
%% =========================================================================
if RUN_ACUTE
    fprintf('\n================ RIPPLES, ACUTE (pre vs %g-%g h) ================\n', LAT/3600, (LAT+POSTLEN)/3600);
    pre  = [zeros(numel(fils),1), inf_t];
    post = [inf_t + LAT, min(durs, inf_t + LAT + POSTLEN)];
    use = hasInf & ismember(cnd, {'ASV','SAL'}) & scOK & ...
          (pre(:,2)-pre(:,1)) >= MIN_EPOCH & (post(:,2)-post(:,1)) >= MIN_EPOCH;
    if FIRST_EXPOSURE_ONLY
        % ASV: only each animal's FIRST ASV infusion session. SAL: only infusions
        % given before that animal's first ASV (all of them if it never had ASV).
        % Order by date, then by recording start clock when both are known.
        rs = cellfun(@recordingStartFromPath, fils);
        tKey = dnm; tKey(isfinite(rs)) = rs(isfinite(rs));
        firstASV = true(numel(fils),1); preASV = true(numel(fils),1);
        for a = unique(ani(:))'
            m = strcmp(ani, a{1});
            mA = find(m & strcmp(cnd, 'ASV') & hasInf);
            if isempty(mA), firstASV(m) = false; continue; end
            [~, o] = min(tKey(mA)); j0 = mA(o);
            firstASV(m) = false; firstASV(j0) = true;
            % before = earlier date, or same date with both start clocks known and earlier
            mm = find(m);
            preASV(mm) = dnm(mm) < dnm(j0) | (dnm(mm) == dnm(j0) & isfinite(rs(mm)) & ...
                isfinite(rs(j0)) & rs(mm) < rs(j0));
        end
        drop = use & ((strcmp(cnd,'ASV') & ~firstASV) | (strcmp(cnd,'SAL') & ~preASV));
        fprintf('first-exposure only: dropping %d later ASV and %d post-ASV SAL session(s)\n', ...
            sum(drop & strcmp(cnd,'ASV')), sum(drop & strcmp(cnd,'SAL')));
        use = use & ~drop;
    end
    idx = find(use);
    rows = {};
    for j = 1:numel(idx)
        i = idx(j);
        W = [pre(i,:); post(i,:)];
        Rr = local_cached(cdir, 'acute', fils{i}, W, W, RATE_RANGE, RF_OPTS, VERSION);
        if ~Rr.ok, fprintf('  skip %s: %s\n', fils{i}, Rr.why); continue; end
        es = local_epochSummary(Rr, MINSPK);
        for e = 1:2
            rows(end+1,:) = [{ani{i}, fils{i}, cnd{i}, ternary(e==1,'pre','post'), NaN, NaN}, es(e,:)]; %#ok<AGROW>
        end
    end
    TA = local_table(rows);
    fprintf('%d sessions (%d ASV / %d SAL), %d animals\n', height(TA)/2, ...
        sum(TA.Drug=='ASV')/2, sum(TA.Drug=='SAL')/2, numel(unique(TA.Animal)));
    RIP.acute = local_acuteStats(TA, METRICS);
    RIP.acute.tbl = TA;
end

%% =========================================================================
%% PERSISTENCE
%% =========================================================================
if RUN_PERSIST
    fprintf('\n================ RIPPLES, PERSISTENCE (drug-naive vs post) ================\n');
    w = [zeros(numel(fils),1), durs];
    w(hasInf,2) = min(inf_t(hasInf), durs(hasInf));
    use = scOK & (w(:,2) - w(:,1)) >= MIN_EPOCH;
    idx = find(use);
    ea = ani(idx); ed = dnm(idx); ec = cnd(idx); eInf = hasInf(idx);
    phase = repmat({'pre'}, numel(idx), 1); dLast = nan(numel(idx),1); days = nan(numel(idx),1);
    uan = unique(ea);
    for a = 1:numel(uan)
        m = find(strcmp(ea, uan{a}));
        days(m) = ed(m) - min(ed(m));
        expD = unique(ed(m(~cellfun(@isempty, regexp(ec(m), EXPOSURE_PAT,'once')))));
        for j = m'
            if eInf(j), prior = expD(expD < ed(j)); else, prior = expD(expD <= ed(j)); end
            if ~isempty(prior), phase{j} = 'post'; dLast(j) = ed(j) - max(prior); end
        end
    end
    sel = ~(strcmp(phase,'post') & dLast < 1);
    rows = {};
    for j = find(sel)'
        i = idx(j);
        W = w(i,:);
        Rr = local_cached(cdir, 'persist', fils{i}, W, W, RATE_RANGE, RF_OPTS, VERSION);
        if ~Rr.ok, fprintf('  skip %s: %s\n', fils{i}, Rr.why); continue; end
        es = local_epochSummary(Rr, MINSPK);
        rows(end+1,:) = [{ea{j}, fils{i}, ec{j}, phase{j}, days(j), NaN}, es(1,:)]; %#ok<AGROW>
    end
    TP = local_table(rows);
    TP.Phase = TP.Period; TP.Days = TP.DaysCol;
    fprintf('%d sessions (%d naive / %d post), %d animals\n', height(TP), ...
        sum(TP.Phase=='pre'), sum(TP.Phase=='post'), numel(unique(TP.Animal)));
    RIP.persist = local_persistStats(TP, METRICS);
    RIP.persist.tbl = TP;
end

save(fullfile(CACHE_ROOT, 'rippleResults.mat'), 'RIP', 'MINSPK', 'LAT', 'POSTLEN', 'RATE_RANGE', 'RF_OPTS');
fprintf('\nsaved -> %s\n', fullfile(CACHE_ROOT, 'rippleResults.mat'));


%% =========================================================================
function Rr = local_cached(cdir, ctx, f, W, UW, rr, rfOpts, ver)
h = typecast(int32(java.lang.String(f).hashCode()), 'uint32');
cf = fullfile(cdir, sprintf('%s_%010u.mat', ctx, h));
if exist(cf, 'file')
    C = load(cf, 'Rr', 'ver');
    % reuse only a successful result that has every field this version needs
    % (a cache written by an older rippleFeatures can carry the right version
    % number but lack the phase histograms); failures are always retried
    need = {'unitHist','unitHistC','phaseEdges','unitN','ev'};
    if isfield(C, 'ver') && C.ver == ver && isfield(C.Rr, 'ok') && C.Rr.ok && ...
            all(isfield(C.Rr, need)) && isequal(C.Rr.wins, W) && strcmp(C.Rr.file, f)
        Rr = C.Rr; return
    end
end
try
    Rr = rippleFeatures(f, W, 'UnitWin', UW, 'RateRange', rr, rfOpts{:});
catch ME
    Rr = struct('ok', false, 'why', ME.message, 'wins', W, 'file', f);
end
save(cf, 'Rr', 'ver');
end


function es = local_epochSummary(Rr, minSpk)
% one row per epoch: {nRip, T, nRipNREM, Tnrem, amp, ampZ, freq, logDur,
% fracActive, unitPPC, popPPC, popPhase, nUnitsPPC, nUnits,
% popHist, uHistC, uPref, uPPC}
%   popHist  1 x 36 counts of all in-ripple spike phases (0 = ripple-band peak)
%   uHistC   1 x 36 mean over units (>= minSpk) of each unit's phase density
%            centred on its own preferred phase: within-unit WIDTH only
%   uPref    1 x N preferred phase per unit (NaN below minSpk)
%   uPPC     1 x N pairwise phase consistency per unit (NaN below minSpk)
nW = size(Rr.wins, 1); es = cell(nW, 18);
for e = 1:nW
    v = Rr.ev(Rr.ev.Epoch == e, :);
    n = Rr.unitN(e,:); C = Rr.unitC(e,:); S = Rr.unitS(e,:);
    ppc = (C.^2 + S.^2 - n) ./ (n .* (n - 1));
    okU = n >= minSpk;
    nP = sum(n); CP = sum(C); SP = sum(S);
    if Rr.hasSleep, nN = sum(v.NREM == 1); else, nN = NaN; end
    es(e,:) = {height(v), Rr.epochT(e), nN, Rr.epochTnrem(e), ...
        mean(v.Amp,'omitnan'), mean(v.AmpZ,'omitnan'), mean(v.Freq,'omitnan'), ...
        mean(log(v.Dur),'omitnan'), mean(v.FracActive,'omitnan'), ...
        mean(ppc(okU)), (CP^2 + SP^2 - nP) / (nP*(nP - 1)), atan2(SP, CP), ...
        sum(okU), Rr.nUnits, ...
        reshape(sum(Rr.unitHist(e,:,:), 2), 1, []), local_centredMean(Rr.unitHistC(e,:,:), okU), ...
        {local_mask(atan2(S, C), okU)}, {local_mask(ppc, okU)}};
end
end


function h = local_centredMean(H, okU)
% mean over units of each unit's normalised, self-centred phase histogram
H = reshape(H, size(H,2), size(H,3));
H = H(okU, :);
if isempty(H), h = nan(1, size(H,2)); return; end
h = mean(H ./ sum(H, 2), 1);
end


function v = local_mask(v, ok)
v(~ok) = NaN;
end


function T = local_table(rows)
vn = {'Animal','Session','Drug','Period','DaysCol','Unused', ...
      'Event','Exposure','EventNREM','ExposureNREM','amp','ampZ','freq','logDur', ...
      'fracActive','unitPPC','popPPC','popPhase','nUnitsPPC','nUnits', ...
      'popHist','uHistC','uPref','uPPC'};
T = cell2table(rows, 'VariableNames', vn);
T.Unused = [];
T.Animal = categorical(T.Animal);
T.Session = categorical(T.Session);
dr = T.Drug; dr(~ismember(dr, {'ASV','SAL'})) = {'SAL'};   % persistence: Drug unused
T.Drug = categorical(dr, {'SAL','ASV'});
T.Period = categorical(T.Period, {'pre','post'});
T.AnimalDrug = categorical(strcat(cellstr(T.Animal), '_', cellstr(T.Drug)));
end


function A = local_acuteStats(T, metrics)
A = struct();
hasN = all(isfinite(T.ExposureNREM) & T.ExposureNREM > 0);
rm = {'rate'}; if hasN, rm{end+1} = 'rateNREM'; end
fprintf('\n%-11s | %-44s | %-44s | %-40s\n', 'metric', 'ASV arm: b (p LME | p RS) k/n up, anim p', ...
    'SAL arm', 'DiD: b (p LME | p RS) anim | paired');
for m = [rm, metrics, {'popPhase'}]
    nm = m{1};
    for arm = {'ASV','SAL'}
        Ta = T(T.Drug == arm{1}, :);
        A.(nm).(arm{1}) = local_one(Ta, nm, 'arm');
    end
    A.(nm).did = local_one(T, nm, 'did');
    A.(nm).did.pAnimal = NaN; A.(nm).did.dPair = []; A.(nm).did.pPair = NaN;
    dA = A.(nm).ASV.dAnimal; dS = A.(nm).SAL.dAnimal;
    if numel(dA) >= 2 && numel(dS) >= 2, A.(nm).did.pAnimal = permtest(dS, dA); end
    [~, ia, is] = intersect(A.(nm).ASV.ids, A.(nm).SAL.ids);
    if strcmp(nm, 'popPhase')
        dp = angle(exp(1i*(dA(ia) - dS(is))));
        A.(nm).did.b = angle(mean(exp(1i*dp)));
    else
        dp = dA(ia) - dS(is);
    end
    A.(nm).did.dPair = dp; A.(nm).did.pPair = local_signflip(dp);
    fprintf('%-11s | %-44s | %-44s | %-40s\n', nm, local_str(A.(nm).ASV, nm), local_str(A.(nm).SAL, nm), ...
        sprintf('%s %.3f | n=%d %.3f', local_bstr(A.(nm).did, nm), A.(nm).did.pAnimal, numel(dp), A.(nm).did.pPair));
end
fprintf(['rate: b = rate ratio post/pre; popPhase: circular mean pre->post shift (rad), animal-level only;\n' ...
         'other metrics: post - pre in their own units.\n']);
if ~hasN, fprintf('rateNREM skipped: not every session has a SleepState file.\n'); end
end


function P = local_persistStats(T, metrics)
P = struct();
hasN = all(isfinite(T.ExposureNREM) & T.ExposureNREM > 0);
rm = {'rate'}; if hasN, rm{end+1} = 'rateNREM'; end
fprintf('\n%-11s | %-44s\n', 'metric', 'post - pre adj. Days: b (p LME | p RS) k/n up, anim p');
for m = [rm, metrics, {'popPhase'}]
    P.(m{1}) = local_one(T, m{1}, 'persist');
    fprintf('%-11s | %-44s\n', m{1}, local_str(P.(m{1}), m{1}));
end
fprintf('animal-level p is not adjusted for days in study.\n');
if ~hasN, fprintf('rateNREM skipped: not every session has a SleepState file.\n'); end
end


function O = local_one(T, nm, ctx)
% one metric, one contrast: model (intercept-only / random slope) + animal level
O = struct('b',NaN,'p',NaN,'bRS',NaN,'pRS',NaN,'dAnimal',[],'ids',{{}},'pAnimal',NaN);
if strcmp(ctx, 'persist'), g = 'Phase'; else, g = 'Period'; end
isRate = any(strcmp(nm, {'rate','rateNREM'}));
isCirc = strcmp(nm, 'popPhase');
if isRate
    t = T;
    if strcmp(nm, 'rateNREM'), t.Event = t.EventNREM; t.Exposure = t.ExposureNREM; end
    t = t(t.Exposure > 0, :);
    t.Exposure = t.Exposure / 3600;              % rate per hour
    switch ctx
        case 'arm',     F = fitHazardGLMM(t, 'Period', 'Period', 'Period_post', 'Period');
        case 'did',     F = fitHazardGLMM(t, 'Drug*Period', 'Drug:Period', 'Drug_ASV:Period_post', 'Period', 'AnimalDrug');
        case 'persist', F = fitHazardGLMM(t, 'Phase + Days', 'Phase', 'Phase_post', 'Phase');
    end
    if F.ok,    O.b = F.HR;    O.p = F.pLRT;     end
    if F.rs.ok, O.bRS = F.rs.HR; O.pRS = F.rs.pLRT; end
    if ~strcmp(ctx, 'did')
        [O.dAnimal, O.ids] = local_animalRate(t, g);
        O.pAnimal = local_signflip(O.dAnimal);
    end
    return
end
t = T; t.y = t.(nm); t = t(isfinite(t.y), :);
if ~isCirc
    switch ctx
        case 'arm'
            [O.b, O.p]     = local_lme(t, 'y ~ Period + (1|Animal) + (1|Session)', {'Period_post'});
            [O.bRS, O.pRS] = local_lme(t, 'y ~ Period + (1+Period|Animal) + (1|Session)', {'Period_post'});
        case 'did'
            [O.b, O.p]     = local_lme(t, 'y ~ Drug*Period + (1|Animal) + (1|Session)', {'Drug_ASV','Period_post'});
            [O.bRS, O.pRS] = local_lme(t, 'y ~ Drug*Period + (1|Animal) + (1+Period|AnimalDrug) + (1|Session)', {'Drug_ASV','Period_post'});
        case 'persist'
            [O.b, O.p]     = local_lme(t, 'y ~ Phase + Days + (1|Animal)', {'Phase_post'});
            [O.bRS, O.pRS] = local_lme(t, 'y ~ Phase + Days + (1+Phase|Animal)', {'Phase_post'});
    end
end
if ~strcmp(ctx, 'did')
    [O.dAnimal, O.ids] = local_animalDiff(t, g, isCirc, ~strcmp(ctx, 'persist'));
    O.pAnimal = local_signflip(O.dAnimal);
    if isCirc, O.b = angle(mean(exp(1i*O.dAnimal))); end
end
end


function [b, p] = local_lme(t, formula, parts)
b = NaN; p = NaN;
if height(t) < 5, return; end
try
    lme = fitlme(t, formula);
    cn = lme.CoefficientNames;
    k = find(cellfun(@(c) all(cellfun(@(q) contains(c, q), parts)) && ...
        numel(strsplit(c, ':')) == numel(parts), cn), 1);
    if isempty(k), return; end
    b = lme.Coefficients.Estimate(k); p = lme.Coefficients.pValue(k);
    tp = cellfun(@(q) strtok(q, '_'), parts, 'UniformOutput', false);
    av = anova(lme, 'DFMethod', 'Satterthwaite');
    j = find(cellfun(@(c) all(cellfun(@(q) contains(c, q), tp)) && ...
        numel(strsplit(c, ':')) == numel(tp), av.Term), 1);
    if ~isempty(j), p = av.pValue(j); end
catch ME
    warning('ripple:lme', '%s: %s', formula, ME.message);
end
end


function [d, ids] = local_animalDiff(t, g, isCirc, paired)
% per animal change pre -> post. paired (acute): per-session post-pre, then
% animal mean; unpaired (persistence): animal mean(post) - mean(pre).
% Circular metrics: wrapped differences and circular means.
ua = unique(t.Animal); d = []; ids = {};
for i = 1:numel(ua)
    m = t.Animal == ua(i);
    if paired
        us = unique(t.Session(m)); x = [];
        for s = us'
            a = t.y(m & t.Session == s & t.(g) == 'pre');
            b = t.y(m & t.Session == s & t.(g) == 'post');
            if numel(a) == 1 && numel(b) == 1
                if isCirc, x(end+1,1) = angle(exp(1i*(b - a))); else, x(end+1,1) = b - a; end %#ok<AGROW>
            end
        end
        if isempty(x), continue; end
        if isCirc, v = angle(mean(exp(1i*x))); else, v = mean(x); end
    else
        a = t.y(m & t.(g) == 'pre'); b = t.y(m & t.(g) == 'post');
        if isempty(a) || isempty(b), continue; end
        if isCirc, v = angle(mean(exp(1i*b)) / mean(exp(1i*a)));
        else, v = mean(b) - mean(a); end
    end
    d(end+1,1) = v; ids{end+1,1} = char(ua(i)); %#ok<AGROW>
end
end


function [d, ids] = local_animalRate(t, g)
% per animal log rate ratio post/pre, events pooled over sessions, +0.5
ua = unique(t.Animal); d = []; ids = {};
for i = 1:numel(ua)
    m = t.Animal == ua(i);
    ma = m & t.(g) == 'pre'; mb = m & t.(g) == 'post';
    if ~any(ma) || ~any(mb), continue; end
    ra = (sum(t.Event(ma)) + 0.5) / sum(t.Exposure(ma));
    rb = (sum(t.Event(mb)) + 0.5) / sum(t.Exposure(mb));
    d(end+1,1) = log(rb / ra); ids{end+1,1} = char(ua(i)); %#ok<AGROW>
end
end


function p = local_signflip(d)
d = d(isfinite(d)); n = numel(d);
if n < 2, p = NaN; return; end
obs = abs(mean(d));
if n <= 14
    S = 2*(dec2bin(0:2^n-1, n) == '1') - 1;
    p = mean(abs(S*d(:)/n) >= obs - 1e-12);
else
    nul = arrayfun(@(k) mean((2*(rand(n,1) > 0.5) - 1) .* d), 1:1e4);
    p = (1 + sum(abs(nul) >= obs)) / (1 + 1e4);
end
end


function s = local_bstr(O, nm)
if any(strcmp(nm, {'rate','rateNREM'})), f = '%.2f'; else, f = '%+.4g'; end
if strcmp(nm, 'popPhase'), s = sprintf('%+.2f rad', O.b); return; end
s = sprintf([f ' (%.3f | %.3f)'], O.b, O.p, O.pRS);
end


function s = local_str(O, nm)
s = sprintf('%s %d/%d p=%.3f', local_bstr(O, nm), sum(O.dAnimal > 0), numel(O.dAnimal), O.pAnimal);
end


function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end