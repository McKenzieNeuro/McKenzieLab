%% printTimescaleSweep.m
% Prints the timescale-sweep tables. Called by runTimescaleSweep after the
% sweep is SAVED, so a printing error cannot lose a finished sweep. Run it on
% its own to reprint: it uses SWEEP from the workspace if present, otherwise
% loads timescaleSweep_N<MatchN>.mat.

if ~exist('CACHE_ROOT','var'),   CACHE_ROOT = 'C:\Users\samckenzie\syncCache'; end
if ~exist('MATCHN_SWEEP','var'), MATCHN_SWEEP = 10; end
if ~exist('SWEEP','var') || ~exist('SCALES','var')
    load(fullfile(CACHE_ROOT, sprintf('timescaleSweep_N%d.mat', MATCHN_SWEEP)), ...
        'SWEEP', 'SCALES', 'RUN_PERSIST_SWEEP');
end
if ~exist('RUN_PERSIST_SWEEP','var')
    RUN_PERSIST_SWEEP = isfield(SWEEP{1}, 'persist');
end

%% ---- ACUTE summary --------------------------------------------------------
fprintf('\n\n================ TIMESCALE SWEEP: ACUTE, MatchN %d ================\n', MATCHN_SWEEP);
fprintf('arm: HR (model p) k/n up, anim p.   DiD: ratio of HR (model p) | anim p | paired\n\n');

for hz = {'entry','exit'}
    fprintf('---------- %s hazard ----------\n', upper(hz{1}));
    fprintf('%-6s %-9s %7s | %-30s | %-30s | %-16s | %8s | %-13s\n', 'width','bout def','min(s)', ...
        'ASV arm', 'SAL arm', 'DiD model', 'DiD anim', 'DiD paired');
    for k = 1:numel(SWEEP)
        R = SWEEP{k}.acute;
        tags = local_boutTags(R);
        for t = 1:numel(tags)
            Q = R.(tags{t});
            fprintf('%-6g %-9s %7.0f | %-30s | %-30s | %-16s | %8s | %-13s\n', SCALES(k,1), tags{t}, Q.minElapsedSec, ...
                local_arm(Q.arm, 'ASV', hz{1}), local_arm(Q.arm, 'SAL', hz{1}), ...
                local_hr(Q.(hz{1})), local_didAnim(Q.arm, hz{1}), local_didPaired(Q.arm, hz{1}));
        end
        fprintf('\n');
    end
end

fprintf('null95, ASV-arm entry HR by width: ');
for k = 1:numel(SWEEP)
    A = SWEEP{k}.acute.null95.arm.ASV;
    fprintf('%g s: %s anim p=%.3f   ', SCALES(k,1), local_hr(A.entry), A.pEntry);
end
fprintf('\n');

%% ---- PERSISTENCE summary --------------------------------------------------
if RUN_PERSIST_SWEEP
    fprintf('\n\n======== TIMESCALE SWEEP: PERSISTENCE (naive vs >=24 h post), MatchN %d ========\n', MATCHN_SWEEP);
    fprintf('HR post/pre adjusted for Days (model p); k/n animals up, sign-flip p.\n\n');
    fprintf('%-6s %-9s %7s | %-32s | %-32s\n', 'width','bout def','min(s)', ...
        'ENTRY', 'EXIT');
    for k = 1:numel(SWEEP)
        P = SWEEP{k}.persist;
        tags = local_boutTags(P);
        for t = 1:numel(tags)
            Q = P.(tags{t});
            fprintf('%-6g %-9s %7.0f | %-32s | %-32s\n', SCALES(k,1), tags{t}, Q.minElapsedSec, ...
                local_cell(Q.entry, Q.animalEntry, Q.pEntry), local_cell(Q.exit, Q.animalExit, Q.pExit));
        end
        fprintf('\n');
    end

    fprintf('---------- threshold-free excess synchrony (no bout definition) ----------\n');
    fprintf('Phase effect adj. Days: b (LME p) k/n animals up, sign-flip p\n');
    flds = {'eMeanPos','eQ90','fracAbove'};
    fprintf('%-6s | %-32s | %-32s | %-32s\n', 'width', flds{:});
    for k = 1:numel(SWEEP)
        TF = SWEEP{k}.persist.threshFreeStats;
        c = cellfun(@(f) sprintf('%+.3f (%.3f) %d/%d p=%.3f', TF.(f).b, TF.(f).p, ...
            sum(TF.(f).dAnimal>0), numel(TF.(f).dAnimal), TF.(f).pAnimal), flds, 'UniformOutput', false);
        fprintf('%-6g | %-32s | %-32s | %-32s\n', SCALES(k,1), c{:});
    end

    fprintf('\n---------- days-in-study checks (drug-naive epochs) ----------\n');
    fprintf(['drift = HR per day within naive epochs; corrected = post/pre HR with the naive\n' ...
             'drift fixed in the offset; placebo = SAL infusion as sham exposure, naive epochs.\n' ...
             'A real persistent effect: corrected ~ uncorrected, placebo ~ 1.\n\n']);
    for hz = {'entry','exit'}
        fprintf('%s\n', upper(hz{1}));
        fprintf('%-6s %-9s | %-30s | %-30s | %-30s\n', 'width','bout def', ...
            'drift HR/day, k/n up p', 'corrected HR, k/n up p', 'PLACEBO HR, k/n up p');
        for k = 1:numel(SWEEP)
            P = SWEEP{k}.persist;
            tags = fieldnames(P.drift);
            for t = 1:numel(tags)
                D = P.drift.(tags{t}).(hz{1}); C = D.corrected; L = P.placebo.(tags{t}).(hz{1});
                fprintf('%-6g %-9s | %-30s | %-30s | %-30s\n', SCALES(k,1), tags{t}, ...
                    local_hrd(D.fit, D.animalSlope, D.pAnimal), ...
                    local_hrd(C.fit, C.dAnimal, C.pAnimal), ...
                    local_hrd(L.fit, L.dAnimal, L.pAnimal));
            end
        end
        fprintf('\n');
    end
    fprintf('threshold-free (eMeanPos | eQ90 | fracAbove): drift/day; corrected k/n up p; placebo b k/n up p\n');
    for k = 1:numel(SWEEP)
        T = SWEEP{k}.persist.threshFreeDrift;
        c = cellfun(@(f) sprintf('%+.4f; %d/%d p=%.3f; %+.3f %d/%d p=%.3f', T.(f).driftB, ...
            sum(T.(f).corrD>0), numel(T.(f).corrD), T.(f).corrPAnimal, ...
            T.(f).placeboB, sum(T.(f).placeboD>0), numel(T.(f).placeboD), T.(f).placeboPAnimal), ...
            {'eMeanPos','eQ90','fracAbove'}, 'UniformOutput', false);
        fprintf('%-6g | %s | %s | %s\n', SCALES(k,1), c{:});
    end

    fprintf('\nnull95 persistence entry HR by width: ');
    for k = 1:numel(SWEEP)
        Q = SWEEP{k}.persist.null95;
        fprintf('%g s: %s anim p=%.3f   ', SCALES(k,1), local_hr(Q.entry), Q.pEntry);
    end
    fprintf('\n');
end

%% ---- WITHIN-BOUT STRENGTH -------------------------------------------------
% mean excess synchrony e = (C-p50)/(p97.5-p50) over each bout's windows
% (e = 1 is the band edge); one row per bout, conditional on bouts existing.
if isfield(SWEEP{1}.acute.null95, 'strength')
    fprintf('\n\n======== WITHIN-BOUT SYNCHRONY STRENGTH (mean e per bout), MatchN %d ========\n', MATCHN_SWEEP);
    fprintf('post-pre in e units: b (p LME | p random slope) k/n animals up, animal p.  Persistence adj. Days (animal p is not).\n\n');
    fprintf('ACUTE\n');
    fprintf('%-6s %-9s %6s | %-40s | %-40s | %-34s\n', 'width','bout def','nBout', 'ASV arm', 'SAL arm', 'DiD b (p LME | p RS) animal | paired');
    for k = 1:numel(SWEEP)
        R = SWEEP{k}.acute; tags = local_boutTags(R);
        for t = 1:numel(tags)
            S = R.(tags{t}).strength;
            if isempty(S.did.b) && isfield(R.(tags{t}), 'strengthTbl')   % sweeps saved before the DiD name fix
                dl = strengthDiD(R.(tags{t}).strengthTbl);
                S.did.b = dl.b; S.did.p = dl.p; S.did.bRS = dl.bRS; S.did.pRS = dl.pRS;
            end
            fprintf('%-6g %-9s %6d | %-40s | %-40s | %-34s\n', SCALES(k,1), tags{t}, S.n, ...
                local_str(S.ASV), local_str(S.SAL), ...
                sprintf('%+.3f (%.3f | %.3f) %.3f | n=%d %.3f', S.did.b, S.did.p, S.did.pRS, ...
                    S.did.pAnimal, numel(S.did.dPair), S.did.pPair));
        end
        fprintf('\n');
    end
    if RUN_PERSIST_SWEEP && isfield(SWEEP{1}.persist.null95, 'strength')
        fprintf('PERSISTENCE (post vs drug-naive)\n');
        fprintf('%-6s %-9s %6s | %-40s\n', 'width','bout def','nBout', 'post-pre b (p LME | p RS) k/n up, animal p');
        for k = 1:numel(SWEEP)
            P = SWEEP{k}.persist; tags = local_boutTags(P);
            for t = 1:numel(tags)
                S = P.(tags{t}).strength;
                fprintf('%-6g %-9s %6d | %-40s\n', SCALES(k,1), tags{t}, S.n, ...
                    sprintf('%+.3f (%.3f | %.3f) %d/%d p=%.3f', S.b, S.p, S.pRS, ...
                        sum(S.dAnimal>0), numel(S.dAnimal), S.pAnimal));
            end
            fprintf('\n');
        end
    end
else
    fprintf('\n(no within-bout strength in this sweep file: re-run runTimescaleSweep)\n');
end

%% =========================================================================
function s = local_str(A)
s = sprintf('%+.3f (%.3f | %.3f) %d/%d p=%.3f', A.b, A.p, A.pRS, ...
    sum(A.dAnimal>0), numel(A.dAnimal), A.pAnimal);
end

function s = local_hr(F)
if ~F.ok, s = 'n/a'; return; end
if F.rs.ok, s = sprintf('%.2f (p=%.3f)', F.rs.HR, F.rs.pLRT);
else,       s = sprintf('%.2f (p=%.3f)*', F.HR, F.pLRT); end   % * = no random-slope fit
end

function s = local_hrd(F, d, p)
% HR (model p) plus animal direction and sign-flip p
s = sprintf('%s %d/%d p=%.3f', local_hr(F), sum(d>0), numel(d), p);
end

function s = local_cell(F, d, p)
s = sprintf('%s %d/%d p=%.3f', local_hr(F), sum(d>0), numel(d), p);
end

function s = local_arm(ARM, name, hz)
if ~isfield(ARM, name), s = 'n/a'; return; end
A = ARM.(name);
if strcmp(hz,'entry'), s = local_cell(A.entry, A.dEntry, A.pEntry);
else,                  s = local_cell(A.exit,  A.dExit,  A.pExit); end
end

function s = local_didAnim(ARM, hz)
f = ['did_' hz];
if ~isfield(ARM, f) || ~isfinite(ARM.(f).p), s = 'n/a'; return; end
s = sprintf('%.3f', ARM.(f).p);
end

function s = local_didPaired(ARM, hz)
f = ['paired_' hz];
if ~isfield(ARM, f), s = 'n/a'; return; end
s = sprintf('n=%d p=%.3f', numel(ARM.(f).d), ARM.(f).p);
end

function tags = local_boutTags(S)
% fields of RES.acute / RES.persist that hold a bout-definition result
f = fieldnames(S);
tags = f(cellfun(@(x) isstruct(S.(x)) && isscalar(S.(x)) && isfield(S.(x), 'minElapsedSec'), f));
end