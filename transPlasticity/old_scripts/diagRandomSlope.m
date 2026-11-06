%% diagRandomSlope.m
% Is the random-slope model actually protecting us? Diagnostics from the saved
% sweep; no refits. For every within-arm and persistence contrast:
%   p plain  intercept-only LRT (random intercepts only)
%   p rs     random-slope LRT -- the number in the sweep tables
%   p rs-t   the rs Wald z referred to t(G-1) instead of the normal, G = animals
%            contributing BOTH conditions: the crude few-clusters correction
%   p anim   exact animal-level sign-flip on raw per-animal log ratios
%            (persistence: NOT adjusted for days in study)
%   tau      estimated SD of the animal-specific slope in the rs fit
%   flag     COLLAPSED  slope variance at the boundary, rs fit == plain fit, so
%                       the "random-slope p" is just the intercept-only p
%            FAILED     rs fit did not converge
%
% Read: p rs systematically << p anim means the rs LRT is still liberal at this
% number of animals. COLLAPSED rows have no protection at all. Uses SWEEP from
% the workspace if present, otherwise the saved timescaleSweep_N<MatchN>.mat.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
MATCHN_SWEEP = 10;
if ~exist('SWEEP','var') || ~exist('SCALES','var')
    load(fullfile(CACHE_ROOT, sprintf('timescaleSweep_N%d.mat', MATCHN_SWEEP)), 'SWEEP', 'SCALES');
end

rows = {};
for k = 1:numel(SWEEP)
    R = SWEEP{k};
    tagsA = local_boutTags(R.acute);
    for t = 1:numel(tagsA)
        Q = R.acute.(tagsA{t});
        for arm = {'ASV','SAL'}
            if ~isfield(Q.arm, arm{1}), continue; end
            A = Q.arm.(arm{1});
            for hz = {'entry','exit'}
                tb = Q.([hz{1} 'Tbl']); tb = tb(tb.Drug == arm{1}, :);
                if strcmp(hz{1},'entry'), pa = A.pEntry; else, pa = A.pExit; end
                rows(end+1,:) = local_row(SCALES(k,1), tagsA{t}, ['acute ' arm{1}], hz{1}, ...
                    A.(hz{1}), pa, tb, 'Period'); %#ok<AGROW>
            end
        end
    end
    if isfield(R, 'persist')
        tagsP = local_boutTags(R.persist);
        for t = 1:numel(tagsP)
            Q = R.persist.(tagsP{t});
            for hz = {'entry','exit'}
                if strcmp(hz{1},'entry'), pa = Q.pEntry; else, pa = Q.pExit; end
                rows(end+1,:) = local_row(SCALES(k,1), tagsP{t}, 'persist', hz{1}, ...
                    Q.(hz{1}), pa, Q.([hz{1} 'Tbl']), 'Phase'); %#ok<AGROW>
            end
        end
    end
end
DIAG = cell2table(rows, 'VariableNames', {'width','bout','context','hazard','G','HR', ...
    'pPlain','pRS','pRSt','pAnim','tau','flag'});

fprintf('\n%-5s %-8s %-11s %-6s %2s | %6s | %-8s %-8s %-8s %-8s | %-5s %-9s\n', ...
    'width','bout','context','hazard','G','HR(rs)','p plain','p rs','p rs-t','p anim','tau','flag');
for i = 1:height(DIAG)
    d = DIAG(i,:);
    fprintf('%-5g %-8s %-11s %-6s %2d | %6.2f | %8.4f %8.4f %8.4f %8.4f | %5.2f %-9s\n', ...
        d.width, d.bout{1}, d.context{1}, d.hazard{1}, d.G, d.HR, d.pPlain, d.pRS, d.pRSt, d.pAnim, d.tau, d.flag{1});
end

use = isfinite(DIAG.pRS) & isfinite(DIAG.pAnim) & isfinite(DIAG.pPlain);
col = strcmp(DIAG.flag, 'COLLAPSED'); bad = strcmp(DIAG.flag, 'FAILED');
fprintf('\n%d contrasts with an rs fit: %d COLLAPSED (rs == plain), %d FAILED\n', sum(use), sum(col & use), sum(bad));
fprintf('median p_anim/p_rs = %.1f; p_anim/p_plain = %.1f; p_anim/p_rs-t = %.1f\n', ...
    median(DIAG.pAnim(use)./DIAG.pRS(use)), median(DIAG.pAnim(use)./DIAG.pPlain(use)), ...
    median(DIAG.pAnim(use & isfinite(DIAG.pRSt))./DIAG.pRSt(use & isfinite(DIAG.pRSt))));
fprintf('p_rs < 0.05 but p_anim >= 0.05: %d cells;  with p_rs-t: %d cells;  p_rs-t < 0.05 but p_anim >= 0.05: %d\n', ...
    sum(use & DIAG.pRS < 0.05 & DIAG.pAnim >= 0.05), ...
    sum(use & DIAG.pRS < 0.05 & DIAG.pAnim >= 0.05 & DIAG.pRSt >= 0.05), ...
    sum(use & DIAG.pRSt < 0.05 & DIAG.pAnim >= 0.05));


%% =========================================================================
function r = local_row(width, bout, ctx, hz, F, pAnim, tb, grp)
G = local_G(tb, grp);
[tau, sdOK] = local_slopeSD(F);
HR = NaN; pRS = NaN; pPlain = NaN; pRSt = NaN; flag = '';
if F.ok, pPlain = F.pLRT; end
if ~F.ok
    flag = 'FAILED';
elseif ~F.rs.ok
    flag = 'FAILED';
else
    HR = F.rs.HR; pRS = F.rs.pLRT;
    se = (log(F.rs.CI(2)) - log(F.rs.CI(1))) / (2*1.96);
    if G >= 2 && isfinite(se) && se > 0
        pRSt = 2 * tcdf(-abs(log(F.rs.HR)/se), G-1);
    end
    same = abs(F.rs.HR - F.HR) < 5e-4 && abs(F.rs.pLRT - F.pLRT) < 5e-3;
    if same || (sdOK && tau < 1e-3), flag = 'COLLAPSED'; end
end
r = {width, bout, ctx, hz, G, HR, pPlain, pRS, pRSt, pAnim, tau, flag};
end

function G = local_G(tb, grp)
% animals with observations in both levels of grp
ua = unique(tb.Animal); G = 0;
for i = 1:numel(ua)
    m = tb.Animal == ua(i);
    if any(tb.(grp)(m) == 'pre') && any(tb.(grp)(m) == 'post'), G = G + 1; end
end
end

function [sd, ok] = local_slopeSD(F)
% SD of the random slope from the stored rs fit (NaN if not stored)
sd = NaN; ok = false;
if ~isfield(F, 'rs') || ~F.rs.ok || ~isfield(F.rs, 'model') || isempty(F.rs.model), return; end
try
    psi = covarianceParameters(F.rs.model);
    for i = 1:numel(psi)
        if size(psi{i}, 1) == 2
            sd = sqrt(max(psi{i}(2,2), 0)); ok = true; return
        end
    end
catch
end
end

function tags = local_boutTags(S)
f = fieldnames(S);
tags = f(cellfun(@(x) isstruct(S.(x)) && isscalar(S.(x)) && isfield(S.(x), 'minElapsedSec'), f));
end
