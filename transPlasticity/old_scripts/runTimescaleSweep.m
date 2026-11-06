%% runTimescaleSweep.m
% Does the ASV effect on bout dynamics depend on the 60 s analysis window --
% acutely (pre vs 4-6 h post infusion, ASV against SAL) and persistently
% (drug-naive vs >=24 h post-exposure, drug-free recordings)?
%
% The bouts in the main analysis last >= ~2 min partly BY CONSTRUCTION: C is
% averaged over 60 s windows and a bout needs a run of them above the
% surrogate band. This reruns the hazard analysis at narrower windows with
% everything else identical -- same units, same shared-unit selection, same
% hazard models, same bout-definition logic.
%
% Reading the result:
%   effect present and similar at every width  -> not a product of the 60 s
%       choice; the change is visible down to the tens-of-seconds scale
%   effect only at 60 s (vanishes at 10-20 s)  -> a slow, state-like effect
%   effect appears or grows at short widths    -> faster structure the 60 s
%       window was averaging over
%
% THE WIDTHS ARE NOT INDEPENDENT TESTS. Same animals, same spikes, four
% smoothings. Agreement across widths is robustness, not replication, and
% the animal-level p at any one width is not strengthened by the others.
%
% HAZARD RATIOS ARE COMPARABLE ACROSS WIDTHS, absolute hazards are not: a
% "bout" at 10 s and at 60 s are different objects (min(s) column).
%
% Bout definitions per scale:
%   null95 / null99  percentile of the surrogate supra-band run lengths at
%                    THAT scale (scale-adaptive by construction)
%   fixed            MINRUN_SEC in ELAPSED seconds, converted to a window
%                    count at each scale (label = window count)
%
% p-VALUES IN THE TABLES
%   model    random-slope GLMM LRT (* = random slope failed, intercept-only
%            LRT shown) -- historically optimistic with ~10 animals
%   anim     animal-level: exact sign-flip on per-animal log rate ratios
%            (within arm / persistence); two-sample permutation for the DiD
%   paired   DiD within animals that had both an ASV and a SAL session
%   "k/n up" counts animals whose hazard ROSE, for entry and exit alike.
%   Quote the animal-level p. The model p is shown so its optimism is visible.
%
% PERSISTENCE carries the days-in-study confound (post epochs are always
% later); Days is in the model, which adjusts for a linear trend only.
%
% COST. Acute: same options as the previous sweep, so syncBatch returns from
% its checkpoints and only the fits rerun. Persistence: one syncBatch pass per
% scale over the WHOLE drug-free periods with 50 surrogates -- the expensive
% part. Windows per second of data scale as 1/step, so 10 s / 2 s costs ~2.5x
% the 60 s / 5 s pass. Checkpointed per scale (ckpt_hzPersist_W..), so an
% interrupted sweep resumes.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
SCALES     = [10 2; 20 4; 30 5; 60 5];     % [width step] in seconds
MINRUN_SEC = [60 175];                     % fixed bout minimum, elapsed seconds
MATCHN_SWEEP = 10;                         % primary subset size (matches the reported runs)
RUN_PERSIST_SWEEP = true;

% The loop variable must not be one runHazardAnalysis uses (k, i, j, ...): it
% runs in this workspace and would overwrite it mid-loop.
SWEEP = cell(size(SCALES,1), 1);
for sweepIdx = 1:size(SCALES,1)
    clear RES
    SW = struct('width',SCALES(sweepIdx,1), 'step',SCALES(sweepIdx,2), 'matchN',MATCHN_SWEEP, ...
        'persist',RUN_PERSIST_SWEEP, 'minRunSec',MINRUN_SEC, ...
        'tag',sprintf('_W%dS%d_N%d', SCALES(sweepIdx,1), SCALES(sweepIdx,2), MATCHN_SWEEP));
    fprintf('\n################ scale %d/%d: width %g s, step %g s ################\n', ...
        sweepIdx, size(SCALES,1), SCALES(sweepIdx,1), SCALES(sweepIdx,2));
    runHazardAnalysis;
    SWEEP{sweepIdx} = RES;
end
clear SW

% save BEFORE printing: a print error must not lose a finished sweep
save(fullfile(CACHE_ROOT, sprintf('timescaleSweep_N%d.mat', MATCHN_SWEEP)), ...
    'SWEEP', 'SCALES', 'MINRUN_SEC', 'MATCHN_SWEEP', 'RUN_PERSIST_SWEEP');
printTimescaleSweep