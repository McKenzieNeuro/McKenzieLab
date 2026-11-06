%% runAllAnalyses.m
% Runs every analysis behind transplasticity_results.tex, in dependency order, with one set of
% shared settings. Each script runs in its own workspace (no variable leaks between scripts),
% its console output is written to a log, and a failure is recorded and the run continues.
%
% SHARED SETTINGS (CFG below) are passed to each script as the struct SHARED. Every script reads
% SHARED when it exists and otherwise keeps its own defaults, so each still runs standalone.
%   cacheRoot, root, controlRoot   data and cache locations
%   minSpk, minUnits               reactivation unit criteria (Methods: >= 10 units, >= 100 spikes)
%   pairExclude                    'none' = every pair of scored units (the published setting)
%   sessionSet                     'all' = every SAL and every ASV session in the reactivation stages
%                                  (4-8); 'first' = the earlier design (SAL before the first ASV, first ASV)
%
% STAGES (manuscript section <- script -> cache it writes)
%   0  exports        runRippleExport, runRippleExportControl -> rippleEvents*.mat next to the data
%                     (slow; off by default: RUN_EXPORTS = true to rebuild)
%   1  rate           compareRateWindowsByCellType   -> compareRateWindowsByCellType.mat
%   1b infusion order compareInfusionOrder          -> compareInfusionOrder.mat (order effect: repeat SAL vs ASV)
%   1c divergence     rateDivergenceTiming          -> rateDivergenceTiming.mat (when the arms part: infusion vs ConB1)
%   2  engagement     compareRippleWindowsByCellType -> compareRippleWindowsByCellType.mat
%   3  ripple events  compareRippleLevels            -> compareRippleLevels.mat
%   4  reactivation   compareReactivation            -> compareReactivation.mat (EV/REV, per-cell fit)
%   5  coupling       reactivationCoupling           -> reactivationCoupling.mat (needs 4 for sec. 4)
%   5b coupling change couplingPlasticity             -> couplingPlasticity.mat (reads 5 for the link)
%   6  covariance     driftInclusionTest             -> driftInclusionTest_out.mat (controls)
%   7  covariance     driftInclusionASV              -> driftInclusionASV.mat (infused arms)
%   7b persistence    pairCoincidenceGLM             -> pairCoincidenceGLM.mat (needs no other stage)
%   7c re-exposure    crossExposureCovariance        -> crossExposureCovariance.mat (needs no other stage)
%   8  day imbalance  checkConA1Rate                 (reads 7)
%   9  figures        plotResultFigures              (reads 1, 2, 4)
%
% USAGE: set ONLY to a subset of stage names to rerun part, e.g. ONLY = {'reactivationCoupling'}.
% Scripts that save to pwd (driftInclusionTest) save into cacheRoot, because every stage runs
% with cacheRoot as the current folder.

CFG = struct( ...
    'cacheRoot',   'C:\Users\samckenzie\syncCache', ...
    'root',        'R:\TransPlasticity', ...
    'controlRoot', 'R:\HuszarR', ...
    'minSpk',      100, ...
    'minUnits',    10, ...
    'pairExclude', 'none', ...
    'sessionSet',  'all');
RUN_EXPORTS = false;
ONLY        = {};          % {} = all stages
CLOSE_FIGS  = true;        % close each stage's figures once its PDFs are written

STAGES = { ...
    'runRippleExport',                'exports'; ...
    'runRippleExportControl',         'exports'; ...
    'compareRateWindowsByCellType',   'rate'; ...
    'compareInfusionOrder',           'infusion order'; ...
    'rateDivergenceTiming',           'divergence timing'; ...
    'compareRippleWindowsByCellType', 'engagement'; ...
    'compareRippleLevels',            'ripple events'; ...
    'compareReactivation',            'reactivation'; ...
    'reactivationCoupling',           'coupling'; ...
    'couplingPlasticity',             'coupling change'; ...
    'driftInclusionTest',             'covariance (controls)'; ...
    'driftInclusionASV',              'covariance (infused)'; ...
    'pairCoincidenceGLM',             'persistence (pair GLM)'; ...
    'crossExposureCovariance',        're-exposure covariance'; ...
    'checkConA1Rate',                 'day imbalance'; ...
    'plotResultFigures',              'figures'};
if ~RUN_EXPORTS, STAGES = STAGES(~strcmp(STAGES(:,2), 'exports'), :); end
if ~isempty(ONLY), STAGES = STAGES(ismember(STAGES(:,1), ONLY), :); end

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
logDir = fullfile(CFG.cacheRoot, 'logs');
if ~exist(logDir, 'dir'), mkdir(logDir); end
stamp = datestr(now, 'yyyymmdd_HHMMSS'); %#ok<TNOW1,DATST>
here = pwd; cleanup = onCleanup(@() cd(here));
cd(CFG.cacheRoot);

fprintf('runAllAnalyses %s | minSpk %d, minUnits %d, pairs %s, sessions %s | %d stages\n', stamp, ...
    CFG.minSpk, CFG.minUnits, CFG.pairExclude, CFG.sessionSet, size(STAGES, 1));
status = cell(size(STAGES, 1), 4);       % name, ok, minutes, message
for s = 1:size(STAGES, 1)
    nm = STAGES{s,1};
    logFile = fullfile(logDir, sprintf('%s_%s.txt', stamp, nm));
    fprintf('\n[%d/%d] %-32s %s ... ', s, size(STAGES, 1), nm, STAGES{s,2});
    t0 = tic;
    [ok, msg] = local_runStage(nm, CFG, logFile);
    status(s,:) = {nm, ok, toc(t0) / 60, msg};
    if ok, fprintf('done (%.1f min)\n', status{s,3}); else, fprintf('FAILED (%.1f min): %s\n', status{s,3}, msg); end
    if CLOSE_FIGS, close all; end
end

fprintf('\n===== summary (%s; logs in %s) =====\n', stamp, logDir);
for s = 1:size(status, 1)
    flag = 'ok'; if ~status{s,2}, flag = 'FAILED'; end
    fprintf('  %-32s %-7s %6.1f min  %s\n', status{s,1}, flag, status{s,3}, status{s,4});
end
save(fullfile(logDir, sprintf('%s_status.mat', stamp)), 'status', 'CFG');


function [ok, msg] = local_runStage(scriptName, CFG, logFile)
% runs one script inside local_exec, whose workspace holds nothing but the script's own variables:
% scripts that define ok / msg / ME (compareReactivation, reactivationCoupling, plotResultFigures
% do) can no longer overwrite this function's status flags.
ok = true; msg = '';
diary(logFile); diary on
try
    local_exec(scriptName, CFG);
catch ME
    ok = false;
    fr = '';
    for q = 1:min(4, numel(ME.stack))
        fr = [fr sprintf(' <- %s:%d', ME.stack(q).name, ME.stack(q).line)]; %#ok<AGROW>
    end
    msg = sprintf('[%s] "%s"%s', ME.identifier, ME.message, fr);
    rep = getReport(ME, 'extended', 'hyperlinks', 'off');
    if isempty(strtrim(ME.message)), rep = sprintf('%s\n(empty message; class %s, %d causes)', rep, class(ME), numel(ME.cause)); end
    fprintf('\n%s\n', rep);                  % stdout, so the diary log and the command window both keep it
end
diary off
end


function local_exec(RUNALL_script__, SHARED) %#ok<INUSD>  SHARED is read by the script
% the script runs in this workspace; it is discarded on return
run(RUNALL_script__);
end