%% =========================================================================
%% CHECK 1: duration-matched truncation sensitivity
%% =========================================================================
% Caps BOTH datasets to the same fixed recording length and re-asks the
% duty-cycle/duration question. If the effect survives, it isn't just "more
% total minutes, more chances to catch a long state-dependent event." If it
% shrinks or vanishes, the earlier result was likely riding on TransPlasticity
% sessions running long enough to reach a state (e.g. sleep/rest) that OptoRAM
% sessions don't.
%
% Cap = OptoRAM median duration: the largest cap that doesn't require
% discarding OptoRAM sessions, since OptoRAM is the shorter dataset.

stepUsed  = cfg{find(strcmp(cfg,'Step'))+1};
widthUsed = cfg{find(strcmp(cfg,'Width'))+1};

capDur = median([b.dur]);
fprintf('truncating both datasets to %.0f s\n', capDur);

tcfg = [cfg, {'MaxDuration', capDur}];
Tt1 = syncBatch(fils1, tcfg{:});
Tt2 = syncBatch(fils2, tcfg{:}, 'Exclude', @optoPulseIntervals);

at = Tt1([Tt1.ok]); bt = Tt2([Tt2.ok]);
%%
fprintf('truncated N sessions: TransPlasticity %d, OptoRAM %d\n', numel(at), numel(bt));
fprintf('truncation actually applied: %d/%d TransPlasticity sessions\n', ...
    sum([at.durTruncated]), numel(at));

[rate1t, occ1t, dur1t] = local_collect(at, stepUsed, widthUsed);
[rate2t, occ2t, dur2t] = local_collect(bt, stepUsed, widthUsed);

fprintf('\n-- truncated to %.0f s --\n', capDur);
fprintf('rate      : %.2f vs %.2f   p = %.4f\n', mean(rate1t), mean(rate2t), permtest(rate1t, rate2t));
fprintf('duty (%%)  : %.2f vs %.2f   p = %.4f\n', mean(occ1t)*100, mean(occ2t)*100, permtest(occ1t, occ2t));
fprintf('duration  : %.1f vs %.1f   p = %.4f\n', mean(dur1t,'omitnan'), mean(dur2t,'omitnan'), ...
    permtest(dur1t(~isnan(dur1t)), dur2t(~isnan(dur2t))));
%%
fprintf('\ncompare to untruncated: rate p=0.064, duty p=0.0006, duration p=0.0108\n');
% Similar p-values -> not a duration artifact. Moving substantially toward
% 1 -> it was.


%% =========================================================================
%% CHECK 2: does the result depend on MatchN?
%% =========================================================================
% More units per session should mean a tighter surrogate band (easier to
% detect a real effect), but MatchN also raises MinUnits' effective bar,
% shrinking the sample and adding subset-draw variance session to session.
% Sweep and look at both costs together. Uses FULL duration -- this is about
% unit count, not time; don't combine with Check 1 in one run.

matchVals = [10 15 20 25];
sweep = struct('MatchN', {}, 'nA', {}, 'nB', {}, 'rateP', {}, 'dutyP', {}, ...
    'rateA', {}, 'rateB', {});

for mi = 1:numel(matchVals)
    mN = matchVals(mi);
    mcfg = cfg;
    mcfg{find(strcmp(mcfg,'MatchN'))+1} = mN;

    Tm1 = syncBatch(fils1, mcfg{:});
    Tm2 = syncBatch(fils2, mcfg{:}, 'Exclude', @optoPulseIntervals);
    am = Tm1([Tm1.ok]); bm = Tm2([Tm2.ok]);

    [r1, o1] = local_collect(am, stepUsed, widthUsed);
    [r2, o2] = local_collect(bm, stepUsed, widthUsed);

    sweep(mi).MatchN = mN;
    sweep(mi).nA = numel(am); sweep(mi).nB = numel(bm);
    sweep(mi).rateA = mean(r1); sweep(mi).rateB = mean(r2);
    sweep(mi).rateP = permtest(r1, r2);
    sweep(mi).dutyP = permtest(o1, o2);

    fprintf('MatchN=%2d  nA=%3d nB=%3d  rate %.2f vs %.2f (p=%.4f)  duty p=%.4f\n', ...
        mN, sweep(mi).nA, sweep(mi).nB, sweep(mi).rateA, sweep(mi).rateB, ...
        sweep(mi).rateP, sweep(mi).dutyP);
end

% What to look for:
%  - nA/nB dropping sharply as MatchN rises: where's the ceiling before
%    you're discarding most of the dataset.
%  - rateP/dutyP stable in size and direction across MatchN: the result
%    isn't an artifact of the specific N=10 choice.
%  - p improving monotonically as MatchN rises (before sessions start
%    dropping): a real effect that N=10 was underpowered to detect cleanly
%    -- report the largest MatchN that doesn't cost too many sessions.
%  - effect flips or vanishes at higher N: distrust the N=10 result --
%    something about which units survive a stricter MatchN is doing work.


%% =========================================================================
function [rate, duty, meanDur] = local_collect(T, step, width)
rate = nan(numel(T),1); duty = nan(numel(T),1); meanDur = nan(numel(T),1);
for i = 1:numel(T)
    [rate(i), duty(i), meanDur(i)] = transientStats(T(i), step, width);
end
end
