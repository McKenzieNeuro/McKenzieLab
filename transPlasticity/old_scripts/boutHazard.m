function H = boutHazard(r, step, minRun)
%BOUTHAZARD  Synchrony bouts and inter-bout intervals as censored survival data.
%
%   H = boutHazard(r, step, minRun)
%
%   r      : one session struct from syncBatch, run with DetectTransients=true
%            (needs .Ct, .band, and for minRun='null95' also .nullRunLen)
%   step   : window step (s) used by syncBatch
%   minRun : minimum run of supra-band windows that counts as a bout.
%            numeric -> fixed number of windows, identical in every session
%            'null95' (default) / 'null99' -> this session's 95th / 99th
%               percentile of ALL surrogate supra-band run lengths, i.e. a
%               fixed per-run false-positive rate
%
%   WHY THIS INSTEAD OF .transients(.significant). The .significant flag uses
%   the percentile of the surrogates' MAXIMUM run, which grows with epoch
%   length, so bout detection got stricter in longer epochs and every epoch
%   had to be clipped to one fixed length. Both rules above are independent
%   of epoch length, so the whole drug-free period of every session can be
%   used.
%
%   STATE MODEL. Each valid window is either IN a bout or OUT. Windows with
%   NaN Ct (excluded, or too few spikes) break the record into separate
%   observed segments; nothing is assumed about what happened inside a gap.
%
%   EXIT observations, one per bout:
%     exposure = time at risk AFTER the bout qualified, (nWin - minRun + 1)*step
%     event    = 1 if the bout ended inside the segment, 0 if the segment ended
%                first (right-censored)
%   ENTRY observations, one per out-interval (including a whole segment with
%   no bout at all -- this is what the per-epoch duration metric dropped):
%     exposure = nWin*step
%     event    = 1 if a qualifying bout started inside the segment, else 0
%
%   These feed a piecewise-exponential Poisson GLMM:
%     event ~ covariates + (1|Animal) + (1|Session), offset log(exposure)
%   whose coefficients are log hazard ratios. A single constant hazard per
%   observation is an approximation (memoryless states); it is what makes
%   left-truncated intervals at segment starts usable without knowing when
%   they began.
%
%   Also returns a threshold-free summary of the trace: excess synchrony
%   e = (Ct - p50) / (p97.5 - p50) per window, so e = 1 is the band edge.
%
%   OUTPUT struct:
%     .exit   [nBout x 2]  [exposure_s event]
%     .entry  [nInt  x 2]  [exposure_s event]
%     .boutWin [nBout x 2]  first/last window index (into r.Ct) of each bout
%     .boutE   [nBout x 2]  mean and peak excess synchrony e inside each bout
%                           (within-bout strength; rows align with .exit/.boutWin)
%     .minRunUsed, .nValidWin, .validSec
%     .eMeanPos   mean(max(e,0)) over valid windows
%     .eQ90       90th percentile of e
%     .fracAbove  fraction of valid windows with e > 1
%
%   See also SYNCBATCH, RUNHAZARDANALYSIS.

if nargin < 3 || isempty(minRun), minRun = 'null95'; end

if ~isfield(r,'band') || isempty(r.band)
    error('boutHazard:noBand', ...
        ['%s has no stored surrogate band. It was loaded from a cache written ' ...
         'before syncBatch stored .band -- rerun syncBatch with a fresh CacheDir.'], r.file);
end

Ct = r.Ct(:);
band = r.band;
if size(band,1) ~= numel(Ct)
    error('boutHazard:size','%s: band (%d) and Ct (%d) lengths differ.', ...
        r.file, size(band,1), numel(Ct));
end

if ischar(minRun)
    if ~isfield(r,'nullRunLen') || isempty(r.nullRunLen) || all(isnan(r.nullRunLen))
        error('boutHazard:noNull','%s: no nullRunLen stored; use a numeric minRun.', r.file);
    end
    switch lower(minRun)
        case 'null95', k = r.nullRunLen(3);
        case 'null99', k = r.nullRunLen(4);
        otherwise, error('boutHazard:minRun','unknown minRun ''%s''', minRun);
    end
    minRun = max(1, ceil(k));
end

valid = isfinite(Ct) & all(isfinite(band),2);
supra = valid & Ct > band(:,3);

H.minRunUsed = minRun;
H.nValidWin  = sum(valid);
H.validSec   = sum(valid) * step;
H.exit  = zeros(0,2);
H.entry = zeros(0,2);
H.boutWin = zeros(0,2);
H.boutE   = zeros(0,2);    % [mean e, peak e] over each bout's windows

% ---- threshold-free summary ----------------------------------------------
den = band(:,3) - band(:,2);
okE = valid & den > 0;
e = (Ct(okE) - band(okE,2)) ./ den(okE);
eFull = nan(numel(Ct),1); eFull(find(okE)) = e; %#ok<FNDSB>   % e per window, NaN where undefined
if isempty(e)
    H.eMeanPos = NaN; H.eQ90 = NaN; H.fracAbove = NaN;
else
    H.eMeanPos  = mean(max(e,0));
    H.eQ90      = local_pct(e, 90);
    H.fracAbove = mean(e > 1);
end

% ---- segment-wise state sequence -----------------------------------------
[sa, sb] = local_runs(valid);
for g = 1:numel(sa)
    seg = supra(sa(g):sb(g));
    Ls  = numel(seg);

    [ra, rb] = local_runs(seg);
    q = (rb - ra + 1) >= minRun;
    ra = ra(q); rb = rb(q);

    inBout = false(Ls,1);
    for k = 1:numel(ra)
        inBout(ra(k):rb(k)) = true;
        atRisk = (rb(k) - ra(k) + 1) - (minRun - 1);
        H.exit(end+1,:) = [atRisk*step, double(rb(k) < Ls)]; %#ok<AGROW>
        H.boutWin(end+1,:) = sa(g) - 1 + [ra(k) rb(k)]; %#ok<AGROW>
        ee = eFull(sa(g) - 1 + (ra(k):rb(k)));
        H.boutE(end+1,:) = [mean(ee,'omitnan') max(ee,[],'omitnan')]; %#ok<AGROW>
    end

    [oa, ob] = local_runs(~inBout);
    for k = 1:numel(oa)
        H.entry(end+1,:) = [(ob(k) - oa(k) + 1)*step, double(ob(k) < Ls)]; %#ok<AGROW>
    end
end
end


function [starts, ends] = local_runs(tf)
tf = [false; tf(:); false];
d = diff(tf);
starts = find(d == 1);
ends = find(d == -1) - 1;
end


function y = local_pct(x, q)
x = sort(x(~isnan(x(:))));
n = numel(x);
if n == 0, y = nan(size(q)); return; end
if n == 1, y = repmat(x, size(q)); return; end
pos = (0.5:n-0.5) / n * 100;
y = interp1(pos, x, q, 'linear');
y(q < pos(1)) = x(1); y(q > pos(end)) = x(end);
end