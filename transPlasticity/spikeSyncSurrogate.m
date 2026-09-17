function [Cnull, Ct, tc, pval] = spikeSyncSurrogate(spikes, varargin)
%SPIKESYNCSURROGATE  Rate-matched null band for a time-resolved C trace.
%
%   [Cnull, Ct, tc, pval] = spikeSyncSurrogate(spikes, ...)
%
%   Why this is not optional: C is not zero under independence. Its chance
%   level depends on the local firing rates, because the adaptive coincidence
%   window tau shrinks with rate (Eq. 15) but the number of candidate partners
%   grows. A raw C(t) trace therefore moves when rate moves, with no change in
%   timing precision anywhere in the population. Any claim that synchrony rose
%   during some epoch has to be a claim about C relative to a null that carries
%   the same rate profile.
%
%   INPUT
%     spikes : cell array as for spikeSync.
%
%   OPTIONS (in addition to everything spikeSync and spikeSyncWindow accept)
%     'NSurr'  : number of surrogates (default 200).
%     'Method' : 'dither' (default) or 'shift'.
%                'dither' displaces every spike independently by U(-D,D).
%                   Destroys coincidence below D, preserves rate above D.
%                   Null = "the observed rate profile with timing scrambled
%                   below D". This is the null you usually want.
%                'shift' circularly rotates each whole train by a random lag.
%                   Preserves each train's own ISI structure and bursts exactly,
%                   destroys only cross-train alignment. Weaker, but the right
%                   null when within-train structure is itself the confound.
%     'Dither' : D for 'dither'. Required for that method. Sets the timescale
%                the test is sensitive to: coincidences tighter than D count as
%                signal, looser than D count as background. Choose it from the
%                conduction/integration window you care about, and report it.
%     'Alpha'  : two-sided band, default 0.05 -> 2.5/97.5 percentiles.
%
%   OUTPUT
%     Cnull : nWin x 3 [lower, median, upper] percentile band across surrogates
%     Ct,tc : observed trace and window centres (from spikeSyncWindow)
%     pval  : nWin x 1 pointwise one-sided p, (1 + #{surr >= obs}) / (1 + NSurr)
%
%   POINTWISE, NOT CORRECTED. With hundreds of overlapping windows these p
%   values are neither independent nor multiplicity-corrected. Use the band for
%   display and, for inference, either a cluster-level statistic (max run of
%   supra-band windows, compared against the same statistic computed on each
%   surrogate) or a single a-priori window.
%
%   See also SPIKESYNC, SPIKESYNCWINDOW.

p = inputParser; p.KeepUnmatched = true;
p.addParameter('NSurr', 200, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('Method', 'dither', @(s) any(strcmpi(s, {'dither','shift'})));
p.addParameter('Dither', [], @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('Interval', [], @(x) isempty(x) || numel(x) == 2);
p.parse(varargin{:});
opt = p.Results;
winArgs = p.Unmatched;

if strcmpi(opt.Method,'dither') && isempty(opt.Dither)
    error('spikeSyncSurrogate:dither','''Dither'' is required for method ''dither''.');
end

syncArgs = {};
if ~isempty(opt.Interval), syncArgs = {'Interval', opt.Interval}; end
winCell = [fieldnames(winArgs), struct2cell(winArgs)]';
winCell = winCell(:)';

% observed
[~, ~, prof] = spikeSync(spikes, syncArgs{:});
[Ct, tc, ~]  = spikeSyncWindow(prof, winCell{:});

% interval for circular shifts
allT = vertcat(spikes{:});
if isempty(opt.Interval), iv = [min(allT) max(allT)]; else, iv = sort(opt.Interval(:))'; end
L = iv(2) - iv(1);

nW = numel(Ct);
Csur = nan(nW, opt.NSurr);
N = numel(spikes);

for s = 1:opt.NSurr
    sur = cell(N,1);
    for n = 1:N
        x = spikes{n}(:);
        switch lower(opt.Method)
            case 'dither'
                x = x + opt.Dither * (2*rand(numel(x),1) - 1);
            case 'shift'
                x = iv(1) + mod(x - iv(1) + rand*L, L);
        end
        sur{n} = sort(x);
    end
    [~, ~, sprof] = spikeSync(sur, syncArgs{:});
    Csur(:,s) = spikeSyncWindow(sprof, winCell{:});
end

q = [100*opt.Alpha/2, 50, 100*(1 - opt.Alpha/2)];
Cnull = zeros(nW, 3);
for w = 1:nW
    Cnull(w,:) = local_pct(Csur(w,:), q);
end
pval  = (1 + sum(Csur >= Ct, 2, 'omitnan')) ./ (1 + opt.NSurr);
pval(isnan(Ct)) = NaN;

end


function y = local_pct(x, q)
% percentiles without the Statistics Toolbox
x = sort(x(~isnan(x(:))));
n = numel(x);
if n == 0, y = nan(size(q)); return; end
if n == 1, y = repmat(x, size(q)); return; end
pos = (0.5:n-0.5) / n * 100;
y = interp1(pos, x, q, 'linear');
y(q < pos(1)) = x(1); y(q > pos(end)) = x(end);

end
