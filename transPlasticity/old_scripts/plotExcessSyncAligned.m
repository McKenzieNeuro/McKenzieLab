%% plotExcessSyncAligned.m
% Excess synchrony e as a function of time around the first ASV infusion
% (ts2) and around the first SAL infusion (ts3), from the per-session
% spikeDistanceProfile.mat files (run addTs3.m first if they predate ts3).
%
% e, per WIN-second block: average the observed SPIKE-distance S and every
% jittered surrogate over the block, then
%     e = (median(null) - S) / (median(null) - p2.5(null))
% S is a distance, so e > 0 = more synchronous than the interval-jitter null;
% e = 1 at the 2.5th-percentile edge. Also stored: block p = (1 + #null <= S)/(nSurr + 1).
% Averaging BEFORE ranking matters: the null narrows with averaging, so a
% smoothed per-bin e understates how far a sustained epoch is outside it.
%
% Figure 1 (hours, infusion day): per-animal traces and the across-animal
%   mean +- SEM (animals averaged first, ZBIN bins), around ASV, around SAL,
%   and the two means overlaid.
% Figure 2 (days): mean e per session against days from the first ASV / SAL.
% Times in date-only mode (meta.ts2Mode / ts3Mode = 'date') are exact on the
% infusion day and +-24 h on other days; the counts are printed.
%
% SAL alignment uses only blocks BEFORE the animal's first ASV infusion
% (ts2 < 0, or the animal has no ASV session): once ASV has been given, time
% from the first SAL infusion no longer describes a drug-naive animal.

CACHE_ROOT = 'C:\Users\samckenzie\syncCache';
ROOT = 'R:\TransPlasticity';
WIN  = 300;            % s per block for e
ZOOM = [-3 8];        % h, figure 1
ZBIN = 10;            % min, bins for the across-animal mean
REFRESH = true;      % t[hrue: recompute e from the session files

cf = fullfile(CACHE_ROOT, sprintf('excessSyncAligned_W%d.mat', WIN));
if ~REFRESH && exist(cf, 'file')
    load(cf, 'E');
else
    ff = getAllExtFiles(ROOT, 'mat', 1);
    ff = ff(endsWith(ff, 'spikeDistanceProfile.mat'));
    E = struct('file',{},'animal',{},'condition',{},'ts2',{},'ts3',{},'e',{},'p',{}, ...
        'mode2',{},'mode3',{},'day2',{},'day3',{},'N',{});
    for i = 1:numel(ff)
        L = load(ff{i}, 'SDIST'); X = L.SDIST;
        if ~isfield(X.bin, 'ts3'), error('%s has no ts3: run addTs3.m first', ff{i}); end
        [e, p, t2, t3] = local_blockE(X, WIN);
        E(end+1) = struct('file',ff{i}, 'animal',X.meta.animal, 'condition',X.meta.condition, ...
            'ts2',t2, 'ts3',t3, 'e',e, 'p',p, 'mode2',X.meta.ts2Mode, 'mode3',X.meta.ts3Mode, ...
            'day2',X.meta.dayFromFirstASV, 'day3',X.meta.dayFromFirstSAL,'N',length(X.meta.keepIdx)); %#ok<SAGROW>
        fprintf('[%3d/%3d] %s  %s  median e %.2f\n', i, numel(ff), X.meta.animal, X.meta.condition, median(e, 'omitnan'));
    end
    save(cf, 'E', 'WIN');
end
fprintf('%d sessions, %d animals\n', numel(E), numel(unique({E.animal})));
fprintf('ts2 modes: %s\n', local_tab({E.mode2}));
fprintf('ts3 modes: %s\n', local_tab({E.mode3}));
nPost = sum(arrayfun(@(q) any(isfinite(q.ts3)) && all(q.ts2 >= 0), E));
fprintf('SAL alignment: %d sessions lie entirely after the first ASV and are excluded from it\n', nPost);

%% ---- figure 1: hours around the infusion ----------------------------------
ax = {'ts2','ASV'; 'ts3','SAL'};
col = {[0.85 0.2 0.1], [0.1 0.35 0.8]};
edges = (ZOOM(1)*60 : ZBIN : ZOOM(2)*60) / 60;            % h
ctr = edges(1:end-1) + diff(edges)/2;
f1 = figure('Color','w','Position',[60 60 1500 420]);
tl = tiledlayout(1, 3, 'TileSpacing','compact');
M = cell(1, 2);
NN = cell(1,2);
for a = 1:2
    nexttile; hold on
    ua = unique({E.animal});
    Y = nan(numel(ua), numel(ctr));
    NCel = nan(numel(ua), 1);
    for k = 1:numel(ua)
        t = []; y = [];
        for s = find(strcmp({E.animal}, ua{k}))
            th = E(s).(ax{a,1}) / 3600;
            in = th >= ZOOM(1) & th <= ZOOM(2) & isfinite(E(s).e) & local_valid(E(s), a);
            t = [t; th(in)]; y = [y; E(s).e(in)]; %#ok<AGROW>
            
            if any(in)
               NCel(k) = E(s).N;
            end
        end
        if isempty(t), continue; end
        [t, o] = sort(t); y = y(o);
        plot(t, y, '-', 'Color', [col{a} 0.25], 'LineWidth', 0.5);
        b = discretize(t, edges);
        ok = ~isnan(b);
        Y(k,:) = accumarray(b(ok), y(ok), [numel(ctr) 1], @mean, NaN)';
     
    end
  
    an{a} = ua(any(isfinite(Y), 2));
    NN{a} = NCel(any(isfinite(Y), 2));
      Y = Y(any(isfinite(Y), 2), :);
    M{a} = Y;
   
    local_band(ctr, Y, col{a});
    xline(0, 'k:'); yline(0, 'k-'); yline(1, 'k--');
    xlabel(sprintf('h from first %s infusion', ax{a,2})); ylabel(sprintf('e (%d s blocks)', WIN));
    title(sprintf('around first %s (%d animals)', ax{a,2}, size(Y,1)));
    xlim(ZOOM);
end
nexttile; hold on
for a = 1:2, local_band(ctr, M{a}, col{a}); end
xline(0, 'k:'); yline(0, 'k-');
xlabel('h from first infusion'); ylabel('e'); xlim(ZOOM);
legend({'', sprintf('ASV (n=%d)', size(M{1},1)), '', sprintf('SAL (n=%d)', size(M{2},1))}, 'Location','best', 'Box','off');
title('ASV vs SAL, aligned to each');
title(tl, sprintf('Excess synchrony around infusion (mean \\pm SEM across animals, %d-min bins)', ZBIN));

%% ---- figure 2: days around the infusion -----------------------------------
f2 = figure('Color','w','Position',[60 520 1100 380]);
tiledlayout(1, 2, 'TileSpacing','compact');
for a = 1:2
    nexttile; hold on
    ua = unique({E.animal});
    for k = 1:numel(ua)
        s = find(strcmp({E.animal}, ua{k}));
        d = arrayfun(@(q) median(E(q).(ax{a,1})(local_valid(E(q), a)), 'omitnan') / 86400, s);
        v = arrayfun(@(q) mean(E(q).e(local_valid(E(q), a)), 'omitnan'), s);
        ok = isfinite(d) & isfinite(v);
        [d, o] = sort(d(ok)); v = v(ok); v = v(o);
        plot(d, v, '-o', 'Color', [col{a} 0.5], 'MarkerFaceColor', col{a}, 'MarkerSize', 3);
        d1{k,a} = d;
        v1{k,a} = v;
    end
    xline(0, 'k:'); yline(0, 'k-');
    xlabel(sprintf('days from first %s infusion', ax{a,2})); ylabel('session mean e');
    title(sprintf('per session, aligned to first %s', ax{a,2}));
end

%%
preASVd= cell2mat(d1(:,1)');

preSAL= cell2mat(v1(:,2)');


preSALd= cell2mat(d1(:,2)');

preASV= cell2mat(v1(:,1)');
dneg =  preASV(preASVd<0 );
dpos = preASV(preASVd>1);

preInf = 1;
kpA = NN{1}>10;
kpS = NN{2}>10;
dof = M{1}(kpA,:);
dof = [nanmean(dof(:,edges(1:end-1)>-preInf&edges(2:end)<0),2) nanmean(dof(:,edges(1:end-1)>0 &edges(2:end)<2),2) nanmean(dof(:,edges(1:end-1)>2 &edges(2:end)<inf),2)];


dnegS =  preSAL(preSALd<0);
dposS = preSAL(preSALd>1);

dofS = M{2}(kpS,:);
dofS = [nanmean(dofS(:,edges(1:end-1)>-preInf&edges(2:end)<0),2) nanmean(dofS(:,edges(1:end-1)>0 &edges(2:end)<2),2) nanmean(dofS(:,edges(1:end-1)>2 &edges(2:end)<inf),2)];
figure
plot([nanmean(dneg) nan nanmean(dof) nan nanmean(dpos)],'-o','color','r')
hold on
errorbar([nanmean(dneg) nan nanmean(dof) nan nanmean(dpos)],[SEM(dneg(:)) nan SEM(dof) nan SEM(dpos(:))],'color','r','linestyle','none')
plot([nanmean(dnegS) nan nanmean(dofS) nan nanmean(dposS)],'-x','color','b')
errorbar([nanmean(dnegS) nan nanmean(dofS) nan nanmean(dposS)],[SEM(dnegS(:)) nan SEM(dofS) nan SEM(dposS(:))],'color','b','linestyle','none')
set(gca,'xtick',1:7,'xticklabel',{'baseline','','pre','2hr','8hr','','>24hr'})
xlim([0 8])
%%
exportgraphics(f1, fullfile(CACHE_ROOT, sprintf('excessSyncAligned_hours_W%d.pdf', WIN)), 'ContentType','vector');
exportgraphics(f2, fullfile(CACHE_ROOT, sprintf('excessSyncAligned_days_W%d.pdf', WIN)), 'ContentType','vector');


%% =========================================================================
function ok = local_valid(Es, a)
% blocks usable on alignment a: all for ASV (a = 1); for SAL (a = 2) only
% those before the first ASV infusion (ts2 < 0) or in animals never given ASV
if a == 1, ok = true(size(Es.e)); return; end
ok = ~(Es.ts2 >= 0);
end


function [e, p, t2, t3] = local_blockE(X, W)
S = X.bin.S(:); Sn = double(X.bin.Snull);
bw = median(diff(X.bin.edges));
nb = max(1, round(W / bw));
nW = floor(numel(S) / nb);
blk = @(v) mean(reshape(v(1:nW*nb), nb, nW), 1)';
Sw = blk(S);
Snw = reshape(mean(reshape(Sn(1:nW*nb, :), nb, nW, []), 1), nW, []);
med = median(Snw, 2); lo = prctile(Snw, 2.5, 2);
e = (med - Sw) ./ (med - lo);
p = (1 + sum(Snw <= Sw, 2)) ./ (size(Snw, 2) + 1);
t2 = blk(X.bin.ts2(:)); t3 = blk(X.bin.ts3(:));
end

function local_band(x, Y, c)
mu = mean(Y, 1, 'omitnan');
n = sum(isfinite(Y), 1);
se = std(Y, 0, 1, 'omitnan') ./ sqrt(max(n, 1));
ok = n >= 2;
if ~any(ok), plot(x, mu, '-', 'Color', c, 'LineWidth', 2); return; end
fill([x(ok) fliplr(x(ok))], [mu(ok) - se(ok), fliplr(mu(ok) + se(ok))], c, 'FaceAlpha',0.2, 'EdgeColor','none');
plot(x, mu, '-', 'Color', c, 'LineWidth', 2);
end

function s = local_tab(c)
[u, ~, j] = unique(c);
s = strjoin(arrayfun(@(k) sprintf('%s %d', u{k}, sum(j == k)), 1:numel(u), 'UniformOutput', false), ', ');
end