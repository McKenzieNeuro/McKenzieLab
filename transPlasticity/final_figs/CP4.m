%% CP4: PYR vs INT firing, WAKE vs NREM, by session
ROOT = 'R:\TransPlasticity'; ANIMAL = 'CP4';
fa = getAllExtFiles(ROOT, 'mat', 1); fa = fa(endsWith(fa, 'spikes.cellinfo.mat'));
fa = fa(strcmp(cellfun(@(f) animalFromPath(f,'TransPlasticity'), fa, 'UniformOutput', false), ANIMAL));
cnd = cellfun(@(f) conditionFromPath(f,'TransPlasticity'), fa, 'UniformOutput', false);
dnm = cellfun(@sessionDateFromPath, fa);  [dnm, o] = sort(dnm); fa = fa(o); cnd = cnd(o);

R = nan(numel(fa), 4);                       % median rate: PYR-WAKE PYR-NREM INT-WAKE INT-NREM
for i = 1:numel(fa)
    d = fileparts(fa{i});
    sp = load(fa{i}); sp = sp.spikes;  ss = dir(fullfile(d,'*.SleepState.states.mat'));
    cm = dir(fullfile(d,'*.cell_metrics.cellinfo.mat'));
    if isempty(ss) || isempty(cm), continue; end
    ss = load(fullfile(d,ss(1).name)); ss = ss.SleepState.ints;
    cm = load(fullfile(d,cm(1).name)); isP = strcmp(cm.cell_metrics.putativeCellType,'Pyramidal Cell');
    st = {ss.WAKEstate, ss.NREMstate};  T = cellfun(@(x) sum(diff(x,1,2)), st);
    r = zeros(numel(sp.times), 2);
    for u = 1:numel(sp.times)
        for s = 1:2
            r(u,s) = sum(InIntervals(sp.times{u}, st{s})) / T(s);   % Hz
        end
    end
    R(i,:) = [median(r(isP,:),1,'omitnan'), median(r(~isP,:),1,'omitnan')];
end

figure; bar(R); set(gca,'XTick',1:numel(fa),'XTickLabel',cellstr(strcat(cnd(:),'_',string(dnm(:)))),'XTickLabelRotation',45);
legend({'PYR wake','PYR NREM','INT wake','INT NREM'}); ylabel('median rate (Hz)'); title(ANIMAL); set(gca,'YScale','log');