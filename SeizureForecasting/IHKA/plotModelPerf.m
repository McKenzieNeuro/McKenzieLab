
% load all predictions for each time point and make long list of TP and FP

topDir = 'R:\IHKA_Scharfman\prediction';
fils = getAllExtFiles(topDir,'mat',1);
kp = contains(fils,'predict') & ~contains(fils,'before');
fils  = fils(kp);

estimateLabel =[];trueLabel=[];time2seizure =[];confLabel =[];ID=[];
for i = 1:length(fils)
    
    v= load(fils{i});
    
    minLen = min(length(v.trueLabel),length(v.estimateLabel));
    trueLabel = [trueLabel;v.trueLabel(1:minLen)'];
    estimateLabel = [estimateLabel;v.estimateLabel(1:minLen,1)];
    confLabel = [confLabel;v.estimateLabel(1:minLen,2:end)];
    time2seizure = [time2seizure;v.time2seizure(1:minLen)'];
    ID = [ID;i*ones(minLen,1)];
    i
end
for i =1:4
    ind(:,i) = estimateLabel==i;
end



%%
close all
figure
bins = -fliplr(logspace(0,log10(5e4),100));
semilogx(bins,avghist(time2seizure,double(estimateLabel(:,1)==1),bins))

hold on
semilogx(bins,avghist(time2seizure,double(estimateLabel(:,1)==2),bins))

semilogx(bins,avghist(time2seizure,double(estimateLabel(:,1)==3),bins))

semilogx(bins,avghist(time2seizure,double(estimateLabel(:,1)==4),bins))
semilogx(bins,avghist(time2seizure,double(estimateLabel(:,1)==5),bins))
