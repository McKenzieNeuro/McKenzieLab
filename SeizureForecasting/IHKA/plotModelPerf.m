
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

%calculate AUC
clear AUC
for i = 1:max(ID)
    kp = ID ==i &trueLabel>0 & ID~=10 ;
    for j = 1:6
        if any(trueLabel(kp)==j)
            [X,Y,T,AUC(i,j)] = perfcurve(trueLabel(kp)==j,confLabel(kp,j),1);
        else
            AUC(i,j) = nan;
        end
        
    end
end

%%
C = confusionmat(trueLabel(trueLabel>0& ID~=10),estimateLabel(trueLabel>0& ID~=10));
%C = confusionmat(trueLabel(trueLabel>0& (ID~=10&ID~=43)),estimateLabel(trueLabel>0& (ID~=10&ID~=43)));

%C = confusionmat(trueLabel(trueLabel>0),estimateLabel(trueLabel>0));

%%
close all
figure
imagesc(100*nanmean(C./sum(C,2),3),[0 100])
set(gca,'ydir','normal')
axis square

set(gca,'fontsize',16)
colormap('hot')
colorbar
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
