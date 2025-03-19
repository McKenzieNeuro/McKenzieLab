% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
ClassifierFileOutputDir =  'R:\IHKA_Scharfman\classification';


FeatureFileOutput = 'E:\data\IHKA\featuresPCP.mat';
load(FeatureFileOutput)

%%

warning off
% loop over files to predict
for i = 1:size(sessions,1)
    ClassifierFileOutput = [ClassifierFileOutputDir filesep 'classification_' num2str(i) '.mat'];
    load(ClassifierFileOutput)
    featureFile =  sessions{i,2};
    seizureFile = sessions{i,1};
    
    %get times used in training
   % trainingTime = sort(cell2mat(cellfun(@(a) a(a(:,1)==i,2),sesID,'UniformOutput',false)'));
    trainingTime = [];
    [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred(featureFile,seizureFile,rusTree,trainingTime,ops);
    outfil = ['R:\IHKA_Scharfman\prediction\predict_' num2str(i) '.mat'];
    save(outfil,'estimateLabel','trueLabel','time2seizure','inTrainingSet','seizure_start')
    disp([' saved: ' outfil])
end


%%

% load all predictions
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
nexthrvuln = nan(size(trueLabel));
for i = 1:length(trueLabel)-600
nexthrvuln(i) = sum(trueLabel(i+1:i+601)==5);
end

%%
XX = ind.*confLabel(:,1:4);
mdl = fitglm(XX,double(nexthrvuln>0),'y ~ x1 +x2+x3+x4','Distribution','binomial');
mdl1 = fitglm(ones(size(nexthrvuln)),double(nexthrvuln>0),'y ~ x1 ','Distribution','poisson');

%%
figure
semilogx(-10000:0,(avghist(time2seizure,predict(mdl,XX),-10000:0))/3600)

hold on
semilogx(-10000:0,avghist(time2seizure,predict(mdl,XX),-10000:0)./avghist(time2seizure,predict(mdl1,ones(size(time2seizure))),-10000:0))

%%
clear AUC
for i = 1:max(ID)
     kp = ID ==i &trueLabel>0 ;
for j = 1:6
   if any(trueLabel(kp)==j)
[X,Y,T,AUC(i,j)] = perfcurve(trueLabel(kp)==j,confLabel(kp,j),1); 
   else
       AUC(i,j) = nan;
   end

end
end

%%
C = confusionmat(trueLabel(trueLabel>0),estimateLabel(trueLabel>0));

figure
imagesc(C./sum(C,2),[0 1])
set(gca,'yticklabel',{'3hrs','1hr','10min','10s','sz','post'})
set(gca,'xticklabel',{'3hrs','1hr','10min','10s','sz','post'})
xlabel('Predicted class')
ylabel('Real class')

set(gca,'ydir','normal')

%%
close all
figure
k  = gaussian2Dfilter([1000 1],0);
for i = 1:5
    d = double(estimateLabel==i);
    %d = nanconvn(d,k);
semilogx(fliplr(avghist(time2seizure,d,-43000:0)))
hold on
end

%%

figure
k  = gaussian2Dfilter([1000 1],3);
for i = 1:5
    d = double(estimateLabel==i);
   
semilogx(fliplr(avghist(time2seizure,d,-86400:0,@numel)))
hold on
end


%%



close all
figure
k = gaussian2Dfilter([1 100],[ 1 50]);
%estimateLabel1 = estimateLabel;
%estimateLabel1(inTrainingSet) = nan;
for i = 1:6
    subplot(6,1,i)
   
    hold on
    plot([seizure_start seizure_start],[0 1],'--','color','r')
     plot(nanconvn(estimateLabel==i,k'),'k')
end



%%


%get ROC''
clear FP TP AUC T AUCtmp lb hb
ix=1;
for i=[1 2 3 4 5 6]
    tmp =  group(ops.N+1:end) ;
%tmp(tmp==2) = 3;
[FP{ix},TP{ix},T{ix},AUC(ix)] = perfcurve(tmp,score(:,i),i);


% get random

for j = 1:100
    tmpt = tmp(randsample(1:length(tmp),length(tmp)));
    [~,~,~,AUCt(j)] = perfcurve(tmpt,score(:,i),i);
end
AUCtmp(ix) = mean(AUCt);
lb(ix) = prctile(AUCt,1);
hb(ix) = prctile(AUCt,99);
ix = ix+1;
end

hb = hb-AUCtmp;
lb = lb-AUCtmp;


%%

    %%
    
     load('E:\data\IHKA\classifier1.mat')
     
    %%
   
   
  
  
    [pred,score] = predict(rusTree,training);

actual_Y = group;
predicted_Y = pred;

%%

%%
close all

figure
x=-([18000,3600 600 10 2 1]);
y= AUC;
semilogx(x,y,'o-')
hold on
xx = [x';flipud(x')];
zz = [AUCtmp'+lb';flipud(AUCtmp')+hb'];



patch(xx,zz,'k','facealpha',.25,'edgecolor','none')
semilogx(x,AUCtmp,'--','color','k')
xlim([-18000 0])
%%


%%

ok = cell2mat(dat');
X  = tsne(nanzscore(ok));
group = cellfun(@(a,b) a*ones(length(b),1),num2cell(1:6),dat,'UniformOutput',false);
group = cell2mat(group');
[~,b] = histc(group,1:6);
kp = ~any(isnan(ok),2);
%%
close all
k  = gaussian2Dfilter([100 100],5);
ax  = tight_subplot(1,6);
for  ii = 1:6
    axes(ax(ii))
bin = histcn([X(b(kp)==ii,1),X(b(kp)==ii,2)],-50:50,-50:50)/sum(b(kp)==ii);
PP(ii,:) = bin(:);
imagesc(nanconvn(bin,k))
end

%%

for i = 1:6
    
    for j = 1:6
        
     kl(i,j) =   KLDiv(PP(i,:),PP(j,:));
        
    end
end

%%