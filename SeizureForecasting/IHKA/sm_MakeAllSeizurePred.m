% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
ClassifierFileOutput =  'E:\data\IHKA\classification_trans.mat';
load(ClassifierFileOutput)

FeatureFileOutput = 'E:\data\IHKA\features_trans.mat';
load(FeatureFileOutput,'sesID')

%%

warning off
% loop over files to predict
for i = 1:size(sessions,1)
    featureFile =  sessions{i,2};
    seizureFile = sessions{i,1};
    
    %get times used in training
    trainingTime = sort(cell2mat(cellfun(@(a) a(a(:,1)==i,2),sesID,'UniformOutput',false)'));
    [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred(featureFile,seizureFile,rusTree,trainingTime,ops);
    [dirOut,fileOut] = fileparts(sessions{i,2});
    save([dirOut filesep fileOut '_predict.mat'],'estimateLabel','trueLabel','time2seizure','inTrainingSet','seizure_start')
    disp([' saved: ' dirOut filesep fileOut '_predict.mat'])
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
    
     load('E:\data\IHKA\classifier1.mat')
     
    %%
    alldat = training;
   
  
    training = training(ops.N+1:end,:);
    
    [pred,score] = predict(rusTree,training);

actual_Y = group(ops.N+1:end);
predicted_Y = pred;
C = confusionmat(actual_Y,predicted_Y);


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
figure
imagesc(C./sum(C,2),[0 1])
set(gca,'yticklabel',{'3hrs','1hr','10min','10s','sz','post'})
set(gca,'xticklabel',{'3hrs','1hr','10min','10s','sz','post'})
xlabel('Predicted class')
ylabel('Real class')

set(gca,'ydir','normal')