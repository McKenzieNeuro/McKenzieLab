% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
ClassifierFileOutput =  'G:\data\IHKA_Haas\classification_red.mat';
load(ClassifierFileOutput)

FeatureFileOutput = 'G:\data\IHKA_Haas\features_red.mat';
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
    
    alldat =  cell2mat(dat');
    group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');
    testix = ops.rix(ops.N+1:end);
    training = alldat(testix,:);
    
    [pred,score] = predict(rusTree,training);

actual_Y = group(testix);
predicted_Y = pred;
C = confusionmat(actual_Y,predicted_Y);

%%
figure
imagesc(C./sum(C,2),[0 1])
set(gca,'yticklabel',{'3hrs','1hr','10min','10s','sz','post'})
set(gca,'xticklabel',{'3hrs','1hr','10min','10s','sz','post'})
xlabel('Predicted class')
ylabel('Real class')

set(gca,'ydir','normal')

%%

Cr = nan(6,6,1000);
for i = 1:1000
    
    % get conf
   actual_Yr =  actual_Y(randsample(1:length(actual_Y),length(actual_Y)));
   tmp = confusionmat(actual_Yr,predicted_Y);
   
   tmp = tmp./sum(tmp,2);
    Cr(:,:,i) =tmp;
    i
end
%%
 

lab = [num2cell(ops.bins),{'sz'},{'post'}];
ok = C./sum(C,2);
ok1 = eye(6)';
okRu = nanmean(Cr,3);
close all
ax  =tight_subplot(6,1);
ixx = 1;
for i = [1:6]
    
    
    axes(ax(i))

semilogx([-ops.bins(1:4) -.1 -.01], 100* ((ok(i,1:6)- okRu(i,1:6))./okRu(i,1:6)),'o-')
    hold on
   % semilogx([-ops.bins(1:4) -.1], 100* ((ok(i,1:5)- ok1(i,1:5))./ok(i,1:5)),'x-')
    x = [-ops.bins(1:4) -.1 -.01]';
    z1 = -prctile(Cr(i,:,:),5,3);
z2 = prctile(Cr(i,:,:),95,3);
    xx = [x;flipud(x)];
zz = 100*[z1';z2'];

col = 'k';
plot([-ops.bins(1:4) -.1 -.01],[0 0 0 0 0 0],'k')
patch(xx,zz,col,'facealpha',.25,'edgecolor','none')

    % semilogx([-ops.bins(1:4) -.1], okRu(i,1:5),'--','color','r')

if i<=4
title(['Model predicts ' num2str(lab{i}) 's to seiz.'])
elseif i==5
title(['Model predicts seiz.'])    
else
    title(['Model predicts post ict.'])    
end
if i ==6
xlabel('Prediction')
end
%xlim([-700 1])
xlim([-1800 0])
if i <5
%ylim([0 .45])

set(gca,'xtick',[-ops.bins(1:4) -.1 -.01],'xticklabel',[])
else
   %ylim([0 .8])
%   ylim(100*[-1 2])
   set(gca,'xtick',[-ops.bins(1:4) -.1 -.01],'xticklabel',lab(1:6))
   xlabel('Real time to seizure')
end
ylim(100*[-1 3.5])
ixx = ixx+1;
end


%%

