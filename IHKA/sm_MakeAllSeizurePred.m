% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
ClassifierFileOutput =  'E:\data\IHKA\classification.mat';
load(ClassifierFileOutput)

FeatureFileOutput = 'E:\data\IHKA\features.mat';
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
    