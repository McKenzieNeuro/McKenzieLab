% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
ClassifierFileOutput =  'F:\data1\IHKA\classification.mat';
load(ClassifierFileOutput)

FeatureFileOutput = 'F:\data1\IHKA\features.mat';
load(FeatureFileOutput,'sesID')

%%


% loop over files to predict
for i = 1:size(sessions,1)
    featureFile =  sessions{i,2};
    seizureFile = sessions{i,1};
    
    %get times used in training
    trainingTime = sort(cell2mat(cellfun(@(a) a(a(:,1)==i,2),sesID,'UniformOutput',false)'));
    [estimateLabel,trueLabel,inTrainingSet,time2seizure] = sm_getSeizurePred(featureFile,seizureFile,rusTree,trainingTime,ops);
    [dirOut,fileOut] = fileparts(sessions{i,2});
    save([dirOut filesep fileOut '_predict.mat'],'estimateLabel','trueLabel','time2seizure','inTrainingSet')
    disp([' saved: ' dirOut filesep fileOut '_predict.mat'])
end


%%

