% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
%ClassifierFileOutputDir =  'C:\Users\AlMaynes\Documents\NeuroLab\models';
ClassifierFileOutputDir =  'R:\IHKA_Scharfman\classification';


%FeatureFileOutput = 'C:\Users\AlMaynes\Documents\NeuroLab\features\features.mat';
FeatureFileOutput = 'R:\IHKA_Scharfman\features\features.mat';
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
   % outfil = ['C:\Users\AlMaynes\Documents\NeuroLab\predictions\predict_' num2str(i) '.mat'];
    outfil = ['R:\McKenzieLab\IHKA_Scharfman\prediction\predict_' num2str(i) '.mat'];
    save(outfil,'estimateLabel','trueLabel','time2seizure','inTrainingSet','seizure_start')
    disp([' saved: ' outfil])
end


%%
