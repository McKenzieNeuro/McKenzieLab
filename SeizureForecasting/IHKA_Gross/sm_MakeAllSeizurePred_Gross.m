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
ClassifierFileOutputDir =  'R:\IHKA_gross\classification';


%FeatureFileOutput = 'C:\Users\AlMaynes\Documents\NeuroLab\features\features.mat';
FeatureFileOutput = 'R:\IHKA_gross\features1.mat';
load(FeatureFileOutput)

%%

warning off
% loop over files to predict
for i = 1:size(sessions,2)
     outfil = ['R:\IHKA_gross\prediction\predict_' num2str(i) '.mat'];
    if ~exist(outfil)
        try
     ClassifierFileOutput = [ClassifierFileOutputDir filesep 'classification_' num2str(i) '.mat'];
    load(ClassifierFileOutput)
    featureFile =  sessions{i};
    seizureFile = sessions{i};
    featureFile = strrep(featureFile,'G:\data','R:');
    seizureFile = strrep(seizureFile,'G:\data','R:');
    %get times used in training
   % trainingTime = sort(cell2mat(cellfun(@(a) a(a(:,1)==i,2),sesID,'UniformOutput',false)'));
    trainingTime = [];
    [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred_gross(featureFile,seizureFile,rusTree,trainingTime,ops);
   % outfil = ['C:\Users\AlMaynes\Documents\NeuroLab\predictions\predict_' num2str(i) '.mat'];
   
    save(outfil,'estimateLabel','trueLabel','time2seizure','inTrainingSet','seizure_start')
    disp([' saved: ' outfil])
        catch
             disp([' failed: ' outfil])
        end
    else
         disp([' cached: ' outfil])
    end
end


%%
