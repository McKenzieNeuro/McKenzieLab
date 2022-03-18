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




%%


% loop over files to predict
for i = 1:size(sessions,1)
    featureFile =  sessions{i,2};
    seizureFile = sessions{i,1};
    [estimateLabel,trueLabel,time2seizure] = sm_SeizurePred(featureFile,seizureFile,rusTree,ops);
    
    save([dirOut filesep d{b(i)} '_predict.mat'],'estimateLabel','trueLabel','time2seizure')
    disp([' saved: E:\Dropbox\UNM\Analysis\IHKA\data\' d{b(i)} '_predict.mat'])
end