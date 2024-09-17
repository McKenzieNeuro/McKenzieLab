% this function will load a feature file to train the RUSboost model to
% predict which time bin the data originated from relative to the nearest
% seizure


% relies on the output of sm_PredictIHKA_getAllFeatures which calculates 
% and saves the feature file 
% features are loaded into the variable 'dat' which is is a cell array
% with each element as a category (time bin) to predict


function sm_PredictIHKA_CARC(idx)

FeatureFileOutput = '/carc/scratch/projects/mckenzie2016183/data/IHKA_Haas/features/features.mat';
%FeatureFileOutput = 'E:\data\IHKA\features_trans.mat';
load(FeatureFileOutput)
ops.ClassifierFileOutput =  ['/carc/scratch/projects/mckenzie2016183/data/IHKA_Haas/classification/classification_' num2str(idx) '.mat'];

%%

% extract the features (each group is an element of the cell array)
ops.nGroup  = length(dat);
all_dat = cell2mat(dat');

%define the groups (1:length(dat))
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');

sesID = cell2mat(sesID');
% randomly sort for cross validation
ops.rix  = sesID(:,1)~=idx;

training = all_dat(ops.rix,:);
group_training = group(ops.rix);

test = all_dat(~ops.rix,:);
group_test = group(~ops.rix);

%set up classifer


ops.N = sum(ops.rix );         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';

%train model
rusTree = fitcensemble(training,group_training,'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate);



estimateLabel = predict(rusTree,test);

C = confusionmat(group_test(group_test>0),estimateLabel(group_test>0));
C_norm = C./repmat(sum(C,2),1,size(C,2));



%save (does not resave features, does save filename)
save(ops.ClassifierFileOutput,'ops','rusTree','dat','estimateLabel','C','C_norm','sessions','-v7.3')

end