% this function will load a feature file to train the RUSboost model to
% predict which time bin the data originated from relative to the nearest
% seizure


% relies on the output of sm_PredictIHKA_getAllFeatures which calculates 
% and saves the feature file 
% features are loaded into the variable 'dat' which is is a cell array
% with each element as a category (time bin) to predict




FeatureFileOutput = 'G:\data\isip\oneTreeData\features.mat';
load(FeatureFileOutput)
ops.ClassifierFileOutput =  'G:\data\isip\oneTreeData\classification.mat';

%%

% extract the features (each group is an element of the cell array)
ops.nGroup  = length(dat);
training = cell2mat(dat');

%define the groups (1:length(dat))
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');


% randomly sort for cross validation
ops.rix  = randsample(1:length(group),length(group));
training = training(ops.rix,:);
group = group(ops.rix);


%set up classifer


ops.N = round(size(training,1)/2);         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';

%train model
rusTree = fitcensemble(training(1:ops.N,:),group(1:ops.N,:),'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate,'ScoreTransform','logit');




%save (does not resave features, does save filename)
save(ops.ClassifierFileOutput,'ops','rusTree','training','group','sessions','-v7.3')

