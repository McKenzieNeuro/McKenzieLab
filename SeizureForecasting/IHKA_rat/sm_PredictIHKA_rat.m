% this function will load a feature file to train the RUSboost model to
% predict which time bin the data originated from relative to the nearest
% seizure


% relies on the output of sm_PredictIHKA_getAllFeatures which calculates 
% and saves the feature file 
% features are loaded into the variable 'dat' which is is a cell array
% with each element as a category (time bin) to predict




FeatureFileOutput = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\features.mat';
%FeatureFileOutput = 'E:\data\IHKA\features_trans.mat';
load(FeatureFileOutput)
ops.ClassifierFileOutput =  'R:\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';

%%
sesID1 = cell2mat(sesID');
% extract the features (each group is an element of the cell array)
ops.nGroup  = length(dat);


training = cell2mat(dat');

%define the groups (1:length(dat))
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');


% randomly sort for cross validation
ops.rix  = (mod(sesID1(:,1),2)==1);



%set up classifer


ops.N = sum(ops.rix );         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';

%train model
rusTree = fitcensemble(training(ops.rix ,:),group(ops.rix ,:),'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate);




%save (does not resave features, does save filename)
save(ops.ClassifierFileOutput,'ops','rusTree','training','group','sessions','-v7.3')
%%

[label,conf] = predict(rusTree,training(~ops.rix ,:));
actual = group(~ops.rix ,:);
C = confusionmat(actual,label);
C1 = C./nansum(C,2);
C2 = C./nansum(C,1);
%%
for i = 1:6
[X,Y,~,auc(i)] = perfcurve(actual==i,nanmean(conf(:,i),2),1);
end
%%
figure
imagesc(C1)