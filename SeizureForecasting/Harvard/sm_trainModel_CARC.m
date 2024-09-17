function sm_trainModel_CARC1(idx)


master_feats = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/bids/allfeat.mat';
outfil = ['/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/model_' num2str(idx) '.mat'];
%load all feauture file (from sm_collectFeats)
load(master_feats)
%%

%get random subset of file
rng(1)

ridx = randsample(1:length(fils),length(fils));
sub_idx = (idx-1)*100+1:(idx-1)*100+100;
ridx = ridx(sub_idx);

fils = fils(ridx);


%
maxi = length(fils);
for i = 1:length(fils)
    
    v(i) = load(fils{i});
    
    ok = whos('v');
    %load <10Gb
    if ok.bytes/1e9>9
        maxi = i;
        break
    end
end
%%
for i = 1:6
    dat{i} = cell2mat(arrayfun(@(a) a.dat{i},v,'uni',0)');
end
for i = 1:6
    output.sesID{i} = cell2mat(arrayfun(@(a) a.sesID{i},v,'uni',0)');
end

%%

for i = 1:length(v)
    sessions_sub(i,:) =  v(i).sessions(i,:);
    sessions_idx(i) = v(i).i;
end

output.feature_fil = fils(1:maxi);
output.raw_fil = sessions_sub;
output.sessions_idx = sessions_idx;
output.feature_ops = rmfield(v(1).ops,'outname');
output.sessions_all = v(1).sessions;
%%


% extract the features (each group is an element of the cell array)
ops.nGroup  = length(dat);
training = cell2mat(dat');

%define the groups (1:length(dat))
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');





%set up classifer


ops.N = round(size(training,1));         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',10);
ops.NumLearningCycles = 1;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';

%train model
rusTree = fitcensemble(training,group,'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate,'ScoreTransform','logit');
output.model_ops = ops;
output.model = rusTree;
save(outfil,'output','sessions_idx')
end
