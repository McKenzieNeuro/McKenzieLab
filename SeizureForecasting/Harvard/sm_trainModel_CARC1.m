function sm_trainModel_CARC1(idx)


outfil = ['/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/model_' num2str(idx) '.mat'];

feat_dir = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/features';
feat_fil = getAllExtFiles(feat_dir,'mat',0);
feat_fil = feat_fil(contains(feat_fil,'_features.mat'));
%%

%get random subset of file
rng(1)

ridx = randsample(1:length(feat_fil),length(feat_fil));
sub_idx = (idx-1)*100+1:(idx-1)*100+100;
ridx = ridx(sub_idx);

feat_fil = feat_fil(ridx);


%
maxi = length(feat_fil);
for i = 1:length(feat_fil)
    try
        v(i) = load(feat_fil{i});
        idx = cellfun(@(a) randsample(1:size(a,1),min(100,size(a,1))),v(i).dat,'uni',0);
        v(i).dat = cellfun(@(a,b) a(b,:),v(i).dat,idx,'uni',0);
        v(i).sesID = cellfun(@(a,b) a(b,:),v(i).sesID,idx,'uni',0);
        ok = whos('v');
        %load <5Gb
        if ok.bytes/1e9>5
            maxi = i;
            break
        end
    catch
       disp('here')
    end
end
%kp = arrayfun(@(a) ~isempty(a.sesID),v);
%v = v(kp);



    %%
    for i = 1:6
        dat{i} = cell2mat(arrayfun(@(a) a.dat{i},v,'uni',0)');
    end
    
    v = rmfield(v,'dat');
    
    for i = 1:6
        output.sesID{i} = cell2mat(arrayfun(@(a) a.sesID{i},v,'uni',0)');
    end
    
    %%
    
    for i = 1:length(v)
        sessions_sub(i,:) =  v(i).sessions(i,:);
        sessions_idx(i) = v(i).i;
    end
    
  %  output.feature_fil = fils(1:maxi);
    output.raw_fil = sessions_sub;
    output.sessions_idx = sessions_idx;
    output.feature_ops = rmfield(v(1).ops,'outname');
    output.sessions_all = v(1).sessions;
    
    clear v
    %%
    
    group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');
    % extract the features (each group is an element of the cell array)
    ops.nGroup  = length(dat);
    training = cell2mat(dat');
    
    clear dat
    
    %define the groups (1:length(dat))
    
    
    
    
    
    
    %set up classifer
    
    
    ops.N = round(size(training,1));         % Number of observations in the training sample
    ops.t = templateTree('MaxNumSplits',ops.N);
    ops.NumLearningCycles = 100;
    ops.Learners = ops.t;
    ops.LearnRate = 0.1;
    ops.Method = 'RUSBoost';
    %train model
    
    rusTree = fitcensemble(training,group,'Method',ops.Method, ...
        'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate,'ScoreTransform','logit');
    
    output.model_ops = ops;
    output.model = rusTree;
    save(outfil,'output','sessions_idx','-v7.3')
    disp(['save: ' outfil])
end

