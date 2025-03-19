function sm_getAllPredictHarvard(idx)
modelDir = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models';

%load global sessions
raw = load('/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/features/raw_sessions.mat');

%get currect file
fname = raw.sessions{idx,2};
[readDir,sesName] = fileparts(raw.sessions{idx,2});

outfil = [readDir filesep sesName '_prediction.mat'];
sz_tim = raw.tims{idx}(:,5:6);
% get all models
models = getAllExtFiles(modelDir,'mat',0);
kp = contains(models,'model_');
models = models(kp);


%choose model that was not trained on session
kp = false(size(models));
for i = 1:length(models)
    try
        v = load(models{i},'sessions_idx');
        train_ses = fileparts(raw.sessions(v.sessions_idx',1));
        kp(i) = ~any(ismember(readDir,train_ses));
    end
end

models = models(kp);

%load relevant models
trainingTime = [];
for i = 1:length(models)
    
    v = load(models{i});
    disp('here')
    [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePredHarvard(fname,sz_tim,v.output.model,trainingTime,v.output.feature_ops);
    pred.estimateLabel{i} = estimateLabel;
    pred.trueLabel = trueLabel;
    pred.time2seizure = time2seizure;
    pred.seizure_start = seizure_start;
    pred.model_fname{i} = models{i};
    
    
    
    clear v
    save(outfil,'pred','-v7.3')
    disp(['save: ' outfil ' model: ' num2str(i)])
end
