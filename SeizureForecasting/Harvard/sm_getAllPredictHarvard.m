function sm_getAllPredictHarvard(idx)

%load global sessions
raw = load('/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/features/raw_sessions.mat');

%get currect file
fname = raw.sessions{idx,2};
readDir = fileparts(raw.sessions{idx,1});


% get all models
models = getAllExtFiles('/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models','mat',0);
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

for i = 1:length(models)
    
      v = load(models{i});
      
       features = sm_PredictHarvard_calcFeatures(fname,tim,v.output.feature_ops)
      clear v
end

%choose file


% loop over file
output.feature_ops.features