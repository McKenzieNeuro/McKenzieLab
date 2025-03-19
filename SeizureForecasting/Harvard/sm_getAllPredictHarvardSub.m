function sm_getAllPredictHarvardSub(idx)

outfil = ['/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/forecast_summary_' num2str(idx) '.mat'];
%raw = load('/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/features/raw_sessions3.mat');

% get all models
modelDir = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models';
models = [modelDir filesep 'model_' num2str(idx) '.mat'];

v = load(models);

output = v.output;
clear v
%get all feature files
feat_dir = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/features';
feat_fil = getAllExtFiles(feat_dir,'mat',0);
feat_fil = feat_fil(contains(feat_fil,'_features.mat'));

%loop over feature files and over models

for i = 1:length(feat_fil)
    try
        prediction(i).fname = feat_fil{i};
        
        
        for k = 1:6
            
            v=load(feat_fil{i},'dat','sesID');
            d = v.dat{k};
            s = v.sesID{k};
            %         %clip d
            %
            %         if size(d,1)>10
            %             idx = randsample(1:size(d,1),10);
            %             d = d(idx,:);
            %             s = s(idx,:);
            %         end
            
            clear v
            if ~isempty(d)
                
                
                prediction(i).model.fname = models;
                
                [prediction(i).model.outpred{k},prediction(i).model.conf{k}] = predict(output.model,d);
                prediction(i).model.sesID{k} = s;
                
                
                
            end
        end
        if mod(i,100)==0
            save(outfil,'prediction','-v7.3')
            i
        end
    end
end
save(outfil,'prediction','-v7.3')
end
