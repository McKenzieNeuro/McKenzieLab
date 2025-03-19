function sm_Asses_pred_Harvard
topDir = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/bids';
outfil = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/forecast_summary.mat';
%get all Pred
fils = getAllExtFiles(topDir,'mat',1);
kp = contains(fils,'_prediction.mat');
fils = fils(kp);

for i = 1:length(fils)
    load(fils{i})
    for j = 1:length(pred.estimateLabel)
        
        minL = min(length(pred.trueLabel),length(pred.estimateLabel{j}));
        C = confusionmat(pred.trueLabel(1:minL),pred.estimateLabel{j}(1:minL,1));
        % output(i).fname = fils{i};
        output(i).model(j).confusion = C;
        output(i).model(j).name = pred.model_fname{j};
        output(i).model(j).estimateLabel = pred.estimateLabel;
        output(i).trueLabel = pred.trueLabel;
    end
end

save(outfil,'output')
end