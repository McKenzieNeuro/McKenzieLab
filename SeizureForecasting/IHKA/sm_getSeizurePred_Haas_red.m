function [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred_Haas_red(fname,seizFil,rusTree,ops)

% this function takes the classifier in rusTree and applease the feature
% space specified in ops and classified ever moment in time for the file
% (fname)
%%

Fs =ops.Fs;
bins = ops.bins;
nCh_featureFil = ops.nCh_featureFile;

%%


powerFil = [fname '_1.dat'];
s = dir(powerFil);
dur = s.bytes/nCh_featureFil/Fs/2;



ts = 0: (dur-ops.durFeat);


if ~isempty(seizFil)
    %get seizure times
    [~,~,ok] = xlsread(seizFil);
    ok = ok(1:2:end,:);
    tmp = (ok(contains(ok(:,1),'load_starts'),2)');
    seizure_start = sort(cell2mat(cellfun(@(a) extractNumFromStr(a),tmp,'uni',0)));
    tmp = (ok(contains(ok(:,1),'load_stop'),2)');
    seizure_end = sort(cell2mat(cellfun(@(a) extractNumFromStr(a),tmp,'uni',0)));
    
    
    
    % get true time to seizure
    
    time2seizure = ts;
    kp = true(size(ts));
    for i = 1:length(seizure_start)
        time2seizure(kp&ts<seizure_start(i)) = ts(kp&ts<seizure_start(i)) - seizure_start(i);
        time2seizure(kp& ts>seizure_start(i) & ts<seizure_end(i)) = .5;
        time2seizure(kp& ts>seizure_end(i) & ts<seizure_end(i)+10) = 1.5;
        
        
        kp(ts<seizure_start(i) | ...
            (ts>seizure_start(i) & ts<seizure_end(i)) | ...
            (ts>seizure_end(i) & ts<seizure_end(i)+10)) = false;
    end
    if ~any(isinf(bins))
        [~,trueLabel] = histc(time2seizure,[-inf -(bins) 0 1 2]);
    else
        [~,trueLabel] = histc(time2seizure,[ -(bins) 0 1 2]);
    end
    trainingTime =[];
    if isempty(trainingTime)
        inTrainingSet = false(size(ts));
    else
        inTrainingSet = histc(trainingTime,ts)>0;
        
        
    end
    
    
else
    
    trueLabel = nan(length(ts),1);
    inTrainingSet = nan(length(ts),1);
    time2seizure = nan(length(ts),1);
    seizure_start = nan;
    
end
%%


%get prediction
estimateLabel =[];
dat1 =[];
for i = ts
    
    
    
    
    tim = i;
    features = ops.features(fname,tim,ops);
    
    
    dat1 = [dat1;features];
    
    if mod(i,100)== 0
        [outpred,conf] = predict(rusTree,dat1);
        estimateLabel = [estimateLabel;outpred conf];
        dat1 =[];
    elseif i > (dur- mod(dur,100))
        [outpred,conf] = predict(rusTree,dat1);
        estimateLabel = [estimateLabel;outpred conf];
        dat1 =[];
    end
    
end




end





