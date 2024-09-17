function [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred_rat(fname,seizFil,rusTree,trainingTime,ops)

% this function takes the classifier in rusTree and applease the feature
% space specified in ops and classified ever moment in time for the file
% (fname)
%%

Fs =ops.Fs;
bins = ops.bins;
nCh_featureFil = ops.nCh_featureFile;
ops.subjName = 'EDS4.0';
%%


ev = LoadEvents(seizFil);
kp_on = contains(ev.description,ops.subjName) & contains(ev.description,'sz_on');
kp_off = contains(ev.description,ops.subjName) & contains(ev.description,'sz_off');


seizure_start = ev.time(kp_on);
seizure_end = ev.time(kp_off);

%%

% get true time to seizure
s = dir(fname);
dur = s.bytes/nCh_featureFil/Fs/2;



ts = 0: (dur-ops.durFeat);
time2seizure = ts;
kp = true(size(ts));
for i = 1:length(seizure_start)
    time2seizure(kp&ts<seizure_start(i)) = ts(kp&ts<seizure_start(i)) - seizure_start(i);
    time2seizure(kp& ts>seizure_start(i) & ts<seizure_end(i)) = .5;
    time2seizure(kp& ts>seizure_end(i) & ts<seizure_end(i)+600) = 1.5;
    
    
    kp(ts<seizure_start(i) | ...
        (ts>seizure_start(i) & ts<seizure_end(i)) | ...
        (ts>seizure_end(i) & ts<seizure_end(i)+600)) = false;
end
if ~any(isinf(bins))
    [~,trueLabel] = histc(time2seizure,[-inf -(bins) 0 1 2]);
else
    [~,trueLabel] = histc(time2seizure,[ -(bins) 0 1 2]);
end

if isempty(trainingTime)
    inTrainingSet = false(size(ts));
else
    inTrainingSet = histc(trainingTime,ts)>0;
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





