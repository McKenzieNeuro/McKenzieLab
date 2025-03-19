function [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePredHarvard(fname,seizFil,rusTree,trainingTime,ops)

% this function takes the classifier in rusTree and applease the feature
% space specified in ops and classified ever moment in time for the file
% (fname)
%%

Fs =ops.Fs_data;
bins = ops.bins;


%%


seizure_start = seizFil(:,1);
seizure_end =  seizFil(:,2);

%%

% get true time to seizure
datfil = matfile(fname);
dur = size(datfil,'data');
dur = dur(2)/datfil.Fs;

ts = 1: (dur-ops.durFeat);
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





