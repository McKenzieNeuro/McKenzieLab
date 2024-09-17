function [estimateLabel,trueLabel,timefromseizure,time2seizure,seizure_start] = sm_getSeizurePred_beforeAfter(fname,seizFil,rusTree,trainingTime,ops)

% this function takes the classifier in rusTree and applease the feature
% space specified in ops and classified ever moment in time for the file
% (fname)
%%

Fs =ops.Fs;
bins = ops.bins;
nCh_featureFil = ops.nCh_featureFile;

%%


warning off
%get all relevant timepoints around the seizure start and end

% TS name in file.
TSname1 = 'ure starts';
TSname2 = 'ends';

     sz = cell(4,5);

    
    
    TSdata = readtable(seizFil);
    TSdata = table2cell(TSdata);
    
    seizure_start = cell2mat(TSdata(cellfun(@any,regexpi(TSdata(:,6),TSname1)),4));
    seizure_end = cell2mat(TSdata(cellfun(@any,regexpi(TSdata(:,6),TSname2)),4));
    
  
    %calc time to seizure start
    
 powerFil = [fname '_1.dat'];
s = dir(powerFil);
dur = s.bytes/nCh_featureFil/Fs/2;



ts = 0: (dur-ops.durFeat);

    clear sz_start sz_end
    for j = 1:length(seizure_start)
          sz_start{j} = nan(size(ts));
        if j==1
            sz_start{j}(ts<seizure_start(j)) =  ts(ts<seizure_start(j)) - seizure_start(j);
        else
            sz_start{j}(ts<seizure_start(j) & ts>seizure_end(j-1)) =  ts(ts<seizure_start(j) & ts>seizure_end(j-1)) - seizure_start(j);
            
        end
    end
    
    
    %calc time from seizure end
    
 
    
    for j = 1:length(seizure_end)
          sz_end{j} = nan(size(ts));
        if j<length(seizure_end)
            sz_end{j}(ts>seizure_end(j) & ts<seizure_start(j+1)) =  ts(ts>seizure_end(j) & ts<seizure_start(j+1)) - seizure_end(j);
        else
            sz_end{j}(ts>seizure_end(j))  =  ts(ts>seizure_end(j)) - seizure_end(j);
            
        end
    end
    
    trueLabel = nan(size(ts));
    for j = 2:length(seizure_end)
        [n,~,~,b] = histcn([sz_start{j}' sz_end{j-1}'],[ops.bins(1:4) 0],[0 ops.bins(5:8)]);
        kp = all(b>0,2);
        ind = sub2ind([4 4],b(kp,1),b(kp,2));
        if all(isnan( trueLabel(kp)))
            trueLabel(kp) = ind;
        else
            error('here')
        end
        
        
    end
    in = InIntervals(ts,[seizure_start seizure_end]);
    trueLabel(in) = 17;
        

%%

% get true time to seizure


time2seizure = nan(size(ts));

kp = true(size(ts));
for i = 1:length(seizure_start)
    time2seizure(kp&ts<seizure_start(i)) = ts(kp&ts<seizure_start(i)) - seizure_start(i);
  
    
    
    
    kp(ts<seizure_start(i) | ...
        (ts>seizure_start(i) & ts<seizure_end(i))) = false;
    
    
    
    
    
end

timefromseizure = nan(size(ts));
kp = true(size(ts));
for i = length(seizure_start):-1:1
   
    
      timefromseizure(kp&ts>seizure_start(i)) = ts(kp&ts>seizure_start(i)) - seizure_start(i);
  
    
    kp(ts>seizure_start(i) | ...
        (ts>seizure_start(i) & ts<seizure_end(i))) = false;
    
    
    
    
    
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
        estimateLabel = [estimateLabel;predict(rusTree,dat1)];
        dat1 =[];
    elseif i > (dur- mod(dur,100))
        estimateLabel = [estimateLabel;predict(rusTree,dat1)];
         dat1 =[];
    end
    i
end




end





