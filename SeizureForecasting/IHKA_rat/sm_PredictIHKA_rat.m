% this function will load a feature file to train the RUSboost model to
% predict which time bin the data originated from relative to the nearest
% seizure


% relies on the output of sm_PredictIHKA_getAllFeatures which calculates
% and saves the feature file
% features are loaded into the variable 'dat' which is is a cell array
% with each element as a category (time bin) to predict




FeatureFileOutput = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\featuresPCP4.mat';
%FeatureFileOutput = 'E:\data\IHKA\features_trans.mat';
load(FeatureFileOutput)
ops.ClassifierFileOutput =  'R:\Analysis\SeizureForecasting\IHKA_rat_RF\modelPCP4.mat';

%%
sesID1 = cell2mat(sesID');
% extract the features (each group is an element of the cell array)
ops.nGroup  = length(dat);


training = cell2mat(dat');

%define the groups (1:length(dat))
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');


% randomly sort for cross validation
%ops.rix  = (mod(sesID1(:,1),2)==1);
ops.rix = true(size(group,1),1);
%ops.rix = sesID1(:,1)~=1;

%set up classifer


ops.N = sum(ops.rix );         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N);
ops.NumLearningCycles = 100;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';

%train model
rusTree = fitcensemble(training(ops.rix ,:),group(ops.rix ,:),'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate);




%save (does not resave features, does save filename)
save(ops.ClassifierFileOutput,'ops','rusTree','training','group','sessions','-v7.3')
%%

[label,conf] = predict(rusTree,training(ops.rix ,:));
actual = group(ops.rix ,:);
C = confusionmat(actual,label);
C1 = C./nansum(C,2);
C2 = C./nansum(C,1);
%%
for i = 1:6
    [X,Y,~,auc(i)] = perfcurve(actual==i,nanmean(conf(:,i),2),1);
end
%%
figure
imagesc(C1)

%%
d = LoadBinary('Labels_4Ch_1Hz.dat','nchannels',4,'channels',1,'frequency',1);
dur = length(d);
%%
ops.ClassifierFileOutput =  'R:\Analysis\SeizureForecasting\IHKA_rat_RF\modelPCP4.mat';
fname = [pwd filesep 'amplifier.lfp'];
load(ops.ClassifierFileOutput)
%%
label = nan(dur,1);
conf = nan(dur,6);
for i = 1:dur
    
    
    feat(i,:) = sm_GetDataFeature_rat(fname,i,ops);
    
    
    [label(i),conf(i,:)] = predict(rusTree,feat(i,:));
    i
end
save('prediction.mat','label','conf')

%%
close all
k = gaussian2Dfilter([100 1],.1);
dirN{1} = 'R:\DGregg\NeuralData\PCP\Recordings\ClosedLoop_model_0uA\9-27-2024(16.44)\RHS_240927_164533';
dirN{2} = 'R:\DGregg\NeuralData\PCP\Recordings\ClosedLoop_model_87uA\9-28-2024(6.8)\RHS_240928_061000';
dirN{3} = 'R:\DGregg\NeuralData\PCP\Recordings\ClosedLoop_model_87uA\10-1-2024(6.8)\RHS_241001_061000';
dirN{4} = 'R:\DGregg\NeuralData\PCP\Recordings\ClosedLoop_model_87uA\10-3-2024(6.8)\RHS_241003_061000';
dirN{5} = 'R:\DGregg\NeuralData\PCP\Recordings\ClosedLoop_model_0uA\10-2-2024(6.8)\RHS_241002_061000';
dirN{6} = 'R:\DGregg\NeuralData\PCP\Recordings\ClosedLoop_model_0-87uA\10-7-2024(6.3)\RHS_241007_060500';
dirN{7} = 'R:\DGregg\NeuralData\PCP\Recordings\ClosedLoop_model_0-87uA\10-8-2024(6.3)\RHS_241008_060500';
dirN{8} = 'R:\DGregg\NeuralData\PCP\Recordings\ClosedLoop_model_0-87uA\10-9-2024(6.3)\RHS_241009_060500';
stimD = [0 1 1 1 0 1 0 1];
figure
for i =8
    cd(dirN{i})
    d = LoadBinary('Labels_4Ch_1Hz.dat','nchannels',4,'channels',1,'frequency',1);
    
    if i<=5
        s = LoadBinary('stim_8Ch_1KHz_.dat','frequency',1000,'nchannels',8,'channels',7);
    else
        s = LoadBinary('stim_8Ch_1KHz_.dat','frequency',1000,'nchannels',40,'channels',7);
    end
    stim = find(diff([0;s]>.5))/1000;
    if any(stim)
        kp = [1;diff(stim)>10];
        stim = stim(kp==1);
        
        %fix timing
        [ix,early,late,ts] = sm_getIndicesAroundEvent(stim,300,300,1,length(d));
        kp = ~early & ~late;
        ok = d(ix(kp,:));
        stim = stim(kp);
        
        if i==1
            
            ok = ok(:,293:315)==400;
            offset = 7;
        elseif i==2
            
            ok = ok(:,294:316)==400;
            offset = 6;
        elseif i==3
            
            ok = ok(:,321:372)==400;
            offset = -21;
        elseif i==4
            ok = ok(:,321:372)==400;
            offset = -21;
        else
            offset = [];
        end
        
        if ~isempty(offset)
            real_ts =[];
            for ii = 1:size(ok)
                
                tmp = ok(ii,:);
                
                ix = [find( diff([-1 tmp]==0)>0)'  find(diff([tmp -1]==0)<0)'];
                ix = ix(diff(ix,[],2)>10);
                ix = ix(1);
                ixx(ii) = ix;
                real_ts(ii) = stim(ii)+ix+offset;
            end
            a=polyfit(stim,real_ts',1);
            fs = a(1);
        else
            fs = 1;
        end
        
    else
        fs = fs;
    end
    
    if ~isempty(offset)
        dur = length(d);
        ts_label = (1:length(d))/fs - 4;
        d1 = interp1(ts_label',double(d),1:dur,'nearest');
        
    else
        d1 = double(d);
    end
    
    [ix,early,late,ts] = sm_getIndicesAroundEvent(stim,600,600,1,length(d1));
    kp = ~early & ~late;
    shg
    hold on
    pr = d1==400;
    pr = nanconvn(pr,k');
    
    ses(i).Stim = pr(ix(kp,:));
    % pr = (pr-nanmean(pr))./nanmean(pr);
    figure
    plot(ts,nanmean(pr(ix(kp,:))),'k')
    
    
    %find surrogate stims
    
    alldur = [0 length(d1)];
    if any(stim)
        noStim = excludeEpochs([0 length(d1)],[stim(1) stim(end)]);
    else
        noStim = [0 length(d1)];
    end
    
    %find 3 in a row
    potentialStim = find(d1(1:end-2)==400 & d1(2:end-1)==400 & d1(3:end)==400);
    
    %5min ISI
    for j = 1:length(potentialStim)
        bd= potentialStim>potentialStim(j) & potentialStim<potentialStim(j)+300;
        potentialStim(bd) = nan;
    end
    potentialStim(isnan(potentialStim)) =[];
    potentialStim = potentialStim+4;
    kp = InIntervals(potentialStim,noStim);
    potentialStim = potentialStim(kp);
    [ix,early,late,ts] = sm_getIndicesAroundEvent(potentialStim,600,600,1,length(d));
    kp = ~early & ~late;
    hold on
    plot(ts,nanmean(pr(ix(kp,:))),'g')
    
    ses(i).stimDay = stimD(i);
    ses(i).mockStim = pr(ix(kp,:));
end
%%

ok1 = [ses(2).mockStim;ses(3).mockStim;ses(4).mockStim;ses(6).mockStim;ses(8).mockStim];
ok2 =[ses(2).Stim;ses(3).Stim;ses(4).Stim;ses(6).Stim;ses(8).Stim];
ok3 = [ses(1).Stim;ses(5).mockStim;ses(7).mockStim];
figure
plotMeanSEM(ts,(ok1),'k')
hold on
plotMeanSEM(ts,(ok2),'r')
plotMeanSEM(ts,(ok3),'b')
xlim([-60 100])
%%
[ix,early,late,ts] = sm_getIndicesAroundEvent(stim,600,600,1,length(d));
for i = 1:6
    tmp = conf(:,i);
    c(i,:) = nanmean(tmp(ix));
end

%%

%%
close all
k = gaussian2Dfilter([1000 1],1);
pr = nanconvn(d1==400,k');

[ix1,early,late,ts] = sm_getIndicesAroundEvent(stim,300,300,1,length(d));
kp = ~early & ~late;
ix1=ix1(kp,:);
[ix,early,late,ts] = sm_getIndicesAroundEvent(potentialStim,300,300,1,length(d));
ix = ix(kp,:);
figure
plotMeanSEM(ts,(pr(ix(kp,:))),'k')
hold on
tmp = pr(ix1);

plotMeanSEM(ts(1:300),tmp(:,1:300),'r')
plotMeanSEM(ts(316:end),tmp(:,316:end),'r')
xlim([-20 200])
plot([0 0],[0 1])
plot([14 14],[0 1])

for i = 1:600
    [p(i)] = ranksum(ok(:,i),tmp(:,i));
end

hold on

%%
figure
h = sm_EventTriggeredSpect('amplifier.lfp',stim,'channel',5,'plotIntervals',[20 200]);

figure
h = sm_EventTriggeredSpect('amplifier.lfp',potentialStim,'channel',5,'plotIntervals',[20 200]);

%%
close all
pr = d1==400;
figure
bar([nanmean(pr(1:floor(stim(1)))) nanmean(pr(floor(stim(1):floor(stim(end)))))  nanmean(pr(floor(stim(end):end)))])
