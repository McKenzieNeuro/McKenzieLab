function output = sm_getIHKARegressors(sessiondata)

bins = 0:max(sessiondata.ts);
ied = sort(cell2mat(sessiondata.IED'));
ied_bin = histc(ied,bins);

vel1 = avghist(sessiondata.poseData.neuralTS, sessiondata.poseData.velocity(1:length(sessiondata.poseData.neuralTS)),bins);
%%
sleepState = sessiondata.sleepData.SleepState.idx.states;
sleepState = avghist(sessiondata.sleepData.SleepState.idx.timestamps, sleepState,bins);
sleepState = round(sleepState);

kp = ismember(sleepState,[1 3 5]);
sleepState1 =sleepState;
sleepState1(sleepState==3) = 2;
sleepState1(sleepState==5) = 3;

epoch_bin = [0 cumsum(sessiondata.seizure_obs)];
[~,ep] = histc(bins, epoch_bin);
TD = avghist(sessiondata.ts,sessiondata.theta_delta,bins);


ep(ep==4) = 3;
ep = categorical(ep(:));
sleepState1 = sleepState1(:);
sleepState1(~kp) = 1;
sleepState1 = categorical(sleepState1);

ied_bin = ied_bin(:);
TD = TD(:);
vel1 = vel1(:);

%cut out seizures

if ~isempty(sessiondata.seizure)
    
    kp1  = ~InIntervals(bins, sessiondata.seizure)';
else
    kp1 = true(size(kp));
end
dat = table(ied_bin,TD,vel1,ep,sleepState1);
dat1 = dat(kp &kp1,:);

output.mdl_full = fitglme(dat1,'ied_bin ~    TD + vel1 + ep +sleepState1 ','link','log','distribution','poisson');
output.no_sleep = fitglme(dat1,'ied_bin ~    TD + vel1 + ep  ','link','log','distribution','poisson');
output.no_vel = fitglme(dat1,'ied_bin ~    TD +  ep +sleepState1 ','link','log','distribution','poisson');
output.no_TD = fitglme(dat1,'ied_bin ~     vel1 + ep +sleepState1 ','link','log','distribution','poisson');
output.no_ep = fitglme(dat1,'ied_bin ~    TD + vel1  +sleepState1 ','link','log','distribution','poisson');
output.subject = sessiondata.subject;
output.StimType = sessiondata.StimType;
output.StimAmp = sessiondata.StimAmp;
output.data  = dat1;
end