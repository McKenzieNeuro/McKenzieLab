function output = sm_getIHKARegressors2(files)


vars =  [...
    {'TD'},...
    {'vel1'},...
    {'sleepState1_2'},...
    {'sleepState1_3'},...
    ];

vars1 =  [...
    
{'sleepState1_1'},...
{'sleepState1_2'},...
];



dat_all = []; ses_id = []; in_Stim =[];
for  ii = 1:length(files)
    ii
    v = load(files{ii});
    sessiondata = v.sessiondata;
    bins = 0:max(sessiondata.ts);
    ied = sort(cell2mat(sessiondata.IED'));
    ied_bin = histc(ied,bins);
    try
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
    
    if ~isempty(sessiondata.stimON)
        kp2  = ~InIntervals(bins, [sessiondata.stimON-.1 sessiondata.stimON+2])';
    else
        kp2 = true(size(kp));
    end
    
    dat = table(ied_bin,TD,vel1,ep,sleepState1);
    dat1 = dat(kp &kp1 & kp2,:);
    dat_all = [dat_all;dat1];
    
    ses_id = [ses_id;ii*ones(size(dat1,1),1)];
    
    
    % not stim epochs
    
    if contains(sessiondata.StimType,'mono') & sessiondata.StimAmp>=100
        stim_ep = [sessiondata.stimON(1) sessiondata.stimON(end)];
        inStimt = InIntervals(bins,stim_ep);
    else
        inStimt = false(length(bins),1);
    end
    in_Stim = [in_Stim;inStimt(kp &kp1 & kp2)];
    
    end
end

% do CV model comparison

for i = 1:length(files)
    
    kp = ses_id~=i;
    dat1 = dat_all(kp,:);
    mdl_full = fitglme(dat1,'ied_bin ~    TD + vel1 + sleepState1  ','link','log','distribution','poisson');
    sleep = fitglme(dat1,'ied_bin ~    sleepState1  ','link','log','distribution','poisson');
    vel = fitglme(dat1,'ied_bin ~    vel1 ','link','log','distribution','poisson');
    TD = fitglme(dat1,'ied_bin ~     TD ','link','log','distribution','poisson');
    base = fitglme(dat1,'ied_bin ~     1 ','link','log','distribution','poisson');
    
    
    y_hat = predict( mdl_full,dat_all(~kp,:));
    
    %  MSE_full(i) = mean((dat_all.ied_bin(~kp)-y_hat).^2);
    output.MSE_full(i) = poisson_log_likelihood(dat_all.ied_bin(~kp),y_hat,1);
    
    y_hat = predict( sleep,dat_all(~kp,:));
    
    output.MSE_Sleep(i) = poisson_log_likelihood(dat_all.ied_bin(~kp),y_hat,1);
    
    y_hat = predict( vel,dat_all(~kp,:));
    output.MSE_Vel(i) =  poisson_log_likelihood(dat_all.ied_bin(~kp),y_hat,1);
    
    y_hat = predict( TD,dat_all(~kp,:));
    
    output.MSE_TD(i) =  poisson_log_likelihood(dat_all.ied_bin(~kp),y_hat,1);
    
    y_hat = predict( base,dat_all(~kp,:));
    
    output.MSE_base(i) =  poisson_log_likelihood(dat_all.ied_bin(~kp),y_hat,1);
    
    
    
    %now just get stim
    
    kp_stim = ~kp & in_Stim;
    
    if any(kp_stim)
        y_hat = predict( mdl_full,dat_all(kp_stim,:));
        
        %  MSE_full(i) = mean((dat_all.ied_bin(~kp)-y_hat).^2);
        output.MSE_full_stim(i) = poisson_log_likelihood(dat_all.ied_bin(kp_stim),y_hat,1);
        
        y_hat = predict( sleep,dat_all(kp_stim,:));
        
        output.MSE_Sleep_stim(i) = poisson_log_likelihood(dat_all.ied_bin(kp_stim),y_hat,1);
        
        y_hat = predict( vel,dat_all(kp_stim,:));
        output.MSE_Vel_stim(i) =  poisson_log_likelihood(dat_all.ied_bin(kp_stim),y_hat,1);
        
        y_hat = predict( TD,dat_all(kp_stim,:));
        
        output.MSE_TD_stim(i) =  poisson_log_likelihood(dat_all.ied_bin(kp_stim),y_hat,1);
        
        y_hat = predict( base,dat_all(kp_stim,:));
        
        output.MSE_base_stim(i) =  poisson_log_likelihood(dat_all.ied_bin(kp_stim),y_hat,1);
    else
        output.MSE_base_stim(i) = nan;
        output.MSE_Vel_stim(i) = nan;
        output.MSE_Sleep_stim(i) = nan;
        output.MSE_TD_stim(i) = nan;
        
        output.MSE_full_stim(i) = nan;
        
        
        
    end
    
    
    
    
    
end


dat1 = dat_all(in_Stim==1,:);

try
    mdl_full = fitglme(dat1,'ied_bin ~  1+ TD + vel1 + sleepState1  ','link','log','distribution','poisson');
    
    [beta,b,c] =   fixedEffects(mdl_full);
    [~,b] = ismember(b.Name,vars);
    
    FE_full(b(b>0)) = beta(b>0);
catch
    FE_full= nan(1,4);
end

try
    sleep = fitglme(dat1,'ied_bin ~    sleepState1  ','link','log','distribution','poisson');
    
    
    [beta,b,c] =   fixedEffects(sleep);
    [~,b] = ismember(b.Name,vars(3:4));
    
    FE_sleep(b(b>0)) = beta(b>0);
catch
    FE_sleep = nan(1,2);
end

output.FE_full = FE_full;
output.FE_sleep = FE_sleep;

dat1 = dat_all(in_Stim==0,:);

try
    mdl_full = fitglme(dat1,'ied_bin ~  1+ TD + vel1 + sleepState1  ','link','log','distribution','poisson');
    
    [beta,b,c] =   fixedEffects(mdl_full);
    [~,b] = ismember(b.Name,vars);
    
    FE_full(b(b>0)) = beta(b>0);
catch
    FE_full= nan(1,4);
end

try
    sleep = fitglme(dat1,'ied_bin ~    sleepState1  ','link','log','distribution','poisson');
    
    
    [beta,b,c] =   fixedEffects(sleep);
    [~,b] = ismember(b.Name,vars(3:4));
    
    FE_sleep(b(b>0)) = beta(b>0);
catch
    FE_sleep = nan(1,2);
end

output.FE_full_baseline = FE_full;
output.FE_sleep_baseline = FE_sleep;


output.subj = sessiondata.subject;
output.StimAmp = sessiondata.StimAmp;
output.StimType = sessiondata.StimType;



end
