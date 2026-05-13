load('pulseInfo.mat')
load('annotations_transitions.mat')
load('amplifier_analogin_auxiliary_int16.spikes.cellinfo.mat')
%%

% to remove the artificate
for i = 1:length(spikes.times)
    
    tims = spikes.times{i};
    kp = abs(tims-(bestmatch(tims,  pulseInfo.time(:))))<.002;
    spikes.times{i}(kp) = [];
    i
end



%%
maxT  =sm_getFileDur('amplifier_analogin_auxiliary_int16.dat');
ConAep = [all_events.Intan_ts(contains(all_events.event,'A_in')) all_events.Intan_ts(contains(all_events.event,'A_out'))];
ConBep = [all_events.Intan_ts(contains(all_events.event,'B_in')) all_events.Intan_ts(contains(all_events.event,'B_out'))];

allEp = [0 maxT];
HC_ep = excludeEpochs(allEp,[ConAep;ConBep]);

%%
kp = InIntervals(pulseInfo.time(:,1),ConAep);

ConA_pulse  = pulseInfo.time(kp,1);

kp = InIntervals(pulseInfo.time(:,1),ConBep);

ConB_pulse  = pulseInfo.time(kp,1);

kp = InIntervals(pulseInfo.time(:,1),HC_ep);

HC_pulse  = pulseInfo.time(kp,1);
%%

[binnedPop,bin_times]=populationMatrix(spikes,0,.095,1,ConA_pulse+.002);
binnedPopA = squeeze(binnedPop);


[binnedPop,bin_times]=populationMatrix(spikes,0,.095,1,ConB_pulse);
binnedPopB = squeeze(binnedPop);

[binnedPop,bin_times]=populationMatrix(spikes,0,.095,1,HC_pulse);
binnedPopHC = squeeze(binnedPop);
%%




opto_mod = nan(length(spikes.times),3);
opto_sign = nan(length(spikes.times),3);

for s = 1:3
    switch s
        
        case 1
            cond1 = binnedPopA*.095;
            cond2 = binnedPopB*.095;
        case 2
            cond1 = binnedPopA*.095;
            cond2 = binnedPopHC*.095;
        case 3
            cond1 = binnedPopB*.095;
            cond2 = binnedPopHC*.095;
    end
    
    
    for i = 1:length(spikes.times)
        
        % multiple rate by lenght of time bin duration to get back to count

        [opto_mod(i,s),opto_sign(i,s)]  = sm_getOptoMod(cond1(i,:)',cond2(i,:)','paired',false);
        
    end
    
    
end


%%
close all

figure

[binnedPop,bin_times]=populationMatrix(spikes,.5,.5,100,ConA_pulse,'zscore',true); % to normalize set zscore to 'true', for raw spike count set zscore to 'false'
imagesc(bin_times,[],nanmean(binnedPop,3),[-.2 .2])
[binnedPop,bin_times]=populationMatrix(spikes,.5,.5,100,ConB_pulse,'zscore',true); % to normalize set zscore to 'true', for raw spike count set zscore to 'false'
figure

imagesc(bin_times,[],nanmean(binnedPop,3),[-.2 .2])

[binnedPop,bin_times]=populationMatrix(spikes,.5,.5,100,HC_pulse,'zscore',true); % to normalize set zscore to 'true', for raw spike count set zscore to 'false'
figure
imagesc(bin_times,[],nanmean(binnedPop,3),[-.2 .2])

