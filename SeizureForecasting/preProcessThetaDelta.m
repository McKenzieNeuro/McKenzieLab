dirN{1} = ('R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023\10-3-2023(8.19)\RHS_231003_082031');
dirN{2} = ('R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023\10-3-2023(8.47)\RHS_231003_084903');
dirN{3} = 'R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023\10-3-2023(13.25)\RHS_231003_132613'
datfils = cellfun(@(a) [a filesep 'amplifier.dat'],dirN,'uni',0)';
stimfils = cellfun(@(a) [a filesep 'stim.dat'],dirN,'uni',0)';


%merge
sm_ConcatDats(datfils,'R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023\amplifier.dat') % merge neural data
sm_ConcatDats(stimfils,'R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023\stim.dat') % merge stim times

%make LFP
bz_LFPfromDat('R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023','basename','amplifier'); % down sample neural

% get stims
ch  = [7 15 23 ]; % these are the LDT stim channels for each rat
on = sm_getStimTimes(['R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023' filesep 'stim.dat'],ch,'Fs',20000,'nchannels',24 );


%%
k = gaussian2Dfilter([100000 1],5000);
ch = [6 14 22];

ok = dir('R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023\amplifier.lfp');
nSamp = ok.bytes/24/2;

TD = nan(nSamp,4);

for i = 1:length(ch)
    d = LoadBinary('R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023\amplifier.lfp','nchannels',24,'frequency',1250,'channels',ch(i));
    
    
    theta = BandpassFilter(double(d),1250,[5 12]);
    power_theta = InstAmplitude(theta);
    ts = (1:length(theta))/1250;
    delta = BandpassFilter(double(d),1250,[1 4]);
    power_delta = InstAmplitude(delta);
    
    
    
    power_theta = nanconvn(power_theta,k);
    power_delta = nanconvn(power_delta,k);
    
    TD(:,i) = power_theta./power_delta;
end

ts = (1:size(TD,1))/1250;
ev = LoadEvents('R:\DGregg\NeuralData\EDS_Cohort2\10-3-2023\pulseTimes.evt.sti');
%%
stimOn = 17955;
stimOFF = 21552;  
    pre = nanmean(TD(ts<stimOn,:));
    stim = nanmean(TD(ts>stimOn & ts<stimOFF,:));
      post = nanmean(TD( ts>stimOFF,:));

 %%
% close all
 figure
 plot([0 0 0 0;1 1 1 1; 2 2 2 2],[pre ;stim; post])
% 
 ylabel('Theta/delta')
 set(gca,'xtick',0:2,'xticklabel',{'pre','stim','post'})
% %%
% ts = (1:size(TD,1))/1250;
% ts0 = (1:size(TD0,1))/1250;
% ts2 = (1:size(TD2,1))/1250;
% %%
% 
% hold on
for i = 1:4
    figure
plot(ts,TD(:,i))
hold on
plot([stimOn stimOn],[0 6],'k')
plot([stimOFF stimOFF],[0 6],'k')
end

