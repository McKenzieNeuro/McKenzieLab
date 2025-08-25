topDir = 'R:\DANEHippocampalResponse\NE2h18\Novel Env\';


ch = 20;% ne2h20
ch = 28;% ne2h19
ch = 29;% ne2h18
fils = getAllExtFiles(topDir,'lfp',1);
%%
clear wavspec
for i = 1:length(fils)
    
    cd(fileparts(fils{i}))
    
    if ~exist('TTL_pulse.mat')
        [ups,dwns]  = sm_getDigitalin(pwd,'digitalin.dat',30000,16);
    else
        load('TTL_pulse.mat')
    end
    
    
    
    stim_on = ups{1};
    
    kp  = diff([0;stim_on])>10;
    stim_on = stim_on(kp);
    
    [h,wavspec{i},ts,freqs] = sm_EventTriggeredSpect('amplifier_analogin_auxiliary_int16.lfp',stim_on,'channel',ch);
    close all
    
end

%%
clear NC HC
for i = 1:length(wavspec)
    NC(:,:,i) = wavspec{i}(:,:,1);
    HC(:,:,i) = wavspec{i}(:,:,2);
    
end


%%

figure
nfreq = length(freqs);
imagesc(ts,[],zscore(nanmedian(NC,3),[],2),[-2 2])
set(gca,'ytick',1:10:nfreq,...
    'yticklabel',round(freqs(1:10:end)*10)/10,'ydir','normal')

xlim([-30 30])