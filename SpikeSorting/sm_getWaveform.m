
function [wfF, shank,channelCorr] = sm_getWaveform(datfil,clu,ts,good)
%%
xml = LoadXml([datfil(1:end-3) 'xml']);
nCh = xml.nChannels;
fs = xml.SampleRate;
nPull = 1000; % number of spikes to pull out
wfWin = 0.008; % Larger size of waveform windows for filterning
filtFreq = 500;
hpFilt = designfilt('highpassiir','FilterOrder',3, 'PassbandFrequency',filtFreq,'PassbandRipple',0.1, 'SampleRate',fs);

wfWin = round((wfWin * fs)/2);
uclu = unique(clu);
nclu = length(uclu);



badCh = cell2mat({xml.AnatGrps.Skip});
shank = nan(nclu,1); channelCorr = nan(nclu,1);
for ii = 1 : nclu
    
    if good(ii)
        spkTmp = ts(clu ==uclu(ii));
        spkTmp = spkTmp(spkTmp< (max(ts) -.005));
        spkTmp = round(spkTmp*fs);
        if length(spkTmp) > nPull
            spkTmp = spkTmp(randperm(length(spkTmp)));
            spkTmp = spkTmp(1:nPull);
        end
        wf = nan((wfWin * 2)+1,nCh,length(spkTmp));
        for jj = 1 : length(spkTmp)
            
            wf = cat(3,wf,bz_LoadBinary(datfil,'offset',spkTmp(jj) - (wfWin),...
                'samples',(wfWin * 2)+1,'frequency',fs,'nChannels',nCh));
        end
        
        tmp = mean(wf,3) - mean(wf(1,:,:),3);
        wf = mean(wf,3);
        
        [~,b] = min(min(tmp,[],1));
        shank(ii) = find(cellfun(@(a) any(a+1==b),{xml.AnatGrps.Channels}));
        
        ch = xml.AnatGrps(shank(ii)).Channels;
        ok = corr(wf(110:130,~badCh));
        
        
        channelCorr(ii) = mean(ok(triu(true(size(ok)),1)));
        wf = wf(:,ch+1);
        
        
        wfF{ii} = [];
        for jj = 1 : size(wf,2)
            wfF{ii}(:,jj) = filtfilt(hpFilt,wf(:,jj) - mean(wf(:,jj)));
        end
        
        wfF{ii} = wfF{ii}(110:130,:);
        
    else
        wfF{ii} = [];
    end
    
end
%%
end