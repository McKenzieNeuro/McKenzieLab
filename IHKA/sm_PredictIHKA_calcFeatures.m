function features = sm_PredictIHKA_calcFeatures(fname,tim,ops)
% function calculates feature space from feature file 
% wavelet phase/zscore power are precalculated and are stored in the feature file
% coherence and phase amplitude coupling are cacluted here



% INPUTs
%  fname = basename for feature file, one per channel (fname_ch.dat)
%  tim = time point in seconds
%  ops = struct with meta data about how to calculate features

% OUTPUT
%  features = 1xN vector of all features for this timepoint

% DEPENDENCIES

% LoadBinary, circ_corrcl, sm_PredictIHKA_getAllFeatures (creates ops),
% getPowerPerChannel (creates feature file)

%%

%get feature meta data

nCh_featureFile = ops.nCh_featureFile;
Fs  = ops.Fs;
durFeat = ops.durFeat;
reScalePower = ops.reScalePower;
ch_phaseAmp = ops.ch_phaseAmp;
nCh_raw = ops.nCh_raw;
reScalePhase = ops.reScalePhase;
amp_idx = ops.amp_idx;
ph_idx = ops.ph_idx;

%%

%loop over channels
features =[];
rD = nan(ops.durFeat*ops.Fs,nCh_raw);
for ch = 1:nCh_raw
    powerFil = [fname '_' num2str(ch) '.dat'];
    
    %load data durFeat(5s) before each relevant start point
    tmp = LoadBinary(powerFil,'nchannels',nCh_featureFile,'frequency',Fs,'channels',1:nCh_featureFile,'duration', durFeat,'start',tim);
    
    %save mean power (every other channel)
    features = [features mean(tmp(:,2:2:40))/reScalePower];
    
    %save time series
    rD(:,ch) = tmp(:,1);
    
    
    %calculate phase/amplitude correlation in time window
    
    if ch==ch_phaseAmp
        comod = nan(length(amp_idx),length(ph_idx));
        ix1 = 1;
      
        %get power and phase for each frquency
        filtered_phase = double(tmp(:,3:2:41))/reScalePhase;
        wavespec_amp = double(tmp(:,2:2:40))/reScalePower;
        
        
        % choose subset of frequencies for power
        for apr = amp_idx
            
            
              ix2 = 1;
            
            % choose subset of frequencies for phase
            for idx = ph_idx
                
                
                comod(ix1,ix2) = circ_corrcl(filtered_phase(:,idx),wavespec_amp(:,apr));
                ix2 = ix2+1;
            end
             ix1 = ix1+1;
        end
    end
end

% calculate coherence
nFreq = length(ops.freqs);
nPairs = nchoosek(ops.nCh_raw,2);
cxy = nan(nFreq,nPairs); ix=1;
for ch1 = 1:4
    for ch2 = ch1+1:4
        [cxy(:,ix),f] = mscohere(rD(:,ch1),rD(:,ch2),[],[],ops.freqs,ops.Fs);
        ix= ix+1;
    end
end


% full feature space
% [Ch1-Ch4 PSD, phase ampl corr., coherence]

features = [features comod(:)' cxy(:)'];
end
