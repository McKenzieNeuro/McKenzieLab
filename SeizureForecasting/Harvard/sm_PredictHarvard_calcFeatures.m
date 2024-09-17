function features = sm_PredictHarvard_calcFeatures(fname,tim,ops)
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
Fs  = ops.Fs_data;
durFeat = ops.durFeat;
ch_phaseAmp = ops.ch_phaseAmp;
nCh_raw = ops.nCh_raw;
amp_idx = ops.amp_idx;
ph_idx = ops.ph_idx;

%%



%loop over channels
features =nan(20,nCh_raw);
comod = nan(length(amp_idx),length(ph_idx));

datfil = matfile(fname);
idx1 = floor(tim*Fs);
idx2 = floor((tim+durFeat)*Fs);
rD = datfil.data(:,idx1:idx2);
for ch = 1:nCh_raw
    
    
    spec = nan(length(ops.freqs)-1,1);
    wavespec_amp = nan(size(rD,2),length(ops.freqs)-1);
    filtered_phase = nan(size(rD,2),length(ops.freqs)-1);
    
    for fr = 1:length(ops.freqs)-1
        wavespec_amp(:,fr) = InstAmplitude(BandpassFilter(rD(ch,:),Fs,[ops.freqs(fr) ops.freqs(fr+1)]));
        filtered_phase(:,fr) = InstPhase(BandpassFilter(rD(ch,:),Fs,[ops.freqs(fr) ops.freqs(fr+1)]));
        spec(fr) = nanmean(  wavespec_amp(:,fr));
    end
    
    
    %save mean power (every other channel)
    features(:,ch) = spec;
    
    
    
    
    %calculate phase/amplitude correlation in time window
    
    if ch==ch_phaseAmp
        
        ix1 = 1;
        
    
        
        
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
nFreq = length(ops.freqs)-1;
nPairs = nchoosek(ops.nCh_raw,2);
cxy = nan(nFreq,nPairs); ix=1;
for ch1 = 1:ops.nCh_raw
    for ch2 = ch1+1:ops.nCh_raw
        
        if ~(any(isnan(rD(:,ch2)))) && ~(any(isnan(rD(:,ch1))))
            [cxy(:,ix),f] = mscohere(rD(:,ch1),rD(:,ch2),[],[],ops.freqs(1:end-1),Fs);
        end
        ix= ix+1;
        
    end
end


% full feature space
% [Ch1-Ch4 PSD, phase ampl corr., coherence]

features = [features(:)' comod(:)' cxy(:)'];

end
