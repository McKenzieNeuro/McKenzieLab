function features = sm_PredictIHKA_calcFeaturesStateTrans(fname,tim,ops)
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
    tmp = LoadBinary(powerFil,'nchannels',nCh_featureFile,'frequency',Fs,'channels',1,'duration', durFeat,'start',tim);
    
    
    
    %save time series
    rD(:,ch) = tmp(:,1);
end
A = ltv_adjacency(rD);

features = A(:)';

end
