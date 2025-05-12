% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


% MakeAll_getPowerPerChannel must have been called on all raw data files to
% prepare feature files





%path where raw data is stored with seizure labels
ops.RawDataPath = 'E:\Dropbox\Scn1a_EEG_for_SM\Chandni';
ops.FeaturePath = {[]};
ops.FeatureFileOutput = 'E:\Dropbox\Scn1a_EEG_for_SM\Chandni\analysis\features.mat';


% define time windows around to predict (s)

ops.bins = [ 10 100 1000 inf];


% define % of time to take within each time bin
%   1-3hrs:     5%
%   10min-1hr:  5%
%   10s-10min:  20%
%   0-10s:      100%
%   seizure:    100%
%   post ictal: 20%

ops.pct = [ 1 1 .2 .2];

ops.nBins = length(ops.pct);

%define postIctal time, 600s after seizure offset
ops.timPost = nan;

% info for full feature space
ops.Fs  = 2000; % sampling rate of lfp files



% frequency bands for coherence. must match sm_getPowerPerChannel
ops.freqs = logspace(log10(2),log10(300),20);

% frequency selection for phase/amplitude (must match index in getPowerPerChannel)


ops.durFeat = 4; % 4s feature bins
ops.ops.art_thres = 5e4;

%information about feature file (getPowerPerChannel)
ops.nCh_featureFile = 16; % if no XML assume 4 rats @ 8ch each
ops.ch_subj = 3:4;
ops.nCh_raw = 2; % N = 2 eeg channels

ops.features = @sm_GetDataFeature_Mattis;




%DEPENDENCIES
% getAllExtFiles, sm_PredictIHKA_calcFeatures,sm_getPowerPerChannel,sm_MakeAll_getPowerPerChannel,
% getXPctTim, circ_corrcl

%%




% find all lfp files and evt.szr  files (with seizure labels)
fils_lfp = getAllExtFiles(ops.RawDataPath,'dat',1);
fils_szr = getAllExtFiles(ops.RawDataPath,'szr',1);

% find all lfp files with annotations
[b_lfp] = cellfun(@fileparts,fils_lfp,'uni',0);
[b_szr] = cellfun(@fileparts,fils_szr,'uni',0);
goodFils = intersect(b_szr,b_lfp);
seizure_fils = fils_szr(ismember(b_szr,goodFils));
lfp_fils = fils_lfp(ismember(b_lfp,goodFils));




%%

sessions = [];

%find feature files that match annotation file
for j = 1:length(seizure_fils)
    sessions{j,1} = seizure_fils{j};
    sessions{j,2} = lfp_fils{ismember(fileparts(lfp_fils),fileparts(seizure_fils{j}))};
end

nSessions = size(sessions,1);


%%

warning off
%get all relevant timepoints around the seizure start and end
%which subject to train
subjName = 'EDS4.0';


% loop over all seizures
for i= 1:size(sessions,1)
    
    %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
    tims{i} = [];
    
    
    
    ev = LoadEvents(sessions{i,1});
    kp_on = contains(ev.description,'sz_on');
    kp_off = contains(ev.description,'sz_off');
    
    
    seizure_start = ev.time(kp_on);
    seizure_end = ev.time(kp_off);
    
    % only deal with the last seizure
    
    seizure_start = seizure_start(end);
    seizure_end = seizure_end(end);
    
        tmp = [seizure_end seizure_end+ops.bins];
        
        
        %get file dur 
        
        ok = dir(sessions{i,2});
        dur = ok.bytes/ops.Fs/ops.nCh_raw/2;
        
        tmp(isinf(tmp)) = dur;
        %save time bins per seizure per subject
        tims{i} = [tims{i};tmp];
        
        
    end
    


%%

% clear feature variable that will be used to aggregate across sesions
clear sz dat sesID


for k = 1:ops.nBins
    
    dat{k} =[];
    sesID{k} = [];
end

%%


%loop over number of sessions
for i = 1:nSessions
    fname = sessions{i,2};
    
    if ~isempty(tims{i})
        %loop of time bins
        for k = 1:ops.nBins
            
            
            %choose random subset (pct) of samples within each session to train
            %model
            
            sz{k} = getXPctTim(tims{i}(:,k:k+1), ops.pct(k),1);
            
            
            
            
            
            
            %loop over all timepoints for each session for each bin
            for ev = 1:numel(sz{k})
                
                %find duration of the file
                
                s = dir(fname);
                dur = s.bytes/ops.nCh_featureFile/ops.Fs/2;
                
                tim = sz{k}(ev)-ops.durFeat;
                %make sure the even does not exceed duration of recording
                if ~isnan(tim) && (dur-tim)>ops.durFeat
                    
                    try
                        %get features (some pre calculated, some calculated on the fly)
                        features = ops.features(fname,tim,ops);
                        
                        
                        %save all features for each time bin
                        dat{k} = [dat{k};features];
                        
                        %keep track of which session matches which feature
                        sesID{k} = [sesID{k}; i tim];
                    catch
                        disp(fname)
                        disp(['file length: ' num2str(dur) ' time of read:' num2str(tim)])
                        
                        
                    end
                    
                end
                
                
            end
            
            
            
        end
        
        %save the full feature space
        save(ops.FeatureFileOutput,'dat','sesID','sessions','ops','-v7.3')
    end
end



