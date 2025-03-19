% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


% MakeAll_getPowerPerChannel must have been called on all raw data files to
% prepare feature files



%DEPENDENCIES
% getAllExtFiles, sm_PredictIHKA_calcFeatures,sm_getPowerPerChannel,sm_MakeAll_getPowerPerChannel,
% getXPctTim, circ_corrcl


warning off
%path where raw data is stored with seizure labels
ops.RawDataPath = 'R:\DGregg\NeuralData\PCP\Recordings';
ops.FeaturePath = {[]};
ops.FeatureFileOutput = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\featuresPCP4.mat';


% define time windows around to predict (s)
ops.bins = [ inf 3600 100 10];


% define % of time to take within each time bin
%   1-3hrs:     5%
%   10min-1hr:  5%
%   10s-10min:  20%
%   0-10s:      100%
%   seizure:    100%
%   post ictal: 20%

ops.pct = [.3 .3 .5 1 1 .5];

ops.nBins = length(ops.pct);

%define postIctal time, 600s after seizure offset
ops.timPost = 300;

% info for full feature space
ops.Fs  = 1250; % sampling rate of lfp files



% frequency bands for coherence. must match sm_getPowerPerChannel
ops.freqs = logspace(log10(2),log10(300),20);

% frequency selection for phase/amplitude (must match index in getPowerPerChannel)


ops.durFeat = 4; % 4s feature bins
ops.art_thres = 5e4; % define amplitude of raw data that is too big to be real (artifact)

%information about feature file (getPowerPerChannel)
ops.nCh_featureFile = 8; % if no XML assume 4 rats @ 8ch each
ops.ch_subj = 1:8;
ops.nCh_raw = length(ops.ch_subj); % N = 8 eeg channels

ops.features = @sm_GetDataFeature_rat;



%%




% find all lfp files and evt.szr  files (with seizure labels)
fils_lfp = getAllExtFiles(ops.RawDataPath,'lfp',1);
topDir = 'R:\DGregg\NeuralData\PCP\Recordings';
fils = getAllExtFiles(topDir,'evt',1);
kp = contains(fils,'BayesOpt');
fils = fils(~kp);

kp = false(length(fils),1);
for  i =1:length(fils)
    dirN = fileparts(fils{i});
    cd(dirN)
    
    ev = LoadEvents(fils{i});
    
    if ~isempty(ev.time)
        kp(i) = true;
    end
    
    
end
fils_szr = fils(kp);
%%
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


% loop over all seizures
for i= 1:size(sessions,1)
    
    %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
    tims{i} = [];
    
    
    
    ev = LoadEvents(sessions{i,1});
    
    
    
    %loop over seizures per subject
    if ~isempty(ev.time)
        
        kp_on = contains(lower(ev.description),'seizurestart');
        kp_off = contains(lower(ev.description),'seizurestop');
        
        
        seizure_start = ev.time(kp_on);
        seizure_end = ev.time(kp_off);
        
        for j = 1:length(seizure_start)
            
            %bins around seizure
            
            if isinf(ops.bins(1))
                
                if j ==1 && seizure_start(j) > ops.bins(2)
                    tmp = [0 seizure_start(j)-ops.bins(2:end) seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
                    
                elseif j ==1 && seizure_end(j) < ops.bins(2)
                    tmp = [nan seizure_start(j)-ops.bins(2:end) seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
                elseif j >1 && seizure_start(j)-ops.bins(2) > seizure_end(j-1) +ops.timPost
                    tmp = [seizure_start(j-1)+ops.timPost seizure_start(j)-ops.bins(2:end) seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
                else
                    tmp = [nan seizure_start(j)-ops.bins(2:end) seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
                end
            else
                tmp = [seizure_start(j)-ops.bins seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
            end
            %if seizure is early, delete time
            tmp(tmp<0) = nan;
            
            
            %if bin overlaps with prior seizure, delete
            if j>1
                tmp(tmp<seizure_start(j-1)) = nan;
            end
            
            % if bin overlaps with next seizure, delete
            if j<length(seizure_start)
                tmp(tmp>seizure_start(j+1)) = nan;
            end
            
            %save time bins per seizure per subject
            tims{i} = [tims{i};tmp];
            
            
        end
    else
        
        in = dir(sessions{i,2});
        xmlf = strrep(sessions{i,2},'lfp','xml');
        xml = LoadXml(xmlf);
        dur = in.bytes/xml.lfpSampleRate/xml.nChannels/2;
        tmp = [0 dur-ops.bins(2) nan nan nan nan nan];
        tims{i} = [tims{i};tmp];
    end
end

%%

% clear feature variable that will be used to aggregate across sesions
clear sz


for k = 1:ops.nBins
    
    dat{k} =[];
    sesID{k} = [];
end

%%


%loop over number of sessions
for i = 7:nSessions
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



