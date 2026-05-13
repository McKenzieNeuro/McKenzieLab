% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


% MakeAll_getPowerPerChannel must have been called on all raw data files to
% prepare feature files





%path where raw data is stored with seizure labels
%ops.RawDataPath = 'R:\DGregg\NeuralData\EDS_Cohort2';
% ops.RawDataPath = 'R:\DGregg\NeuralData\PCP\Recordings\';
ops.RawDataPath = {'R:\DGregg\NeuralData\PTP\Recordings','R:\DGregg\NeuralData\PCP\Recordings'};
ops.FeaturePath = {[]};
ops.FeatureFileOutput = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\features_pilo2.mat';
ops.bin_label = ['inter-ictal','ictal','noise'];
ops.pct = [.01 1 .2];

ops.nBins = length(ops.pct);

% info for full feature space
ops.Fs  = 1250; % sampling rate of lfp files

% frequency bands for coherence. must match sm_getPowerPerChannel
ops.freqs = logspace(log10(2),log10(300),20);

% frequency selection for phase/amplitude (must match index in getPowerPerChannel)

ops.durFeat = 4; % 4s feature bins
ops.art_thres = 5e4; % NOTE probably needs to be decreased 

%information about feature file (getPowerPerChannel)
ops.nCh_featureFile = 32; % if no XML assume 4 rats @ 8ch each
ops.ch_subj = 1:8; % this needs to be updated for >1 rat per dat file
ops.nCh_raw = length(ops.ch_subj); % N = 8 eeg channels

ops.features = @sm_GetDataFeature_rat;




%DEPENDENCIES
% getAllExtFiles, sm_getPowerPerChannel,sm_MakeAll_getPowerPerChannel,
% getXPctTim, circ_corrcl

%%




% find all lfp files and evt.szr  files (with seizure labels)
fils_szr =[]; fils_lfp =[];

for i = 1:length(ops.RawDataPath)
fils_lfp = [fils_lfp;getAllExtFiles(ops.RawDataPath{i},'lfp',1)];
fils_szr = [fils_szr;getAllExtFiles(ops.RawDataPath{i},'szr',1)];
end
% find all sessiondata files

% find all lfp files with annotations
[b_lfp] = cellfun(@fileparts,fils_lfp,'uni',0);
[b_szr] = cellfun(@fileparts,fils_szr,'uni',0);

% find all sessiondata files (ALEX)
fils = [];
for i = 1:length(b_lfp)
    
    tmp = getAllExtFiles(b_lfp{i},'mat',1);
    
    kp = contains(tmp,'sessiondata');
    fils = [fils;tmp(kp)];
end
% subjID =  'PTP 4.2';
% 
% kp = contains(fils,subjID);
% 
% fils1 = fils(kp);
sessionfils = cellfun(@fileparts,fils,'uni',0);


goodFils = intersect(b_szr,b_lfp);
goodFils = intersect(goodFils, sessionfils); % keep only files with subject sessiondata (ALEX)
seizure_fils = fils_szr(ismember(b_szr,goodFils));
lfp_fils = fils_lfp(ismember(b_lfp,goodFils));
session_fils = fils(ismember(sessionfils,goodFils)); % keep only files with subject sessiondata (ALEX)





%%

sessions = [];

%find feature files that match annotation file
for j = 1:length(seizure_fils)
    sessions{j,1} = seizure_fils{j};
    sessions{j,2} = lfp_fils{ismember(fileparts(lfp_fils),fileparts(seizure_fils{j}) )};
    sessions{j,3} = session_fils{ismember(fileparts(session_fils), fileparts(seizure_fils{j}))}; % Add sessiondata to the fileparts (ALEX)
end

nSessions = size(sessions,1);


%%

warning off
%get all relevant timepoints around the seizure start and end
%which subject to train
% subjName = 'PTP4.2';
clear tims

% loop over all seizures
for i= 1:size(sessions,1)
    
    %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
    tims{i} = [];
    
    
    
    ev = LoadEvents(sessions{i,1});
    %  kp_on = contains(ev.description,subjName) & contains(ev.description,'sz_on');
    %  kp_off = contains(ev.description,subjName) & contains(ev.description,'sz_off');
    if ~isempty(ev.description)
        kp_sz_on = contains(ev.description,'on') & contains(ev.description,'sz');
        kp_sz_off =  contains(ev.description,'off')& contains(ev.description,'sz');
        kp_noise_on = contains(ev.description,'on') & contains(ev.description,'noise');
        kp_noise_off =  contains(ev.description,'off')& contains(ev.description,'noise');
        
        seizure_start = ev.time(kp_sz_on);
        seizure_end = ev.time(kp_sz_off);
        noise_start = ev.time(kp_noise_on);
        noise_end = ev.time(kp_noise_off);
        
        
        if ~isempty(seizure_start)
            tims{i} = [tims{i};seizure_start seizure_end ones(length(seizure_end),1)];
        end
        
        if ~isempty(noise_start)
            tims{i} = [tims{i};noise_start noise_end 2*ones(length(noise_start),1)];
        end
        
    else
        tims{i} = [];
    end
    % ops.ch_subj = v.sessiondata.channelID; % Save the channel ID numbers for the subject (ALEX)
end
%%

% clear feature variable that will be used to aggregate across sesions


for k = 1:ops.nBins
    
    dat{k} =[];
    sesID{k} = [];
end

%%
%loop over number of sessions
for i = 1:nSessions
% for i = nSessions
    clear sz
    sz = cell(1,3);
    % for i = nSessions
    fname = sessions{i,2};
    fxml = strrep(fname,'lfp','xml');
    xml = LoadXml(fxml);
    v = load(sessions{i,3});
    % ops.ch_subj = v.sessiondata.channelID;
    ops.ch_subj = v.sessiondata.channelID(5:6); %only grab hippocampal
    %loop of time bins
    for k = 1:ops.nBins
    % for k = 2
        
        
        %choose random subset (pct) of samples within each session to train
        %model
        if k ==1
            fileDur = sm_getFileDur(sessions{i,2});
            interictal = excludeEpochs([0 fileDur],tims{1}(:,1:2));
            sz{k} = getXPctTim(interictal, ops.pct(k),1);
        else
            
            if ~isempty(tims{i})
                sz{k}  = getXPctTim(tims{i}(tims{i}(:,3)==(k-1),1:2), ops.pct(k),1);
            end
        end
        
        
        
        
        %loop over all timepoints for each session for each bin
        for ev = 1:numel(sz{k})
            
            %find duration of the file
            
            s = dir(fname);
            
            % GET number of channels
            
            dur = s.bytes/xml.nChannels/ops.Fs/2;
            
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
        
        % save the full feature space
        save(ops.FeatureFileOutput,'dat','sesID','sessions','ops','-v7.3')
    end
end



