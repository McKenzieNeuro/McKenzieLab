% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


% MakeAll_getPowerPerChannel must have been called on all raw data files to
% prepare feature files





%path where raw data is stored with seizure labels
ops.RawDataPath = 'E:\Dropbox\Data SamMcKenzie\Data_raw';
ops.FeaturePath = {'G:\data\IHKA_Haas'};
ops.FeatureFileOutput = 'G:\data\IHKA_Haas\features_red.mat';


% define time windows around to predict (s)

ops.bins = [900 600 300 10];


% define % of time to take within each time bin
%   1-3hrs:     5%
%   10min-1hr:  5%
%   10s-10min:  20%
%   0-10s:      100%
%   seizure:    100%
%   post ictal: 20%

ops.pct = [.5 .5 .5 1 1 1];

ops.nBins = length(ops.pct);

%define postIctal time, 600s after seizure offset
ops.timPost = 60;

% info for full feature space
ops.Fs  = 500; % sampling rate of edf files
ops.reScalePhase = 1000;  %
ops.reScalePower = 1000;  % 32767/maxPower < 1,000


% frequency bands for coherence. must match sm_getPowerPerChannel
ops.freqs = logspace(log10(.5),log10(200),20);

% frequency selection for phase/amplitude (must match index in getPowerPerChannel)
ops.amp_idx = 15:20;
ops.ph_idx = 1:10;

ops.durFeat = 2; % 2s feature bins

%channel for phase amplitude coupling
%ch1 = ipsi HPC
%ch2 = ipsi cortex
%ch3 = contra HPC
%ch4 = contral cortex

ops.ch_phaseAmp =2;

%information about feature file (getPowerPerChannel)
ops.nCh_featureFile = 41; %20 power, 20 phase, 1 time series
ops.nCh_raw = 2; % N = 4 eeg channels

ops.features = @sm_PredictIHKA_calcFeatures;




%DEPENDENCIES
% getAllExtFiles, sm_PredictIHKA_calcFeatures,sm_getPowerPerChannel,sm_MakeAll_getPowerPerChannel,
% getXPctTim, circ_corrcl

%%




% find all edf files and txt files (with seizure labels)
fils_h5 = getAllExtFiles(ops.RawDataPath,'.h5',1);


% get unique session
fils_h5 = unique(cellfun(@(a,b) a(1:b(3)-1),fils_h5,regexp(fils_h5,'_'),'UniformOutput',false));

%%

sessions = [];

%find feature files that match annotation file
for j = 1:length(ops.FeaturePath)
    d = dir(ops.FeaturePath{j});
    d = {d(cell2mat({d.isdir})).name}';
    [a,seizure_fils] = fileparts(fils_h5);
    
    [~,b] = ismember(seizure_fils,d);
    kp = b>0;
    b(~kp) =[];
    seizure_filst = fils_h5(kp);
    
    
    % save list of file names for the paired annotation/feature files
    
    sessions = [sessions;seizure_filst cellfun(@(a) [ops.FeaturePath{j} filesep a filesep a],d(b),'uni',0)];
    
    
    
end

nSessions = size(sessions,1);


%%

warning off
%get all relevant timepoints around the seizure start and end

% TS name in file.
TSname1 = 'Seizure starts';
TSname2 = 'ends';

T = readtable('E:\Dropbox\Data SamMcKenzie\Focal_seizures_start_stop_times.xlsx');
% loop over all seizures
for i= 1:size(sessions,1)
    
    %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
    tims{i} = [];
    
    [a,b] = fileparts(sessions{i,1});
    
    
    
    
    seizure_start = table2array(T(contains(table2cell(T(:,1)),b),2));
    seizure_end = table2array(T(contains(table2cell(T(:,1)),b),3));
    
    
    %loop over seizures per subject
    for j = 1:length(seizure_start)
        
        %bins around seizure
        tmp = [seizure_start(j)-ops.bins seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
        
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
    sz_on{i} = seizure_start;
    
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
for i = 1:nSessions
    fname = sessions{i,2};
    
    %loop of time bins
    for k = 1:ops.nBins
        
        
        %choose random subset (pct) of samples within each session to train
        %model
        
        sz{k} = getXPctTim(tims{i}(:,k:k+1), ops.pct(k),1);
        
        
        
        
        
        
        %loop over all timepoints for each session for each bin
        for ev = 1:numel(sz{k})
            
            %find duration of the file
            powerFil = [fname '_1.dat'];
            s = dir(powerFil);
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
    i
end


%%
filenameps = ['G:\data\IHKA_Haas\all_seiz.ps'];
filenamepdf= ['G:\data\IHKA_Haas\all_seiz.pdf'];
for i = 1:nSessions
    for j = 1:length(sz_on{i})
       h= plot_sz_gross(fileparts(sessions{i,2}),sz_on{i}(j),'fs',ops.Fs,'nCh_rec',3);
        print(h, '-dpsc2',filenameps ,'-append')
        close all
    end
end
sm_ps2pdf(filenameps,filenamepdf,[])