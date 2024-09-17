% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


% MakeAll_getPowerPerChannel must have been called on all raw data files to
% prepare feature files





%path where raw data is stored with seizure labels
ops.RawDataPath = 'R:\IHKA_gross\KA Spike2 + Videos';
ops.FeaturePath = {'G:\data\IHKA_gross'};
ops.FeatureFileOutput = 'G:\data\IHKA_gross\features1.mat';%8/21/2024


% define time windows around to predict (s)

ops.bins = [ inf 3600 100 10];


% define % of time to take within each time bin
%   1-3hrs:     5%
%   10min-1hr:  5%
%   10s-10min:  20%
%   0-10s:      100%
%   seizure:    100%
%   post ictal: 20%

ops.pct = [.05 .05 .2 1 1 .2];

ops.nBins = length(ops.pct);

%define postIctal time, 600s after seizure offset
ops.timPost = 600;

% info for full feature space
ops.Fs  = 2000; % sampling rate of edf files
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
%ch2 = contra HPC
%ch3 = EEG


ops.ch_phaseAmp = 2;

%information about feature file (getPowerPerChannel)
ops.nCh_featureFile = 41; %20 power, 20 phase, 1 time series
ops.nCh_raw = 3; % N = 4 eeg channels

ops.features = @sm_PredictIHKA_calcFeatures;




%DEPENDENCIES
% getAllExtFiles, sm_PredictIHKA_calcFeatures,sm_getPowerPerChannel,sm_MakeAll_getPowerPerChannel,
% getXPctTim, circ_corrcl

%%




subject_IDs = {'KA5','KA6','KA8','KA12','KA13','KA14','KA11'};
%these are the subject IDs for all the animals on file
%subject_IDs = {'KA11'};  %for testing, only read in one animals data at a time


sz_files = true;                   %whether to load only the files with seizures (true) or all files (false)



cd(ops.RawDataPath)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add paths to necessary folders

% Here's a url to the CED website where you can download the file (after
% giving them your name/email): https://ced.co.uk/upgrades/spike2matson
% I also included the files on the hard drive (CEDMATLAB folder), but you
% may have to download them separately to get them to work - I followed the
% instructions on the s64mat file and didn't have much issue. Let me know
% if you have trouble.

%load CED MATSON library for reading in .smrx files
cedpath = getenv('R:\Analysis\CEDMATLAB\CEDS64ML');        %this part may be where you need to edit;
cedpath = 'R:\Analysis\CEDMATLAB\CEDS64ML';
%you might be able to switch to the manual path location rather than use
%the 'getenv' function; see pg. 1-6 of the s64mat pdf for more info.
addpath(cedpath)
CEDS64LoadLib(cedpath)

%other libraries - again, may have to adjust when you start
%addpath(genpath('D:\file_for_Sam\'))
%%
clear sessions sz
ix=1;
subject_IDs = {'KA5','KA6','KA8','KA12','KA13','KA14','KA11'};

for k1 =1:length(subject_IDs) %loops through all the animal IDs you provide at the top of the code
    
    
    [rows_to_extract, stim_metatable] = read_metatable_KAsz(subject_IDs{k1}, sz_files);
    
    sz_o = cellfun(@str2num,stim_metatable.Seizures_sec_(rows_to_extract),'UniformOutput',false);
    sz_off = cellfun(@str2num,stim_metatable.SeizuresEnd_sec_(rows_to_extract),'UniformOutput',false);
    ID = [stim_metatable(rows_to_extract,:).AnimalIdentification stim_metatable(rows_to_extract,:).TrialNumber];
    for k2 = 1:size(ID,1)
        
        str = [ID{k2,1} filesep ID{k2,1} '_' ID{k2,2}];
        sessions{ix}= [ops.FeaturePath{1} filesep str];
        sz{ix} = [sz_o{k2} sz_off{k2}];
        ix   = ix+1;
    end
    
    
end




%%

% loop over all seizures
clear tims
for i= 1:length(sessions)
    
    %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
    tims{i} = [];
    
    
    seizure_start = sz{i}(:,1);
    seizure_end = sz{i}(:,2);
    
    
    %loop over seizures per subject
    for j = 1:length(seizure_start)
        
        %bins around seizure
        
        if isinf(ops.bins(1))
            
            if j ==1 && seizure_start(j) > ops.bins(2)
                tmp = [0 seizure_start(j)-ops.bins(2:end) seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
                
            elseif j ==1 && seizure_start(j) < ops.bins(2)
                tmp = [nan seizure_start(j)-ops.bins(2:end) seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
            elseif j >1 && seizure_start(j)-ops.bins(2) > seizure_start(j-1) +ops.timPost
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
    
end

%%

% clear feature variable that will be used to aggregate across sesions
clear sz dat


for k = 1:ops.nBins
    
    dat{k} =[];
    sesID{k} = [];
end

%%

nSessions = length(sessions);
%loop over number of sessions
% start at 103 next time
for i = 1:nSessions
    
    [~,basename] = fileparts(sessions{i});
    fname = [sessions{i} filesep basename] ;
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
                    
                    if isempty(features)
                        features = nan(1,180);
                    end
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


%%
ok = cell2mat(dat');
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');
d = tsne(nanzscore(ok));

%%
close all
k  = gaussian2Dfilter([100 100],2);
lab = [{'-inf to -3600s'},{'-3600s to -100s'},{'-100 to -10'},{'-10 to seiz.'},{'seiz'},{'post'}];
for ii = 1:6
    subplot(3,2,ii)
    imagesc(nanconvn(histcn(d(group==ii,:),-50:50,-50:50),k))
    axis off
    title(lab{ii})
end

%%
close all
filenameps = ['G:\data\IHKA_gross\all_seiz4.ps'];
filenamepdf= ['G:\data\IHKA_gross\all_seiz.pdf'];
for i = 1:length(sz)
    
    for j = 1:size(sz{i},1)
        
        h=plot_sz_gross(sessions{i},sz{i}(j,1));
        title([sessions{i} ' time: ' num2str(sz{i}(j,1)) 's'])
        
        print(h, '-dpsc2',filenameps ,'-append')
        close all
    end
end
sm_ps2pdf(filenameps,filenamepdf,[])