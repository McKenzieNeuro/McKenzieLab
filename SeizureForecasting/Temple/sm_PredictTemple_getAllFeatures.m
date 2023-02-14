% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


% MakeAll_getPowerPerChannel must have been called on all raw data files to
% prepare feature files





%path where raw data is stored with seizure labels
ops.RawDataPath = 'G:\data\isip\oneTreeData\';
ops.FeaturePath = {'G:\data\isip\oneTreeData\binaries'};
ops.FeatureFileOutput = 'G:\data\isip\oneTreeData\features1.mat';


% define time windows around to predict (s)

ops.bins = [ 600 300 100 10];


% define % of time to take within each time bin
%   1-3hrs:     5%
%   10min-1hr:  5%
%   10s-10min:  20%
%   0-10s:      100%
%   seizure:    100%
%   post ictal: 20%

ops.pct = [.25 .25 .25 .5 .25 .125];

ops.nBins = length(ops.pct);

%define postIctal time, 600s after seizure offset
ops.timPost = 600;

% info for full feature space

ops.reScalePhase = 1000;  %
ops.reScalePower = 1000;  % 32767/maxPower < 1,000


% frequency bands for coherence. must match sm_getPowerPerChannel
ops.freqs = logspace(log10(.5),log10(200),20);

% frequency selection for phase/amplitude (must match index in getPowerPerChannel)
ops.amp_idx = 15:20;
ops.ph_idx = 1:10;

ops.durFeat = 1; % 5s feature bins

%channel for phase amplitude coupling
%ch1 = ipsi HPC
%ch2 = ipsi cortex
%ch3 = contra HPC
%ch4 = contral cortex

ops.ch_phaseAmp =2;

%information about feature file (getPowerPerChannel)
ops.nCh_featureFile = 41; %20 power, 20 phase, 1 time series
ops.nCh_raw = 30; % N = 30 eeg channels

ops.features = @sm_PredictTemple_calcFeatures;




%DEPENDENCIES
% getAllExtFiles, sm_PredictIHKA_calcFeatures,sm_getPowerPerChannel,sm_MakeAll_getPowerPerChannel,
% getXPctTim, circ_corrcl

%%




% find all edf files and txt files (with seizure labels)
fils_edf = getAllExtFiles(ops.RawDataPath,'edf',1);
fils_txt = getAllExtFiles('G:\data\isip\data\tuh_eeg_seizure\v1.5.4\edf','_bi',1);

% find all edf files with annotations
[~,b_edf] = cellfun(@fileparts,fils_edf,'uni',0);
[~,b_txt] = cellfun(@fileparts,fils_txt,'uni',0);
goodFils = intersect(b_txt,b_edf);
seizure_fils = fils_txt(ismember(b_txt,goodFils));
edf_fils = fils_edf(ismember(b_edf,goodFils));




%%

sessions = [];

%find feature files that match annotation file
for j = 1:length(ops.FeaturePath)
d = dir(ops.FeaturePath{j});
d = {d.name}';
d1 = cellfun(@(a) a(1:end-11),d,'uni',0);
[a,seizure_fils_csv] = fileparts(seizure_fils);
[~,b] = ismember(seizure_fils_csv,d1);
kp = b>0;
b(~kp) =[];
seizure_filst = seizure_fils(kp);


% save list of file names for the paired annotation/feature files

sessions = [sessions;seizure_filst cellfun(@(a) [ops.FeaturePath{j} filesep a],d(b),'uni',0)];



end

nSessions = size(sessions,1);


%%

warning off
%get all relevant timepoints around the seizure start and end

% TS name in file.
TSname1 = 'Seizure starts';
TSname2 = 'ends';


% loop over all seizures
for i= 1:size(sessions,1)
    
    %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
    tims{i} = [];
    
    
    TSdata = readtable(sessions{i,1},'FileType','text');
    TSdata = table2cell(TSdata);
    
    if size(TSdata,2)==2
        
        C = strsplit(TSdata{end,1},',');
        seizure_start = str2double(C{2});
        seizure_end = str2double(C{3});
    else
    seiz = find(cell2mat(regexp(TSdata(:,4),'seiz')));
    seizure_start = cell2mat(TSdata(seiz,2));
    seizure_end =  cell2mat(TSdata(seiz,3));
    
    end
    %loop over seizures per subject
    for j = 1:length(seizure_start)
        
        %bins around seizure
        tmp = [seizure_start(j)-ops.bins seizure_start(j) min(seizure_end(j),seizure_start(j)+20) seizure_end(j)+ops.timPost];
        
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
clear sz


for k = 1:ops.nBins
    
    dat{k} =[];
    sesID{k} = [];
end

%%


%loop over number of sessions
for i = 1:nSessions
    fname = sessions{i,2};
    fname = fname(1:end-11);
    %loop of time bins
    for k = 1:ops.nBins
        
        
        %choose random subset (pct) of samples within each session to train
        %model
        
        sz{k} = getXPctTim(tims{i}(:,k:k+1), ops.pct(k),1);
        
        
        
        
        
        
        %loop over all timepoints for each session for each bin
        for ev = 1:numel(sz{k})
            
            %find duration of the file
            powerFil = [fname '_ch_000.dat'];
            s = dir(powerFil);
            
            % get fs
            txtFil = [fname '_ch_000.txt'];
            fsInfo = readtable(txtFil,'FileType','text');
            ops.Fs = cell2mat(table2cell(fsInfo(1,2)));
            dur = s.bytes/ops.nCh_featureFile/ops.Fs/2;
            
            tim = sz{k}(ev)-ops.durFeat;
            %make sure the even does not exceed duration of recording
            if ~isnan(tim) && (dur-tim)>ops.durFeat && tim<dur
                
              %  try
                %get features (some pre calculated, some calculated on the fly)
                features = ops.features(fname,tim,ops);
                
                
                %save all features for each time bin
                dat{k} = [dat{k};features];
                
                %keep track of which session matches which feature
                sesID{k} = [sesID{k}; i tim];
               % catch
               %     disp(fname)
               %     disp(['file length: ' num2str(dur) ' time of read:' num2str(tim)])
                    
                    
               % end
                
            end
            
           
        end
        
        
        
    end
    
     %save the full feature space
        save(ops.FeatureFileOutput,'dat','sesID','sessions','ops','-v7.3')
        i
end



