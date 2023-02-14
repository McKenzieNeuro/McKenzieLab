% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


% MakeAll_getPowerPerChannel must have been called on all raw data files to
% prepare feature files





%path where raw data is stored with seizure labels
ops.RawDataPath = 'R:\IHKA_Scharfman\IHKA data';
ops.FeaturePath = {'F:\data1\IHKA','E:\data\IHKA'};
ops.FeatureFileOutput = 'E:\data\IHKA\features_dur.mat';


% define time windows around to predict (s)

ops.bins = 10:10:70;


% define % of time to take within each time bin


ops.pct = ones(length(ops.bins)-1,1);

ops.nBins = length(ops.pct);



% info for full feature space
ops.Fs  = 2000; % sampling rate of edf files
ops.reScalePhase = 1000;  %
ops.reScalePower = 1000;  % 32767/maxPower < 1,000


% frequency bands for coherence. must match sm_getPowerPerChannel
ops.freqs = logspace(log10(.5),log10(200),20);

% frequency selection for phase/amplitude (must match index in getPowerPerChannel)
ops.amp_idx = 15:20;
ops.ph_idx = 1:10;

ops.durFeat = 5; % 5s feature bins

%channel for phase amplitude coupling
%ch1 = ipsi HPC
%ch2 = ipsi cortex
%ch3 = contra HPC
%ch4 = contral cortex

ops.ch_phaseAmp =2;

%information about feature file (getPowerPerChannel)
ops.nCh_featureFile = 41; %20 power, 20 phase, 1 time series
ops.nCh_raw = 4; % N = 4 eeg channels

ops.features = @sm_PredictIHKA_calcFeatures__durSeizure;




%DEPENDENCIES
% getAllExtFiles, sm_PredictIHKA_calcFeatures,sm_getPowerPerChannel,sm_MakeAll_getPowerPerChannel,
% getXPctTim, circ_corrcl

%%




% find all edf files and txt files (with seizure labels)
fils_edf = getAllExtFiles(ops.RawDataPath,'edf',1);
fils_txt = getAllExtFiles(ops.RawDataPath,'txt',1);

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
    d = {d(cell2mat({d.isdir})).name}';
    [a,b] = fileparts(seizure_fils);
    seizure_fils_txt = regexprep(b,' ','_');
    [~,b] = ismember(seizure_fils_txt,d);
    kp = b>0;
    b(~kp) =[];
    seizure_filst = seizure_fils(kp);
    
    
    % save list of file names for the paired annotation/feature files
    
    sessions = [sessions;seizure_filst cellfun(@(a) [ops.FeaturePath{j} filesep a filesep a],d(b),'uni',0)];
    
    
    
end

nSessions = size(sessions,1);


%%

warning off
%get all relevant timepoints around the seizure start and end

% TS name in file.
TSname1 = 're starts';
TSname2 = 're ends';


% loop over all seizures
clear tims
for i= 1:size(sessions,1)
    
    %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
    tims{i} = [];
    
    
    TSdata = readtable(sessions{i,1});
    TSdata = table2cell(TSdata);
    
    seizure_start = cell2mat(TSdata(cellfun(@any,regexpi(TSdata(:,6),TSname1)),4));
    seizure_end = cell2mat(TSdata(cellfun(@any,regexpi(TSdata(:,6),TSname2)),4));
    
    if any( (seizure_end-seizure_start)>100)
        disp('ok')
    end
    
    
    tims{i} = [tims{i};seizure_start seizure_end];
    
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
    
    %loop of seizures
    for k = 1:size(tims{i},1)
        
        
        
        try
            features = ops.features(fname,tims{i}(k,:),ops);
            
            
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


%save the full feature space
save(ops.FeatureFileOutput,'dat','sesID','sessions','ops','-v7.3')




