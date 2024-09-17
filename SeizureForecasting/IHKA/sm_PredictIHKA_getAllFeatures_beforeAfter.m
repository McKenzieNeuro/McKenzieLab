% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


% MakeAll_getPowerPerChannel must have been called on all raw data files to
% prepare feature files





%path where raw data is stored with seizure labels
ops.RawDataPath = 'R:\IHKA_Scharfman\IHKA data';
ops.FeaturePath = {'F:\data1\IHKA','E:\data\IHKA'};
ops.FeatureFileOutput = 'E:\data\IHKA\features_before_after.mat';


% define time windows around to predict (s)

ops.bins = [ -inf -3600 -100 -10 10 100 3600 inf];


% define % of time to take within each time bin
%   1-3hrs:     5%
%   10min-1hr:  5%
%   10s-10min:  20%
%   0-10s:      100%
%   seizure:    100%
%   post ictal: 20%

ops.pct = [.05 .5 .1 1 1 1 .5 .05];




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
%ch2 = ipsi cortex
%ch3 = contra HPC
%ch4 = contral cortex

ops.ch_phaseAmp =2;

%information about feature file (getPowerPerChannel)
ops.nCh_featureFile = 41; %20 power, 20 phase, 1 time series
ops.nCh_raw = 4; % N = 4 eeg channels

ops.features = @sm_PredictIHKA_calcFeatures;




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
TSname1 = 'ure starts';
TSname2 = 'ends';

clear sz
% loop over all seizures
for i= 1:size(sessions,1)
     sz(:,:,i) = cell(4,5);
    %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
    tims{i} = [];
    
    
    TSdata = readtable(sessions{i,1});
    TSdata = table2cell(TSdata);
    
    seizure_start = cell2mat(TSdata(cellfun(@any,regexpi(TSdata(:,6),TSname1)),4));
    seizure_end = cell2mat(TSdata(cellfun(@any,regexpi(TSdata(:,6),TSname2)),4));
    
  
    %calc time to seizure start
    
    ts = [0:24*3600]';

    clear sz_start sz_end
    for j = 1:length(seizure_start)
          sz_start{j} = nan(size(ts));
        if j==1
            sz_start{j}(ts<seizure_start(j)) =  ts(ts<seizure_start(j)) - seizure_start(j);
        else
            sz_start{j}(ts<seizure_start(j) & ts>seizure_end(j-1)) =  ts(ts<seizure_start(j) & ts>seizure_end(j-1)) - seizure_start(j);
            
        end
    end
    
    
    %calc time from seizure end
    
 
    
    for j = 1:length(seizure_end)
          sz_end{j} = nan(size(ts));
        if j<length(seizure_end)
            sz_end{j}(ts>seizure_end(j) & ts<seizure_start(j+1)) =  ts(ts>seizure_end(j) & ts<seizure_start(j+1)) - seizure_end(j);
        else
            sz_end{j}(ts>seizure_end(j))  =  ts(ts>seizure_end(j)) - seizure_end(j);
            
        end
    end
    
    
    for j = 2:length(seizure_end)
        [n,~,~,b] = histcn([sz_start{j} sz_end{j-1}],[ops.bins(1:4) 0],[0 ops.bins(5:8)]);
      
        % get 1000 or max of every category
       ub = unique(b,'rows');
       ub = ub(~any(ub==0,2),:);
        for k = 1:size(ub,1)
            
            tst = find(ismember(b,ub(k,:),'rows'));
            
            if length(tst)<1000
                ix = 1:length(tst);
            else
                ix = randsample(1:length(tst),1000);
            end
            tst = tst(ix);
            
            tst = tst;
            tst = tst(tst>0);
            sz{ub(k,1),ub(k,2),i} = [ sz{ub(k,1),ub(k,2),i};tst];
        end
        
    end
    
    for k = 1:length(seizure_start)
        sz{1,5,i} = [ sz{1,5,i} seizure_start(k):seizure_end(k)];
    end
end

    


%%

% clear feature variable that will be used to aggregate across sesions



for k = 1:17
    
    dat{k} =[];
    sesID{k} = [];
end

%%


%loop over number of sessions
for i = 1:nSessions
    fname = sessions{i,2};
    szt = sz(:,:,i);
    %loop of time bins
    for k = 1:17
        
        
        %choose random subset (pct) of samples within each session to train
        %model
        
      
        
        
        
        
        
        
        %loop over all timepoints for each session for each bin
        for ev = 1:numel(szt{k})
            
            %find duration of the file
            powerFil = [fname '_1.dat'];
            s = dir(powerFil);
            dur = s.bytes/ops.nCh_featureFile/ops.Fs/2;
            
            tim = szt{k}(ev)-ops.durFeat;
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



