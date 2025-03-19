function sm_PredictHarvard_getAllFeatures_CARC(idx)

% this is the top level function for calculating feature space for all
% sessions. For computational efficacy, features are only calculated for a
% subset of times, a fixed percentage per seizure for each time bin


%%
%path where raw data is stored with seizure labels

sz_idx = 1;

%%
basefil = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/features/raw_sessions3.mat';
outdir = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/models/features/';
if ~exist(basefil)
    
    
    ops.RawDataPath = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/bids';
    
    
    ops.bins = [ inf 3600 100 10];
    
    % define % of time to take within each time bin
    %   1-3hrs:     5%
    %   10min-1hr:  5%
    %   10s-10min:  20%
    %   0-10s:      100%
    %   seizure:    100%
    %   post ictal: 20%
    
    ops.pct = [.01 .05 .2 .5 .5 .2];
    
    ops.nBins = length(ops.pct);
    
    %define postIctal time, 600s after seizure offset
    ops.timPost = 600;
    
    % info for full feature space
    
    % frequency bands for coherence. must match sm_getPowerPerChannel
    ops.freqs = logspace(log10(.5),log10(99),6);
    
    
    
    ops.durFeat = 2; % 2s feature bins
    
    
    
    
    
    %channels={'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};
    ops.nCh_raw = 20; % N = 20 eeg channels
    
    ops.features = @sm_PredictHarvard_calcFeatures;
    ops.score_thres = .9;
    ops.Fs_data = 200;
    ops.Fs_score = .5;
    %%
    clear sessions
    
    fils = getAllExtFiles(ops.RawDataPath,'mat',1);
    
    fils = fils(contains(fils,'eeg_score'));
    
    sessions(:,1) = fils;
    
    rawData = strrep(fils,'_score','');
    sessions(:,2) = rawData;
    kp = all(cellfun(@exist,sessions),2);
    sessions = sessions(kp,:);
    nSessions = size(sessions,1);
    
    %%
    % loop over all seizures
    clear tims sz_conf
    for i= 1:nSessions
        
        
        v=load(sessions{i,1});
        
        
        %  tims{i} = [NxM] , N = seizure number, M = timestamp for bin of interest length(bins)+2
        tims{i} = [];
        sz_conf{i} =[];
        seizure_start = ((find(diff([0;v.pred_plus]==sz_idx)>0))/ops.Fs_score );
        seizure_end = ((find(diff([0;v.pred_plus]==sz_idx)<0))/ops.Fs_score );
        
        
        if ~isempty(seizure_start)
            if length(seizure_end)<length(seizure_start)
                % ends high, so clip at last point
                seizure_end = [seizure_end;length(v.Y)/ops.Fs_score];
                
            end
            
            if length(seizure_start)<length(seizure_end)
                
                %starts high, start at first point
                seizure_start = [0;seizure_start];
                
            end
            
            
            %loop over seizures per subject
            
            %loop over seizures per subject
            for j = 1:length(seizure_start)
                
                %bins around seizure
                
                if isinf(ops.bins(1))
                    
                    if j ==1 && seizure_start(j) > ops.bins(2)
                        tmp = [0 seizure_start(j)-ops.bins(2:end) seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
                        
                    elseif j ==1 && seizure_start(j) < ops.bins(2)
                        tmp = [nan seizure_start(j)-ops.bins(2:end) seizure_start(j) seizure_end(j) seizure_end(j)+ops.timPost];
                    elseif  j >1 && seizure_start(j)-ops.bins(2) > seizure_end(j-1) +ops.timPost
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
                
                
                tims{i} = [tims{i};tmp];
                iix = round(seizure_start(j))*ops.Fs_score:round(seizure_end(j))*ops.Fs_score;
                sz_conf{i} = [ sz_conf{i};nanmean(v.Y(iix,sz_idx+1))];
            end
        else
            
            
            %deal with sessions that don't have any seizures.
            % grab a random time >1 hr from end of sesssion
            if length(v.Y)/ops.Fs_score>3600+ops.durFeat
                tmp = [0 size(v.Y,1)-3600*ops.Fs_score-ops.durFeat nan(1,5)];
                tims{i} = [tims{i};tmp];
                sz_conf{i} = [ sz_conf{i};nan];
            end
            
        end
        if ~isempty(tims{i})
            kp = diff([0;tims{i}(:,length(ops.bins)+1)])>ops.timPost;
            tims{i} = tims{i}(kp,:);
            sz_conf{i} = sz_conf{i}(kp,:);
        end
    end
    
    kp = ~cellfun(@isempty,tims);
    sessions = sessions(kp,:);
    tims = tims(kp);
    sz_conf = sz_conf(kp);
    
    save(basefil,'tims','sessions','basefil','ops','sz_conf')
else
    
    v = load(basefil);
    tims = v.tims;
    sessions = v.sessions;
    ops = v.ops;
    sz_conf = v.sz_conf;
end

%%

% clear feature variable that will be used to aggregate across sesions
clear sz dat sesID







%loop over number of sessions
idx = (idx-1)*16+1:(idx-1)*16+16;
%%
for i = idx
    for k = 1:ops.nBins
        
        dat{k} =[];
        sesID{k} = [];
    end
    
    [~,fname]= fileparts(sessions{i,2});
    
    ops.outname = [outdir fname '_features.mat'];
    
    if ~exist(ops.outname )
        %loop of time bins
        for k = 1:ops.nBins
            
            v=load(sessions{i,1});
            dur = length(v.Y)/ops.Fs_score;
            %choose random subset (pct) of samples within each session to train
            %model
            
            [sz{k},bad] = getXPctTim(tims{i}(:,k:k+1), ops.pct(k),1);
            
            sz_conf1 =  repmat(sz_conf{i}(~bad),1,size(sz{k} ,2));
            
            
            
            
            %loop over all timepoints for each session for each bin
            for ev = 1:numel(sz{k})
                
                
                tim = sz{k}(ev)-ops.durFeat;
                c = sz_conf1(ev);
                %make sure the even does not exceed duration of recording
                try
                    if ~isnan(tim) && (dur-tim)>ops.durFeat && tim<dur
                        
                        %  try
                        %get features (some pre calculated, some calculated on the fly)
                        features = ops.features(sessions{i,2},tim,ops);
                        if ~isempty(features)
                            
                            %save all features for each time bin
                            dat{k} = [dat{k};features];
                            
                            %keep track of which session matches which feature
                            sesID{k} = [sesID{k}; i tim c];
                            
                        end
                        
                    end
                catch
                    disp(fname)
                    disp(['file length: ' num2str(dur) ' time of read:' num2str(tim)])
                    
                    
                end
                
                
            end
            
            
            
        end
        
        %save the full feature space
        
        
        save(ops.outname,'dat','sesID','sessions','ops','i','tims','-v7.3')
        disp(['saved: ' ops.outname])
        
    end
end
end

