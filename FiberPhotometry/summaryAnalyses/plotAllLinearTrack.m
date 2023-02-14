masterDir = 'R:\McKenzieLab\DANEHippocampalResponse';

masterDir = 'Y:\DANEHippocampalResponse';
masterDir = 'R:\DANEHippocampalResponse';
fils = getAllExtFiles(masterDir,'mat',1);

kp = cellfun(@any,regexpi(fils,'Linear'));

fils = fils(kp);

[dirs] = cellfun(@fileparts,fils,'uni',0);

dirs =  unique(dirs);
%%
trackThres = 325; % position in pixels
k  = gaussian2Dfilter([1 100000],100);
% first detect if mouse is on track
% makee this struct sessiondata.contextEntry
for i = 1:length(dirs)
    i
    cd(dirs{i})
    
    if exist('sessiondata.mat')
        load('sessiondata.mat');
        %  if ~isfield(sessiondata,'contextEntry')
        % if ycoord < 325 then mouse is on track
        
        nSamples = length(sessiondata.behavior.position.left_ear);
        sessiondata.behavior.position.context = cell(nSamples,1);
        YCoord = nanmean([sessiondata.behavior.position.left_ear(:,2) sessiondata.behavior.position.right_ear(:,2)],2);
        ts = sessiondata.behavior.ts_video(:)';
        YCoord = nanconvn(YCoord,k');
        trackON = YCoord<trackThres;
        trackON = ts([find(diff(trackON)>0) find(diff(trackON)<0)]);
        homeCage = excludeEpochs([0 ts(end)],trackON);
        
        epochs = [trackON;homeCage];
        labels = [repmat({'track'},size(trackON,1),1); repmat({'homecage'},size(homeCage,1),1)];
        
        [epochs,b] = sortrows(epochs);
        labels = labels(b);
        data = [labels num2cell(epochs(:,1),2)];
        sessiondata.contextEntry = data;
        
        
        %define entry times
        epochs_on  = cell2mat(sessiondata.contextEntry(:,2));
        
        %define exit times
        epochs_off = [epochs_on(2:end); sessiondata.behavior.ts_video(end)];
        
        %build matrix of onsets and offsets
        epochs = [epochs_on epochs_off];
        
        for j = 1:size(epochs,1)
            
            
            
            %select the subset of times (in the video) between the onset
            %and offset of the jth epoch
            kp = sessiondata.behavior.ts_video >= epochs(j,1) & sessiondata.behavior.ts_video<epochs(j,2);
            sessiondata.behavior.position.context(kp) = sessiondata.contextEntry(j,1);
            
            save('sessiondata.mat','sessiondata','-v7.3');
            
        end
    end
end
%%
k  = gaussian2Dfilter([100 100],[5 5]);
ixx = 1;
clear ratemap
for i = 1:length(dirs)
    i
    cd(dirs{i})
    
    
    
    if exist('sessiondata.mat')
        %load tracking data
        load('sessiondata.mat');
        
        
        
        
        %load context edges
        %   load('contextTransition1.mat');
        
        kp = cellfun(@any,regexp(sessiondata.behavior.position.context,'track'));
        
        [a,LE,c] = pca(sessiondata.behavior.position.left_ear(kp,:));
        
        
        
        ts = sessiondata.behavior.ts_video;
        ts_neural  = (1:length(sessiondata.neural.signal_DFoF))/sessiondata.neural.fs_neural;
        dsNeural = interp1(ts_neural,sessiondata.neural.signal_DFoF,ts);
        
        
        ts = ts(kp);
        dsNeural = double(dsNeural(kp));
        
        bins1 = -260:1:140;
        bins2 = [-20 20];
        [occ,~,~,ix] = histcn(LE,bins1,bins2);
        inbin = all(ix>0,2);
        vel = sessiondata.behavior.vel(kp);
        tmp = accumarray(ix(inbin,:),vel(inbin),[length(bins1),length(bins2)],@nanmean,nan);
        %  tmp(occ<10) = nan;
        ratemap(:,ixx) = nanconvn(tmp(:,1),k,'nanout',false);
        sess{ixx} = dirs{i};
        ixx = ixx+1;
        
        
    end
end

%%

ix=1;
for i = 1:length(dirs)
    i
    cd(dirs{i})
    
    
    
    if exist('sessiondata.mat')
        %load tracking data
        load('sessiondata.mat');
        
        
        
        
        %load context edges
        %   load('contextTransition1.mat');
        
        kp = cellfun(@any,regexp(sessiondata.behavior.position.context,'track'));
        
        %downsample neural
        
        ts_neural = (1:length(sessiondata.neural.signal_DFoF))/sessiondata.neural.fs_neural;
        neural1 = interp1(ts_neural,sessiondata.neural.signal_DFoF,sessiondata.behavior.ts_video);
        
        vel_u(ix,:) = avghist(sessiondata.behavior.vel,double(neural1),0:15:150);
        ix = ix+1;
    end
end

%%

