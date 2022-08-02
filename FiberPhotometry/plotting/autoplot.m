function autoplot(experiment)

topDir = uigetdir();
files = dir(topDir);
dirFlags = [files.isdir];
subdirs = files(dirFlags);
dirs = {subdirs(3:end).name};


%%


switch experiment
    
    case 'drug'
        labels = {'newContext','homeCage', 'newContext_drug'};
        
    case 'NEE'
        labels = {'newContext','homeCage'};
end

%%
newCon = [];
homeCage =[];
checkStability = false;
k = gaussian2Dfilter([1000 1],[100 1]);
saveAllDays = true;
k  = gaussian2Dfilter([1000 1],[200 1]);
newCon1 =[];
homeCage1= [];
newConDrug1 =[];
for i = 1:length(dirs)
    currentDir = [topDir filesep dirs{i} filesep];
    cd(currentDir)
    
    if checkStability
        [signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(pwd,'photoBleachCorrection','exp2');
        plot(ts_data,nanconvn(signal_DFoF,k'))
        ylim([-5 5])
        waitforbuttonpress
    end
    
    if exist('newContext.mat') && strcmp(experiment,'drug')
        [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(currentDir,'newContext.mat',labels,'photoBleachCorrection','exp2','plotIntervals',[300 300],'returnedDataType','corrected');
        % average all new context and home cage
        newCon(i,:) = nanmean(signal_DFoF(ix{1}),1); %first element of ix is the first event type which is newContext
        homeCage(i,:) = nanmean(signal_DFoF(ix{2}),1);
        newConDrug(i,:) = nanmean(signal_DFoF(ix{3}),1);
        
        
        newCon1(i,:) = nanconvn(newCon(i,:),k');
        homeCage1(i,:) = nanconvn(homeCage(i,:),k');
        newConDrug1(i,:) = nanconvn(newConDrug(i,:),k');
    elseif exist('newContext.mat') && strcmp(experiment,'NEE')
        [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(currentDir,'newContext.mat',labels,'photoBleachCorrection','exp2','plotIntervals',[300 300],'returnedDataType','corrected');
        % average all new context and home cage
        newCon(i,:) = nanmean(signal_DFoF(ix{1}),1); %first element of ix is the first event type which is newContext
        homeCage(i,:) = nanmean(signal_DFoF(ix{2}),1);
        
        newCon1(i,:) = nanconvn(newCon(i,:),k');
        homeCage1(i,:) = nanconvn(homeCage(i,:),k');
        
    end
    
    cd(topDir)
    close all
end


%%
if strcmp(experiment,'drug')
    close all
    figure
    ax = tight_subplot(1,2);
    
  
    axes(ax(1))
    imagesc(ts_PETH,[],newCon1,[-1.5 2])
    title('Home cage')
    xlabel('Time from context entry')%%
    xlim([-300 300])
    title('New context')
    
    axes(ax(2))
    imagesc(ts_PETH,[],newConDrug1,[-1.5 2])
    title('Home cage')
    xlabel('Time from context entry')%%
    xlim([-300 300])
    title('New context drug')
    %caxis([-1.5 2.5])
    
     figure
    plot(ts_PETH,nanmean(newCon1,1),'b')
    hold on
    plot(ts_PETH,nanmean(newConDrug1,1), 'k')
    xlim([-300 300])
    ylim([-0.5 2.5])
    
    xlabel ('Time (s)')
    ylabel ('Z-Score')
    legend ('Control','Drug')
    set(gca,'box','off')
    
    
    
elseif strcmp(experiment,'NEE')
    close all
    figure
    ax = tight_subplot(1,3);
    
    axes(ax(1))
    imagesc(ts_PETH,[],homeCage1,[-1.5 2])
    title('Novel Context')
    ylabel('Day')
    xlabel('Time From Context Entry (s)')
    xlim([-300 300])
    
    title('Home cage')
    
    axes(ax(2))
    imagesc(ts_PETH,[],newCon1,[-1.5 2])
    title('Home cage')
    xlabel('Time from context entry')%%
    xlim([-300 300])
    title('New context')
    
    figure
    plot(ts_PETH,nanmean(newCon1,1),'b')
    hold on
    plot(ts_PETH,nanmean(homeCage1,1), 'k')
    xlim([-300 300])
    ylim([-0.5 2.5])
    
    xlabel ('Time (s)')
    ylabel ('Z-Score')
    legend ('Novel','Home')
    set(gca,'box','off')
    
    
end
%%
%close all
