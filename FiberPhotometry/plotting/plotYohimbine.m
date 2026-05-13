 dirs =  [...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_251022'};... %yoh
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_251022'};... %yoh
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_251022'};... %yoh
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_251022'};... %yoh
    ];
%%



for i = 1:length(dirs)
    
    cd(dirs{i})
    
    %find TDT directory
    dirN = fileparts(getAllExtFiles(dirs{i},'Tbk',1));
    if iscell(dirN)
        dirN = dirN{1};
    end
    %find context annotations
    fils = getAllExtFiles(dirs{i},'mat',1);
    clear data
    for j = 1:length(fils)
        
        load(fils{j},'data')
        
        if exist('data') && ~any(contains(data(:,1),'yohimbine'))
            ses(i).context_transitions = data;
        elseif exist('data')
             ses(i).yohimbine_injection = data{1,2};
        end
    end
    [signal_DFoF2,ts_data2,fs,data] = sm_getSignal_DFoF(dirN,'baseline',[10 600]);
    
    % set up smoothing kernel
    
    %get all TDT info
   
    x465A = data.streams.x465A.data;
    x405A = data.streams.x405A.data;
    %ts_data2 = (1:length(signal_DFoF2))/fs;
    
    k = gaussian2Dfilter([10*fs 1],fs);
    if isfield(data.epocs,'Wi2_')
    stim_times = data.epocs.Wi2_.onset;
    elseif isfield(data.epocs,'PC1_')
       stim_times =  data.epocs.PC1_.onset;
    else
        stim_times = [];
    end
    ses(i).stim_times = stim_times;
    %get indices around stim
    [ix,early,late,ts] = sm_getIndicesAroundEvent(stim_times,20,50,data.streams.x405A.fs,length(signal_DFoF2));
    
    
    kp = ~(early | late);
    ses(i).stim_times = ses(i).stim_times(kp);
    ix = ix(kp,:);
    ses(i).FP_smooth = nanconvn(signal_DFoF2,k');
    ses(i).FP_smooth465 = nanconvn(x465A,k');
    ses(i).FP_smooth405 = nanconvn(x405A,k');
    
    ses(i).PETH_FP = ses(i).FP_smooth (ix);
    ses(i).PETH_FP465 = ses(i).FP_smooth465 (ix);
    ses(i).PETH_FP405 = ses(i).FP_smooth405 (ix);
    
    
    [~,ix] = bestmatch(0,ts);
    ses(i).fs = fs;
    baseline = mean(ses(i).PETH_FP465(:,1:ix),2);
    ses(i).PETH_FP_norm465 = ses(i).PETH_FP465 - repmat(baseline,1,size(ses(i).PETH_FP465,2));
    
    
       baseline = mean(ses(i).PETH_FP405(:,1:ix),2);
    ses(i).PETH_FP_norm405 = ses(i).PETH_FP405 - repmat(baseline,1,size(ses(i).PETH_FP405,2));
    
    
     baseline = mean(ses(i).PETH_FP(:,1:ix),2);
    ses(i).PETH_FP_norm = ses(i).PETH_FP - repmat(baseline,1,size(ses(i).PETH_FP,2));
    
    
    
    ses(i).directory = dirN;
    i
end

%%

PETH_FP_norm = cell2mat({ses.PETH_FP_norm}');


%%
clear HCStims NCStims
for i = 1:length(ses)
    try
    maxT = length(ses(i).FP_smooth)/fs;
    
      %get yoh epoch
      yoh = [ses(i).yohimbine_injection maxT];
      baseline = [0 ses(i).yohimbine_injection];
      
    %get HC epoch
    HC = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'HC'),2));
    NCA = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'NCA'),2));
    NCB = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'NCB'),2));
    other_con = [NCA NCB];
    clear HC_ep
    for j = 1:length(HC)
        if j<length(HC)
            HC_ep(j,:) = [HC(j) other_con(find(other_con>HC(j),1,'first'))];
        else
            HC_ep(j,:) = [HC(j) maxT];
        end
        
    end
    
    other_con_ep = excludeEpochs([HC_ep(1,1) maxT],HC_ep);
    
    % find stims that are in each block
    kp = InIntervals(ses(i).stim_times,HC_ep);
    
    if length(kp)>size(ses(i).PETH_FP_norm,1)
        kp = kp(1:end-1);
    end
    
    kp_yoh = InIntervals(ses(i).stim_times,yoh);
    
    HCStims_baseline{i} =  ses(i).PETH_FP_norm(kp & ~kp_yoh,:);
    HCStims_yohimbine{i} =  ses(i).PETH_FP_norm(kp & kp_yoh,:);
    NCStims{i} =  ses(i).PETH_FP_norm(~kp,:);
    
    
    HCStims_baseline465{i} =  ses(i).PETH_FP_norm465(kp & ~kp_yoh,:);
    HCStims_yohimbine465{i} =  ses(i).PETH_FP_norm465(kp & kp_yoh,:);
    NCStims465{i} =  ses(i).PETH_FP_norm465(~kp,:);
    
    HCStims_baseline405{i} =  ses(i).PETH_FP_norm405(kp & ~kp_yoh,:);
    HCStims_yohimbine405{i} =  ses(i).PETH_FP_norm405(kp & kp_yoh,:);
    NCStims405{i} =  ses(i).PETH_FP_norm405(~kp,:);
    end
end

%%

HCStims_baseline0 = cell2mat(HCStims_baseline465');
HCStims_yohimbine0 = cell2mat(HCStims_yohimbine465');

figure
plotMeanSEM(ts(1:1000:end),HCStims_baseline0(:,1:1000:end),'k')

plotMeanSEM(ts(1:1000:end),HCStims_yohimbine0(:,1:1000:end),'r')

%%
clear conTran

for i = 1:length(ses)
    
    
    kp = find(contains(ses(i).context_transitions(:,1),'NCA'));
    trans = ses(i).context_transitions{kp,2};
     [ix,early,late,ts1] = sm_getIndicesAroundEvent(trans,300,600,ses(i).fs,length(ses(i).FP_smooth));
     
            [ixx,early,late,ts] = sm_getIndicesAroundEvent(ses(i).stim_times,20,50,fs,length(ses(i).FP_smooth));
        
        tmp = ses(i).FP_smooth;
        
        tmp(ixx) = tmp(ixx) - repmat(nanmean(HCStims{i}),size(ixx,1),1);
        
        sig_corrected = tmp;
        
        
    conTran(i,:) = sig_corrected(ix);
end

%%


clear inj

for i = 1:length(ses)
    
    

    trans = ses(i).yohimbine_injection;
     [ix,early,late,ts1] = sm_getIndicesAroundEvent(trans,300,600,ses(i).fs,length(ses(i).FP_smooth));
     
            [ixx,early,late,ts] = sm_getIndicesAroundEvent(ses(i).stim_times,20,50,fs,length(ses(i).FP_smooth));
        
        tmp = ses(i).FP_smooth;
        
        tmp(ixx) = tmp(ixx) - repmat(nanmean(HCStims_baseline{i}),size(ixx,1),1);
        
        sig_corrected = tmp;
        
        
    inj(i,:) = sig_corrected(ix);
end

%%
figure
plotMeanSEM(ts1(1:1000:end),inj(:,1:1000:end),'r')
%%
close all
figure
hold on

ts_FP =   (1:length(ses(2).FP_smooth))/ses(2).fs;
plot([ses(2).yohimbine_injection ses(2).yohimbine_injection],[-3.75 7.5],'k')
plot([ses(2).context_transitions{3,2} ses(2).context_transitions{3,2}],[-3.75 7.5],'k')

plot(ts_FP,ses(2).FP_smooth,'r')
plot(ses(2).stim_times,-3.5,'.','color','r')
ylim([-3.75 7.5])
xlim([10 5378])
%%

[~,ii] = bestmatch(0,ts1);

figure

hold on
plotMeanSEM(ts1(1:1000:end),conTran(:,1:1000:end)-repmat(nanmean(conTran(:,1:ii),2),1,size(conTran(:,1:1000:end),2)),'k')
plotMeanSEM(ts1(1:1000:end),conTranYOH(:,1:1000:end)-repmat(nanmean(conTranYOH(:,1:ii),2),1,size(conTranYOH(:,1:1000:end),2)),'r')
xlim([-300 600])

%%
ix=1;
clear conTranA conTranB s conTranHCA conTranHCB
for i = 1:length(ses)
    %find sessions with both contexts
    
    if any(contains(ses(i).context_transitions(:,1),'NCA')) && any(contains(ses(i).context_transitions(:,1),'NCB'))
        % correct for stim
        [~,s(ix).mouse] = fileparts(fileparts(fileparts(ses(i).directory)));
        [ixx,early,late,ts] = sm_getIndicesAroundEvent(ses(i).stim_times,20,50,fs,length(ses(i).FP_smooth));
        
        tmp = ses(i).FP_smooth;
        
        tmp(ixx) = tmp(ixx) - repmat(nanmean(HCStims{i}),size(ixx,1),1);
        
        sig_corrected = tmp;
        % get transitions times
        NCA = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'NCA'),2));
        NCB = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'NCB'),2));
        
        % getHC
        
        HC = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'HC'),2));
        HCA = HC(find(HC>NCA,1,'first'));
        HCB = HC(find(HC>NCB,1,'first'));
        
        
        [ixA,early,late,ts] = sm_getIndicesAroundEvent(NCA,300,1000,fs,length(ses(i).FP_smooth));
        [ixB,early,late,ts] = sm_getIndicesAroundEvent(NCB,300,1000,fs,length(ses(i).FP_smooth));
        
        tmpA = sig_corrected(ixA);
        conTranA(ix,:) = tmpA;
        
        
        
        tmpB =  sig_corrected(ixB);
        conTranB(ix,:) = tmpB;
        
        % get predicted decay
        ok = ses(i).stim_times-NCB;
        
        kp = ok>-300 & ok<1000;
        ok = ok(kp);
        [~,iix] = bestmatch(ok,ts);
        iix = iix+repmat(round(-60*fs): round(fs*90),length(iix),1);
        kp = all(iix<=length(tmpA) & iix>0,2);
        iix = iix(kp,:);
        
        tmpB = tmpB(iix);
        tmpBc = tmpB - repmat(tmpB(:,1),1,size(tmpB,2));
        tmpA = tmpA(iix);
        tmpAc = tmpA - repmat(tmpA(:,1),1,size(tmpA,2));
        s(ix).stimA = tmpA;
        s(ix).stimA_bl = tmpAc;
        s(ix).stimB = tmpB;
        s(ix).stimB_bl = tmpBc;
        s(ix).stim_ts_relative = ok(kp);
        
        %now for HC
        
        
        
        [ixA,early,late,ts] = sm_getIndicesAroundEvent(HCA,300,500,fs,length(ses(i).FP_smooth));
        [ixB,early,late,ts] = sm_getIndicesAroundEvent(HCB,300,500,fs,length(ses(i).FP_smooth));
        
        tmpA = sig_corrected(ixA);
        conTranHCA(ix,:) = tmpA;
        
        
        
        tmpB =  sig_corrected(ixB);
        conTranHCB(ix,:) = tmpB;
        
       
        
        
        
        if ~isempty(ses(i).stim_times)
            s(ix).stim_ts = ses(i).stim_times;
        else
            s(ix).stim_ts  = nan;
        end
        ix = ix+1;
        
        
    elseif  any(contains(ses(i).context_transitions(:,1),'NCA'))
        
        [~,s(ix).mouse] = fileparts(fileparts(fileparts(ses(i).directory)));

        [ixx,early,late,ts] = sm_getIndicesAroundEvent(ses(i).stim_times,20,50,fs,length(ses(i).FP_smooth));
        
        tmp = ses(i).FP_smooth;
        
        tmp(ixx) = tmp(ixx) - repmat(nanmean(HCStims{i}),size(ixx,1),1);
        
        sig_corrected = tmp;
        % get transitions times
        NCB = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'NCA'),2));
         % getHC
        
        HC = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'HC'),2));
      
        HCB = HC(find(HC>NCB,1,'first'));
        [ixB,early,late,ts] = sm_getIndicesAroundEvent(NCB,300,1000,fs,length(ses(i).FP_smooth));
        
        tmpB =  sig_corrected(ixB);
        conTranB(ix,:) = tmpB;
        conTranA(ix,:) = nan(size(tmpB));
        % get predicted decay
        ok = ses(i).stim_times-NCB;
        
        kp = ok>-300 & ok<1000;
        ok = ok(kp);
        [~,iix] = bestmatch(ok,ts);
        iix = iix+repmat(round(-60*fs): round(fs*90),length(iix),1);
        kp = all(iix<=length(tmpB) & iix>0,2);
        iix = iix(kp,:);
        
        tmpB = tmpB(iix);
        tmpBc = tmpB - repmat(tmpB(:,1),1,size(tmpB,2));
      
        s(ix).stimB = tmpB;
        s(ix).stimB_bl = tmpBc;
        s(ix).stim_ts_relative = ok(kp);
        
        
         %now for HC
        
        
        
     
        [ixB,early,late,ts_HC] = sm_getIndicesAroundEvent(HCB,300,500,fs,length(ses(i).FP_smooth));
        
       
        
        tmpB =  sig_corrected(ixB);
        conTranHCB(ix,:) = tmpB;
        
       
        
        
        if ~isempty(ses(i).stim_times)
            s(ix).stim_ts = ses(i).stim_times;
        else
            s(ix).stim_ts  = nan;
        end
        ix = ix+1;
    end
    
    
    
end



%%
stimB_bl = cell2mat({s.stimB_bl}');
stimA_bl = cell2mat({s.stimA_bl}');

stimB = cell2mat({s.stimB}');
stimA = cell2mat({s.stimA}');

stim_ts = cell2mat({s.stim_ts}');

%%

close all
figure
ax  = tight_subplot(1,7);
for i = 1:7
    axes(ax(i))
    kp = stim_ts>(i-1)*100 & stim_ts<i*100;

   plotMeanSEM(-60:1/(fs/10):90,(stimA(kp,1:10:end)),'k')
hold on
plotMeanSEM(-60:1/(fs/10):90,(stimB(kp,1:10:end)),'r')
ylim([-3 4])
if i>1
    xlim([-10 80])
else
    ylabel('Signal_N_E')
end
title([num2str((i-1)*100) '-' num2str((i)*100) 's'])
xlabel('time from stim. (s)')
end
%%


stimB_bl = cell2mat({sy.stimB_bl}');


stimB = cell2mat({sy.stimB}');  


stim_ts = cell2mat({sy.stim_ts}');

%%

%close all
%figure
%ax  = tight_subplot(1,5);
for i = 1:7
    axes(ax(i))
    kp = stim_ts>(i-1)*100 & stim_ts<i*100;

   %plotMeanSEM(-60:1/(fs/10):90,(stimA(kp,1:10:end)),'k')
hold on
plotMeanSEM(-60:1/(fs/10):90,(stimB(kp,1:10:end)),'b')
ylim([-3 4])
if i>1
    xlim([-10 80])
else
    ylabel('Signal_N_E')
end
title([num2str((i-1)*100) '-' num2str((i)*100) 's'])
xlabel('time from stim. (s)')
end

%%


close all
figure
plotMeanSEM((round(-60*fs):round(fs*90))/fs,(stimB_bl-stimA_bl),'k')


%%

stim1 = arrayfun(@(a) min(a.stim_ts),s)'<1000;
stim2 = arrayfun(@(a) min(a.stim_ts),s)'>1000;
nostim = isnan(arrayfun(@(a) min(a.stim_ts),s)');

baselineA = repmat(nanmedian(conTranA(:,1:200000),2),1,size(conTranB(:,1:100:end),2));
baselineB = repmat(nanmedian(conTranB(:,1:200000),2),1,size(conTranB(:,1:100:end),2));

%%
baselineA = zeros(size(baselineA));
baselineB = zeros(size(baselineB));
 close all
 figure
plotMeanSEM(ts(1:100:end),conTranA(stim1,1:100:end)-baselineA(stim1,:),'b') % no stim 2nd
hold on
plotMeanSEM(ts(1:100:end),conTranB(stim2,1:100:end)-baselineB(stim2,:),'k') % no stim 2nd
plotMeanSEM(ts(1:100:end),conTranB(nostim,1:100:end)-baselineB(nostim,:),'r')% stim 2nd

figure
plotMeanSEM(ts(1:100:end),conTranA(nostim | stim2,1:100:end)-baselineA(nostim | stim2,:),'r') % no stim 1st
hold on
plotMeanSEM(ts(1:100:end),conTranB(stim1,1:100:end)-baselineB(stim1,:),'k') % stim 1st


%%
[~,~,mid] = unique({s.mouse}');
stims = cell2mat({s.stim_ts_relative}');
stims1 = arrayfun(@(a) a.stim_ts_relative(find(a.stim_ts_relative>0,1,'first')),s,'UniformOutput',false)'; 


baselineA = zeros(size(baselineA));
baselineB = zeros(size(baselineB));
close all

allStim = conTranB((stim1 | stim2),1:100:end);
allNOStim = [conTranB(nostim,1:100:end);conTranA(:,1:100:end)];

A = sum(allStim(:,3100:8300),2);
B = sum(allNOStim(:,3100:8300),2);

mid_stim = mid(stim1 | stim2);
mid_nostim = [mid(nostim );mid];



figure
plotMeanSEM(ts(1:100:end),allStim,'r') % no stim 1st
hold on
plotMeanSEM(ts(1:100:end),allNOStim,'k') % stim 1st


for i = 1:size(allStim,2)
[~,hh(i)] =  ttest2(allStim(:,i),allNOStim(:,i));
end

hh(hh>.05) = nan;
hh(~isnan(hh)) = 5;
plot(stims,4,'.','color','r')
plot(ts(1:100:end),hh,'k','linewidth',6)
xlim([-300 1000])
ylim([-1 5.5])

%%

s_mouse = accumarray(mid_stim,double(A),[],@nanmean,nan);
ns_mouse = accumarray(mid_nostim,double(B),[],@nanmean,nan);
s_mouse_S = accumarray(mid_stim,double(A),[],@SEM,nan);
ns_mouse_S = accumarray(mid_nostim,double(B),[],@SEM,nan);


%%


baselineA = zeros(size(baselineA));
baselineB = zeros(size(baselineB));
close all

allStim = conTranHCB(stim1 | stim2,1:100:end);
allNOStim = [conTranHCB(nostim,1:100:end);conTranHCA(:,1:100:end)];


figure
plotMeanSEM(ts_HC(1:100:end),allStim,'k') % no stim 1st
hold on
plotMeanSEM(ts_HC(1:100:end),allNOStim,'r') % stim 1st
clear hh

for i = 1:size(allStim,2)
[~,hh(i)] =  ttest2(allStim(:,i),allNOStim(:,i));
end

hh(hh>.05) = nan;
hh(~isnan(hh)) = 5;
plot(cell2mat({s.stim_ts_relative}'),4,'.','color','k')
plot(ts_HC(1:100:end),hh,'k','linewidth',6)
xlim([-300 1000])
ylim([-1 5.5])
%%
 close all
figure
plotMeanSEM(ts(1:100:end),conTranA(end-7:end,1:100:end),'r')
hold on
plotMeanSEM(ts(1:100:end),conTranB(end-7:end,1:100:end),'k')
xlim([-300 1000])
ylim([-1 4])

%%
figure
plotMeanSEM(ts(1:100:end),conTranB_yoh(:,1:100:end),'b')
xlim([-300 550])
ylim([-1 4])
%%
cd('R:\DANEHippocampalResponse\NE2h23\Pos_Control_Homecage_23\NE2h23_250811_152342')
[ups,dwns]  = sm_getDigitalin(pwd,'digitalin.dat',30000,16);
[spikes] = sm_LoadPhy('basepath',pwd,'forceReload',true,'getWaveforms',false,'basename','amplifier_analogin_auxiliary_int16');

%%
[binnedPop_post,bin_times]=populationMatrix(spikes,10,10,1000, ups{1});

%%

plot([data.epocs.Wi2_.onset data.epocs.Wi2_.onset],[-1 3])

%%
figure

plotMeanSEM(ts,PETH_FP_norm(kp1,:),'k')
hold on
%plot(bin_times,nanmean(nanmean(binnedPop_post,3)))
xlim([-20 50])

%%
[~,ix1]  = bestmatch(0,ts);
[~,ix2]  = bestmatch(20,ts);
pre = nanmean(PETH_FP_norm(:,1:ix1),2);
post = nanmean(PETH_FP_norm(:,ix1:ix2),2);


%%
figure
plot([0 1],[pre post],'k')

set(gca,'xtick',[0 1],'xticklabel',{'pre','post'})
ylabel('GRABNE (baseline norm)')