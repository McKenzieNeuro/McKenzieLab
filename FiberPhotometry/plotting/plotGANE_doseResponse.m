%%

dirs = [...
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_250812'};...
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_250815'};...
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_250812'};...
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_250813'};...
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_250814'};...
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_250815'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_250812'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_250813'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_250814'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_250815'};...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_250812'};...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_250813'};...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_250814'};...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_250815'};...
    ];

%%


for i = 1:length(dirs)
    
    cd(dirs{i})
    
    %find TDT directory
    dirN = fileparts(getAllExtFiles(dirs{i},'Tbk',1));
    
    %find context annotations
    fils = getAllExtFiles(dirs{i},'mat',0);
    clear data
    for j = 1:length(fils)
        
        load(fils{j},'data')
        
        if exist('data')
            ses(i).context_transitions = data;
            break
        end
    end
    [signal_DFoF2,ts_data2,fs] = sm_getSignal_DFoF(dirN);
    
    % set up smoothing kernel
    
    %get all TDT info
    data = TDTbin2mat(dirN);
    
    k = gaussian2Dfilter([10*fs 1],fs);
    stim_times = data.epocs.Wi2_.onset;
    ses(i).stim_times = stim_times;
    %get indices around stim
    [ix,early,late,ts] = sm_getIndicesAroundEvent(stim_times,20,50,data.streams.x405A.fs,length(signal_DFoF2));
    
    
    kp = ~(early | late);
    ses(i).stim_times = ses(i).stim_times(kp);
    ix = ix(kp,:);
    ses(i).FP_smooth = nanconvn(signal_DFoF2,k');
    ses(i).sig465 = data.streams.x465A.data(ix);
    ses(i).sig405 = data.streams.x405A.data(ix);

    ses(i).PETH_FP = ses(i).FP_smooth (ix);
    [~,ix] = bestmatch(0,ts);
    ses(i).fs = fs;
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
    
    HCStims{i} =  ses(i).PETH_FP_norm(kp,:);
    NCStims{i} =  ses(i).PETH_FP_norm(~kp,:);
    end
end

%%
figure
plotMeanSEM(ts,cell2mat(HCStims'),'k')

%%
ix=1;
clear conTranA conTranB s
for i = 1:length(ses)
    %find sessions with both contexts
    
    if any(contains(ses(i).context_transitions(:,1),'NCA')) && any(contains(ses(i).context_transitions(:,1),'NCB'))
        
        
        % get transitions times
        NCA = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'NCA'),2));
        NCB = cell2mat(ses(i).context_transitions(contains(ses(i).context_transitions(:,1),'NCB'),2));
        [ixA,early,late,ts] = sm_getIndicesAroundEvent(NCA,300,600,fs,length(ses(i).FP_smooth));
        [ixB,early,late,ts] = sm_getIndicesAroundEvent(NCB,300,600,fs,length(ses(i).FP_smooth));
        
   tmpA = ses(i).FP_smooth(ixA);
    conTranA(ix,:) = tmpA;
    
    
    
    tmpB =  ses(i).FP_smooth(ixB);
    conTranB(ix,:) = tmpB;
    
    % get predicted decay
  ok = ses(i).stim_times-NCB;
  
  kp = ok>-300 & ok<600;
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
  s(ix).stim_ts = ok(kp);
  ix=ix+1;
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
ax  = tight_subplot(1,5);
for i = 1:5
    axes(ax(i))
    kp = stim_ts>(i-1)*100 & stim_ts<i*100;

   plotMeanSEM(-60:1/(fs/10):90,(stimA(kp,1:10:end)),'k')
hold on
plotMeanSEM(-60:1/(fs/10):90,(stimB(kp,1:10:end)),'r')
ylim([-1 2.5])
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
 close all
plotMeanSEM(ts,(conTranB),'k')
hold on
plotMeanSEM(ts,(conTranA),'r')

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
