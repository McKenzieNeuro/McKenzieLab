%%

dirs = [...
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_250812'};... %g
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_250815'};...  %ng
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_250812'};...  %ng
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_250813'};... %g
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_250814'};... %ng
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_250815'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_250812'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_250813'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_250814'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_250815'};...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_250812'};...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_250813'};...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_250814'};...
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_250815'};...
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_251024'};...
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_251024'};...
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_251024'};...
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_251027'};... % stim first
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_251027'};... % stim first
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_251027'};... % stim first
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_251027'};... % stim first
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_251028'};... % stim first
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_251028'};... % stim first
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_251028'};... % stim first
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_251028'};... % stim first
   
    {'R:\DANEHippocampalResponse\NE2h21\NE2h21_251029'};... % no stim
    {'R:\DANEHippocampalResponse\NE2h22\NE2h22_251029'};... % no stim
    {'R:\DANEHippocampalResponse\NE2h23\NE2h23_251029'};... % no stim
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_251029'};... % no stim
    {'R:\DANEHippocampalResponse\NE2h24\NE2h24_251024'};... % no stim
    

   
    
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
    fils = getAllExtFiles(dirs{i},'mat',0);
    clear data
    for j = 1:length(fils)
        
        load(fils{j},'data')
        
        if exist('data')
            ses(i).context_transitions = data;
            break
        end
    end
    [signal_DFoF2,ts_data2,fs,data] = sm_getSignal_DFoF(dirN,'baseline',[10 600]);
    
    % set up smoothing kernel
    
    %get all TDT info
    %data = TDTbin2mat(dirN);
    %fs = data.streams.x465A.fs;
    %signal_DFoF2 = data.streams.x465A.data;
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
plotMeanSEM(ts,cell2mat(HCStims(1:end-2)'),'k')
plotMeanSEM(ts,cell2mat(HCStims(end-1:end)'),'b')
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
        
        
        [ixA,early,late,ts_tran] = sm_getIndicesAroundEvent(NCA,300,1000,fs,length(ses(i).FP_smooth));
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



        s(ix).stimA_corrected = tmpA;
        s(ix).stimA_bl = tmpAc;
        s(ix).stimB_corrected = tmpB;
        s(ix).stimB_bl = tmpBc;
        s(ix).stim_ts_relative = ok(kp);
        


        tmpA =  ses(i).FP_smooth(ixA);
     
        
        
        
        tmpB =   ses(i).FP_smooth(ixB);
      
        
        % get predicted decay
        ok = ses(i).stim_times-NCB;
        
        kp = ok>-300 & ok<1000;
        ok = ok(kp);
        [~,iix] = bestmatch(ok,ts);
        iix = iix+repmat(round(-60*fs): round(fs*90),length(iix),1);
        kp = all(iix<=length(tmpA) & iix>0,2);
        iix = iix(kp,:);
        
        tmpB = tmpB(iix);
       
        tmpA = tmpA(iix);

        s(ix).stimA_uncorrected = tmpA;

        s(ix).stimB_uncorrected = tmpB;

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
      
        s(ix).stimB_corrected = tmpB;
        s(ix).stimB_bl = tmpBc;
        s(ix).stim_ts_relative = ok(kp);
        

      
     
        
        
        
        tmpB =   ses(i).FP_smooth(ixB);
      
        
        % get predicted decay
        ok = ses(i).stim_times-NCB;
        
        kp = ok>-300 & ok<1000;
        ok = ok(kp);
        [~,iix] = bestmatch(ok,ts);
        iix = iix+repmat(round(-60*fs): round(fs*90),length(iix),1);
        kp = all(iix<=length(tmpA) & iix>0,2);
        iix = iix(kp,:);
        
        tmpB = tmpB(iix);
       

        s(ix).stimB_uncorrected = tmpB;

        
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

stim1 = arrayfun(@(a) min(a.stim_ts),s)'<1000;
stim2 = arrayfun(@(a) min(a.stim_ts),s)'>1000;
nostim = isnan(arrayfun(@(a) min(a.stim_ts),s)');

kp = ~arrayfun(@(a) isempty(a.stimB_bl),s) & ~arrayfun(@(a) isempty(a.stimA_bl),s);
stimB_bl = cell2mat({s(kp).stimB_bl}');

stimA_bl = cell2mat({s(kp).stimA_bl}');

stimB_uncorrected = cell2mat({s(kp).stimB_uncorrected}');

stimA_uncorrected = cell2mat({s(kp).stimA_uncorrected}');
stimB_corrected = cell2mat({s(kp).stimB_corrected}');

stimA_corrected = cell2mat({s(kp).stimA_corrected}');
stim_ts = cell2mat({s(kp).stim_ts_relative}');

%%

close all
figure
ax  = tight_subplot(1,7);
for i = 1:7
    axes(ax(i))
    kp = stim_ts>(i-1)*100 & stim_ts<i*100;
 
   plotMeanSEM(-60:1/(fs/10):90,(stimA_corrected(kp,1:10:end)),'k')
hold on
plotMeanSEM(-60:1/(fs/10):90,(stimB_corrected(kp,1:10:end)),'r')
ylim([-3 4])
if i>1
  
else
    ylabel('Signal_N_E')
end
  xlim([-45 45])
title([num2str((i-1)*100) '-' num2str((i)*100) 's'])
xlabel('time from stim. (s)')
end

figure
ax  = tight_subplot(1,7);
for i = 1:7
    axes(ax(i))
    kp = stim_ts>(i-1)*100 & stim_ts<i*100;
 
   plotMeanSEM(-60:1/(fs/10):90,(stimA_uncorrected(kp,1:10:end)),'k')
hold on
plotMeanSEM(-60:1/(fs/10):90,(stimB_uncorrected(kp,1:10:end)),'r')
ylim([-3 4])
if i>1
  
else
    ylabel('Signal_N_E')
end
  xlim([-45 45])
title([num2str((i-1)*100) '-' num2str((i)*100) 's'])
xlabel('time from stim. (s)')
end


%%
% 
% %close all
% %figure
% %ax  = tight_subplot(1,5);
% for i = 1:7
%     axes(ax(i))
%     kp = stim_ts>(i-1)*100 & stim_ts<i*100;
% 
%    %plotMeanSEM(-60:1/(fs/10):90,(stimA(kp,1:10:end)),'k')
% hold on
% plotMeanSEM(-60:1/(fs/10):90,(stimB(kp,1:10:end)),'b')
% ylim([-3 4])
% if i>1
%     xlim([-10 80])
% else
%     ylabel('Signal_N_E')
% end
% title([num2str((i-1)*100) '-' num2str((i)*100) 's'])
% xlabel('time from stim. (s)')
% end

%%

% 
% close all
% figure
% plotMeanSEM((round(-60*fs):round(fs*90))/fs,(stimB_bl-stimA_bl),'k')


%%

stim1 = arrayfun(@(a) min(a.stim_ts),s)'<1000;
stim2 = arrayfun(@(a) min(a.stim_ts),s)'>1000;
nostim = isnan(arrayfun(@(a) min(a.stim_ts),s)');

baselineA = repmat(nanmedian(conTranA(:,1:200000),2),1,size(conTranB(:,1:100:end),2));
baselineB = repmat(nanmedian(conTranB(:,1:200000),2),1,size(conTranB(:,1:100:end),2));


%%
ok = ses(nostim);
clear ConA_double ConB_double
for i = 1:length(ok)

kpA = find(contains(ok(i).context_transitions(:,1),'NCA'));
kpB = find(contains(ok(i).context_transitions(:,1),'NCB'));

conA = ok(i).context_transitions{kpA,2};
conB = ok(i).context_transitions{kpB,2};
[ix] = sm_getIndicesAroundEvent(conA,100,600,ok(i).fs,length(ok(i).FP_smooth));
ConA_double(i,:) = ok(i).FP_smooth(ix);
[ix] = sm_getIndicesAroundEvent(conB,100,600,ok(i).fs,length(ok(i).FP_smooth));

ConB_double(i,:) = ok(i).FP_smooth(ix);

end


%%
%baselineA = zeros(size(baselineA));
%baselineB = zeros(size(baselineB));
 close all
 figure
plotMeanSEM(ts_tran(1:100:end),conTranA(stim1,1:100:end)-baselineA(stim1,:),'b') % no stim 2nd
hold on
plotMeanSEM(ts_tran(1:100:end),conTranB(stim2,1:100:end)-baselineB(stim2,:),'k') % no stim 2nd
plotMeanSEM(ts_tran(1:100:end),conTranB(nostim,1:100:end)-baselineB(nostim,:),'r')% stim 2nd

figure
plotMeanSEM(ts_tran(1:100:end),conTranA(nostim | stim2,1:100:end)-baselineA(nostim | stim2,:),'r') % no stim 1st
hold on
plotMeanSEM(ts_tran(1:100:end),conTranB(stim1,1:100:end)-baselineB(stim1,:),'k') % stim 1st




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
plotMeanSEM(ts_tran(1:100:end),allStim,'r') % no stim 1st
hold on
plotMeanSEM(ts_tran(1:100:end),allNOStim,'k') % stim 1st


for i = 1:size(allStim,2)
[~,hh(i)] =  ttest2(allStim(:,i),allNOStim(:,i));
end

hh(hh>.05) = nan;
hh(~isnan(hh)) = 5;
plot(stims,4,'.','color','r')
plot(ts_tran(1:100:end),hh,'k','linewidth',6)
xlim([-300 1000])
ylim([-1 5.5])

%%

s_mouse = accumarray(mid_stim,double(A),[],@nanmean,nan);
ns_mouse = accumarray(mid_nostim,double(B),[],@nanmean,nan);
s_mouse_S = accumarray(mid_stim,double(A),[],@SEM,nan);
ns_mouse_S = accumarray(mid_nostim,double(B),[],@SEM,nan);
figure
plot([ns_mouse s_mouse]')
set(gca,'xtick',[1 2],'xticklabel',{'Control','CA1 stim'},'fontsize',16)
ylabel('Total Signal_N_E')
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
 