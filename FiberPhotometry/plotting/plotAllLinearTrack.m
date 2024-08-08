%plot NE response around track on/off
%plot NE response around reward


fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);

kp = cellfun(@any,regexp(fils,'Linear')) & cellfun(@any,regexp(fils,'NE2'));

fils = fils(kp);




[dirs] = fileparts(fils);
% [dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);
%%
subj =[];
for i = 1:length(dirs)
    [~,subjt] = fileparts(dirs{i});
    subj_split = split(subjt, "-");
    subjt = subj_split{1};
    subj = [subj;{subjt}];
end
usubjs = unique(subj);
%%
ix=  1;
clear home_resp track_resp rew_resp sesID subjID acc_avg track_resp_pred home_resp_pred
for i = 1:length(dirs)   
    
    cd(dirs{i})
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        if size(sessiondata.contextEntry,1) ==3
            
            fs = sessiondata.neural.fs_neural;
            
            k  = gaussian2Dfilter([fs*10 1],fs/2);
            signal_DFoF = nanconvn(sessiondata.neural.signal_DFoF, k');
            
            
            ts_home  = sessiondata.contextEntry{3,2};
            ts_lt  = sessiondata.contextEntry{2,2};
            ts_rew = sort([sessiondata.behavior.rewardTimeLeft;sessiondata.behavior.rewardTimeRight]);
            
            [ix_hc,early1,late1,ts_PETH1] = sm_getIndicesAroundEvent(ts_home,250,500,fs,length(signal_DFoF));
            [ix_lt,early2,late2,ts_PETH2] = sm_getIndicesAroundEvent(ts_lt,250,250,fs,length(signal_DFoF));
            [ix_rew,early,late,ts_PETH] = sm_getIndicesAroundEvent(ts_rew,10,10,fs,length(signal_DFoF));
          
            
           
            if (~( any(late1) | any(late2) | any(early1)| any(early2))) && isfield(sessiondata.neural,'pred') && (ts_home- ts_lt) > 500
                ts_neural = (1:length(signal_DFoF))/fs;
                pred1 = sessiondata.neural.pred.yhat;
                pred1 = interp1(sessiondata.behavior.ts_video,pred1,ts_neural);
                kp = sessiondata.behavior.ts_video>ts_lt+3 &sessiondata.behavior.ts_video<ts_home-3;
                
                signal_DFoF_b = avghist(ts_neural,double(signal_DFoF),sessiondata.behavior.ts_video);
                ttt = sessiondata.behavior.ts_video(kp)+5;
              
                % acc_avg(ix,:) = avghist(sessiondata.behavior.acc(kp,:),signal_DFoF_b(kp),0:2:60);
                
                [~,~,~,b] = histcn([sessiondata.behavior.acc(kp) ttt],0:1:30,ttt(1):10:ttt(1)+600);
                kp = all(b>0,2);
                acc_avg(:,:,ix) = accumarray(b(kp,:),signal_DFoF_b(kp),[31 61],@nanmean,nan);
                home_resp(ix,:) = signal_DFoF(ix_hc);
                track_resp(ix,:) = signal_DFoF(ix_lt);
                track_resp_pred(ix,:) = pred1(ix_lt);
                home_resp_pred(ix,:) = pred1(ix_hc);
                rew_resp(ix,:) = nanmean(signal_DFoF(ix_rew));
            
            sesID(ix) = i;
            subjID(ix) = subj(i);
            ix = ix+1
            end
        end
        
    end
end

%%
close all
figure
 plotMeanSEM(ts_PETH2(1:100:end),(track_resp(:,1:100:end)),'k')
hold on
plotMeanSEM(ts_PETH1(1:100:end)+700,(home_resp(:,1:100:end)),'k')



 plotMeanSEM(ts_PETH2(1:100:end),(track_resp_pred(:,1:100:end)),'r')
hold on
plotMeanSEM(ts_PETH1(1:100:end)+700,(home_resp_pred(:,1:100:end)),'r')
plot([-300 -300],[.5 1],'k')
plot([-300 -100],[.5 .5],'k')

plot([0 0],[-1 2],'k')
plot([-400 1000],[0 0],'k')
set(gca,'fontsize',16)

%%
figure
plotMeanSEM(ts_PETH(1:100:end),(rew_resp(:,1:100:end)),'k')
ylim([-1 2])

%%
figure
col = flipud(linspecer(10,'jet'));
for i = 2:11
   idx = ((i-1)*5+1):((i-1)*5+5) ;
    plotMeanSEM(0:1:30,squeeze(nanmean(acc_avg(:,idx,:),2))',col{i-1})
    hold on
    
    
end
%%
% from plotAllNoverCon
[ix_hc1,early1,late1,ts_PETH1] = sm_getIndicesAroundEvent(ts_home,250,500,fs,length(signal_DFoF));

[ix_hc2,early1,late1,ts_PETH2] = sm_getIndicesAroundEvent(ts_home,300,500,fs,length(signal_DFoF));

figure 
hold on

plotMeanSEM(ts_PETH2(1:100:end),(home_resp_novel(:,1:100:end)),'r')
plotMeanSEM(ts_PETH1(1:100:end),(home_resp(:,1:100:end)),'k')

xlim([-250 500])

plot([-250 500],[0 0],'k')
plot([0 0],[-1 2],'k')
plot([-200 -200],[1 1.5],'k')

plot([-200 -100],[1 1],'k')

%%

% example
load('R:\DANEHippocampalResponse\NE2h4\LinearTrack\NE2h4-220615-143617\sessiondata.mat')
%%
ts = (1:length(sessiondata.neural.signal_DFoF))/sessiondata.neural.fs_neural;
rew = sort([sessiondata.behavior.rewardTimeLeft;sessiondata.behavior.rewardTimeRight]);
close all
plot(ts,nanconvn(sessiondata.neural.signal_DFoF,k'))
hold on
plot([rew],3,'.','color','k')
xlim([500 1700])
ylim([-1 4])
figure
plot(sessiondata.behavior.ts_video,sessiondata.behavior.acc_cor,'r')

xlim([500 1700])


