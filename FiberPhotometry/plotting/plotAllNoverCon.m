%plot NE response around track on/off
%plot NE response around reward


fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);

kp = contains(fils,'Novel') & contains(fils,'NE2') & ~contains(fils,'Drug');

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

home_resp_pred_novel =[];
track_resp_pred =[];


home_resp_novel =[];
track_resp =[];
home_resp_predred_novel =[];
track_resp_predred =[];
for i = 1:length(dirs)
    
    cd(dirs{i})
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        if   isfield(sessiondata.neural,'pred')%&& isfield(sessiondata.neural.pred,'yhat_red')  
            
            fs = sessiondata.neural.fs_neural;
            
            k  = gaussian2Dfilter([fs*10 1],fs);
            
            signal_DFoF = nanconvn(sessiondata.neural.signal_DFoF, k');
            ts_neural = (1:length(signal_DFoF))/fs;
            pred1 = sessiondata.neural.pred.yhat;
            pred1 = interp1(sessiondata.behavior.ts_video,pred1,ts_neural);
             pred2 = sessiondata.neural.pred.yhat_red;
            pred2 = interp1(sessiondata.behavior.ts_video,pred2,ts_neural);
            
            ts_NO  = cell2mat(sessiondata.contextEntry(2:2:end,2));
            ts_home  = cell2mat(sessiondata.contextEntry(3:2:end,2));
            
            
            
            [ix_hc1,early1,late1,ts_PETH1] = sm_getIndicesAroundEvent(ts_home,300,500,fs,length(signal_DFoF));
            [ix_lt1,early2,late2,ts_PETH2] = sm_getIndicesAroundEvent(ts_NO,300,300,fs,length(signal_DFoF));
            
            
            
            
            
            home_resp_novel = [home_resp_novel;signal_DFoF(ix_hc1)];
            track_resp = [track_resp;signal_DFoF(ix_lt1)];
            
            
            home_resp_pred_novel = [home_resp_pred_novel;pred1(ix_hc1)];
            track_resp_pred = [track_resp_pred;pred1(ix_lt1)];
          
            home_resp_predred_novel = [home_resp_predred_novel;pred2(ix_hc1)];
            track_resp_predred = [track_resp_predred;pred2(ix_lt1)];
          
          i  
        end
        
    end
end

%%
close all
figure
hold on
plotMeanSEM(ts_PETH2(1:100:end),[track_resp(:,1:100:end)],'k')
plotMeanSEM(ts_PETH1(1:100:end)+600,[home_resp_novel(:,1:100:end)],'k')

plotMeanSEM(ts_PETH2(1:100:end),[track_resp_pred(:,1:100:end)],'r')
plotMeanSEM(ts_PETH1(1:100:end)+600,[home_resp_pred_novel(:,1:100:end)],'r')

plotMeanSEM(ts_PETH2(1:100:end),[track_resp_predred(:,1:100:end)],'b')
%plotMeanSEM(ts_PETH1(1:100:end)+600,[home_resp_predred_novel(:,1:100:end)],'b')


plot([-300 1000],[0 0],'k')
plot([0 0;600 600 ]',[-1.5 3],'k')

plot([-300 -300],[.5 1],'k')
plot([-300 -100],[.5 .5],'k')

set(gca,'fontsize',16)