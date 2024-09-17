%clear all
close all

fils1 = getAllExtFiles('R:\McKenzieLab\ASommer\FP experiments DA-NE\Running Wheel\DA2h9\Single Obj No Delay','mat',1);
fils2 =  getAllExtFiles('R:\McKenzieLab\ASommer\FP experiments DA-NE\Running Wheel\DA2h9\Single Obj Long Delay','mat',1);
fils = [fils1;fils2];
kp = contains(fils,'ContextTransitionRevised');
fils = fileparts(fils(kp));

%%
clear wheel obj_nowheel obj_wheel HC_post HC_pre date HC_pre_nowheel HC_pre_wheel HC_post_nowheel HC_post_wheel
iix1 =1;
iix = 1;
date =[];
for i = 1:length(fils)
    cd(fils{i})
    
    if exist('signal_data.mat')
        load('signal_data.mat')
        load('ContextTransitionRevised.mat')
        ix = regexp(fils{i},'-');
        tmp = fils{i}(ix(2)+1:ix(3)-1);
        if length(regexp(tmp,'[0-9]')) == length(tmp)
            tmp = [tmp(5:6) '/' tmp(3:4) '/' '20' tmp(1:2)];
        end
        
        datet = datenum(tmp);
        date = [date;datet];
        fs  = 1000;
        newts = (1:length(analyzed.filtered))/fs;
        filtered = interp1(analyzed.time*60,analyzed.filtered,newts);
        
        inArena = cell2mat(data(contains(data(:,1),'large'),2));
        inArena = inArena(1);
        obj = cell2mat(data(contains(data(:,1),'object'),2));
        obj = obj(obj>inArena);
        
        HC = cell2mat(data(contains(data(:,1),'homecage'),2));
        HC1 = HC(HC<inArena);
        HC2 = HC(HC>inArena);
        [ix_HC1,early,late,ts_PETH_HC1] = sm_getIndicesAroundEvent(HC1,-2,600,fs,length(filtered));
        [ix_HC2,early,late,ts_PETH_HC2] = sm_getIndicesAroundEvent(HC2,30,600,fs,length(filtered));
        [ix,early,late,ts_PETH1] = sm_getIndicesAroundEvent(obj,30,660,fs,length(filtered));
        
        
        kp_wheel = contains(data(:,1),'wheel');
        
        if any(kp_wheel)
            wheelT =   cell2mat(data(kp_wheel,2));
            [ix1,early,late,ts_PETH] = sm_getIndicesAroundEvent(wheelT,30,600,fs,length(filtered));
            obj_wheel(iix,:) = filtered(ix);
            wheel(iix,:) = filtered(ix1);
            HC_pre_wheel(iix,:) = filtered(ix_HC1);
            HC_post_wheel(iix,:) = filtered(ix_HC2);
            
            iix = iix+1;
        else
            obj_nowheel(iix1,:) = filtered(ix);
            HC_pre_nowheel(iix,:) = filtered(ix_HC1);
            HC_post_nowheel(iix,:) = filtered(ix_HC2);
            iix1 = iix1+1
        end
        
    end
    i
end

%%
figure
plotMeanSEM(ts_PETH1,obj_wheel,'r')
hold on
plotMeanSEM(ts_PETH1,obj_nowheel,'k')


%%
% p = [];ix=1;
% for i = 1:10:60000
%     p(ix) = ranksum(obj_nowheel(:,i),obj_wheel(:,i));
%     ix = ix+1;
% end
%%
% p1 =  p;
% p1(p>.05) = nan;
% p1(~isnan(p1)) = 1.5;
% plot(ts_PETH1(1:10:60000),p1,'k','linewidth',6)
% ylim([-1,2])
xlim([-100,300])
title('Object Response')
legend('Wheel','', 'No Wheel')
xlabel('Time from object introduction (s)')
ylabel('DA Response (z-score)')

%%
figure
plotMeanSEM(ts_PETH_HC1,HC_pre_wheel,'r')
hold on
plotMeanSEM(ts_PETH_HC1,HC_pre_nowheel,'k')
legend('Wheel','', 'No Wheel')
title('Pre-Experiment Baseline')
xlabel('Time from context transition (s)')
ylabel('DA Response (z-score)')

%%
% p = [];ix=1;
% for i = 1:10:60000
%     p(ix) = ranksum(HC_pre_nowheel(:,i),HC_pre_wheel(:,i));
%     ix = ix+1;
% end
% %%
% p1 =  p;
% p1(p>.05) = nan;
% p1(~isnan(p1)) = 1.5;
% plot(ts_PETH_HC1(1:10:60000),p1,'k','linewidth',6)


%%
figure
plotMeanSEM(ts_PETH_HC2,HC_post_wheel,'r')
hold on
plotMeanSEM(ts_PETH_HC2,HC_post_nowheel,'k')
title('Post-Experiment Baseline')
legend('Wheel','', 'No Wheel')
xlabel('Time from context transition (s)')
ylabel('DA Response (z-score)')

%%
% p = [];ix=1;
% for i = 1:10:60000
%     p(ix) = ranksum(HC_post_nowheel(:,i),HC_post_wheel(:,i));
%     ix = ix+1;
% end
% %%
% p1 =  p;
% p1(p>.05) = nan;
% p1(~isnan(p1)) = 2.5;
% plot(ts_PETH_HC2(1:10:60000),p1,'k','linewidth',6)
% ylim([-1,3])
xlim([-100,300])

