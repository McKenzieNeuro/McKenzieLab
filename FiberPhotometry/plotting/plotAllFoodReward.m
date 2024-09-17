clear all
fils = getAllExtFiles('R:\ASommer\GNemer\Data\Food Reward','mat',1);
kp = contains(fils,'ContextTransitionRevised');
fils = fileparts(fils(kp));

%%
clear wheel obj_nowheel obj_wheel
iix1 =1;
iix = 1;
for i = 1:length(fils)
    cd(fils{i})
    
    if exist('signal_data.mat')
        load('signal_data.mat')
        load('ContextTransitionRevised.mat')
        fs  = 1000;
        newts = (1:length(analyzed.filtered))/fs;
        filtered = interp1(analyzed.time*60,analyzed.filtered,newts);
        
        inArena = cell2mat(data(contains(data(:,1),'large'),2));
        inArena = inArena(1);
        obj = cell2mat(data(contains(data(:,1),'object'),2));
        obj = obj(obj>inArena);
        
        [ix,early,late,ts_PETH1] = sm_getIndicesAroundEvent(obj,30,660,fs,length(filtered));
        
        
        kp_wheel = contains(data(:,1),'wheel');
        
        if any(kp_wheel)
            wheelT =   cell2mat(data(kp_wheel,2));
            [ix1,early,late,ts_PETH] = sm_getIndicesAroundEvent(wheelT,30,600,fs,length(filtered));
            obj_wheel(iix,:) = filtered(ix);
            wheel(iix,:) = filtered(ix1);
            iix = iix+1;
        else
            obj_nowheel(iix1,:) = filtered(ix);
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

