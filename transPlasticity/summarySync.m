fils = getAllExtFiles('R:\TransPlasticity','mat',1);
kp = contains(fils,'spikes.cellinfo.mat');
fils= fils(kp);
kp_ASV = contains(fils,'ASV');
%%

%%
clear Ct1
idx = 1;
for i  = 1:length(fils)
    load(fils{i})
    maxT = max(cellfun(@max,spikes.times));
    FR = cellfun(@length,spikes.times)/maxT;
    
    kp = FR>.1 & FR<10;
    
    if sum(kp)>1
        [C1, Cmat1, prof1, S1] = spikeSync(spikes.times(kp));
        [Ct1{idx},tc1]    = spikeSyncWindow(prof1,'Width',60,'Step',1);
        sesID{idx} = fils{i};
        idx= idx+1;
    end
    i
end

%%

fils = getAllExtFiles('R:\WSun\ePhys\OptoRAM','mat',1);
kp = contains(fils,'spikes.cellinfo.mat');
fils= fils(kp);

%%
clear Ct2
%%
for i  = 1:length(fils)
    load(fils{i})
    maxT = max(cellfun(@max,spikes.times));
    FR = cellfun(@length,spikes.times)/maxT;
    
    kp = FR>.1 & FR<10;
    
    if sum(kp)>1
        [C2, Cmat2, prof2, S2] = spikeSync(spikes.times(kp));
        [tmp,tc2]    = spikeSyncWindow(prof2,'Width',60,'Step',1);
        [kp,b] = histc(pulseInfo.time(:,1),tc2);
        tmp(b(b>0)) = nan;
         Ct2{i} = tmp;
        
    end
    i
end

%%
figure
PKR_all = cell2mat(Ct1(kp_ASV)');
ctr_all = cell2mat(Ct2');
figure
plot(0:.001:.4,histc(ctr_all,0:.001:.4)/sum(histc(ctr_all,0:.001:.4)),'k')
hold on
plot(0:.001:.4,histc(PKR_all,0:.001:.4)/sum(histc(PKR_all,0:.001:.4)),'r')