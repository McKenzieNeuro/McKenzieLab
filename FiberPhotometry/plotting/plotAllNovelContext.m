
close all

fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'contextTransition')) & cellfun(@any,regexp(fils,'Novel Environment'));
fils = fils(kp);
[dirs,bi] = fileparts(fils);
% sl  =regexp(fils,'\');
 da  =regexp(fils,'-');
% subjs = cellfun(@(a,b) a(b(3)+1:b(4)-1),fils,sl,'uni',0);
dates = cellfun(@(a,b) a(b(1)+1:b(2)-1),fils,da,'uni',0);
%%

% for i = 1:length(dirs)
% 
% 
%     cd(dirs{i})
% 
%     if ~exist('sessiondata.mat')
% 
%         disp(dirs{i})
%     end
% 
% end
%%




%%
k  = gaussian2Dfilter([10000 1],[ 1017.3 1]);
clear context  sessionDate subj PETH
for i = 1:length(dirs)
    
    
    cd(dirs{i})
    
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        subj{i} = sessiondata.subject;
       % sessionDate{i} = dates{i};
       
        ts = cell2mat(sessiondata.contextEntry(:,2));
        fs =sessiondata.neural.fs_neural;
        %get PETH for each conext
        
        
        ok = nanconvn(sessiondata.neural.signal_DFoF,k');
        [ix,early,late,ts_PETH] = sm_getIndicesAroundEvent(ts,300,300,fs,length(sessiondata.neural.signal_DFoF));
        
         kp = ~early & ~late;
         
         ix = ix(kp,:);
         contexts{i,1} = sessiondata.contextEntry(kp,[1 ]);
          contexts{i,2} = sessiondata.contextEntry(kp,[5 ]);
         PETH{i} = ok(ix);
         
    end

i
end
subj = cellfun(@(a) a(1:5),subj,'uni',0)';

contexts(:,2) = cellfun(@(a) cellfun(@(b) nanPad(b,1),a,'uni',0) , contexts(:,2),'uni',0);
%%

kp_DA = cellfun(@any,regexp(subj,'DA2')) |  cellfun(@any,regexp(subj,'DA3'));

kp_NE = cellfun(@any,regexp(subj,'NE2')) |  cellfun(@any,regexp(subj,'NE3'));

homeCage = cellfun(@(a) find(cellfun(@any,regexp(a,'home'))),contexts(:,1),'uni',0);
other = cellfun(@(a) find(~cellfun(@any,regexp(a,'home'))),contexts(:,1),'uni',0);

DA_homeCage = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp_DA),homeCage(kp_DA)','uni',0)');
NE_homeCage = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp_NE),homeCage(kp_NE)','uni',0)');

DA_otherContext = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp_DA),other(kp_DA)','uni',0)');
NE_otherContext = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp_NE),other(kp_NE)','uni',0)');


%%
figure
plotMeanSEM(ts_PETH,NE_otherContext,'k')

hold on
plotMeanSEM(ts_PETH,NE_homeCage,'r')
xlim([-100,300])
p = [];
ix = 1;
% for jj = 1:1:size()
%     p(ix) = ranksum(NE_homeCage(:,jj),NE_otherContext(:,jj));
%     ix = ix+1;
% end
% p1 =  p;
% p1(p>.05) = nan;
% p1(~isnan(p1)) = 3;
% plot(ts_PETH,p1,'k','linewidth',6)


%%

figure
plotMeanSEM(ts_PETH,DA_otherContext,'k')

hold on
plotMeanSEM(ts_PETH,DA_homeCage,'r')
title('DA')
xlim([-100,300])


%%



kp  = kp_DA;
figure
ax  = tight_subplot(10,1);
for i = 0:9
   


axes(ax(i+1))
homeCage = cellfun(@(a) find(cellfun(@any,regexp(a(:,1),'ome'))),contexts(:,1),'uni',0);
other = cellfun(@(a,b) find(~cellfun(@any,regexp(a(:,1),'ome')) & cell2mat(b(:,1))==i),contexts(:,1),contexts(:,2),'uni',0);

homeCageResp = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp),homeCage(kp)','uni',0)');
otherContextResp = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp),other(kp)','uni',0)');




if i ==0
    plotMeanSEM(ts_PETH(1:100:end),otherContextResp(:,1:100:end),'r')
    tmp = otherContextResp(:,1:100:end);
    xlim([-100 300])
    ylim([-1 4])
    
else
    plotMeanSEM(ts_PETH(1:100:end),otherContextResp(:,1:100:end),'k')
    hold on
    plotMeanSEM(ts_PETH(1:100:end),tmp,'r')
    xlim([-100 300])
    ylim([-1 4])
    p = [];
    ix = 1;
    downsampled = otherContextResp(:,1:100:end);
    downsampledts = ts_PETH(1:100:end);
    for jj = 1:1:length(tmp)
        p(ix) = ranksum(tmp(:,jj),downsampled(:,jj));
        ix = ix+1;
    end
    p1 =  p;
    p1(p>.05) = nan;
    p1(~isnan(p1)) = 3;
    plot(downsampledts,p1,'k','linewidth',6)
end
end
%%
kp  = kp_NE;
figure
ax  = tight_subplot(10,1);
for i = 0:9
   


axes(ax(i+1))
homeCage = cellfun(@(a) find(cellfun(@any,regexp(a(:,1),'ome'))),contexts(:,1),'uni',0);
other = cellfun(@(a,b) find(~cellfun(@any,regexp(a(:,1),'ome')) & cell2mat(b(:,1))==i),contexts(:,1),contexts(:,2),'uni',0);

homeCageResp = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp),homeCage(kp)','uni',0)');
otherContextResp = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp),other(kp)','uni',0)');




if i ==0
    plotMeanSEM(ts_PETH(1:100:end),otherContextResp(:,1:100:end),'r')
    tmp = otherContextResp(:,1:100:end);
    xlim([-100 300])
    ylim([-1 4])
    
else
    plotMeanSEM(ts_PETH(1:100:end),otherContextResp(:,1:100:end),'k')
    hold on
    plotMeanSEM(ts_PETH(1:100:end),tmp,'r')
    xlim([-100 300])
    ylim([-1 4])
    p = [];
    ix = 1;
    downsampled = otherContextResp(:,1:100:end);
    downsampledts = ts_PETH(1:100:end);
    for jj = 1:1:length(tmp)
        p(ix) = ranksum(tmp(:,jj),downsampled(:,jj));
        ix = ix+1;
    end
    p1 =  p;
    p1(p>.05) = nan;
    p1(~isnan(p1)) = 3;
    plot(downsampledts,p1,'k','linewidth',6)

end
end



