dirs{1} = ['R:\DANEHippocampalResponse\DA3h7\CPP\DA3h7-231218-124808']
dirs{2} = ['R:\DANEHippocampalResponse\DA3h4\CPP\DA3h4-231218-122315']
dirs{3} = ['R:\DANEHippocampalResponse\DA3h8\CPP\DA3h8-231218-133742']
dirs{4} = ['R:\DANEHippocampalResponse\DA2h9\CPP\DA2h9-231218-131205']

makeField = {'MergeDLCwNeural'};
for i = 2:4
    
    sm_MakeSessionStruct(dirs{i},'makeField',makeField)
end
%%
clear DA_map

k  = gaussian2Dfilter([100 100],[2 2]);
for i = 1:4
    
    cd(dirs{i})
    load('sessiondata.mat')
    N_rec =length(sessiondata.behavior.ts_video);
    clear pos
    pos(:,1) = nanmean([sessiondata.behavior.position.left_ear(:,1), sessiondata.behavior.position.right_ear(:,1)],2);
    pos(:,2) = nanmean([sessiondata.behavior.position.left_ear(:,2), sessiondata.behavior.position.right_ear(:,2)],2);
    
    pos = pos(1:N_rec,:);
    ts = (1:length(sessiondata.neural.signal_DFoF))/sessiondata.neural.fs_neural;
    
    
    
    [occt,~,~,b] = histcn(pos,200:5:450,175:5:255);
   
    neural_beh = interp1(ts',sessiondata.neural.signal_DFoF,sessiondata.behavior.ts_video);
    ts_vid = sessiondata.behavior.ts_video;
    ts_vid = ts_vid(sessiondata.behavior.ts_video>600);
    cross_in = find(diff(pos(sessiondata.behavior.ts_video>600,1)<300)>0);
    cross_out = find(diff(pos(sessiondata.behavior.ts_video>600,1)<300)<0);
    
    
    if cross_out(1)< cross_in(1)
        cross_out = cross_out(2:end);
    end
       if length(cross_in)> length(cross_out)
        cross_in = cross_in(1:end-1);
      end
      
    cr_in = ts_vid([cross_in cross_out]);
    kp = diff(cr_in,[],2)>1;
    
    cr_in = cr_in(kp,:);
     
    cross_in = find(diff(pos(sessiondata.behavior.ts_video>600,1)>300)>0);
    cross_out = find(diff(pos(sessiondata.behavior.ts_video>600,1)>300)<0);
    
    
    if cross_out(1)< cross_in(1)
        cross_out = cross_out(2:end);
    end
    
      if length(cross_in)> length(cross_out)
        cross_in = cross_in(1:end-1);
      end
      
      
    
    
    cr_out = ts_vid([cross_in cross_out]);
    kp = diff(cr_out,[],2)>1;
    
    cr_out = cr_out(kp,:);
    
     [idx_in,early,late,ts]  = sm_getIndicesAroundEvent(cr_in(:,1),30,30,sessiondata.neural.fs_neural,length(sessiondata.neural.signal_DFoF));
     idx_in = idx_in(~early & ~late,:);
     [idx_out,early,late,ts] = sm_getIndicesAroundEvent(cr_out(:,1),30,30,sessiondata.neural.fs_neural,length(sessiondata.neural.signal_DFoF));
   idx_out = idx_out(~early & ~late,:);
    
  %  k  = gaussian2Dfilter([10000 1],1250);
  %  signal_DFoF = nanconvn(sessiondata.neural.signal_DFoF,k');
    
    kp = all(b>0,2) & sessiondata.behavior.ts_video>600;
    neural_beh = neural_beh(kp,:);
    b = b(kp,:);
    tmp = accumarray(b,double(neural_beh),[],@nanmean,nan);
    tmp(occt<1) = nan;
    DA_map(:,:,i) = nanconvn(tmp,k,'nanout',true);
     [occt,~,~,b] = histcn(pos(kp,:),200:3:450,175:3:255);
     occ(:,:,i) = occt;
end
ok = nanmean(DA_map,3)';
close all
figure
plot(nanmean(ok(:,2:end-1)))

%%

plot([1 2],[squeeze(nansum(nansum(occ(1:30,:,:)))) squeeze(nansum(nansum(occ(50:80,:,:))))]'/10)

set(gca,'xtick',1:2,'xticklabel',{'METH','SAL'})
ylabel('time spent (s)')
