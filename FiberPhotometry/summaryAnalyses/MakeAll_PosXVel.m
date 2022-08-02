fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
fils = fils(cellfun(@any,regexp(fils,'sessiondata')));

%%
clear task vel_avg
for i = 1:length(fils)
    load(fils{i})
    [a,subj] = fileparts(sessiondata.path);
    sl = regexp(sessiondata.path,'\');
    taskt = sessiondata.path(sl(3)+1:sl(4)-1);
    ts_neural = (1:length(sessiondata.neural.signal_DFoF))/sessiondata.neural.fs_neural;
    
    good_ix = 1:length(sessiondata.behavior.ts_video);
    
    vel = sessiondata.behavior.vel( good_ix);
    vel_upsample = interp1(sessiondata.behavior.ts_video,vel,ts_neural);
    vel_avg(i,:) =  avghist(vel_upsample,double(sessiondata.neural.signal_DFoF),0:20:200);
   % r(i) = corr(vel_upsample',double(sessiondata.neural.signal_DFoF)','rows','pairwise','type','spearman');
    task{i} = taskt;
    d = regexp(subj,'-');
    subject{i} = subj(1:d(1)-1);
    i
end


%%

for i = 1:size(vel_avg,1)
    rr(i) = corr(vel_avg(i,:)',[1:11]','type','spearman','rows','pairwise');
end

%%
kp_DA = cellfun(@any,regexp(subject,'DA'));
kp_NE = cellfun(@any,regexp(subject,'NE'));
figure
plot(nanmean(vel_avg(kp_DA,:)))


hold on
plot(nanmean(vel_avg(kp_NE,:)))'

%%
figure
t  =unique(task);
ax = tight_subplot(1,4); 
ix = 1;
for i = [1 3 4 5]
    
  axes(ax(ix))

kp_t = cellfun(@any,regexp(task,t{i})) & kp_DA;

plot(0:20:200,nanmean(vel_avg(kp_t,:)),'k')

hold on


kp_t = cellfun(@any,regexp(task,t{i})) & kp_NE;

plot(0:20:200,nanmean(vel_avg(kp_t,:)))
title(t{i})
ix = ix+1;
end
%%
close all
i =3;

kp_t = cellfun(@any,regexp(task,t{i})) & kp_DA;


figure
plotMeanSEM(0:20:200,vel_avg(kp_t,:),'r')
ylim([-.5 1.5])
xlim([0 140])
ylabel('Response')
xlabel('Velocity (pixel/s)')

figure

kp_t = cellfun(@any,regexp(task,t{i})) & kp_NE;

plotMeanSEM(0:20:200,vel_avg(kp_t,:),'k')
ylim([-.5 1.5])
xlim([0 140])
ylabel('Response')
xlabel('Velocity (pixel/s)')

%%
i=4;
kp_t = find(cellfun(@any,regexp(task,t{i})) & kp_NE);
k  = gaussian2Dfilter([1000 1000],[1 1]);
clear cat_map
for i = 1:length(fils)
    
    load(fils{i})
    %   cd(fils{kp_t(i)})
    [x1] = sessiondata.behavior.position.left_ear;
    [x2] = sessiondata.behavior.position.right_ear;
    ts_neural = (1:length(sessiondata.neural.signal_DFoF))/sessiondata.neural.fs_neural;
    x = nanmean([x1(:,1) x2(:,1)],2);
    y = nanmean([x1(:,2) x2(:,2)],2);
    good_ix = 1:length(sessiondata.behavior.ts_video);
    xbin = min(x_upsample):range(x_upsample)/100:max(x_upsample);
    ybin = min(y_upsample):range(y_upsample)/100:max(y_upsample);
    x_upsample = interp1(sessiondata.behavior.ts_video,x(good_ix),ts_neural);
    y_upsample = interp1(sessiondata.behavior.ts_video,y(good_ix),ts_neural);
    [occ,~,~,idx] = histcn([x_upsample' y_upsample'],xbin,ybin);
    sz = size(occ);
    kp  = all(idx>0,2);
    tmp = accumarray(idx(kp,:),double(sessiondata.neural.signal_DFoF(kp)),sz,@nanmean,nan)./occ;
    tmp = nanconvn(tmp,k,'nanout',true);
    figure
    imagesc(tmp,[-.5 .5])
    title(sessiondata.task)
    waitforbuttonpress
    close all
   
end


