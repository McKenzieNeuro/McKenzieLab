
files = getAllExtFiles('R:\DGregg\NeuralData','mat',1);
keepFiles = contains(files,'sessiondata');

files = files(keepFiles);

%%

% get pct REM
missing  = true(size(files));
clear wavspec
state_pct = nan(length(files),3,3);

state_tim = nan(length(files),3,3);

obs_tim = nan(length(files),3,3);

wake_pct = nan(length(files),3);
stimSes = nan(length(files),1);

for i = 1:length(files)
    v = load(files{i});
    sessiondata = v.sessiondata;
    if isfield(sessiondata,'sleepData') 
        try
        neuralTS = sessiondata.poseData.neuralTS;
        vel = sessiondata.poseData.velocity;
        dt = .03;
        newTS = dt:dt:max(neuralTS);
        
        minTS = min(length(vel),length(neuralTS));
        vel = interp1(neuralTS(1:minTS),vel(1:minTS),newTS);
        
        idx=1;
        for j = [1 3 5]
            REMON = find(diff([0;sessiondata.sleepData.SleepState.idx.states==j])>0);
            REMOFF = find(diff([0;sessiondata.sleepData.SleepState.idx.states==j])<0);
            bins = cumsum([0 sessiondata.seizure_obs]);
            binep = [bins(1:end-1)'+.001 bins(2:end)'];
            if length(REMON)>length(REMOFF) && max(REMON)>max(REMOFF)
                REMOFF = [REMOFF;bins(end)];
            end
            
          
         
            rem_ep = [REMON REMOFF];
            
            for k = 1:3
                
                goodepoch = excludeEpochs(binep(k,:),[sessiondata.stimON-.1 sessiondata.stimON+2]);
                totdur = sum(diff(IntersectEpochs2(goodepoch,rem_ep),[],2));
                
                if isempty(goodepoch)
                    error('here')
                end
                state_pct(i,k,idx) = totdur/sum(diff(goodepoch,[],2));
                state_tim(i,k,idx) = totdur;
                obs_tim(i,k,idx) = sum(diff(goodepoch,[],2));
                if j==5
                    evs = rem_ep(InIntervals(REMON(:,1),goodepoch),:);
                    ch = sessiondata.channelID(contains(sessiondata.channel,'CA1(R)'));
                    
                    if k==1
%                         [h,wavspec1,ts,freqs] = sm_EventTriggeredSpect(sessiondata.lfpFile,evs(:,1), ...
%                             'channel',ch,'plotIntervals',[10 11],'plotIt',false);
%                         
%                         [h,wavspec2,ts,freqs] = sm_EventTriggeredSpect(sessiondata.lfpFile,evs(:,2), ...
%                             'channel',ch,'plotIntervals',[11 10],'plotIt',false);
%                         
%                         
%                         base_REM{i} =    nanmean(zscore([wavspec1(:,1:end-1250,:) wavspec2(:,1251:end,:)],[],2),3);
%                         
                        
                        [ix1,early1,late1,ts] = sm_getIndicesAroundEvent(evs(:,1),10,10,1/dt,length(vel));
                        [ix2,early2,late2,ts] = sm_getIndicesAroundEvent(evs(:,2),10,10,1/dt,length(vel));
                       ix1 = ix1(~(early1 | late1 | early2 |late2),:);
                     
                        ix2 = ix2(~(early1 | late1 | early2 |late2),:);
                        base_vel{i} = nanmean(vel([ix1 ix2]),1);
                    elseif k==2
%                          [h,wavspec1,ts,freqs] = sm_EventTriggeredSpect(sessiondata.lfpFile,evs(:,1), ...
%                             'channel',ch,'plotIntervals',[10 11],'plotIt',false);
%                         
%                         [h,wavspec2,ts,freqs] = sm_EventTriggeredSpect(sessiondata.lfpFile,evs(:,2), ...
%                             'channel',ch,'plotIntervals',[11 10],'plotIt',false);
%                         
%                         
%                         stim_REM{i} =    nanmean(zscore([wavspec1(:,1:end-1250,:) wavspec2(:,1251:end,:)],[],2),3);
%                         
                        
                        [ix1,early1,late1,ts] = sm_getIndicesAroundEvent(evs(:,1),10,10,1/dt,length(vel));
                        
                     
                        
                        [ix2,early2,late2,ts] = sm_getIndicesAroundEvent(evs(:,2),10,10,1/dt,length(vel));
                        ix1 = ix1(~(early1 | late1 | early2 |late2),:);
                     
                        ix2 = ix2(~(early1 | late1 | early2 |late2),:);
                        stim_vel{i} = nanmean(vel([ix1 ix2]),1);
                    end
                end
                
            end
            
            
            idx = idx+1;
        end
    end
    
    if  sessiondata.StimAmp>=0 && contains(sessiondata.StimType,'mono')
        
        stimSes(i) = true;
    elseif contains(sessiondata.StimType,'none')
        stimSes(i) = false;
    end
  
    i
    end
      subj{i} = sessiondata.subject;
end


%%
figure
%plot REM spectrogram and vel
ts = 1/1250:1/1250:40;
ts1 = dt:dt:40;
kp = ~cellfun(@isempty,stim_REM);
ok = stim_REM(kp);
ok = cellArrayTo3D(ok);
vel = cell2mat(stim_vel');
imagesc(ts,[],nanmean(ok,3),[-.5 .5])
 set(gca,'ytick',1:10:length(freqs),...
        'yticklabel',round(freqs(1:10:end)*10)/10,'ydir','normal')
    
    hold on
plot(ts1(1:end-1),nanmean(vel)*10,'w')
ylim([0 75])
%%

%plot REM for each subject
[usub,~,subj_id] = unique(subj);

%close all
figure
ax = tight_subplot(4,4);
u = unique(subj_id)';
clear p_TD
for i = 1:length(usub)
    axes(ax(i))
    kp  =subj_id==(i)&stimSes==1;
    errorbar(1:3,squeeze(nanmean(state_pct(kp,:,3))),SEM(squeeze(state_pct(kp,:,3))))
    title(usub(i))
    
    xlim([0 4])
   % ylim([.5 1.5])
    
    if ismember(i,[4 8 12])
        set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
    else
        set(gca,'xtick',1:3,'xticklabel','')
    end
    uRem_st(i,:) = squeeze(nanmean(state_pct(kp,:,3)));
%    [~,p_TD(i)] = ttest(TD_stim(subj_id==u(i),1),TD_stim(subj_id==u(i),2));
    
    
end

%%

%plot percent time in REM for control
[usub,~,subj_id] = unique(subj);

%close all
figure
ax = tight_subplot(4,4);
u = unique(subj_id)';
clear p_TD
for i = 1:length(usub)
    axes(ax(i))
    kp  =subj_id==(i)&stimSes==0;
    errorbar(1:3,squeeze(nanmean(state_pct(kp,:,3))),SEM(squeeze(state_pct(kp,:,3))))
    title(usub(i))
    
    xlim([0 4])
   % ylim([.5 1.5])
    
    if ismember(i,[4 8 12])
        set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
    else
        set(gca,'xtick',1:3,'xticklabel','')
    end
    uRem_nost(i,:) = squeeze(nanmean(state_pct(kp,:,3)));
%    [~,p_TD(i)] = ttest(TD_stim(subj_id==u(i),1),TD_stim(subj_id==u(i),2));
    
    
end
%%

%plot percent time in REM for stim
close all
figure
kp_subj = ~ismember(usub,{'EDS 2.2','EDS 2.0','PCP 4.0'});
errorbar(1:3,100*nanmean(uRem_st(kp_subj,:)),100*SEM(uRem_st(kp_subj,:)))
set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
ylabel('% time in REM')
figure
hold on
errorbar(1:3,100*nanmean(uRem_nost(kp_subj,:)),100*SEM(uRem_nost(kp_subj,:)))
%set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
ylabel('% time in REM')
xlim([0 4])

%%

kp = (stimSes==0 ) & ~contains(subj',{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
%kp = (stimSes==1 ) & ~contains(subj',{'EDS 2.2','EDS 2.0'}) ;
uRem_pre = squeeze(state_pct(kp,1,3));
uRem_stim = squeeze(state_pct(kp,2,3));
uRem_post = squeeze(state_pct(kp,3,3));

%uRem_pre = squeeze(state_tim(kp,1,3));
%uRem_stim = squeeze(state_tim(kp,2,3));
%uRem_post = squeeze(state_tim(kp,3,3));

sleep_effort = log([squeeze(obs_tim(kp,1,3));squeeze(obs_tim(kp,2,3));squeeze(obs_tim(kp,3,3))]);





REMpct = ([uRem_pre;uRem_stim;uRem_post]);
%epsilon = 1e-6;  % small constant to avoid 0 or 1
%REMpct_adj = max(min(REMpct, 1 - epsilon), epsilon);
%REMpct = log(REMpct_adj ./ (1 - REMpct_adj));

%REMpct(REMpct==0) = .00001;
%REMpct = log10(REMpct);
%REMpct = round(REMpct);

nes = size(uRem_pre,1);

block = (categorical([ones(nes,1);2*ones(nes,1);3*ones(nes,1)]));
sesID = categorical(repmat([1:nes]',3,1));


subjIDt = repmat(subj(kp)',3,1);
[~,~,subjIDt] = unique(subjIDt);
dat1 = table(REMpct,block,sesID,subjIDt);
%%
y2 = fitglme(dat1,'REMpct ~   block   + (1|sesID:subjIDt) + (1|subjIDt)');
%y3 = fitglme(dat1,'REMpct ~   block   + (1|sesID:subjIDt) + (1|subjIDt)','link','log','distribution','poisson','Offset',sleep_effort)

figure

d1 = y2.coefCI;
d1(2:3,:) = d1(2:3,:)+y2.Coefficients.Estimate(1);
d2 = y2.Coefficients.Estimate;
d2(2:3) = d2(2:3)+d2(1);
errorbar(1:3,d2,d1(:,1)-d2,d1(:,2)-d2)
%%
clear vel
for i = 1:length(files)
    v = load(files{i});
    sessiondata = v.sessiondata;
    if isfield(sessiondata,'sleepData')
        
        
        if ~isempty(sessiondata.seizure)
            evs = sessiondata.seizure(:,1);
            [ix,early,late,ts] = sm_getIndicesAroundEvent(evs,3,3,fs,length(sessiondata.poseData.velocity));
            ix = ix(~early & ~late,:);
            vel{i} = nanmean(sessiondata.poseData.velocity(ix),1);
        end
    end
    i
end



%%
chReg = [...
    {'PrL(L)' };...
    {'PrL(R)' };...
    {'AVT(L)' };...
    {'BLA(R)' };...
    {'CA1(L)' };...
    {'CA1(R)' };...
    {'gRSC (L)'};...
    {'LDT(L1)'}];


clear datInj
datInj(:,1) = [...
    {'EDS 1.0'} ; ...
    {'EDS 1.3'}; ...
    {'EDS 2.0'}; ...
    {'EDS 2.1'}; ...
    {'EDS 2.2'}; ...
    {'EDS 2.3'}; ...
    {'LDT 10.0'}; ...
    {'EDS 3.0'}; ...
    {'EDS 3.2'}; ...
    {'EDS 4.0'}; ...
    {'EDS 4.1'}; ...
    {'EDS 4.2'}; ...
    {'EDS 5.1'}; ...
    {'EDS 1.1'} ; ...
    {'PCP 4.0'} ; ...
    ];


datInj{1,2} = datenum(2022,6,17);
datInj{2,2} = datenum(2022,6,18);
datInj{3,2} = datenum(2022,7,7);
datInj{4,2} = datenum(2022,7,7);
datInj{5,2} = datenum(2022,7,8);
datInj{6,2} = datenum(2022,7,8);
datInj{7,2} = datenum(2022,3,25);
datInj{8,2} = datenum(2023,3,21);
datInj{9,2} = datenum(2023,3,23);
datInj{10,2} = datenum(2023,3,27);
datInj{11,2} = datenum(2023,3,24);
datInj{12,2} = datenum(2023,3,27);
datInj{13,2} = datenum(2023,4,20);
datInj{14,2} = datenum(2022,6,17);
datInj{15,2} = datenum(2024,3,18);

datInj{1,3} = '6/17/2022';
datInj{2,3} = '6/18/2022';
datInj{3,3} = '7/7/2022';
datInj{4,3} = '7/7/2022';
datInj{5,3} = '7/8/2022';
datInj{6,3} = '7/8/2022';
datInj{7,3} = '3/24/2022';
datInj{8,3} = '3/21/2023';
datInj{9,3} = '3/23/2023';
datInj{10,3} = '3/27/2023';
datInj{11,3} = '3/24/2023';
datInj{12,3} = '3/27/2023';
datInj{13,3} = '4/20/2023';
datInj{14,3} = '6/17/2022';
datInj{15,3} = '3/18/2024';


IED_rate_control =[];
TD_control =[];
IED_syn_stim =[];
IED_syn_con =[];
IED_rate_stim =[];
TD_stim =[];
stim_sub =[];
con_sub =[];
binnedPopStart =[];
binnedPopEnd =[];
ISI_stim =[];
sz_stim =[];
sz_stim_num =[];

daysPost_stim = [];
daysPost_con =[];
f_stim= [];
IED_rate_pre =[];
binnedTDEnd =[];
binnedTDStart =[];
sz_control =[];
subj_id = [];
filID_stim = [];
IED_stim_dur = [];
IED_stim_num =[];
sz_stim_obs = [];
condition = [];
sz_stim_dur =[];
for i = 1:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    [~,chIDX] = ismember(chReg,sessiondata.channel);
    
    
    
    tmp = nan(8,3);
    tmp(chIDX>0,:) = sessiondata.IED_rate(chIDX(chIDX>0),:);
    IED_rate_stim = cat(3,IED_rate_stim,tmp);
    
    tmp = nan(8,3);
    tmp(chIDX>0,:) = sessiondata.IED_num(chIDX(chIDX>0),:);
    IED_stim_num = cat(3,IED_stim_num,tmp);
    
    f_stim= [f_stim; {files{i}}];
    IED_syn_stim = cat(4,IED_syn_stim,sessiondata.IED_syn);
    
    TD_stim = [TD_stim;sessiondata.TD];
    sz_stim = [sz_stim;sessiondata.seizure_rate(:)'];
    sz_stim_num = [sz_stim_num;sessiondata.seizure_num(:)'];
    sz_stim_dur = [sz_stim_dur;sessiondata.seizure_duration];
    IED_stim_dur = [IED_stim_dur;sessiondata.IED_dur];
    sz_stim_obs = [sz_stim_obs;sessiondata.seizure_obs(:)'];
    stim_sub = [stim_sub;{sessiondata.subject}];
    binnedPopEnd = [binnedPopEnd;sessiondata.binnedPopEnd];
    binnedPopStart = [binnedPopStart;sessiondata.binnedPopStart];
    
    if ~isempty(sessiondata.stimON)
        tim(1)  = sessiondata.stimON(1);
         tim(2)  = sessiondata.stimON(end);
    else
        tim(1) = sessiondata.seizure_obs(1);
         tim(2) = sum(sessiondata.seizure_obs(1:2));
    end
    
    tmp = avghist(sessiondata.ts-tim(1),sessiondata.theta_delta',-3600:1800);
    binnedTDStart = [binnedTDStart; tmp];
    
    tmp = avghist(sessiondata.ts-tim(2),sessiondata.theta_delta',-1800:3600);
    binnedTDEnd = [binnedTDEnd; tmp];
    
    if isfield(sessiondata,'days_post_IHKA')
       
    daysPost_stim = [daysPost_stim;sessiondata.days_post_IHKA];
    else
            daysPost_stim = [daysPost_stim;sessiondata.days_post_pilo];

    end
    ISI_stim = [ISI_stim; median(diff(sessiondata.stimON)) sessiondata.StimAmp];
    
    
%    [~,b] = histc(sessiondata.theta_delta(sessiondata.ts<floor(sessiondata.stimON(1))),0:.05:3);
%    ok = sum(cell2mat(cellfun(@(a) histoc(a,sessiondata.ts(sessiondata.ts<floor(sessiondata.stimON(1))))',sessiondata.IED,'uni',0)'));
    
 %   IED_rate_pret =  accumarray(b(b>0),ok(b>0)',[length(0:.05:3) 1],@nanmean,nan);
    
    
    
  %  IED_rate_pre = [IED_rate_pre; IED_rate_pret'];
    
    
    [~,ID] = ismember(sessiondata.subject,datInj(:,1));
    filID_stim = [filID_stim;i];
    subj_id = [subj_id;ID];
    if contains(sessiondata.StimType,'mono') & sessiondata.StimAmp >0
        condition = [condition;1];
    elseif strmatch(sessiondata.StimType,'none')
        condition = [condition;0];
    else
        condition = [condition;-1];
        
    end
    i
end



%%


close all
figure
ax = tight_subplot(4,4);
u = unique(subj_id)';
clear p_sz
for i = 1:length(u)
    axes(ax(i))
    kp = condition==1& subj_id==u(i);
    errorbar(1:3,nanmean(sz_stim(kp,:)),SEM(sz_stim(kp,:)))
    title(datInj{u(i),1})
    
    xlim([0 4])
    %ylim([.5 1.5])
    
    if ismember(i,[4 8 12])
        set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
    else
        set(gca,'xtick',1:3,'xticklabel','')
    end
    
    [~,p_sz(i)] = ttest(sz_stim(kp,1),sz_stim(kp,2));
    
    
end

%%

%plot change in sz rate for stim
kp_subj = ~ismember(datInj(:,1),{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
kp = condition==1;
pre = accumarray(subj_id(kp),sz_stim(kp,1),[15 1],@nanmean,nan)*3600;
stim = accumarray(subj_id(kp),sz_stim(kp,2),[15 1],@nanmean,nan)*3600;
post = accumarray(subj_id(kp),sz_stim(kp,3),[15 1],@nanmean,nan)*3600;

tmp = 100*([stim post]-pre)./pre;

figure
errorbar(1:3,[0 nanmean(tmp(kp_subj,:))],[0 SEM([tmp(kp_subj,:)])])
        set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
        ylabel('% decrease in sz rate')
 xlim([0 4])
 ylim([-50 40])
 set(gca,'fontsize',16)
 
 %%
 %plot change in sz rate for control
kp_subj = ~ismember(datInj(:,1),{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
kp = condition==0;
pre = accumarray(subj_id(kp),sz_stim(kp,1),[15 1],@nanmean,nan)*3600;
stim = accumarray(subj_id(kp),sz_stim(kp,2),[15 1],@nanmean,nan)*3600;
post = accumarray(subj_id(kp),sz_stim(kp,3),[15 1],@nanmean,nan)*3600;

tmp = 100*([stim post]-pre)./pre;

figure
errorbar(1:3,[0 nanmean(tmp(kp_subj,:))],[0 SEM([tmp(kp_subj,:)])])
        set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
        ylabel('% decrease in sz rate')
 xlim([0 4])
 
 set(gca,'fontsize',16)
  ylim([-50 40])

 
%%
%close all

%plot T/D for each subject
figure
ax = tight_subplot(3,4);
u = unique(subj_id)';
clear p_TD
for i = 1:length(u)
    axes(ax(i))
     kp = condition==1& subj_id==u(i);
    errorbar(1:3,nanmean(TD_stim(kp,:)),SEM(TD_stim(kp,:)))
    title(datInj{u(i),1})
    
    xlim([0 4])
    ylim([.5 1.5])
    
    if ismember(i,[4 8 12])
        set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
    else
        set(gca,'xtick',1:3,'xticklabel','')
    end
    
    [~,p_TD(i)] = ttest(TD_stim(kp,1),TD_stim(kp,2));
    
    
end
%%
%plot IED for each subject
close all
figure
ax = tight_subplot(3,4);

tt=[];ix=1;
for i = 1:length(u)
    axes(ax(i))
    tmp = squeeze(nansum(IED_rate_stim(:,:,subj_id==u(i))))';
    
    tt(ix,:) = nanmean(tmp);
    errorbar(1:3,nanmean(tmp),SEM(tmp))
    title(datInj{u(i),1})
    
    xlim([0 4])
    %  ylim([.5 1.5])
    
    if ismember(i,[4 8 12])
        set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'})
    else
        set(gca,'xtick',1:3,'xticklabel','')
    end
    ix = ix+1;
end




%%
%plot sz rate pre vs stim for each subj (stim)
close all

any_sz = accumarray(subj_id,sz_stim(:,1),[],@sum,nan);
any_sz = ismember(subj_id,find(any_sz>0));
gd_target = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1'});

 kp = condition==1;

kp_subj = ~ismember(datInj(:,1),{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0','PCP 4.0'});
pre = accumarray(subj_id(kp),sz_stim(kp,1),[15 1],@mean,nan)*3600;
stim = accumarray(subj_id(kp),sz_stim(kp,2),[15 1],@mean,nan)*3600;
post = accumarray(subj_id(kp),sz_stim(kp,3),[15 1],@mean,nan)*3600;

figure
 loglog(pre(kp_subj),stim(kp_subj),'o')
 hold on
 plot([0:.01e-3:1.6e-3]*3600,[0:.01e-3:1.6e-3]*3600)
 xlabel('Pre stim rate')
 ylabel('Stim rate')
 
 %%
 
 tmp = pre(kp_subj);
 tmp(isnan(tmp)) = 0;
 nsub = length(tmp);
 close all
 figure
 hold on
 
 plot(1+ones(nsub,1).*(.5-rand(nsub,1))/5,tmp,'.','color','k','markersize',10)
  boxplot(tmp)
set(gca,'xtick',[],'fontsize',16)
ylabel('# sz/hr')

%%

kp_subj = ~ismember(datInj(:,1),{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0','PCP 4.0'});
pre = accumarray(subj_id(kp),sz_stim_dur(kp,1),[15 1],@nanmean,nan);
stim = accumarray(subj_id(kp),sz_stim_dur(kp,2),[15 1],@nanmean,nan);
post = accumarray(subj_id(kp),sz_stim_dur(kp,3),[15 1],@nanmean,nan);

figure
 loglog(pre(kp_subj),stim(kp_subj),'o')
 hold on
 plot(0:1:200,0:1:200)
 xlabel('Pre stim rate')
 ylabel('Stim rate')
 %%
 
 tmp = pre(kp_subj);
 %tmp(isnan(tmp)) = 0;
 nsub = length(tmp);
 close all
 figure
 hold on
 
 plot(1+ones(nsub,1).*(.5-rand(nsub,1))/5,tmp,'.','color','k','markersize',10)
  boxplot(tmp)
set(gca,'xtick',[],'fontsize',16)
ylabel('Mean sz duration (s)')
 
 
 %%
 %plot sz rate pre vs stim for each subj (control)

 
close all


 kp = condition==0;

kp_subj = ~ismember(datInj(:,1),{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
pre = accumarray(subj_id(kp),sz_stim(kp,1),[15 1],@mean,nan)*3600;
stim = accumarray(subj_id(kp),sz_stim(kp,2),[15 1],@mean,nan)*3600;
post = accumarray(subj_id(kp),sz_stim(kp,3),[15 1],@mean,nan)*3600;
figure
 loglog(pre(kp_subj),stim(kp_subj),'o')
 hold on
 plot([0:.01e-3:1.6e-3]*3600,[0:.01e-3:1.6e-3]*3600)
 xlabel('Pre stim rate')
 ylabel('Stim rate')
 
 
 %%
 
 %plot T/D block for stim
  kp = condition==1;
 figure
pre = accumarray(subj_id(kp),TD_stim(kp,1),[15 1],@mean,nan);
stim = accumarray(subj_id(kp),TD_stim(kp,2),[15 1],@mean,nan);
post = accumarray(subj_id(kp),TD_stim(kp,3),[15 1],@mean,nan);
 

d = nanmean([pre(kp_subj) stim(kp_subj) post(kp_subj)]);
s = SEM([pre(kp_subj) stim(kp_subj) post(kp_subj)]);
errorbar(1:3,d,s)
set(gca,'fontsize',16)
set(gca,'xtick',1:3,'xticklabel',{'Pre','Stim','Post'})
xlim([0 4])

%%
%plot T/D block for control
 
 kp = condition==0;
 figure
pre = accumarray(subj_id(kp),TD_stim(kp,1),[15 1],@mean,nan);
stim = accumarray(subj_id(kp),TD_stim(kp,2),[15 1],@mean,nan);
post = accumarray(subj_id(kp),TD_stim(kp,3),[15 1],@mean,nan);
 

d = nanmean([pre(kp_subj) stim(kp_subj) post(kp_subj)]);
s = SEM([pre(kp_subj) stim(kp_subj) post(kp_subj)]);
errorbar(1:3,d,s)
set(gca,'fontsize',16)
set(gca,'xtick',1:3,'xticklabel',{'Pre','Stim','Post'})
xlim([0 4])
ylim([0.75 1.1])
%%
%plotIED block for stim
kp_subj = ~ismember(datInj(:,1),{'EDS 3.0','EDS 5.1'});
kp_subj = ~ismember(datInj(:,1),{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
 kp = condition==1;
figure
tmp = squeeze(nansum1(IED_rate_stim(4:6,:,:),1))';
pre = accumarray(subj_id(kp),tmp(kp,1),[15 1],@nanmean,nan);
stim = accumarray(subj_id(kp),tmp(kp,2),[15 1],@nanmean,nan);
post = accumarray(subj_id(kp),tmp(kp,3),[15 1],@nanmean,nan);
 

d = nanmean([pre(kp_subj) stim(kp_subj) post(kp_subj)],1);
s = SEM([pre(kp_subj) stim(kp_subj) post(kp_subj)],1);
errorbar(1:3,d,s)
set(gca,'fontsize',16)
set(gca,'xtick',1:3,'xticklabel',{'Pre','Stim','Post'})
xlim([0 4])

%%

%plot prcent change
figure
d = ([stim post]-pre)./pre;
d = d(kp_subj,:);
s = SEM([ d]);
d = nanmean([ d]);

d = [0 d];
s = [0 s];
errorbar(1:3,d,s)
set(gca,'fontsize',16)
set(gca,'xtick',1:3,'xticklabel',{'Pre','Stim','Post'})
xlim([0 4])

%%


%plotIED block for stim
 kp = condition==0;
figure
tmp = squeeze(nansum1(IED_rate_stim(4:6,:,:),1))';
pre = accumarray(subj_id(kp),tmp(kp,1),[15 1],@nanmean,nan);
stim = accumarray(subj_id(kp),tmp(kp,2),[15 1],@nanmean,nan);
post = accumarray(subj_id(kp),tmp(kp,3),[15 1],@nanmean,nan);
 

d = nanmean([pre(kp_subj) stim(kp_subj) post(kp_subj)]);
s = SEM([pre(kp_subj) stim(kp_subj) post(kp_subj)]);
errorbar(1:3,d,s)
set(gca,'fontsize',16)
set(gca,'xtick',1:3,'xticklabel',{'Pre','Stim','Post'})
xlim([0 4])

%%

%plot prcent change
figure
d = ([stim post]-pre)./pre;
d = d(kp_subj,:);
s = SEM([ d]);
d = nanmean([ d]);

d = [0 d];
s = [0 s];
errorbar(1:3,d,s)
set(gca,'fontsize',16)
set(gca,'xtick',1:3,'xticklabel',{'Pre','Stim','Post'})
xlim([0 4])

 %%
 close all
clear p

 kp = condition==1;
for r = 1:8
kp_subj = ~ismember(datInj(:,1),{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
IED_stim = squeeze(nansum1(IED_rate_stim(r,:,:),1))';
preIED = accumarray(subj_id(kp),IED_stim(kp,1),[15 1],@mean,nan);
stimIED = accumarray(subj_id(kp),IED_stim(kp,2),[15 1],@mean,nan);
postIED = accumarray(subj_id(kp),IED_stim(kp,3),[15 1],@mean,nan);
[p(r)] = signtest(preIED(kp_subj),stimIED(kp_subj));
end
IED_stim = squeeze(nansum1(IED_rate_stim(4:6,:,:),1))';
preIED = accumarray(subj_id(kp),IED_stim(kp,1),[15 1],@mean,nan);
stimIED = accumarray(subj_id(kp),IED_stim(kp,2),[15 1],@mean,nan);
postIED = accumarray(subj_id(kp),IED_stim(kp,3),[15 1],@mean,nan);
figure
 loglog(preIED(kp_subj),stimIED(kp_subj),'o')
 hold on
 plot([0:.01e-3:10],[0:.01e-3:10])
 xlabel('Pre stim rate')
 ylabel('Stim rate')
 
%%
kp = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0','PCP 4.0'})& condition==1;
kp1 = repmat(kp,3,1);
nes = size(sz_stim,1);
sz_num =[sz_stim_num(:,1);sz_stim_num(:,2);sz_stim_num(:,3)];

tmp = squeeze(nansum1(IED_stim_num(4:7,:,:),1))';
ied_num = [tmp(:,1);tmp(:,2);tmp(:,3)];
ied_effort = log([IED_stim_dur(:,1);IED_stim_dur(:,2);IED_stim_dur(:,3)]);
effort = log([sz_stim_obs(:,1);sz_stim_obs(:,2);sz_stim_obs(:,3)]);
block = (categorical([ones(nes,1);2*ones(nes,1);3*ones(nes,1)]));


sesID = categorical(repmat([1:nes]',3,1));
subjID = categorical(repmat(subj_id,3,1));
theta = TD_stim(:);
dat = table(sz_num,block,sesID,subjID,theta,ied_num,ied_effort);
dat1 = dat(kp1,:);
y = fitglme(dat1,'sz_num ~   block   + (1|sesID:subjID) + (1|subjID)','link','log','distribution','poisson','Offset',effort(kp1));
y_hat = predict(y,dat1,'Offset',effort(kp1));


kp = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1'});
kp = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'}) & condition==1;

kp1 = repmat(kp,3,1);
dat1 = dat(kp1,:);
y2 = fitglme(dat1,'ied_num ~   block   + (1|sesID:subjID) + (1|subjID)','link','log','distribution','poisson','Offset',ied_effort(kp1));
%y_hat = predict(y,dat1,'Offset',ied_effort(kp1));

[a,~,b] =fixedEffects(y);
fits_all_allsz  = [a(2) b.pValue(2) b.Lower(2) b.Upper(2);a(3) b.pValue(3) b.Lower(3) b.Upper(3)];
[a,~,b] =fixedEffects(y2);
fits_all_allIED =  [a(2) b.pValue(2) b.Lower(2) b.Upper(2);a(3) b.pValue(3) b.Lower(3) b.Upper(3)];


%%
close all
%plot beta estimates for limbic IED
figure
 y = (exp(fits_all_allIED(:,1))-1)*100;
        neg= (exp(fits_all_allIED(:,3))-1)*100 - y;
        pos  =(exp(fits_all_allIED(:,4))-1)*100 - y;
        errorbar(1:3,[0; y],[0 ;neg],[0 ;pos])
        hold on
        xlim([0 4])
        set(gca,'xtick',1:3,'xticklabel',{'pre','stim','post'},'fontsize',16)
      %%
      
      figure
 y = (exp(fits_all_allsz(:,1))-1)*100;
        neg= (exp(fits_all_allsz(:,3))-1)*100 - y;
        pos  =(exp(fits_all_allsz(:,4))-1)*100 - y;
        errorbar(1:3,[0; y],[0 ;neg],[0 ;pos])
        hold on
        
%%

kp = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0','PCP 4.0'})& condition==0;
kp1 = repmat(kp,3,1);
nes = size(sz_stim,1);
sz_num =[sz_stim_num(:,1);sz_stim_num(:,2);sz_stim_num(:,3)];

tmp = squeeze(nansum1(IED_stim_num(4:7,:,:),1))';
ied_num = [tmp(:,1);tmp(:,2);tmp(:,3)];
ied_effort = log([IED_stim_dur(:,1);IED_stim_dur(:,2);IED_stim_dur(:,3)]);
effort = log([sz_stim_obs(:,1);sz_stim_obs(:,2);sz_stim_obs(:,3)]);
block = (categorical([ones(nes,1);2*ones(nes,1);3*ones(nes,1)]));


sesID = categorical(repmat([1:nes]',3,1));
subjID = categorical(repmat(subj_id,3,1));
theta = TD_stim(:);
dat = table(sz_num,block,sesID,subjID,theta,ied_num,ied_effort);
dat1 = dat(kp1,:);
y = fitglme(dat1,'sz_num ~   block   + (1|sesID:subjID) + (1|subjID)','link','log','distribution','poisson','Offset',effort(kp1));
y_hat = predict(y,dat1,'Offset',effort(kp1));


kp = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1'});
kp = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0','PCP 4.0'}) & condition==0;

kp1 = repmat(kp,3,1);
dat1 = dat(kp1,:);
y2 = fitglme(dat1,'ied_num ~   block   + (1|sesID:subjID) + (1|subjID)','link','log','distribution','poisson','Offset',ied_effort(kp1));
%y_hat = predict(y,dat1,'Offset',ied_effort(kp1));

[a,~,b] =fixedEffects(y);
fits_all_allsz  = [a(2) b.pValue(2) b.Lower(2) b.Upper(2);a(3) b.pValue(3) b.Lower(3) b.Upper(3)];
[a,~,b] =fixedEffects(y2);
fits_all_allIED =  [a(2) b.pValue(2) b.Lower(2) b.Upper(2);a(3) b.pValue(3) b.Lower(3) b.Upper(3)];


%%

%plot beta estimates for limbic IED
figure
 y = (exp(fits_all_allIED(:,1))-1)*100;
        neg= (exp(fits_all_allIED(:,3))-1)*100 - y;
        pos  =(exp(fits_all_allIED(:,4))-1)*100 - y;
        errorbar(1:3,[0; y],[0 ;neg],[0 ;pos])
        hold on
      %%
      
      figure
 y = (exp(fits_all_allsz(:,1))-1)*100;
        neg= (exp(fits_all_allsz(:,3))-1)*100 - y;
        pos  =(exp(fits_all_allsz(:,4))-1)*100 - y;
        errorbar(1:3,[0; y],[0 ;neg],[0 ;pos])
        hold on
        


%%
clear y2 
fits_all = nan(8,4);

for r = 1:8
IED_stim_all = squeeze(nansum1(IED_stim_num([r],:,:),1))';
ied_num =[IED_stim_all(:,1);IED_stim_all(:,2);IED_stim_all(:,3)];
%kp = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'})& condition==1;
kp = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1','PCP 4.0'})& condition==1;
dat = table(sz_num,block,sesID,subjID,theta,ied_num,ied_effort);

kp1 = repmat(kp,3,1);
dat1 = dat(kp1,:);
y2{r} = fitglme(dat1,'ied_num ~   block   + (1|sesID:subjID) + (1|subjID)','link','log','distribution','poisson','Offset',ied_effort(kp1));

[a,~,b] =fixedEffects(y2{r});
fits_all(r,:) =  [a(2) b.pValue(2) b.Lower(2) b.Upper(2)];
%y1_hat = predict(y,dat1,'Offset',ied_effort(kp1));
end

%%
close all
figure
offset = (0:9)/10-.5;
for i = 1:size(fits_all,1)
   
   
        y = (exp(fits_all(i,1))-1)*100;
        neg= (exp(fits_all(i,3))-1)*100 - y;
        pos  =(exp(fits_all(i,4))-1)*100 - y;
        errorbar(i*2,y,neg,pos)
        hold on
        plot(i*2,y,'o')
        
    end

plot([0 20],[0 0])
set(gca,'xtick',2:2:16,'xticklabel',chReg,'fontsize',16)


%%
clear y2 
fits = nan(8,9,4);
gd_sub = setdiff(stim_sub,{'EDS 3.0','EDS 5.1','PCP 4.0'});
%gd_sub = setdiff(stim_sub,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
for s = 1:length(gd_sub)
    kp = ismember(stim_sub,gd_sub{s});
    for r = 1:8
        IED_stim_all = squeeze(nansum1(IED_stim_num([r],:,:),1))';
        ied_num =[IED_stim_all(:,1);IED_stim_all(:,2);IED_stim_all(:,3)];
        
        dat = table(sz_num,block,sesID,subjID,theta,ied_num,ied_effort);
        
        kp1 = repmat(kp,3,1);
        dat1 = dat(kp1,:);
        if ~all(isnan(dat1.ied_num))
        y2 = fitglme(dat1,'ied_num ~   block   + (1|sesID)','link','log','distribution','poisson','Offset',ied_effort(kp1));
        
        [a,~,b] =fixedEffects(y2);
        fits(r,s,:) = [a(2) b.pValue(2) b.Lower(2) b.Upper(2)];
        %y1_hat = predict(y2,dat1,'Offset',ied_effort(kp1));
        end
    end
    
end
%%
figure
offset = (0:9)/10-.5;
for i = 1:size(fits,1)
    for j = 1:size(fits,2)
   
        y = (exp(fits(i,j,1))-1)*100;
        neg= (exp(fits(i,j,3))-1)*100-y;
        pos  =(exp(fits(i,j,4))-1)*100-y;
        errorbar(i*2+offset(j),y,neg,pos)
        hold on
        
        
    end
end
plot([0 20],[0 0])
set(gca,'xtick',2:2:16,'xticklabel',chReg,'fontsize',16)
%nanmean(sz_num(double(block)==1 & kp1))-nanmean(sz_num(double(block)==2 & kp1))
%nanmean(y_hat(double(block(kp1))==1)) - nanmean(y_hat(double(block(kp1))==2))
%%
close all
kp_IED = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1','PCP 4.0'});%,'EDS 2.2','EDS 2.0'});
td1 = binnedTDStart(kp_IED,1:3600);
IED1 = binnedPopStart(kp_IED,1:3600);
bins = 0:.1:3;
pre =[];
for i = 1:size(td1,1)
    pre(i,:) =  avghist(td1(i,:),IED1(i,:),bins) ;
end



td2 = [binnedTDStart(:,3600:end-1) binnedTDEnd(:,1:1800)];
IED2 = [binnedPopStart(:,3600:end) binnedPopEnd(:,1:1800)];


%td2 = [binnedTDEnd(:,1:1800)];
%IED2 = [binnedPopEnd(:,1:1800)];


stim =[];
for i = 1:size(td1,1)
    stim(i,:) =  avghist(td2(i,:),IED2(i,:),bins) ;
end

p= [];
for i = 1:size(stim,2)
    
    [p(i)] = signtest(pre(:,i),stim(:,i));
end
p(p>.05) = nan;

p(~isnan(p)) = .2;
figure
plotMeanSEM(bins,pre,'k','yAxisFunction',@nanmean)
plotMeanSEM(bins,stim,'r','yAxisFunction',@nanmean)

plot(bins,p,'o')
xlim([.6 1.3])
%%
close all
k = gaussian2Dfilter([10000 1],1250/16);
% typical example
144
154
for i = 144:length(files)
    load(files{i})
    if isfield(sessiondata,'IED') & ~isempty( sessiondata.stimON)
        %load(  'R:\DGregg\NeuralData\EDS\OL3\2-7-2023(12.59)\RHS_230207_130000\EDS 2.2.sessiondata.mat')
        figure
        ts1 = (0: sessiondata.ts(end)) - sessiondata.stimON(1);
        set(gcf,'position',[697    57   560   479]);
        plot(sessiondata.ts-sessiondata.stimON(1),sessiondata.theta_delta)
        hold on
        plot(ts1,nanconvn(histc( sessiondata.IED{2},0: sessiondata.ts(end)),k)*50)
        plot([sessiondata.stimON(1) sessiondata.stimON(1)]- sessiondata.stimON(1),[0 7])
        plot([sessiondata.stimON(end) sessiondata.stimON(end)]- sessiondata.stimON(1),[0 7])
        xlim([-2100 7200])
        waitforbuttonpress
        close all
        
    end
end
%%

ISI_stim =[];CA1_AMYG =[];

for i = 1:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    
    if sessiondata.StimAmp>0 && isfield(sessiondata,'IED') && strcmp(sessiondata.StimType,'mono')
        
        
        preEpoch = [0 sessiondata.stimON(1)];
        postEpoch = [sessiondata.stimON(end) sessiondata.ts(end)];
        stimEpoch =  [sessiondata.stimON(1:end-1)+3 sessiondata.stimON(2:end)-.25];
        ok =[];
        for j = 1:3
            switch j
                case 1
                    ep = preEpoch;
                    
                case 2
                    ep = stimEpoch;
                    
                case 3
                    ep = postEpoch;
                    
            end
            
            totDur = sum(diff(ep,[],2));
            
            
            
            [chIDX2] = ismember(sessiondata.channel,'BLA(R)');
            [chIDX1] = ismember(sessiondata.channel,'CA1(R)');
            
            if any(chIDX2)
                
                kp1 = InIntervals(sessiondata.IED{chIDX1},ep);
                kp2 = InIntervals(sessiondata.IED{chIDX2},ep);
                ok(j,:) = CrossCorr(sessiondata.IED{chIDX1}(kp1),sessiondata.IED{chIDX2}(kp2),1,101)/sum(kp1);
                
            end
        end
        CA1_AMYG = cat(3,CA1_AMYG,ok);
        
    end
    i
end

%%
idx=1;
IED_syn_stim_s = nan(8,8,3,13);
for i = 1:1:13
    
    IED_syn_stim_s(:,:,:,idx) = nanmean(IED_syn_stim(:,:,:,subj_id==i),4);

        IED_syn_con_s(:,:,:,idx) = nanmean(IED_syn_con(:,:,:,subj_id_c==i),4);
idx = idx+1;
end


%%
IED_diff_C = IED_syn_con_s(:,:,2,:) - IED_syn_con_s(:,:,1,:);
IED_diff_S = IED_syn_stim_s(:,:,2,:)- IED_syn_stim_s(:,:,1,:);
figure
imagesc(nanmean(IED_diff_S,4))
figure
imagesc(nanmean(IED_diff_C,4))
%%
close all
kp_IED = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1','EDS 2.2','EDS 2.0'});
p = nan(8); d = nan(8);
for i = 1:8
    for j =1:8
        try
            d(i,j) = nanmean(squeeze(IED_syn_stim(i,j,2,kp_IED))- squeeze(IED_syn_stim(i,j,1,kp_IED)));
   [~,p(i,j)]   = ttest(squeeze(IED_syn_stim(i,j,2,kp_IED)),squeeze(IED_syn_stim(i,j,1,kp_IED)));
        end
    end
end
[a,b] =  find(p<.05);
figure
imagesc(d)
colormap('bluewhitered')
hold on
plot(b,a,'o')
%%
close all
ax  = tight_subplot(2,2)

for i = 1:4
    axes(ax(i))
    if i<4
        imagesc(nanmean(IED_syn_stim(:,:,i,:) ,4),[0 .5])
       %   colormap('jet')
    else
       imagesc(nanmean(IED_syn_stim(:,:,2,:),4) - nanmean(IED_syn_con(:,:,2,:),4))
       colormap('bluewhitered')
    end
    
end
%%
k  = gaussian2Dfilter([ 1 10000 ],400);
k1 =k;
k2 = k;
k1(1:5000) = 0;

k1 = k1*2;
k2(5001:end) = 0;

k2 = k2*2;

binnedPopStart1 =[];
binnedPopEnd1 =[];
close all
for i = 1:size(binnedPopStart,1)
    
    binnedPopStart1(i,:) = (nanconvn(binnedPopStart(i,:) - nanmean(binnedPopStart(i,1800:3600)),k1));%./ nanmedian(binnedPopStart(i,1:3600));
    binnedPopEnd1(i,:) = nanconvn(binnedPopEnd(i,:) - nanmean(binnedPopStart(i,1800:3600)),k1);%./ nanmedian(binnedPopStart(i,1:3600));
    binnedTDStart1(i,:) = nanconvn(binnedTDStart(i,:) - nanmean(binnedTDStart(i,1800:3600)),k1);
    binnedTDEnd1(i,:) = nanconvn(binnedTDEnd(i,:) - nanmean(binnedTDEnd(i,1800:3600)),k1);
    
end
%%
binnedPopStart1(isinf(binnedPopStart1)) = nan;
binnedPopEnd1(isinf(binnedPopEnd1)) = nan;

kp = ismember(subj_id,u(p_TD<.025));
binnedPopStart1 = binnedPopStart1(kp,:);
binnedPopEnd1 = binnedPopEnd1(kp,:);
binnedTDStart1 = binnedTDStart1(kp,:);
binnedTDEnd1 = binnedTDEnd1(kp,:);
close all
figure
plotMeanSEM(-3600:1799-60,binnedPopStart1(:,1:end-60),'k')
hold on
plotMeanSEM(1801+100:7200,(binnedPopEnd1(:,101:end)),'k')


plotMeanSEM(-3600:1799-60,binnedTDStart1(:,1:end-61),'r')
hold on
plotMeanSEM(1801+100:7200,(binnedTDEnd1(:,102:end)),'r')


set(gca,'xtick',[-3600:600:1730 1800:600:7200],'xticklabel',[-3600:600:1730 -1800:600:3600])
xlim([-1800 5400])
pp=[];
for i = 1:5400
    pp(i) = signtest(binnedPopStart1(:,i));
end

pp(pp>.05) = nan;
pp(~isnan(pp)) = .27;

plot(-3600:1799-60,pp(1:end-60),'k','linewidth',6)

pp =[];
for i = 1:5400
    pp(i) = signtest(binnedPopEnd1(:,i));
end

pp(pp>.05) = nan;
pp(~isnan(pp)) = .27;

plot(1801+100:7200,pp(101:end),'k','linewidth',6)



pp=[];
for i = 1:5400
    pp(i) = signtest(binnedTDStart1(:,i));
end

pp(pp>.05) = nan;
pp(~isnan(pp)) = .25;

plot(-3600:1799-60,pp(1:end-60),'r','linewidth',6)

pp =[];
for i = 1:5400
    pp(i) = signtest(binnedTDEnd1(:,i));
end

pp(pp>.05) = nan;
pp(~isnan(pp)) = .25;

plot(1801+100:7200,pp(101:end),'r','linewidth',6)


plot([0 0],[-.1 .05],'k')
plot([3600 3600],[-.1 .05],'k')
%%
IED_mean_stim_12 =[];
for i = 1:6
    
    if i<7
        IED_mean_stim_12(i,:) = (IED_rate_stim(i,2,:) - IED_rate_stim(i,1,:))./(IED_rate_stim(i,1,:));
    else
        IED_mean_stim_12(i,:) = (sum(IED_rate_stim(7:8,2,:),1) - sum(IED_rate_stim(7:8,1,:),1))./( sum(IED_rate_stim(7:8,1,:),1));
    end
end


IED_mean_stim_23 =[];
for i = 1:6
    
    
    IED_mean_stim_23(i,:) = (IED_rate_stim(i,2,:) - IED_rate_stim(i,3,:))./(IED_rate_stim(i,2,:) + IED_rate_stim(i,3,:));
end

IED_mean_stim_13 =[];
for i = 1:6
    
    if i<7
        IED_mean_stim_13(i,:) = (IED_rate_stim(i,3,:) - IED_rate_stim(i,1,:))./(IED_rate_stim(i,3,:) + IED_rate_stim(i,1,:));
    else
        IED_mean_stim_13(i,:) = (sum(IED_rate_stim(7:8,3,:),1) - sum(IED_rate_stim(7:8,1,:),1))./(sum(IED_rate_stim(7:8,3,:),1) + sum(IED_rate_stim(7:8,1,:),1));
    end
end

TD_stim_mean = (TD_stim(:,2) - TD_stim(:,1))./( TD_stim(:,1));
IED_mean_stim_12(isinf(IED_mean_stim_12)) = nan;


%%
usub = unique(stim_sub);

[~,b] = ismember(stim_sub,usub);

close all
for i = 1:6
    try
    figure
    
    violin(num2cell(IED_mean_stim_12(:,b==i),2)')
    subj_mean(i,:) = mean(IED_mean_stim_12(:,b==i),2);
    title(usub{i})
    end
end

%%
figure
kp_IED = ~ismember(stim_sub,{'EDS 3.0','EDS 5.1'});


x = 100*[IED_mean_stim_12(:,kp_IED)];
nes = size(x,2);
x= x(:);

g = [linearize(repmat([1:6]',1,nes))];


boxplot(x,g,'notch','on')
hold on
plot([0 8],[0 0],'k')



%%
ok = squeeze((nansum(IED_rate_stim(:,:,:))- repmat(nansum(IED_rate_stim(:,1,:)),1,3,1))./repmat(nansum(IED_rate_stim(:,1,:)),1,3,1));
ok1 = squeeze((nansum(IED_rate_control(:,:,:))- repmat(nansum(IED_rate_control(:,1,:)),1,3,1))./repmat(nansum(IED_rate_control(:,1,:)),1,3,1));


ok = squeeze((nansum(IED_rate_stim(:,:,:),1)));
%ok1 = squeeze((nansum(IED_rate_control(:,:,:),1)));


%ok = ok(:,daysPost_stim>190);
%ok1 = ok1(:,daysPost_con>190);
x = [ok(:);ok1(:)];
g = [linearize(repmat([1:3]',1,size(ok,2))); linearize(repmat([4:6]',1,size(ok1,2)))];
figure
boxplot(x,g,'notch','on')
%%
close all
figure
ok = squeeze((nansum(IED_rate_stim(:,:,:),1)));
plot(0:.01:1,cumsum(histc(ok(1,:),0:.01:1))/size(ok,2))
hold on
plot(0:.01:1,cumsum(histc(ok(2,:),0:.01:1))/size(ok,2),'k')

figure
plot(.6:.01:1.6,cumsum(histc(TD_stim(:,1),.6:.01:1.6))/57)
hold on
plot(.6:.01:1.6,cumsum(histc(TD_stim(:,2),.6:.01:1.6))/57)

%%
usub = unique(stim_sub);

[~,b] = ismember(stim_sub,usub);

IED_mean_stim_12_all = squeeze(nansum(IED_rate_stim));
IED_mean_stim_12_all = (IED_mean_stim_12_all(2,:)-IED_mean_stim_12_all(1,:))./(IED_mean_stim_12_all(1,:));
close all
figure
col = linspecer(14,'jet');
for i = 1:14
    
    
    plot(TD_stim_mean(b==i),IED_mean_stim_12_all(b==i),'.','color',col{i},'markersize',30)
    hold on
    if i ==1
        plot([-.15 .25],[0 0],'k')
        plot([0 0],[-1 1],'k')
    end
end

%%
[~,b] = ismember(stim_sub,usub);
[~,b1] = ismember(con_sub,usub);
close all
col = flipud(linspecer(6,'jet'))
for i = 1:7
    
    plot(daysPost_stim(b==i),squeeze(sum(IED_rate_stim(:,2,b==i),1)) ,'x','color',col{i})
    hold on
    plot(daysPost_con(b1==i),squeeze(sum(IED_rate_control(:,2,b1==i),1)) ,'o','color',col{i})
end
plot([160 240],[0 0],'k')

figure

col = flipud(linspecer(6,'jet'))
for i = 1:6
    
    plot(daysPost_stim(b==i),squeeze(sum(IED_rate_stim(:,1,b==i),1)) ,'x','color',col{i})
    hold on
    plot(daysPost_con(b1==i),squeeze(sum(IED_rate_control(:,1,b1==i),1)) ,'o','color',col{i})
end
plot([160 240],[0 0],'k')

%%
[~,b_stim] = ismember(stim_sub,usub);
[~,b_con] = ismember(con_sub,usub);

N_stim = size(TD_stim,1);
N_con = size(TD_control,1);

IEDtot_stim = squeeze((nansum(IED_rate_stim(:,:,:),1)))';
IEDtot_con = squeeze((nansum(IED_rate_control(:,:,:),1)))';

%b_stim = repmat(b_stim,1,3);
condition_stim = repmat(1:3,N_stim,1);
session_stim = repmat([1:N_stim]',1,3);



%b_con = repmat(b_con,1,3);
condition_con = repmat(1:3,N_con,1);
session_con = repmat([1:N_con]',1,3);

IED = [IEDtot_stim(:,2) - IEDtot_stim(:,1);IEDtot_con(:,2)-IEDtot_con(:,1)];
TD = [TD_stim(:,2)-TD_stim(:,1);TD_control(:,2)-TD_control(:,1)];
%condition = categorical([condition_stim(:);condition_con(:)]);
subject = usub([b_stim(:);b_con(:)])';
%session = [session_stim(:);session_con(:)];
stim_on = (categorical([ones(numel(b_stim(:)),1);2*ones(numel(b_con(:)),1)]));


tbl = table(IED,TD,stim_on,subject);
lme = fitlme(tbl,'IED ~  -1+ stim_on*TD+(1/subject) ')


%%

%get mean seizure rate
topdir = 'R:\DGregg\NeuralData\EDS';
files2 = getAllExtFiles(topdir,'szr',1);
keepFiles = contains(files2,'sessiondata');
files2 = files2(keepFiles);
files = [files1;files2];

