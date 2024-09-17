
% get the sessions we want


%get all mat files
fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);


%only keep files with NE2 and Novel Env
kp = cellfun(@any,regexp(fils,'Novel Env')) & (cellfun(@any,regexp(fils,'DA2')) | cellfun(@any,regexp(fils,'DA3')));

fils = fils(kp);



%find directories of those files
[dirs] = fileparts(fils);
dirs =  unique(dirs);
%%
subj =[];
for i = 1:length(dirs)
[~,subjt] = fileparts(dirs{i});
subj_split = split(subjt, "-");
subjt = subj_split{1};
subj = [subj;{subjt}];
end
subj = strrep(subj,'NE2m3' ,'NE2h2');

usubjs = unique(subj);

%% initializing the variables of interest


neural =[]; %neural data
vel =[]; % velocity
acc = []; % acceleration
kp_Novel =[]; % in novel env?
sesID = []; % session ID info
subjID =[]; % subject ID
distEdg = []; % distance from edge of envir
timeFrmEntry = [];% time from transition
dayNum = []; % day number
inRear =[]; % is animal rearing?
faceEdge =[]; % is the animal facing the edge of the env?
time_rear =[];
ts_tot =[];
pos = [];
conID =[];
% loop over all directories
for i = 1:length(dirs)
    
    
    %change directory into each
    cd(dirs{i})
    
    %check if we have the sessiondata struct
    if exist('sessiondata.mat')
        
        
      
        
        %load struct with neural/behavioral data
        load('sessiondata.mat')
          %get times in home cage and in novel context
        ts_video = sessiondata.behavior.ts_video;
        
        %make sure the struct has been organize by context
        isfield(sessiondata,'contextEntry')
        
        
          %find which subject ID
        [~,ID] = ismember(subj{i},usubjs);
        
        
        %check if rearing has been scored
        if exist('rear.mat')
            vv = load('rear.mat');
            
            %define every moment as rear or not rear
            rears = [cell2mat(vv.data(contains(vv.data(:,1),'start','IgnoreCase',true),2)) cell2mat(vv.data(contains(vv.data(:,1),'end','IgnoreCase',true),2))];
            rears = MergeEpochs2(rears);
            inReart = double(InIntervals(sessiondata.behavior.ts_video,[rears(:,1) rears(:,2)]));
             time_reart =  sm_timeFromEvent(ts_video,rears(:,1));
        else
            
            inReart = nan(length(sessiondata.behavior.ts_video),1);
            time_reart =nan(length(sessiondata.behavior.ts_video),1);
        end
        
        
        %load distance to edge
        totalEdgDist = sessiondata.behavior.totalEdgDist;
        
        %get context transitions
        data = sessiondata.contextEntry;
        
        %get times for each context entry
        context_entry = cell2mat(data(:,2));
        
        %sort and get which index goes in which order (b)
        [context_entry,b] = sort(context_entry);
        
        %sort the edge/transition by time
        data = data(b,:);
        
        %define entry times
        epochs_on  = context_entry;
        
        %define exit times
        epochs_off = [context_entry(2:end); sessiondata.behavior.ts_video(end)];
        
        %build matrix of onsets and offsets
        epochs = [epochs_on epochs_off];
        
        %find epochs in home cage
        home = cellfun(@any,regexp(data(:,1),'home'));
        inNovel = find(~home);
        
        
      
       
        
        kp_homet = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(home,:),2),'uni',0);
        kp_homet = cell2mat(kp_homet');
        kp_homet = any(kp_homet,2);
        
        kp_Novelt = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(inNovel,:),2),'uni',0);
        kp_Novelt = cell2mat(kp_Novelt');
        kp_Novelt = any(kp_Novelt,2);
        
        %build dayNum
        novEp  = epochs(inNovel,:);
        dayNumt = 100*ones(size(ts_video));
        for jj = 1:size(novEp)
            in = InIntervals(ts_video,novEp(jj,:));
            dayNumt(in) = sessiondata.contextEntry{inNovel(jj),5};
        end
        dayNumt(dayNumt==100) = min(dayNumt);
        %get velocity
        velt = sessiondata.behavior.vel_cor;
        acct = sessiondata.behavior.acc_cor;
        distEdgt = mean(totalEdgDist,2);
        faceEdget = nanmean(totalEdgDist(:,1:2),2)-totalEdgDist(:,3);
        %                dayNumt = sessiondata.behavior.dayNum;
        timeFrmEntryt = sessiondata.behavior.timeFrmEntry;
        
        %getRear
        
        %downsample neural data to video clock
        neuralt = sessiondata.neural.signal_DFoF;
        fs = sessiondata.neural.fs_neural;
        
        
        k = gaussian2Dfilter([1 20*fs],3*fs);
        neuralt = nanconvn(neuralt,k);
        ts_neural = (1:length(neuralt))/fs;
        
        neuralt = double(interp1(ts_neural,neuralt,ts_video));
       
        
        
        % get position and context
        post = nanmean(cat(3,sessiondata.behavior.position.left_ear_cor,sessiondata.behavior.position.right_ear_cor),3);
        conIDt = sessiondata.behavior.position.context;
        conIDt = strrep(conIDt,' ' ,'');
        nSamples = length(neuralt);
        post = post(1:nSamples,:);
        conIDt = conIDt(1:nSamples,:);
        
        
        % concatenate each session
        neural = [neural;neuralt];
        acc = [acc;acct];
        vel = [vel;velt];
        distEdg = [distEdg; distEdgt];
        dayNum = [dayNum;dayNumt];
        pos = [pos;post];
        conID = [conID;conIDt];
        timeFrmEntry = [timeFrmEntry; timeFrmEntryt];
        %   dayNum = [dayNum; dayNumt];
        inRear = [inRear;inReart];
        time_rear =[time_rear;time_reart];
        kp_Novel = [kp_Novel;kp_Novelt];
        ts_tot = [ts_tot;ts_video(:)];
        faceEdge = [faceEdge;faceEdget];
        sesID = [sesID; i *ones(nSamples,1)];
        subjID = [subjID;ID*ones(nSamples,1)];
        i
        pwd
    end
    
end


%%
dayNum1 = dayNum+1;
acc1 = acc;
acc1(acc>10) = 10;
vel1 = vel;
vel1(vel>50) = 50;
ix = 1;
clear XX YY

intercept =true;

maxDistt = accumarray(sesID,distEdg,[],@max);
maxDist = distEdg;

for i = 1:max(sesID)
    
   maxDist(sesID ==i) = maxDistt(i);
end
distEdg1 = distEdg./maxDist;

%%


velDay = nan(10,8);
accDay = nan(10,8);
edgeDay = nan(10,8);
rearDay = nan(10,8);


for i = 1:8
   kp = kp_Novel==1 & subjID==i;%  & ts_tot>1200;
   rearDay(:,i) = avghist(dayNum1(kp),inRear(kp),1:10);
   velDay(:,i) = avghist(dayNum1(kp),vel1(kp),1:10);
   accDay(:,i) = avghist(dayNum1(kp),acc1(kp),1:10);
   edgeDay(:,i) = avghist(dayNum1(kp),distEdg1(kp),1:10);
  
   
   
end


%%

veltim = nan(8,31);
acctim = nan(8,31);
edgetim = nan(8,31);
reartim = nan(8,31);


for i = 1:8
   kp = kp_Novel==1 & subjID==i & dayNum1==1;%  & ts_tot>1200;
   reartim(i,:) = avghist(timeFrmEntry(kp),inRear(kp),1:10:310);
   veltim(i,:) = avghist(timeFrmEntry(kp),vel1(kp),1:10:310);
   acctim(i,:) = avghist(timeFrmEntry(kp),acc1(kp),1:10:310);
   edgetim(i,:) = avghist(timeFrmEntry(kp),distEdg1(kp),1:10:310);
  
   
   
end
edgetim = edgetim(:,1:end-1);
acctim = acctim(:,1:end-1);
veltim = veltim(:,1:end-1);
reartim = reartim(:,1:end-1);


%%

velNE = nan(8,31);
accNE = nan(8,31);
edgeNE = nan(8,37);

timNE = nan(8,31);


velBin = min(vel1):range(vel1)/30:max(vel1);
accBin = min(acc1):range(acc1)/30:max(acc1);
distBin = -.1:.5/30:.5;
rearBin = min(inRear):range(inRear)/30:max(inRear);
timBin =1:10:310;


for i = 1:8
   kp = kp_Novel==1 & subjID==i & dayNum1>1;%  & ts_tot>1200;
   velNE(i,:) = avghist(vel1(kp),neural(kp),velBin);
   accNE(i,:) = avghist(acc1(kp),neural(kp),accBin);
   edgeNE(i,:) = avghist(distEdg1(kp),neural(kp),distBin);
  
   timNE(i,:) = avghist(timeFrmEntry(kp),neural(kp),timBin);
   
   
end
velNE = velNE(:,1:end-1);
accNE = accNE(:,1:end-1);
edgeNE = edgeNE(:,1:end-1);
timNE = timNE(:,1:end-1);

clear rearNE rearNu
for i = 1:8
    ix = find(diff(inRear==1 &kp_Novel==1 & subjID==i )>0);
rearNE{i} = (neural(repmat(ix,1,1001) + repmat(-500:500,size(ix,1),1)));
rearNEu(i,:) = nanmean(neural(repmat(ix,1,1001) + repmat(-500:500,size(ix,1),1)));
end

%%



clear velNEr accNEr edgeNEr
for i = 1:8
   kp = kp_Novel==1 & subjID==i;%  & ts_tot>1200;
   velNEr(i,:) = polyfit(neural(kp),vel1(kp),1);
   accNEr(i,:) = polyfit(neural(kp),acc1(kp),1);
   
   kp = kp&~isnan(distEdg1);
   edgeNEr(i,:) = polyfit(neural(kp),distEdg1(kp),1);
  
   
   
   
end



%%

% get rate maps
uCon = unique(conID);
[~,~,c] = unique(conID);
clear NEbin

k = gaussian2Dfilter([100 100],[1 1]);
posBin = -20:.5:60;
for i = 1:length(uCon)
    kp = c==i;
    
    
    [occ,~,~,b] = histcn(pos(kp,:),posBin,posBin);
    
   % occ = nanconvn(occ,k,'nanout',true);
    NEbin{i} = accumarray(b+1,neural(kp),[length(posBin) length(posBin)],@nanmean,nan);
     NEbin{i} =  NEbin{i}(1:end-1,1:end-1);
    NEbin{i} = nanconvn(NEbin{i},k,'nanout',true);
    NEbin{i}(occ<.01) = nan;
end

%%

% stat1



subj = linearize(repmat(1:8,10,1));
day  = linearize(repmat((1:10)',1,8));
kp = day==1 | day ==2;
tbl = table(rearDay(:),accDay(:),velDay(:),edgeDay(:),subj,day,'VariableNames',{'rearDay','accDay','velDay','edgeDay','subj','day'});
lme = fitlme(tbl(kp,:),'edgeDay ~ day +(day|subj)');
%%
set(0, 'DefaultFigureRenderer', 'painters');

filename = 'E:\Dropbox\UNM\Papers\NE_timeConstants\Figures\FigureS1A_DA.eps';
%Figure S1a
h1=figure;
subplot(1,4,1)
boxplot(linearize(accDay(1:2,:)'),[ones(8,1);2*ones(8,1)])
hold on
plot(accDay(1:2,:))
title('Acc')
subplot(1,4,2)
boxplot(linearize(velDay(1:2,:)'),[ones(8,1);2*ones(8,1)])
hold on
plot(velDay(1:2,:))
title('Vel')
subplot(1,4,3)
boxplot(linearize(rearDay(1:2,:)'),[ones(8,1);2*ones(8,1)])
hold on
plot(rearDay(1:2,:))
title('Rear')
subplot(1,4,4)
boxplot(linearize(edgeDay(1:2,:)'),[ones(8,1);2*ones(8,1)])
hold on
plot(edgeDay(1:2,:))
title('Edge')
saveas(h1,filename,'epsc')


%%

%Figure S1b

filename = 'E:\Dropbox\UNM\Papers\NE_timeConstants\Figures\FigureS1B_DA.eps';
h2 = figure;
subplot(1,4,1)

plotMeanSEM(1:10:300,acctim,'k')
title('Acc')


subplot(1,4,2)
plotMeanSEM(1:10:300,veltim,'k')
title('Vel')

subplot(1,4,3)
plotMeanSEM(1:10:300,reartim,'k')
title('Rear')

subplot(1,4,4)
plotMeanSEM(1:10:300,edgetim,'k')
title('Edge')

saveas(h2,filename,'epsc')

%%
%stat2

[h,p,~,stat] = ttest(accNEr(:,1));
[h,p,~,stat] = ttest(velNEr(:,1));
[h,p,~,stat] = ttest(edgeNEr(:,1));


%%
%Figure 1A


close all

filename = 'E:\Dropbox\UNM\Papers\NE_timeConstants\Figures\Figure1A_DA.eps';
h3 = figure;



subplot(1,5,1)


plotMeanSEM(timBin(1:end-1),timNE,'k')

title('Time')


subplot(1,5,2)



plotMeanSEM(accBin(1:end-1),accNE,'k')
title('Acc')


subplot(1,5,3)
plotMeanSEM(velBin(1:end-1),velNE,'k')
title('Vel')

subplot(1,5,4)
Fs = 1./mode(diff(ts_tot));
plotMeanSEM((-500:500)/ Fs,rearNEu,'k')
title('Rear')

subplot(1,5,5)
plotMeanSEM(distBin(1:end-1),edgeNE,'k')
title('Edge')

saveas(h3,filename,'epsc')

%%

%Figure S2
close all

filename = 'E:\Dropbox\UNM\Papers\NE_timeConstants\Figures\FigureS2-DA.eps';

close all
h4 = figure;
for i = 1:9
   subplot(3,3,i)
    imagesc(NEbin{i},[0 1.5])
    
end
saveas(h4,filename,'epsc')

%%

%Figure S3, stat3


for j = 1:8
    
    kp = kp_Novel==1 & subjID==j ;
    
 [rr(:,:,i),p(:,:,i)] =    corr( [vel1(kp) acc1(kp) distEdg1(kp) time_rear(kp) timeFrmEntry(kp)],'rows','pairwise');
end
%%


b_all_NE =[];MSE_full =[];MSE_red=[];b_full=[];pred_tim=[];

    for j = 1:max(subjID)
        
        kp = kp_Novel==1 & subjID==j ;%& ts_tot>1000;
        if any(kp) 
        %cross val
        dayN = dayNum(kp);
        
        kp_train = mod(dayN,2)==0;
       % kp_train = true(sum(kp),1);
        inter = ones(sum(kp),1);
        x1 = timeFrmEntry(kp);
        x2 =  [vel1(kp) acc1(kp) distEdg1(kp) ];
        x3 = time_rear(kp);
        y = neural(kp);
        IV = [inter x2 x1 ];
      %  IV =[inter x1];
        DV = neural(kp);
        isTime = false(5,1);
        isTime([5 ]) = true;
        
        %isTime = [false true];
        [MSE_full(j), MSE_red(j,:),b_all_NE(j,:),pred] = sm_model_comp_MSE(IV,DV,isTime,kp_train,true);
        
      %  pred_tim(i,:,j) = avghist(x1,pred,0:600);
        else
             MSE_full(j) = nan;
             MSE_red(j,:) = nan;
             b_all_NE_day(j,:) = nan;
            
    end
    end


%%
nVar = length(isTime);
nSubj = length(MSE_full);

ok=((MSE_red - MSE_full'))./MSE_full';
close all
figure
bar(100*nanmean(ok),'facecolor','w')
hold on

plot(repmat(1:nVar,nSubj,1)'+(rand(nVar,nSubj)-.5)/5,100*ok','o')
set(gca,'xtick',1:5,'xticklabel',{'Intercept','velocity','acceleration','distance2dge','time from transfer'})

%%


b_all_HC =[];MSE_full_HC =[];MSE_red_HC=[];pred_tim=[];

    for j = 1:6
        
        kp = kp_Novel==0 & subjID==j &ts_tot>1000;
        if any(kp) 
        %cross val
        dayN = dayNum(kp);
        
        kp_train = mod(dayN,2)==1;
       % kp_train = true(sum(kp),1);
        inter = ones(sum(kp),1);
        x1 = timeFrmEntry(kp);
        x2 =  [vel1(kp) acc1(kp) distEdg1(kp) ];
        x3 = time_rear(kp);
        y = neural(kp);
        IV = [inter x2 x1];
      %  IV =[inter x1];
        DV = neural(kp);
        isTime = false(5,1);
        isTime([5 ]) = true;
        
        %isTime = [false true];
        [MSE_full_HC(j), MSE_red_HC(j,:),b_all_HC(j,:),pred] = sm_model_comp_MSE(IV,DV,isTime,kp_train,true);
        
      %  pred_tim(i,:,j) = avghist(x1,pred,0:600);
        else
             MSE_full_HC(j) = nan;
             MSE_red_HC(j,:) = nan;
             b_all_HC(j,:) = nan;
            
    end
    end


%%
nVar = length(isTime);
nSubj = length(MSE_full);
ok=((MSE_red_HC - MSE_full_HC'))./MSE_full_HC';
close all
figure
bar(100*nanmean(ok),'facecolor','w')
hold on

plot(repmat(1:nVar,nSubj,1)'+(rand(nVar,nSubj)-.5)/5,100*ok','o')


%%
b_all_NE_day =[];pred_tim_day=[];
for i = 1:10
    for j = 1:max(subjID)
        
        kp = kp_Novel==1 & dayNum1==i & subjID==j ;%& ts_tot>1000;
        if any(kp) 
        %cross val
        dayN = dayNum1(kp);
        
        kp_train = mod(dayN,2)==1;
        kp_train = true(sum(kp),1);
        inter = ones(sum(kp),1);
        x1 = timeFrmEntry(kp);
        x2 =  [vel1(kp) acc1(kp) distEdg(kp) ];
        x3 = time_rear(kp);
        y = neural(kp);
        IV = [inter x2 x1 ];
        %IV =[ x1];
        DV = neural(kp);
        isTime = false(1,5);
        isTime([5 ]) = true;
        
        %isTime = [false true];
        [~,~,b_all_NE_day(i,:,j),pred,MSE_full_tim(i,j)] = sm_model_comp_MSE(IV,DV,isTime,kp_train,false);
        
        pred_tim_day(i,:,j) = avghist(x1,pred,0:600);
        else
             MSE_full(i) = nan;
             MSE_red(i,:) = nan;
             b_all_NE_day(i,:,j) = nan;
            
    end
    end
end
%%
amp = linearize(squeeze(b_all_NE_day(:,5,:)));
tau = linearize(squeeze(b_all_NE_day(:,6,:)));
amp1 = linearize(squeeze(b_all_NE_day(:,7,:)));
tau1 = linearize(squeeze(b_all_NE_day(:,8,:)));

subj = linearize(repmat(1:nSubj,10,1));
day  = linearize(repmat((1:10)',1,nSubj));
kp = day <=10;
tbl = table(amp,log10(-tau),amp1,(-tau1),subj,day,'VariableNames',{'amp','tau','amp1','tau1','subj','day'});
lme = fitlme(tbl(kp,:),'amp ~ day + (1+day|subj)')

%%

%%
%close all
figure
tau = linearize(squeeze(b_all_NE_day(:,6,:)));
ok = [squeeze(((b_all_NE_day(:,6,:))))' b_all_HC(:,6)];
taut = [tau;  b_all_HC(:,6)];
day1 = [day; 11*ones(nSubj,1)];
semilogy(1,1)
hold on
boxplot((-taut),day1)
hold on
ylim([.001 .12])

figure
tau = linearize(squeeze(b_all_NE_day(:,7,:)));
ok1 = [squeeze(((b_all_NE_day(:,7,:))))' b_all_HC(:,7)];
taut = [tau;  b_all_HC(:,7)];
day1 = [day; 11*ones(nSubj,1)];

hold on
boxplot((taut),day1)
hold on


%ylim([.0001 .012])

%%






b_all_HC =[];MSE_full =[];MSE_red=[];b_full=[];pred_tim_HC=[];

    for j = 1:nSubj
        
        kp = kp_Novel==0 & subjID==j & ts_tot>1000;
        if any(kp) 
        %cross val
        dayN = dayNum(kp);
        
        kp_train = mod(dayN,2)==1;
        kp_train = true(sum(kp),1);
        inter = ones(sum(kp),1);
        x1 = timeFrmEntry(kp);
        x2 =  [vel1(kp) acc1(kp) distEdg(kp) ];
        x3 = time_rear(kp);
        y = neural(kp);
        IV = [inter x2 x1 x3];
      %  IV =[inter x1];
        DV = neural(kp);
        isTime = false(5,1);
        isTime([5 ]) = true;
        
        %isTime = [false true];
        [~,~,b_all_HC(j,:),pred] = sm_model_comp_MSE(IV,DV,isTime,kp_train,false);
        
        pred_tim_HC(j,:) = avghist(x1,pred,0:600);
        else
             MSE_full(j) = nan;
             MSE_red(j) = nan;
             b_all_HC(j,:) = nan;
            
    end
    end

%%

col = flipud(linspecer(10,'jet'));
col = [col ;{[0 0 0]}];
close all
figure
for i = 1:10
plot(nanmean(pred_tim_day(i,:,:),3),'color',col{i})
hold on
end

plot(nanmean(pred_tim_HC,1),'color','k')
colormap(cell2mat(col))
colorbar

xlim([0 500])
xlabel('Seconds after transition')
ylabel('NE')
figure

for i = 1:10
    kp = kp_Novel==1 & dayNum1==i ;
plot(1:5:600,avghist(timeFrmEntry(kp),neural(kp),1:5:600,@nanmean),'color',col{i})
hold on
end
    kp = kp_Novel==0 & ts_tot>1000 ;
plot(1:5:600,avghist(timeFrmEntry(kp),neural(kp),1:5:600),'color','k')
colormap(cell2mat(col))
colorbar

xlim([0 500])
xlabel('Seconds after transition')
ylabel('NE')

%%

figure
for i = 1:10
    u1 = nanmean(squeeze(b_all_NE_day(i,5,:)));
       u2 = nanmean(squeeze(b_all_NE_day(i,6,:)));
        s1 = SEM(squeeze(b_all_NE_day(i,5,:)));
       s2 = SEM(squeeze(b_all_NE_day(i,6,:)));
       
       
        errorbar(u1,u2,s2,s2,s1,s1,'color',col{i})
hold on
        
end

    u1 = nanmean(b_all_HC(:,5,:));
       u2 = nanmean(b_all_HC(:,6,:));
        s1 = SEM(b_all_HC(:,5,:));
       s2 = SEM(b_all_HC(:,6,:));
       
       
        errorbar(u1,u2,s2,s2,s1,s1,'color','k')
        
colormap(cell2mat(col))
colorbar
%%
% %%
% %%
% pred = nan(size(neural));
% for i = 1:8
%     kp = kp_Novel==1 & subjID==i;
% 
% inter = ones(sum(kp),1);
% x1 = timeFrmEntry(kp);
% x2 =  [inter vel1(kp) acc1(kp) distEdg(kp) ];
% x3 = time_rear(kp);
% y = neural(kp);
% IV = [x2 x1 x3];
% %IV = [inter x1];
% DV = neural(kp);
% 
%     pred(kp) = b_full(i,5)*exp(b_full(i,6)*IV(:,5)) + ...
%     b_full(i,7)*exp(b_full(i,8)*IV(:,6));
% end
% %%
% 
% 
% 
% clear b_all_NE_day
% 
% for j=1:10
% for i = 1:8
% 
% kp = kp_Novel==1 & subjID==i & dayNum==j;
% if any(kp)
% inter = ones(sum(kp),1);
% x1 = timeFrmEntry(kp);
% x2 =  [inter vel1(kp) acc1(kp) distEdg(kp) ];
% y = neural(kp);
% fitfun = @(b) b(1).*exp(b(2).*x1)+  b(3).*exp(b(4).*x1) + x2*b(5:end)';
% f = @(b) nanmean((y - fitfun(b)).^2 ) + 0.001*sqrt(b(1)^2 + b(3)^2 + sum(b(5:end).^2)); 
% 
% x0 = [0 -.05 0 -.05 0 0 0 0 ];
% Aeq =[];
% beq =[];
% A =[];
% b = [];
% lb = [0 -.1  -20 -.1 -2*ones(1,4)];
% ub = [20 -.001 0 -.001 2*ones(1,4)];
% b = fmincon(f,x0,A,b,Aeq,beq,lb,ub)% Objective Function
% b_all_NE_day{j}(i,:) = b;
% else
%     b_all_NE_day{j}(i,:) = nan(1,8);
% end
% end
% end
% %%
% 
% 
% b_all_HC =[];
% for i = 1:8
% 
% kp = kp_Novel==0 & subjID==i;
% inter = ones(sum(kp),1);
% x1 = timeFrmEntry(kp);
% x2 =  [inter vel1(kp) acc1(kp) distEdg(kp) ];
% y = neural(kp);
% fitfun = @(b) b(1).*exp(b(2).*x1)+  b(3).*exp(b(4).*x1) + x2*b(5:end)';
% f = @(b) nanmean((y - fitfun(b)).^2 ) + 0.001*sqrt(b(1)^2 + b(3)^2 + sum(b(5:end).^2)); 
% 
% x0 = [0 -.05 0 -.05 0 0 0 0 ];
% Aeq =[];
% beq =[];
% A =[];
% b = [];
% lb = [0 -.1  -20 -.1 -2*ones(1,4)];
% ub = [20 -.001 0 -.001 2*ones(1,4)];
% b = fmincon(f,x0,A,b,Aeq,beq,lb,ub)% Objective Function
% b_all_HC(i,:) = b;
% end
% 
% 
% 
% %%
% idx = 1;
% clear beta_NE
% 
% %loop over time constants
% for t = fliplr(-logspace(log10(.0001),log10(.1),30))
%     ix = 1;
%     clear YY XX
%     
%     %loop over subjects
%     for i = 1:max(subjID)
%         
%         % only keep data for each subj in novel env
%         kp = subjID==i & kp_Novel==1;
%         if any(kp)
%             YY{ix} = neural(kp);
%             
%             % usually large negative value are a sign of mechanical instability
%             YY{ix}( YY{ix}<-5) = nan;
%            
%             %build our design matrix of IVs: time vel acc and distance
%             XX{ix} = [  exp(t.*timeFrmEntry(kp)) vel(kp)  acc1(kp)  distEdg(kp)];
%            
%             ix  = ix+1;
%             
%         end
%     end
%     
%     %fit the model
%     
%     if intercept
%         stats_NE = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc','dist'});
%         beta_NE(:,idx) = stats_NE.first_level.beta(2,:);
%         for i =1:length(XX)
%             
%             % find the MSE for each time constant
%             mse_NE(i,idx) = nanmean(sqrt((YY{i} - [ones(size(XX{i},1),1) XX{i}]* stats_NE.first_level.beta(:,i)).^2));
%         end
%     else
%         
%         stats_NE = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc','dist'},'noint');
%         beta_NE(:,idx) = stats_NE.first_level.beta(1,:);
%         for i =1:length(XX)
%             
%             % find the MSE for each time constant
%             mse_NE(i,idx) = nanmean(sqrt((YY{i} - [ XX{i}]* stats_NE.first_level.beta(:,i)).^2));
%         end
%     end
%     idx=  idx+1
% end
% 
% 
% %%
% ok =  fliplr(-logspace(log10(.0001),log10(.1),30));
% 
% % find the lowest MSE per subject
% [a,b] = min(sum(mse_NE));
% [~,ix] = min(mse_NE,[],2);
% 
% %the best model for each subject  = the best timeconstant
% tau_hat_NE = ok(ix);
% beta_NE = beta_NE(sub2ind(size(beta_NE),1:size(beta_NE,1),ix'));
% 
% all_subj = usubjs(unique(subjID));
% 
% %%
% 
% idx = 1;
% clear beta_HC
% for t = fliplr(-logspace(log10(.0001),log10(.1),30))
%     ix = 1;
%     clear YY XX
% for i = 1:10
%     kp = subjID==i & kp_Novel==0; 
%    if any(kp) && ~all(isnan(inRear(kp)))
%     YY{ix} = neural(kp);
%    YY{ix}( YY{ix}<-5) = nan;
%     %XX{ix} = [  exp(t.*timeFrmEntry(kp)) ];
%     
%     XX{ix} = [  exp(t.*timeFrmEntry(kp)) vel(kp)  acc1(kp)  distEdg(kp) inRear(kp)];
%  %  D = dummyvar(kp_Novel(kp)+1);
%   %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
%     ix  = ix+1;
%    
%    end
% end
% 
% if intercept
%     
%     stats_HC = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc','dist','rear'});
%     beta_HC(:,idx) = stats_HC.first_level.beta(2,:);
%     for i =1:length(XX)
%         mse_HC(i,idx) = nanmean(sqrt((YY{i} - [ones(size(XX{i},1),1) XX{i}]* stats_HC.first_level.beta(:,i)).^2));
%     end
%     
% else
%     
%     stats_HC = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc','dist','rear'},'noint');
%     beta_HC(:,idx) = stats_HC.first_level.beta(1,:);
%     for i =1:length(XX)
%         mse_HC(i,idx) = nanmean(sqrt((YY{i} - [ XX{i}]* stats_HC.first_level.beta(:,i)).^2));
%     end
%     
% end
% idx=  idx+1
% end
% 
% 
% %%
% ok =  fliplr(-logspace(log10(.0001),log10(.1),30));
% [a,b] = min(sum(mse_HC));
% [~,ix] = min(mse_HC,[],2);
% tau_hat_HC = ok(ix);
% beta_HC = beta_HC(sub2ind(size(beta_HC),1:size(beta_HC,1),ix'));
% 
% all_subj = usubjs(unique(subjID));
% 
% %%
% % do per day
% tau_hat_NE_days = nan(9,10);
% 
% 
% beta_NE_days = nan(9,10);
% 
% 
% for day = 0:9
% 
% 
% idx = 1;
% clear beta_NEt mse_NE
% for t = fliplr(-logspace(log10(.0001),log10(.1),30))
%     ix = 1;
%     clear YY XX
% for i = 1:10
%     kp = subjID==i & kp_Novel==1 & dayNum==day; 
%    if any(kp) 
%     YY{ix} = neural(kp);
%    YY{ix}( YY{ix}<-5) = nan;
%     %XX{ix} = [  exp(t.*timeFrmEntry(kp)) ];
%     
%     XX{ix} = [  exp(t.*timeFrmEntry(kp)) vel(kp)  acc1(kp)  distEdg(kp) ];
%  %  D = dummyvar(kp_Novel(kp)+1);
%   %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
%     ix  = ix+1;
%    
%    end
% end
% stats_NE = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc','dist'},'noint');
% beta_NEt(:,idx) = stats_NE.first_level.beta(1,:);
% for i =1:length(XX)
% mse_NE(i,idx) = nanmean(sqrt((YY{i} - [ XX{i}]* stats_NE.first_level.beta(:,i)).^2));
% end
% idx=  idx+1
% end
% 
% 
% 
% ok =  fliplr(-logspace(log10(.0001),log10(.1),30));
% [a,b] = min(sum(mse_NE));
% [~,ix] = min(mse_NE,[],2);
% tau_hat_NE_days(unique(subjID(dayNum==day)),day+1) = ok(ix);
% beta_NE_days(unique(subjID(dayNum==day)),day+1) = beta_NEt(sub2ind(size(beta_NEt),1:size(beta_NEt,1),ix'));
% end
% 
% 



% 
% 
% %%
% ix=1;
% clear YY XX
% 
% 
% for i = 1:10
%     kp = subjID==i & kp_Novel==1; 
%    if any(kp)
%     YY{ix} = neural(kp);
%   YY{ix}( YY{ix}<-5) = nan;
%     XX{ix} = [  exp(tau_hat_NE(ix).*timeFrmEntry(kp)) ];
%    
%  %  D = dummyvar(kp_Novel(kp)+1);
%   %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
%     ix  = ix+1;
%    
%    end
% end
% stats = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time'},'noint');
% neural1 = neural;
% timeEst = neural;
% ix=1;
% for i = 1:10
%     kp = subjID==i & kp_Novel==1; 
%    if any(kp) && ~all(isnan(inRear(kp)))
%        
%        Y_hat = stats.first_level.beta(1,ix).*exp(tau_hat_NE(ix).*timeFrmEntry(kp));
%        neural1(kp) = neural(kp) - Y_hat;
%        timeEst(kp) = Y_hat;
%        ix=1+ix;
%    end
% end
% 
% %%
% neural1 = neural;
% neural1(abs(neural1)>5) = nan;
% clear XX YY
% ix=1;
% for i =  1:10
%     kp = subjID==i & kp_Novel==1 & timeFrmEntry>10; 
%    if any(kp) && ~all(isnan(inRear(kp))) 
%     YY{ix} = neural1(kp);
%       
%     XX{ix} = [ vel1(kp) acc1(kp) distEdg1(kp) exp(tau_hat_NE(ix).*timeFrmEntry(kp)) inRear(kp) ];
%       %  XX{ix} = [ vel(kp) acc(kp) distEdg(kp)  inRear(kp)];
%   
%    D = dummyvar(kp_Novel(kp)+1);
%   %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
%     ix  = ix+1;
%    end
% end
% stats = glmfit_multilevel(YY, XX, [], 'verbose', 'weighted','names',{'vel','Acc','dist','time','intRear'},'noint');
% %stats = glmfit_multilevel(YY, XX, [], 'verbose', 'weighted','names',{'intAcc','intRear','intDist','intTime'});
% 
% for i =1:length(XX)
% mse(i) = nanmean(sqrt((YY{i} - [ XX{i}]* stats.first_level.beta(:,i)).^2));
% end
% Est = neural;
% ix=1;
% 
% %%
% for i = 1:10
%     kp = subjID==i & kp_Novel==1 & timeFrmEntry>10; 
%    if any(kp) && ~all(isnan(inRear(kp)))
%        
%        Y_hat = stats.first_level.beta(1,ix) + ...
%           stats.first_level.beta(2,ix).*vel(kp) + ...
%           stats.first_level.beta(3,ix).*acc(kp) + ...
%           stats.first_level.beta(4,ix).*distEdg(kp);
% 
%      
%       Est(kp) = Y_hat;
%        ix=1+ix;
%    end
% end
% %%
% 
% clear ru con
% for i = 1:10
%     kp = kp_Novel==1 & subjID==i;
% idx=  find(diff(inRear(kp))>0);
% ok = neural(kp);
% idx = repmat(idx,1,2001)+repmat(-1000:1000,length(idx),1);
% kp = all(idx>0 &idx<=numel(ok),2);
% idx = idx(kp,:);
% 
% ru(i,:) = nanmean(ok(idx));
%  kp = kp_Novel==1 & subjID==i;
% ok = Est(kp);
% con(i,:) = nanmean(ok(idx));
% end
