fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);

kp = cellfun(@any,regexp(fils,'SOR')) & cellfun(@any,regexp(fils,'NE2'));

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

ts_decay = -.3;
neural =[];
vel =[];
acc = [];
kp_Novel =[];
sesID = [];
ts_tot =[];
ts_trans =[];
subjID =[];
distEdg = [];
timeFrmEntry = [];
dayNum = [];
inSample =[];
faceEdge =[];
time2rew = [];
posout =[];
inObjBlock = [];
posin =[];
for i = 1:length(dirs)
        
    cd(dirs{i})
 
    
    if exist('sessiondata.mat') && exist('novelObject.mat')
        load('sessiondata.mat')
        v=load('novelObjectIntro.mat');
        vv= load('novelObject.mat');
        obj = cell2mat(v.data(:,2));
        [~,ID] = ismember(subj{i},usubjs);
        
        
        if isfield(sessiondata,'contextEntry') 
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
            inNovel = ~home;
            
            
            %get times in home cage and in novel context
            ts_video = sessiondata.behavior.ts_video;
            
            %get whether sampling Object
            in = InIntervals(ts_video, MergeEpochs2([cell2mat(vv.data(:,2)) cell2mat(vv.data(:,2))+1]));
            
            
            kp_homet = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(home,:),2),'uni',0);
            kp_homet = cell2mat(kp_homet');
            kp_homet = any(kp_homet,2);
            
            kp_Novelt = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(inNovel,:),2),'uni',0);
            kp_Novelt = cell2mat(kp_Novelt');
            kp_Novelt = any(kp_Novelt,2);
            
            %get velocity
            velt = sessiondata.behavior.vel_cor;
            acct = sessiondata.behavior.acc_cor;
            
            if length(in) ~=length(velt)
                error('here')
            end
          
            %                dayNumt = sessiondata.behavior.dayNum;
            timeFrmEntryt = sessiondata.behavior.timeFrmObj;
            distEdgt = nanmean(sessiondata.behavior.totalEdgDist(:,1:2),2);
            %getRear
            
            %downsample neural data to video clock
            neuralt = sessiondata.neural.signal_DFoF;
            fs = sessiondata.neural.fs_neural;
            
            
            k = gaussian2Dfilter([1 10*fs],fs);
            neuralt = nanconvn(neuralt,k);
            ts_neural = (1:length(neuralt))/fs;
            
            neuralt = double(interp1(ts_neural,neuralt,ts_video));
            [~,inObjBlockt] = histc(ts_video,[obj ;obj(end)+300]);
            nSamples = length(neuralt);
            % concatenate each session
            neural = [neural;neuralt];
            acc = [acc;acct];
            vel = [vel;velt];
            distEdg = [distEdg; distEdgt];
            inSample = [inSample;in];
            inObjBlock = [inObjBlock;inObjBlockt];
            timeFrmEntry = [timeFrmEntry; timeFrmEntryt];
            %   dayNum = [dayNum; dayNumt];
          ts_tot = [ts_tot;ts_video(:)];
            kp_Novel = [kp_Novel;kp_Novelt];
    
          
            sesID = [sesID; i *ones(nSamples,1)];
            subjID = [subjID;ID*ones(nSamples,1)];
            
            i
            pwd
        end
        
        
  end
end
%%
acc1 = acc;
acc1(acc>10) = 10;
vel1 = vel;
vel1(vel>50) = 50;
ix = 1;
clear XX YY
distEdg1 = distEdg;
distEdg1(distEdg1<0) = -1;
intercept = true;
%%



b_all_SOR =[];MSE_full =[];MSE_red=[];pred_tim_SOR=[];real_tim_SOR =[];

    for j = 1:7
        
        kp = kp_Novel==1 & subjID==j;% & ts_tot>1000;
        if any(kp) 
        %cross val
       
        
        kp_train = ts_tot(kp)>1200;
       % kp_train = true(sum(kp),1);
        inter = ones(sum(kp),1);
        x1 = timeFrmEntry(kp);
        x2 =  [inter vel(kp)  acc1(kp)  distEdg(kp) inSample(kp) ];
       
        y = neural(kp);
        IV = [x2 x1 ];
      %  IV =[inter x1];
        DV = neural(kp);
        isTime = false(6,1);
        isTime([6]) = true;
        
        %isTime = [false true];
        [MSE_full(j), MSE_red(j,:)] = sm_model_comp_MSE(IV,DV,isTime,kp_train,true);
        
     
        else
            
             MSE_full(j) = nan;
             MSE_red(j,:) = nan(1,6);
            
    end
    end

    
    %%
    
ok=((MSE_red - MSE_full'))./MSE_full';
close all
figure
bar(100*nanmean(ok),'facecolor','w')
hold on

plot(repmat(1:6,7,1)'+(rand(6,7)-.5)/5,100*ok','o')

    set(gca,'xtick',1:6,'xticklabel',{'Intercept','velocity','acceleration','distance2dge','sampling object','time from object'})

%%


b_all_SOR =[];MSE_full =[];MSE_red=[];pred_tim_SOR=[];real_tim_SOR =[];

    for j = 1:7
        
        kp = kp_Novel==1 & subjID==j;% & ts_tot>1000;
        if any(kp) 
        %cross val
       
        
       % kp_train = mod(sesID(kp),2)==0;
        kp_train = true(sum(kp),1);
        inter = ones(sum(kp),1);
        x1 = timeFrmEntry(kp);
        x2 =  [inter vel(kp)  acc1(kp)  distEdg(kp) inSample(kp) ];
       
        y = neural(kp);
        IV = [x2 x1 ];
      %  IV =[inter x1];
        DV = neural(kp);
        isTime = false(6,1);
        isTime([6]) = true;
        
        %isTime = [false true];
        [~, ~,b_all_SOR(j,:),pred] = sm_model_comp_MSE(IV,DV,isTime,kp_train,false);
        
        pred_tim_SOR(j,:) = avghist(x1,pred,0:300);
         real_tim_SOR(j,:) = avghist(x1,y,0:300);
        else
             pred_tim_SOR(j,:) = nan(1,301);
             real_tim_SOR(j,:) = nan(1,301);
             b_all_SOR(j,:) = nan(1,9);
            
    end
    end

    
    

%%
figure
plotMeanSEM(0:300,pred_tim_SOR,'k')
plotMeanSEM(0:300,real_tim_SOR,'r')


%%
uiopen('E:\Dropbox\UNM\Papers\NE_timeConstants\Figures\delta_Tau_LT_NE.fig',1)
hold on
      u1 = nanmean(b_all_SOR(:,6,:));
       u2 = nanmean(b_all_SOR(:,7,:));
        s1 = SEM(b_all_SOR(:,6,:));
       s2 = SEM(b_all_SOR(:,7,:));
       
       
        errorbar(u1,u2,s2,s2,s1,s1,'color','k')
        ylim([-.1 0])


%%

b_all_SOR_obj =[];pred_tim_obj=[];real_tim_obj=[];MSE_red =[];
for i = 1:5
    for j = 1:7
        
        kp = kp_Novel==1 & inObjBlock==i & subjID==j ;%& ts_tot>1000;
        if any(kp)
            
            kp_train = true(sum(kp),1);
            inter = ones(sum(kp),1);
            x1 = timeFrmEntry(kp);
            x2 =  [inter vel(kp)  acc1(kp)  distEdg(kp) inSample(kp) ];
            
            y = neural(kp);
            IV = [x2 x1 ];
            %  IV =[inter x1];
            DV = neural(kp);
            isTime = false(6,1);
            isTime([6]) = true;
            %isTime = [false true];
            [~,~,b_all_SOR_obj(i,:,j),pred] = sm_model_comp_MSE(IV,DV,isTime,kp_train,false);
            
            pred_tim_obj(i,:,j) = avghist(x1,pred,0:300);
            real_tim_obj(i,:,j) = avghist(x1,y,0:300);
        else
            MSE_full(i) = nan;
            MSE_red(i,:) = nan(1,9);
             b_all_SOR_obj(i,:,j) = nan(1,9);
            real_tim_obj(i,:,j) = nan(1,301);
            pred_tim_obj(i,:,j) = nan(1,301);
    end
    end
end

amp = linearize(squeeze(b_all_SOR_obj(:,6,:)));
tau = linearize(squeeze(b_all_SOR_obj(:,7,:)));
subj = linearize(repmat(1:7,5,1));
obj  = linearize(repmat((1:5)',1,7));
kp = obj<5;
tbl = table(amp,log10(-tau),subj,obj,'VariableNames',{'amp','tau','subj','day'});
lme = fitlme(tbl(kp,:),'tau ~ day +(1/subj)+ (-1+day|subj)');


%%
close all
col = flipud(linspecer(5,'jet'));

for i = 1:5
    plotMeanSEM(0:300,squeeze(pred_tim_obj(i,:,:))',col{i})
    hold on
end
ylim([-1 2])
colormap(cell2mat(col))
cbh = colorbar;
xlabel('Time from object intro (s)')
ylabel('NE')

 cbh.Ticks = linspace(0, 1, 5) ; %Create 8 ticks from zero to 1
 cbh.TickLabels = num2cell(1:5) ;    
%%

close all
figure
for i = 1:5
    u1 = nanmean(squeeze(b_all_SOR_obj(i,6,:)));
       u2 = nanmean(squeeze(b_all_SOR_obj(i,7,:)));
        s1 = SEM(squeeze(b_all_SOR_obj(i,6,:)));
       s2 = SEM(squeeze(b_all_SOR_obj(i,7,:)));
       
       
        errorbar(u1,u2,s2,s2,s1,s1,'color',col{i})
hold on
        
end
colormap(cell2mat(col))
cbh = colorbar;
xlabel('Amp')
ylabel('decay')

 cbh.Ticks = linspace(0, 1, 5) ; %Create 8 ticks from zero to 1
 cbh.TickLabels = num2cell(1:5) ;   

%%
ok =  fliplr(-logspace(log10(.0001),log10(.1),30));
[a,b] = min(sum(mse));

[~,ix] = min(mse,[],2);
tau_hat = ok(ix);
beta_SOR = beta_SOR(sub2ind(size(beta_SOR),1:size(beta_SOR,1),ix'));

%%
all_subj = usubjs(unique(subjID));


tau_hat = ok(b);
ix=1;
clear YY XX


for i = 1:10
    kp = subjID==i & kp_Novel==1; 
   if any(kp) 
    YY{ix} = neural(kp);
  YY{ix}( YY{ix}<-5) = nan;
    XX{ix} = [  exp(tau_hat.*timeFrmEntry(kp)) ];
   
 %  D = dummyvar(kp_Novel(kp)+1);
  %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
    ix  = ix+1;
   
   end
end
stats = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time'});
neural1 = neural;
timeEst = neural;
ix=1;
for i = 1:10
    kp = subjID==i & kp_Novel==1; 
   if any(kp) 
       
       Y_hat = stats.first_level.beta(1,ix) +stats.first_level.beta(2,ix).*exp(tau_hat.*timeFrmEntry(kp));
       neural1(kp) = neural(kp) - Y_hat;
       timeEst(kp) = Y_hat;
       ix=1+ix;
   end
end

%%
neural1(neural1<-5) = nan;
clear XX YY
ix=1;
for i =  1:10
    kp = subjID==i & kp_Novel==1 & timeFrmEntry>10; 
   if any(kp) 
    YY{ix} = neural1(kp);
    time = exp(tau_hat.*timeFrmEntry);
   
   % XX{ix} = [ vel(kp) acc(kp) distEdg1(kp) exp(tau_hat.*timeFrmEntry(kp)) inRear(kp) ];
        XX{ix} = [ vel(kp) acc(kp) distEdg(kp) time2rew(kp) ];
  
   D = dummyvar(kp_Novel(kp)+1);
  %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
    ix  = ix+1;
   end
end
stats = glmfit_multilevel(YY, XX, [], 'verbose', 'weighted','names',{'vel','Acc','dist','time2rew'});
%stats = glmfit_multilevel(YY, XX, [], 'verbose', 'weighted','names',{'intAcc','intRear','intDist','intTime'});
clear mse
for i =1:length(XX)
mse(i) = nanmean(sqrt((YY{i} - [ones(size(XX{i},1),1) XX{i}]* stats.first_level.beta(:,i)).^2));
end
Est = neural;
ix=1;
for i = 1:6
    kp = subjID==i & kp_Novel==1 & timeFrmEntry>10; 
   if any(kp) 
       
       Y_hat = stats.first_level.beta(1,ix) + ...
          stats.first_level.beta(2,ix).*vel(kp) + ...
          stats.first_level.beta(3,ix).*acc(kp) + ...
          stats.first_level.beta(4,ix).*distEdg(kp);

     
      Est(kp) = Y_hat;
       ix=1+ix;
   end
end
%%

clear ru con
for i = 1:10
    kp = kp_Novel==1 & subjID==i;
idx=  find(diff(inSample(kp))>0);
ok = neural(kp);
idx = repmat(idx,1,2001)+repmat(-1000:1000,length(idx),1);
kp = all(idx>0 &idx<=numel(ok),2);
idx = idx(kp,:);

ru(i,:) = nanmean(ok(idx));
 kp = kp_Novel==1 & subjID==i;
ok = Est(kp);
con(i,:) = nanmean(ok(idx));
end
