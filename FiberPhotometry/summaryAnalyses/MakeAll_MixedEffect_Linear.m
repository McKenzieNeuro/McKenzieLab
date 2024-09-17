fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);

kp = cellfun(@any,regexp(fils,'Linear')) & cellfun(@any,regexp(fils,'NE2'));

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
ts_trans =[];
subjID =[];
distEdg = [];
timeFrmEntry = [];
dayNum = [];
inSample =[];
faceEdge =[];
time2rew = [];
timefromrew = [];
posout =[];
posin =[];
posH =[];
posB =[];
ts_tot =[];

for i = 1:length(dirs)
        
    cd(dirs{i})
    
    
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        if size(sessiondata.contextEntry,1)==3
        [~,ID] = ismember(subj{i},usubjs);
        nsamples = length(sessiondata.behavior.ts_video);
        %load struct with neural/behavioral data
        linpos  =nanmean([sessiondata.behavior.position.left_ear_cor(:,1) sessiondata.behavior.position.right_ear_cor(:,1)],2);
        linpos = linpos(1:nsamples);
        
        k  = gaussian2Dfilter([10000 1],10);
        inbound  = nanconvn(diff(linpos),k)<0;
        outbound  = nanconvn(diff(linpos),k)>0;
        distEdgt = nan(length(linpos),1);
        posint = nan(length(linpos),1);
        posoutt = nan(length(linpos),1);
        distEdgt(outbound) =  121-(linpos(nanconvn(diff(linpos),k)>0));
        distEdgt(inbound) =  linpos(nanconvn(diff(linpos),k)<0) ;
        posint(inbound) = linpos(inbound);
        posoutt(outbound) = linpos(outbound);
        
        xt = nanmean([sessiondata.behavior.position.left_ear_cor(:,1) sessiondata.behavior.position.right_ear_cor(:,1)],2);
        yt = nanmean([sessiondata.behavior.position.left_ear_cor(:,2) sessiondata.behavior.position.right_ear_cor(:,2)],2);
        posHt = [xt yt];
        posHt = posHt(1:nsamples,:); 
        
        
        xt = sessiondata.behavior.position.tail_base_cor(:,1);
        yt = sessiondata.behavior.position.tail_base_cor(:,2);
        posBt = [xt yt];
        posBt = posBt(1:nsamples,:); 
        
        
        % get time to reward
        time2rewt = nan(length(sessiondata.behavior.ts_video),1);
        rewts = sort([sessiondata.behavior.rewardTimeLeft;sessiondata.behavior.rewardTimeRight]);
        for j = 1:length(rewts)
            if j==1
                kp = sessiondata.behavior.ts_video<rewts(j);
            else
                kp = sessiondata.behavior.ts_video<rewts(j) & sessiondata.behavior.ts_video>rewts(j-1);
            end
            time2rewt(kp) =   rewts(j) - sessiondata.behavior.ts_video(kp);
            
        end
        
          % get time to reward
        timefromrewt = nan(length(sessiondata.behavior.ts_video),1);
       
        for j = 1:length(rewts)
            if j==length(rewts)
                kp = sessiondata.behavior.ts_video>rewts(j);
            else
                kp = sessiondata.behavior.ts_video>rewts(j) & sessiondata.behavior.ts_video<rewts(j+1);
            end
            timefromrewt(kp) =   rewts(j) - sessiondata.behavior.ts_video(kp);
            
        end
        
        
        
        
        %             totalDistanceFrmEdg = struct2cell(totalDistanceFrmEdg);
        %             totalDistanceFrmEdg = cell2mat(totalDistanceFrmEdg);
        
        if isfield(sessiondata,'contextEntry') && (size(sessiondata.contextEntry,1)==3)
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
            
            %get time from context transition
            
            
            
            kp_homet = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(home,:),2),'uni',0);
            kp_homet = cell2mat(kp_homet');
            kp_homet = any(kp_homet,2);
            
            kp_Novelt = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(inNovel,:),2),'uni',0);
            kp_Novelt = cell2mat(kp_Novelt');
            kp_Novelt = any(kp_Novelt,2);
            
            %get velocity
            velt = sessiondata.behavior.vel_cor;
            acct = sessiondata.behavior.acc_cor;
          
            %                dayNumt = sessiondata.behavior.dayNum;
            timeFrmEntryt = sessiondata.behavior.timeFrmEntry;
            
            %getRear
            
            %downsample neural data to video clock
            neuralt = sessiondata.neural.signal_DFoF;
            fs = sessiondata.neural.fs_neural;
            
            
            k = gaussian2Dfilter([1 10*fs],fs);
            neuralt = nanconvn(neuralt,k);
            ts_neural = (1:length(neuralt))/fs;
            
            neuralt = double(interp1(ts_neural,neuralt,ts_video));
            
            nSamples = length(neuralt);
            % concatenate each session
            neural = [neural;neuralt];
            acc = [acc;acct];
            vel = [vel;velt];
            distEdg = [distEdg; distEdgt];
            posin = [posin;posint];
            posout = [posout;posoutt];
            posH = [posH;posHt];
            posB = [posB;posBt];
            timeFrmEntry = [timeFrmEntry; timeFrmEntryt];
            timefromrew = [timefromrew;timefromrewt];
            
            %   dayNum = [dayNum; dayNumt];
        
            kp_Novel = [kp_Novel;kp_Novelt];
            time2rew = [time2rew;time2rewt];
            ts_tot = [ts_tot;ts_video(:)];
            sesID = [sesID; i *ones(nSamples,1)];
            subjID = [subjID;ID*ones(nSamples,1)];
            
            i

            pwd
        end
        
        end
    end
end

%%
%load('R:\DANEHippocampalResponse\Workspaces\LinearTrack_doubleExp.mat')
%%
%%
[occ,~,~,bb] = histcn(posH,-20:.2:140,-20:.2:20);
figure
imagesc(-20:2:140,-20:2:20,accumarray(bb+1,abs(posH(:,2)-posB(:,2)),[801 201],@nanmean,nan)')

%%
faceEdge = abs(posH(:,2)-posB(:,2));
acc1 = acc;
acc1(acc>30) = nan;
vel1 = vel;
vel1(vel>40) = nan;
ix = 1;
clear XX YY
distEdg1 = distEdg;
distEdg1(distEdg1<0) = -1;

k  = gaussian2Dfilter([10000 1],10);
idx = find( diff(abs(posH(:,2)-5)>5 & nanconvn(faceEdge,k)>3 & nanconvn(faceEdge,k)<10)>0);
idx = repmat(idx,1,151)+repmat(50:200,length(idx),1);
inRear = false(size(acc1));
inRear(idx(:)) = true;
intercept = true;
%%
idx = find(diff(nanconvn(faceEdge,k)>3)>0);

for i = 1:length(idx)
    figure
plot(pos(idx(i)-1000:idx(i)+1000,1),pos(idx(i)-1000:idx(i)+1000,2),'.')
hold on
plot(pos(idx(i)-30:idx(i)+30,1),pos(idx(i)-30:idx(i)+30,2))

waitforbuttonpress
close all
end

%%
load('R:\DANEHippocampalResponse\Workspaces\LinearTrack_doubleExp.mat')
%%
MSE_full =[];MSE_red=[];b_all_LT=[];pred_tim=[];

    for j = 1:6
        
        kp = kp_Novel==1 & subjID==j ;%& ts_tot>1000;
        if any(kp) 
        %cross val
       
        
        kp_train = mod(sesID(kp),2)==0;
       % kp_train = true(sum(kp),1);
        inter = ones(sum(kp),1);
        x1 = timeFrmEntry(kp);
        x2 =  [inter vel1(kp) acc1(kp) distEdg(kp) ];
        x3 = time2rew(kp);
        y = neural(kp);
        IV = [x2 x1 x3];
      %  IV =[inter x1];
        DV = neural(kp);
        isTime = false(6,1);
        isTime([5 ,6]) = true;
        
        %isTime = [false true];
        [MSE_full(j), MSE_red(j,:),b_all_LT(j,:),pred] = sm_model_comp_MSE(IV,DV,isTime,kp_train,true);
        
      %  pred_tim(i,:,j) = avghist(x1,pred,0:600);
        else
             MSE_full(j) = nan;
             MSE_red(j,:) = nan;
             b_all_NE_day(j,:) = nan;
            
        end
    
           %resave
%            usesID = unique(sesID(kp))';
%          for k  = usesID
%              saveID = dirs{k};
%              cd(saveID)
%              kp_sav = sesID(kp);
%              kp_sav = kp_sav(~kp_train);
%              kp_sav = kp_sav==k;
%              
%              if any(kp_sav)
%                  ts_sav = ts_tot(kp);
%                  ts_sav = ts_sav(~kp_train);
%                  ts_sav = ts_sav(kp_sav);
%                  
%                  pred_sav = pred(kp_sav);
%                  load('sessiondata.mat')
%                  
%                  pred_ses = nan(length(sessiondata.behavior.ts_video),1);
%                  pred_ses(ismember(sessiondata.behavior.ts_video,ts_sav)) = pred_sav;
%                  sessiondata.neural.pred.yhat = pred_ses;
%                  save('sessiondata.mat','sessiondata','-v7.3')
%              end
%          end
           
    end


%%

nsubj = size(MSE_red,1);
nvar = size(MSE_red,2);
ok=((MSE_red - MSE_full'))./MSE_full';
close all
figure
bar(100*nanmean(ok),'facecolor','w')
hold on

plot(repmat(1:nvar,nsubj,1)'+(rand(nvar,nsubj)-.5)/5,100*ok','o')
    %%
    
    b_all_HC =[];MSE_full =[];MSE_red=[];pred_tim_HC=[];real_tim_HC =[];
    
    for j = 1:6
        
        kp = kp_Novel==0 & subjID==j & ts_tot>1000;
        if any(kp)
            %cross val
            
            
            kp_train = mod(sesID(kp),2)==0;
            %  kp_train = true(sum(kp),1);
            inter = ones(sum(kp),1);
            x1 = timeFrmEntry(kp);
            x2 =  [inter vel1(kp) acc1(kp)  ];
            x3 = time2rew(kp);
            y = neural(kp);
            IV = [x2 x1 ];
            %  IV =[inter x1];
            DV = neural(kp);
            isTime = false(4,1);
            isTime([4]) = true;
            
            %isTime = [false true];
            [MSE_full(j), MSE_red(j,:),b_all_HC(j,:),pred] = sm_model_comp_MSE(IV,DV,isTime,kp_train,true);
            
          %  pred_tim_HC(j,:) = avghist(x1,pred,0:600);
          %  real_tim_HC(j,:) = avghist(x1,y,0:600);
            %
            %resave
%             usesID = unique(sesID(kp))';
%             for k  = usesID
%                 saveID = dirs{k};
%                 cd(saveID)
%                 kp_sav = sesID(kp);
%                 kp_sav = kp_sav(~kp_train);
%                 kp_sav = kp_sav==k;
%                 
%                 if any(kp_sav)
%                     ts_sav = ts_tot(kp);
%                     ts_sav = ts_sav(~kp_train);
%                     ts_sav = ts_sav(kp_sav);
%                     
%                     pred_sav = pred(kp_sav);
%                     load('sessiondata.mat')
%                     
%                     pred_ses = sessiondata.neural.pred.yhat;
%                     pred_ses(ismember(sessiondata.behavior.ts_video,ts_sav)) = pred_sav;
%                     sessiondata.neural.pred.yhat = pred_ses;
%                     save('sessiondata.mat','sessiondata','-v7.3')
%                 end
%             end
            
        else
            MSE_full(j) = nan;
            MSE_red(j,:) = nan;
            b_all_NE_day(j,:) = nan;
            
        end
    end
    

    %%
    

b_all_LT =[];MSE_full =[];MSE_red=[];pred_tim_LT=[];

    for j = 1:6
        
        kp = kp_Novel==1 & subjID==j ;%& ts_tot>1000;
        if any(kp) 
        %cross val
       
        
       % kp_train = mod(sesID(kp),2)==0;
        kp_train = true(sum(kp),1);
        inter = ones(sum(kp),1);
        x1 = timeFrmEntry(kp);
        x2 =  [inter vel1(kp) acc1(kp) distEdg(kp) ];
        x3 = time2rew(kp);
        y = neural(kp);
        IV = [x2 x1 x3];
      %  IV =[inter x1];
        DV = neural(kp);
        isTime = false(6,1);
        isTime([5 ,6]) = true;
        
        %isTime = [false true];
        [~, ~,b_all_LT(j,:),pred] = sm_model_comp_MSE(IV,DV,isTime,kp_train,false);
        
        pred_tim_LT(j,:) = avghist(x1,pred,0:600);
         real_tim_LT(j,:) = avghist(x1,y,0:600);
        else
             MSE_full(j) = nan;
             MSE_red(j,:) = nan;
             b_all_NE_day(j,:) = nan;
            
        end
    
     
        
    end

    
    
    
%%
figure
plotMeanSEM(0:600,pred_tim_LT,'k')
plotMeanSEM(0:600,real_tim_LT,'r')
%%
figure
plotMeanSEM(0:600,pred_tim_HC,'k')
plotMeanSEM(0:600,real_tim_HC,'r')


%%
close all
figure
uiopen('E:\Dropbox\UNM\Papers\NE_timeConstants\Figures\beta_tau_NE_days.fig',1)


    
hold on

    u1 = nanmean(b_all_LT(:,5,:));
       u2 = nanmean(b_all_LT(:,6,:));
        s1 = SEM(b_all_LT(:,5,:));
       s2 = SEM(b_all_LT(:,6,:));
       
       
        errorbar(u1,u2,s2,s2,s1,s1,'color',[.7 .7 .7])
        
          u1 = nanmean(b_all_HC(:,4,:));
       u2 = nanmean(b_all_HC(:,5,:));
        s1 = SEM(b_all_HC(:,4,:));
       s2 = SEM(b_all_HC(:,5,:));
       
       
        errorbar(u1,u2,s2,s2,s1,s1,'color','k')
        ylim([-.1 0])
        

%%
nsubj = size(MSE_red,1);
nvar = size(MSE_red,2);
ok=((MSE_red - MSE_full'))./MSE_full';
close all
figure
bar(100*nanmean(ok),'facecolor','w')
hold on

plot(repmat(1:nvar,nsubj,1)'+(rand(nvar,nsubj)-.5)/5,100*ok','o')


%%


b_all_LT =[];
for i = 1:10

kp = kp_Novel==1 & subjID==i;
if ~all(isnan(inRear(kp)))
inter = ones(sum(kp),1);
x1 = timeFrmEntry(kp);
x2 =  [inter vel1(kp) acc1(kp) distEdg(kp) time2rew(kp) inRear(kp)];
y = neural(kp);
fitfun = @(b) b(1).*exp(b(2).*x1)+  b(3).*exp(b(4).*x1) + x2*b(5:10)';
f = @(b) nanmean((y - fitfun(b)).^2 ) + 0.001*sqrt(b(1)^2 + b(3)^2 + sum(b(5:end).^2)); 

x0 = [0 -.05 0 -.05 0 0 0 0 0 0];
Aeq =[];
beq =[];
A =[];
b = [];
lb = [0 -.1  -20 -.1 -2*ones(1,6)];
ub = [20 -.001 0 -.001 2*ones(1,6)];
b = fmincon(f,x0,A,b,Aeq,beq,lb,ub)% Objective Function
b_all_LT(i,:) = b;
end
end

%%
b_all_hc=[];
for i = 1:10

kp = kp_Novel==0 & subjID==i;
if ~all(isnan(inRear(kp)))
inter = ones(sum(kp),1);
x1 = timeFrmEntry(kp);
x2 =  [inter vel1(kp) acc1(kp)];
y = neural(kp);
fitfun = @(b) b(1).*exp(b(2).*x1)+  b(3).*exp(b(4).*x1) + x2*b(5:7)';
f = @(b) nanmean((y - fitfun(b)).^2 ) + 0.001*sqrt(b(1)^2 + b(3)^2 + sum(b(5:end).^2)); 

x0 = [0 -.05 0 -.05 0 0 0];
Aeq =[];
beq =[];
A =[];
b = [];
lb = [0 -.1  -20 -.1 -2*ones(1,3)];
ub = [20 -.001 0 -.001 2*ones(1,3)];
b = fmincon(f,x0,A,b,Aeq,beq,lb,ub)% Objective Function
b_all_hc(i,:) = b;
end
end

%%
neural(abs(neural)>10) = nan;
idx = 1;
clear beta_LT mse_LT
for t = fliplr(-logspace(log10(.0001),log10(.1),30))
    ix = 1;
    clear YY XX
for i = 1:10
    kp = subjID==i & kp_Novel==1; 
   if any(kp) 
    YY{ix} = neural(kp);
   YY{ix}( YY{ix}<-5) = nan;
    XX{ix} = [ exp(t.*timeFrmEntry(kp))  vel1(kp) acc1(kp) distEdg(kp) time2rew(kp) inRear(kp)];
 %  D = dummyvar(kp_Novel(kp)+1);
  %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
    ix  = ix+1;
   
   end
end

if intercept
    
      
    stats_NE = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc','dist','time2rew','inRear'});
    beta_LT(:,idx) = stats_NE.first_level.beta(2,:);
    for i =1:length(XX)
        mse_LT(i,idx) = nanmean(sqrt((YY{i} - [ones(size(XX{i},1),1) XX{i}]* stats_NE.first_level.beta(:,i)).^2));
    end
    
else
    
    stats_NE = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc','dist','time2rew','inRear'},'noint');
    beta_LT(:,idx) = stats_NE.first_level.beta(1,:);
    for i =1:length(XX)
        mse_LT(i,idx) = nanmean(sqrt((YY{i} - [ XX{i}]* stats_NE.first_level.beta(:,i)).^2));
    end
    
end


idx=  idx+1
end


%%
ok =  fliplr(-logspace(log10(.0001),log10(.1),30));
[a,b] = min(sum(mse_LT));

[~,ix] = min(mse_LT,[],2);
tau_hat_LT = ok(ix);
beta_LT = beta_LT(sub2ind(size(beta_LT),1:size(beta_LT,1),ix'));
tau_hat = ok(b);

%%

neural(abs(neural)>10) = nan;
idx = 1;
clear mse_HC beta_HC
for t = fliplr(-logspace(log10(.0001),log10(.1),30))
    ix = 1;
    clear YY XX
for i = 1:10
    kp = subjID==i & kp_Novel==0; 
   if any(kp) 
    YY{ix} = neural(kp);
   YY{ix}( YY{ix}<-5) = nan;
    XX{ix} = [ exp(t.*timeFrmEntry(kp))  vel1(kp) acc1(kp) ];
 %  D = dummyvar(kp_Novel(kp)+1);
  %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
    ix  = ix+1;
   
   end
end

if intercept
      stats_HC = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc'});
    beta_HC(:,idx) = stats_HC.first_level.beta(2,:);
    for i =1:length(XX)
        mse_HC(i,idx) = nanmean(sqrt((YY{i} - [ones(size(XX{i},1),1) XX{i}]* stats_HC.first_level.beta(:,i)).^2));
    end
else
    
    
    stats_HC = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time','vel','acc'},'noint');
    beta_HC(:,idx) = stats_HC.first_level.beta(1,:);
    for i =1:length(XX)
        mse_HC(i,idx) = nanmean(sqrt((YY{i} - [ XX{i}]* stats_HC.first_level.beta(:,i)).^2));
    end
end
idx=  idx+1
end


%%
ok =  fliplr(-logspace(log10(.0001),log10(.1),30));
[a,b] = min(sum(mse_HC));

[~,ix] = min(mse_HC,[],2);
tau_hat_HC = ok(ix);
beta_HC = beta_HC(sub2ind(size(beta_HC),1:size(beta_HC,1),ix'));


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
neural1(neural1<-5 | neural1>5) = nan;

clear XX YY
ix=1;
for i =  1:10
    kp = subjID==i & kp_Novel==1 & time2rew>5; 
   if any(kp) 
    YY{ix} = neural1(kp);
     
    time = exp(tau_hat.*timeFrmEntry);
   
   % XX{ix} = [ vel(kp) acc(kp) distEdg1(kp) exp(tau_hat.*timeFrmEntry(kp)) inRear(kp) ];
        XX{ix} = [ vel1(kp) acc1(kp) distEdg(kp) time2rew(kp) inRear(kp) timefromrew(kp)];
  
   D = dummyvar(kp_Novel(kp)+1);
  %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
    ix  = ix+1;
   end
end
stats = glmfit_multilevel(YY, XX, [], 'verbose', 'weighted','names',{'vel','Acc','dist','time2rew','inRear','timefromrew'});
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
