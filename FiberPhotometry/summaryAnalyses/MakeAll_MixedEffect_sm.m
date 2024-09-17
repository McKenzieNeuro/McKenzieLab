%GLME Plots
%These are the GLME 1st, 2nd, and 3rd best fit models
%Without HomeCage (HC) data, and without dayNum data


%NE 1st Best: glmebfPacd = fitglme(dat,'neural ~ 1   + acceleration*timeFrmEntry  + velocity + kp_novel  + distEdg  +  (1|sesID:subjID) + (1|subjID)');
%NE 2nd Best: glmebcPaf

%GLME Plots

fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);

kp = cellfun(@any,regexp(fils,'Novel Env')) & cellfun(@any,regexp(fils,'NE2'));

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
inRear =[];
faceEdge =[];

for i = 1:length(dirs)
        
    cd(dirs{i})
 

    if exist('sessiondata.mat')
    
            [~,ID] = ismember(subj{i},usubjs);
            
            %load struct with neural/behavioral data
            load('sessiondata.mat')
            if exist('rear.mat')
            vv = load('rear.mat');
             rears = [cell2mat(vv.data(contains(vv.data(:,1),'start','IgnoreCase',true),2)) cell2mat(vv.data(contains(vv.data(:,1),'end','IgnoreCase',true),2))];
             rears = MergeEpochs2(rears);   
             inReart = double(InIntervals(sessiondata.behavior.ts_video,[rears(:,1) rears(:,2)]));
            else
                
                inReart = nan(length(sessiondata.behavior.ts_video),1);
            end
           
            totalEdgDist = sessiondata.behavior.totalEdgDist;
            
            
%             totalDistanceFrmEdg = struct2cell(totalDistanceFrmEdg);
%             totalDistanceFrmEdg = cell2mat(totalDistanceFrmEdg);
           
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
                
                %get time from context transition
                
                ts_tran = sessiondata.behavior.timeFrmEntry;
                
                kp_homet = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(home,:),2),'uni',0);
                kp_homet = cell2mat(kp_homet');
                kp_homet = any(kp_homet,2);
                
                kp_Novelt = cellfun(@(a)  ts_video>a(1) & ts_video<a(2),num2cell(epochs(inNovel,:),2),'uni',0);
                kp_Novelt = cell2mat(kp_Novelt');
                kp_Novelt = any(kp_Novelt,2);
                
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
                
            
                timeFrmEntry = [timeFrmEntry; timeFrmEntryt];
             %   dayNum = [dayNum; dayNumt];
                inRear = [inRear;inReart]; 
                kp_Novel = [kp_Novel;kp_Novelt];
                ts_trans = [ts_trans;ts_tran];
                faceEdge = [faceEdge;faceEdget];
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
%%
idx = 1;

for t = fliplr(-logspace(log10(.0001),log10(.1),30))
    ix = 1;
    clear YY XX
for i = 1:10
    kp = subjID==i & kp_Novel==1; 
   if any(kp) && ~all(isnan(inRear(kp)))
    YY{ix} = neural(kp);
   YY{ix}( YY{ix}<-5) = nan;
    XX{ix} = [  exp(t.*timeFrmEntry(kp))];
 %  D = dummyvar(kp_Novel(kp)+1);
  %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
    ix  = ix+1;
   
   end
end
stats = glmfit_multilevel(YY, XX, [], 'verbose','names',{'time'});

for i =1:length(XX)
mse(i,idx) = nanmean(sqrt((YY{i} - [ones(size(XX{i},1),1) XX{i}]* stats.first_level.beta(:,i)).^2));
end
idx=  idx+1
end


%%
ok =  fliplr(-logspace(log10(.0001),log10(.1),30));
[a,b] = min(sum(mse));
[~,ix] = min(mse,[],2);
tau_hat_all = ok(ix);
tau_hat = ok(b);

all_subj = usubjs(unique(subjID(~(isnan(inRear)))));
ix=1;
clear YY XX


for i = 1:10
    kp = subjID==i & kp_Novel==1; 
   if any(kp) && ~all(isnan(inRear(kp)))
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
   if any(kp) && ~all(isnan(inRear(kp)))
       
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
   if any(kp) && ~all(isnan(inRear(kp))) 
    YY{ix} = neural1(kp);
    time = exp(tau_hat.*timeFrmEntry);
   
   % XX{ix} = [ vel(kp) acc(kp) distEdg1(kp) exp(tau_hat.*timeFrmEntry(kp)) inRear(kp) ];
        XX{ix} = [ vel(kp) acc(kp) distEdg(kp)];
  
   D = dummyvar(kp_Novel(kp)+1);
  %  XX{ix} = [D.*vel(kp)  D.*acc1(kp)  D.*distEdg(kp) D.*timeFrmEntry(kp) D.*inRear(kp)  ];
    ix  = ix+1;
   end
end
stats = glmfit_multilevel(YY, XX, [], 'verbose', 'weighted','names',{'vel','Acc','dist'});
%stats = glmfit_multilevel(YY, XX, [], 'verbose', 'weighted','names',{'intAcc','intRear','intDist','intTime'});

for i =1:length(XX)
mse(i) = nanmean(sqrt((YY{i} - [ones(size(XX{i},1),1) XX{i}]* stats.first_level.beta(:,i)).^2));
end
Est = neural;
ix=1;
for i = 1:10
    kp = subjID==i & kp_Novel==1 & timeFrmEntry>10; 
   if any(kp) && ~all(isnan(inRear(kp)))
       
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
idx=  find(diff(inRear(kp))>0);
ok = neural(kp);
idx = repmat(idx,1,2001)+repmat(-1000:1000,length(idx),1);
kp = all(idx>0 &idx<=numel(ok),2);
idx = idx(kp,:);

ru(i,:) = nanmean(ok(idx));
 kp = kp_Novel==1 & subjID==i;
ok = Est(kp);
con(i,:) = nanmean(ok(idx));
end
%%
save('R:\DANEHippocampalResponse\Workspaces\novelContext.mat')