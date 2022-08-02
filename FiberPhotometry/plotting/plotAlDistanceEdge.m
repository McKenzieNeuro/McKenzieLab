

fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'contextTransition')) & cellfun(@any,regexp(fils,'Novel Environment'));
fils = fils(kp);
[dirs,bi] = fileparts(fils);
% sl  =regexp(fils,'\');
 da  =regexp(fils,'-');
% subjs = cellfun(@(a,b) a(b(3)+1:b(4)-1),fils,sl,'uni',0);
dates = cellfun(@(a,b) a(b(1)+1:b(2)-1),fils,da,'uni',0);
%%
tot = 0;
for i = 50:length(dirs)
    
    
    cd(dirs{i})
    
    if ~exist('arena_edges.mat')
        tot = tot+1;
        %cd(dirs{i})
        %error('work here')
    end
    
end

%%








for i = 1:length(dirs)
    
    
    cd(dirs{i})
    
    if exist('sessiondata.mat')
        load('sessiondata.mat')
    end
    
    
    v = load(fils{i});
   
    
    sessiondata.contextEntry = v.data;
    save('sessiondata.mat','sessiondata','-v7.3')
    i
end

%%
k  = gaussian2Dfilter([10000 1],[ 1017.3 1]);
clear context  sessionDate subj PETH
for i = 1:length(dirs)
    
    
    cd(dirs{i})
    
    if exist('sessiondata.mat')
        load('sessiondata.mat')
        
        subj{i} = sessiondata.subject;
        sessionDate{i} = dates{i};
       
        ts = cell2mat(sessiondata.contextEntry(:,2));
        fs =sessiondata.neural.fs_neural;
        %get PETH for each conext
        
        
        ok = nanconvn(sessiondata.neural.signal_DFoF,k');
        [ix,early,late,ts_PETH] = sm_getIndicesAroundEvent(ts,300,300,fs,length(sessiondata.neural.signal_DFoF));
        
         kp = ~early & ~late;
         
         ix = ix(kp,:);
         contexts{i} = sessiondata.contextEntry(kp,1);
         PETH{i} = ok(ix);
         
    end

i
end
subj = cellfun(@(a) a(1:5),subj,'uni',0)';
%%

kp_DA = cellfun(@any,regexp(subj,'DA2'));

kp_NE = cellfun(@any,regexp(subj,'NE'));


homeCage = cellfun(@(a) find(cellfun(@any,regexp(a,'home'))),contexts,'uni',0);
other = cellfun(@(a) find(~cellfun(@any,regexp(a,'home'))),contexts,'uni',0);

DA_homeCage = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp_DA),homeCage(kp_DA),'uni',0)');
NE_homeCage = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp_NE),homeCage(kp_NE),'uni',0)');

DA_otherContext = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp_DA),other(kp_DA),'uni',0)');
NE_otherContext = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp_NE),other(kp_NE),'uni',0)');


%%
figure
plotMeanSEM(ts_PETH,NE_otherContext,'k')

hold on
plotMeanSEM(ts_PETH+700,NE_homeCage,'k')

%%

figure
plotMeanSEM(ts_PETH,DA_otherContext,'k')

hold on
plotMeanSEM(ts_PETH+700,DA_homeCage,'k')
title('DA')

%%
context1=  contexts;
% get context day # for each subject
all_subj = unique(subj)';
for i =1:length(all_subj)
    
    %find all sessions for a subject
   kp = find(ismember(subj,all_subj{i}));
   cnttmp = contexts(kp);
   
   %sort by date
   [~,b] = sort(dates(kp));
   [~,bb] = sort(b);
   cnttmp  = cnttmp(b);
   
   
   % get all unique labels
   c =[];
   for j = 1:length(kp)
      
        c = [c ;cnttmp{j}(~cellfun(@any,regexp(cnttmp{j},'home')))];
       
        
   end
   c = unique(c);
   
   %now loop over dates to count how many times previously that label
   %appeared
   
   cnt = zeros(length(kp),length(c));
   for j = 1:length(kp)
       expo = 100*ones(length(cnttmp{j}),1);
       c1 = cnttmp{j}(~cellfun(@any,regexp(cnttmp{j},'home')));
       kpo =  ~cellfun(@any,regexp(cnttmp{j},'home'));
       if j==1
           
       [cnt(j,:),~] = ismember(c,c1);
         [~,bx] = ismember(c1,c);
     
       else
          [tt,~]  = ismember(c,c1);
           [~,bx]  = ismember(c1,c);
          cnt(j,:) = tt'+cnt(j-1,:);
       
       end
        expo(kpo) = cnt(j,bx);
       cnttmp{j} = [cnttmp{j} num2cell(expo,2)];
   end
cnttmp = cnttmp(bb);
context1(kp) = cnttmp;
end
%%


kp_NE = cellfun(@any,regexp(subj,'NE'));


kp_DA = cellfun(@any,regexp(subj,'DA'));
kp  = kp_DA;
figure
ax  = tight_subplot(10,1);
for i = 1:10
   


axes(ax(i))
homeCage = cellfun(@(a) find(cellfun(@any,regexp(a(:,1),'home'))),context1,'uni',0);
other = cellfun(@(a) find(~cellfun(@any,regexp(a(:,1),'home')) & cell2mat(a(:,2))==i),context1,'uni',0);

homeCageResp = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp),homeCage(kp),'uni',0)');
otherContextResp = cell2mat(cellfun(@(a,b) a(b,:),PETH(kp),other(kp),'uni',0)');

plotMeanSEM(ts_PETH(1:100:end),otherContextResp(:,1:100:end),'k')
    
end
