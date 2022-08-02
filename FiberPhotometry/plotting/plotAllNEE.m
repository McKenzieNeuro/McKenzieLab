
masterDir = 'R:\DonaldsonT\DA_animals_NEE-210817-124922\';
dirs = [...
    {[masterDir 'NE2m3']} ; ...
    {[masterDir 'NE2h4']} ; ...
    {[masterDir 'NE2h6']} ; ...
    {[masterDir 'NE2h7']} ; ...
    {[masterDir 'NE2h8']} ; ...
    {[masterDir 'NE2h9']} ; ...
  ];
    

%%


%clear newContext homeCage
for i = 5:length(dirs)
    
    
    cd(dirs{i})
    subDir = dir(dirs{i});
    subDir = {subDir(cell2mat({subDir.isdir})).name}';
    subDir = subDir(3:end);
    for j = 1:length(subDir)
        cd(subDir{j})
        
        if exist('newContext.mat')
            [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(pwd,'newContext.mat',{'newContext','homeCage'},'photoBleachCorrection','exp2','plotIntervals',[300 300],'returnedDataType','corrected');
            newContext{i}(j,:) = signal_DFoF(ix{1}(1,:));
            if ~isempty(ix{2})
            homeCage{i}(j,:) = signal_DFoF(ix{2}(1,:));
            else
                 homeCage{i} = nan(1,length(ix{1}(1,:)));
            end
           [i j] 
        end
           cd(dirs{i})
    end
    close all
    
end

%%

newContext = cellArrayTo3D(cellfun(@(a) doubble(a(1:11,:)),newContext,'uni',0));
homeCage = cellArrayTo3D(cellfun(@(a) double(a(1:11,:)),homeCage,'uni',0));


%%
k  = gaussian2Dfilter([1000 1],[500 1]);
homeCage1 =homeCage;
newContext1 =newContext;
for i = 1:size(homeCage,1)
    
    for j = 1:size(homeCage,3)
        homeCage1(i,:,j) = nanconvn(homeCage(i,:,j),k');
        newContext1(i,:,j) = nanconvn(newContext(i,:,j),k');
    end
end
%%
close all

figure
imagesc(ts_PETH,[],nanmean(newContext1,3),[-1 1])
ylabel('Days of exposure')
xlabel('Time to context entrance (s)')

figure
imagesc(ts_PETH,[],nanmean(homeCage1,3),[-1 1])
ylabel('Days of exposure')
xlabel('Time to context entrance (s)')
%%


figure


plotMeanSEM(ts_PETH,squeeze(homeCage1(1,:,:))','r')
shg
hold on
plotMeanSEM(ts_PETH,squeeze(newContext1(1,:,:))','k')
plot([-300 -200],[1 1],'k')
plot([-300 300],[0 0],'k')
plot([0 0],[-1 2.5],'k')
plot([-300 -300],[1 1.5],'k')
axis off

ylim([-1 3])
%%
figure


plotMeanSEM(ts_PETH,squeeze(homeCage1(end,:,:))','r')
shg
hold on
plotMeanSEM(ts_PETH,squeeze(newContext1(end,:,:))','k')
plot([-300 -200],[1 1],'k')
plot([-300 300],[0 0],'k')
plot([0 0],[-1 2.5],'k')
plot([-300 -300],[1 1.5],'k')
axis off
ylim([-1 3])
%%

%plot phase over days

close all
% get mean response for the day (AUC from 0-60s )

%get index at time 0
[a,ix_0] = bestmatch(0,ts_PETH);

%get index at time 100
[a,ix_10] = bestmatch(10,ts_PETH);
[a,ix_100] = bestmatch(30,ts_PETH);

u_response_HC_phasic = squeeze(nanmean(homeCage1(:,ix_0:ix_10,:),2));
u_response_HC_tonic = squeeze(nanmean(homeCage1(:,ix_10:ix_100,:),2));

u_response_novel_phasic = squeeze(nanmean(newContext1(:,ix_0:ix_10,:),2));
u_response_novel_tonic = squeeze(nanmean(newContext1(:,ix_10:ix_100,:),2));


figure
 plotMeanSEM(1:11,u_response_HC_phasic','k')
hold on
plotMeanSEM(1:11,u_response_HC_tonic','r')
ylim([-.5 2])


figure
 plotMeanSEM(1:11,u_response_novel_phasic','k')
hold on
plotMeanSEM(1:11,u_response_novel_tonic','r')
ylim([-.5 2])

%%

u_response_HC = mean(homeCage1(:,ix_0:ix_10),2);

m_response_novel = max(newCon1(:,ix_0:ix_10),[],2);
m_response_HC = max(homeCage1(:,ix_0:ix_10),[],2);


%%


newCon1 =[];
homeCage1= [];
newConDrug1 = [];
for i = 1:size(newCon,1)
    
    
    newCon1(i,:) = nanconvn(newCon(i,:),k');
    homeCage1(i,:) = nanconvn(homeCage(i,:),k');
    newConDrug1(i,:) = nanconvn(newConDrug(i,:),k');
end

if saveAllDays
    
    outDir = 'R:\DonaldsonT\DA_animals_NEE-210817-124922\DA2h1';
    outfil = [outDir filesep 'DA2h1_allNovel.mat'];
    
    save(outfil,'-v7.3')
end

%%

close all
figure
ax = tight_subplot(1,3);

axes(ax(1))
imagesc(ts_PETH,[],newCon1,[-1.5 2])
title('Novel Context')
ylabel('Day')
xlabel('Time From Context Entry (s)')
xlim([-300 300])
caxis([-1.5 2.5])
axes(ax(2))
imagesc(ts_PETH,[],homeCage1,[-1.5 2])
title('Home cage')
axes(ax(3))
imagesc(ts_PETH,[],newConDrug1,[-1.5 2])
title('Novel Context Drug')

xlabel('Time from context entry')%%


%%
%close all
figure
plot(ts_PETH,nanmean(newCon1,1),'b')
hold on
plot(ts_PETH,nanmean(homeCage1,1), 'k')
hold on 
plot(ts_PETH,nanmean(newConDrug1,1), 'r')
xlim([-300 300])
ylim([-0.5 2.5])

xlabel ('Time (s)')
ylabel ('Z-Score')
legend ('Novel','Home')
set(gca,'box','off')
%%




%%
figure
hold on

plot(m_response_novel,'b')
plot(m_response_HC,'k')

xlim([0 15])
ylim([0 4])

xlabel ('Day')
ylabel ('Z-Score')
legend ('Novel','Home')
set(gca,'box','off')
%%


% make scatter plot of homecage vs novel

figure
plot(mean(homeCage1(:,ix_0:ix_100),2),mean(newCon1(:,ix_0:ix_100),2),'.')
hold on
plot(-.3:.1:1,-.3:.1:1)
xlim([-.3 1])
ylim([-.3 1])




xlabel ('Familiar')
ylabel ('Novel')
set(gca,'box','off')



%%

clear all
%load all sessions
ok = getAllExtFiles('R:\DonaldsonT\DA_animals_NEE-210817-124922','mat',1);
fils = ok(cellfun(@any,regexp(ok,'allNovel')));
clear u_response_novel u_response_HC


u_response_novel = nan(length(fils),15);

u_response_HC = nan(length(fils),15);

for  i = 1:length(fils)
    
    v= load(fils{i});
    
    nDays = size(v.newCon,1);
    %get index at time 0
    [a,ix_0] = bestmatch(0,v.ts_PETH);
    
    %get index at time 100
    [a,ix_100] = bestmatch(100,v.ts_PETH);
    
    
    u_response_novel(i,1:nDays) = max(v.newCon1(:,ix_0:ix_100),[],2);
    u_response_HC(i,1:nDays) = max(v.homeCage1(:,ix_0:ix_100),[],2);
    i
  
end
%%


nTr = 2;
EL_novel = [nanmean(u_response_novel(:,1:nTr),2) nanmean(u_response_novel(:,end-nTr+1:end),2)];
EL_HC = [nanmean(u_response_HC(:,1:nTr),2) nanmean(u_response_HC(:,end-nTr+1:end),2)];


close all
figure
plot([0,1], EL_novel')
ylim([-1 2.5])
set(gca,'xtick',[0 1],'xticklabel',{'Early','Late'})
legend({'DA2h3','DA2h4','NE2m3','NE2h4','NE2h7'})
title('Novel environment')
figure
plot([0,1],EL_HC')
ylim([-1 2.5])
set(gca,'xtick',[0 1],'xticklabel',{'Early','Late'})
legend({'DA2h3','DA2h4','NE2m3','NE2h4','NE2h7'})
title('Home Cage')