%% change for how the R drive is mounted on your computer
masterDir = 'R:\';
masterDir =  'R:\McKenzieLab\';

%% DA2h3 (KEEP)
dirs = [...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210817-161523']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210818-135519']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210819-115556']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210820-145938']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210821-084618']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210822-091253']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210823-161754']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210824-091830']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210825-094340']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210826-091103']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210827-112030']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210830-134951']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210831-142633']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210901-111232']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h3-210902-122150']}; ...
    ];


%% NE2m3 (KEEP)
dirs = [...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210819-123801']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210820-140801']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210821-081259']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210822-094611']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210823-165050']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210824-124626']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210825-101604']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210826-105416']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210827-104741']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210830-130702']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210831-135215']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210901-103622']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210902-114718']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210903-092753']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2m3-210907-135954']} ; ...
    ];

%% NO RESPONSE
% 
% dirs = [...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210817-130415']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210819-140301']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210823-150453']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210824-132044']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210825-110018']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210826-145040']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210827-135245']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210830-142332']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210831-145937']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210901-114549']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210902-125554']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210903-104015']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h1\DA2h1-210907-132533']} ; ...
%     ];
%
%% DA2h4 (KEEP)
dirs = [...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-210927-152919']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-210928-135333']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-210929-141310']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-210930-152412']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211001-122043']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211002-133705']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211003-121913']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211004-103754']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211005-115604']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211006-125255']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211007-112940']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211008-092022']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211009-125249']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211010-093531']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\DA2h4-211011-135337']} ; ...
    ];

%% NE2h4 (KEEP)
dirs = [...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220119-104713']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220131-131223']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220201-135617']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220202-091338']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220203-142033']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220204-135854']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220207-104742']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220208-134340']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220209-114835']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220210-135040']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h4\NE2h4-220211-122642']} ; ...
    ];

%% NE2h5 (??)
dirs = [...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h5-220119-120430']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h5-220120-120308']} ; ...
];





%% NE2h6 (unstable, signal unknown)

dirs = [...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h6-220207-124942']} ; ...

];



%% NE2h7 (KEEP) first 3 exclude had to start sequence over
dirs = [...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220119-135318']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220120-123703']} ; ...
%     {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220121-105533']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220131-120800']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220201-163937']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220202-084149']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220203-145222']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220204-132714']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220207-111933']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220208-130905']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220209-122221']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220210-131745']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h7\NE2h7-220211-125841']} ; ...
    ];



%% NE2h8 (unstable but signal)
dirs = [...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h8-220207-115255']} ; ...
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h8-220207-121417']} ; ...
    
    ];

%% NE2h9 (stable and signal)
dirs = [...
    
    {[masterDir 'DonaldsonT\DA_animals_NEE-210817-124922\NE2h9-220207-141915']} ; ...
    
    ];



%%

newCon = [];homeCage =[];

checkStability = false;
k  = gaussian2Dfilter([1000 1],[100 1]);

saveAllDays = true;


for i = 1:length(dirs)
    
    
    cd(dirs{i})
    
    if checkStability
        
        [signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(pwd,'photoBleachCorrection','exp2');
        plot(ts_data,nanconvn(signal_DFoF,k'))
        ylim([-5 5])
        waitforbuttonpress
    end
    
    
    if exist('newContext.mat')
        [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(pwd,'newContext.mat',{'newContext','homeCage'},'photoBleachCorrection','exp2','plotIntervals',[300 300],'returnedDataType','corrected');
        
        
        % average all new context and home cage
        newCon(i,:) = nanmean(signal_DFoF(ix{1}),1); %first element of ix is the first event type which is newContext
        homeCage(i,:) = nanmean(signal_DFoF(ix{2}),1);
        
    else
        newCon(i,:) = nan(1,size(ts_PETH,2));
        homeCage(i,:) =  nan(1,size(ts_PETH,2));
        
    end
    close all
    
end





k  = gaussian2Dfilter([1000 1],[200 1]);
newCon1 =[];
homeCage1= [];
for i = 1:size(newCon,1)
    
    
    newCon1(i,:) = nanconvn(newCon(i,:),k');
    homeCage1(i,:) = nanconvn(homeCage(i,:),k');
    
end

if saveAllDays
    
    outDir = 'R:\DonaldsonT\DA_animals_NEE-210817-124922\DA2h1';
    outfil = [outDir filesep 'DA2h1_allNovel.mat'];
    
    save(outfil,'-v7.3')
end

%%

close all
figure
ax = tight_subplot(1,2);

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

xlabel('Time from context entry')%%


%%
%close all
figure
plot(ts_PETH,nanmean(newCon1,1),'b')
hold on
plot(ts_PETH,nanmean(homeCage1,1), 'k')
xlim([-300 300])
ylim([-0.5 2.5])

xlabel ('Time (s)')
ylabel ('Z-Score')
legend ('Novel','Home')
set(gca,'box','off')
%%

% get mean response for the day (AUC from 0-60s )

%get index at time 0
[a,ix_0] = bestmatch(0,ts_PETH);

%get index at time 100
[a,ix_100] = bestmatch(60,ts_PETH);


u_response_novel = mean(newCon1(:,ix_0:ix_100),2);
u_response_HC = mean(homeCage1(:,ix_0:ix_100),2);

m_response_novel = max(newCon1(:,ix_0:ix_100),[],2);
m_response_HC = max(homeCage1(:,ix_0:ix_100),[],2);


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