%% change for how the R drive is mounted on your computer
masterDir = 'R:\';
masterDir =  'R:\McKenzieLab\';

%% DA2h1

% is there a positive control where we do see a DA response?

dirs = [...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211112-133652']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211115-090942']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211116-155846']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211117-111858']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211118-111628']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211119-114356']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211122-142745']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211123-112439']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211124-105724']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\DA2h1-211129-111811']} ; ...
];

%% NE2h2
dirs = [...
 {[masterDir 'DonaldsonT\FR-211025-121822\NE2h2-211110-114609']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\NE2h2-211112-124932']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\NE2h2-211115-095513']} ; ...
 {[masterDir 'DonaldsonT\FR-211025-121822\NE2h2-211116-150808']} ; ...
];
%%
loopOn = [];loopOff =[];
for i = 1:length(dirs)
    
    
    cd(dirs{i})
 [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(pwd,'foodReward.mat',{'loopOn','loopOff'},'photoBleachCorrection','quadratic','plotIntervals',[300 300],'returnedDataType','corrected');
   
 for j = 1:2
        newObj{i,j} = signal_DFoF(ix{j});
    end
    close all 

end
%%
figure
kernel  = gaussian2Dfilter([1000 1],[200 1]);
col = linspecer(2,'jet');
%close all
plot(ts_data,nanconvn(signal_DFoF,kernel'))
hold on
ylim([-1 5])

for i = 1:2
   plot([ev_tims{i} ev_tims{i}]',[ones(length(ev_tims{i}),1) 2*ones(length(ev_tims{i}),1)]','color',col{i})
end


%%

%get index at time 0
[a,ix_0] = bestmatch(-5,ts_PETH);

%get index at time 30
[a,ix_30] = bestmatch(2.5,ts_PETH);

clear u_obj


newObj1 = newObj;
for i = 1:size(newObj,1)
    
    for j = 1:size(newObj,2)
        
        for k = 1:size(newObj1{i,j},1)
        newObj1{i,j}(k,:) = nanconvn( newObj1{i,j}(k,:) ,kernel');
        
        
        u_obj{i,j}(k) = nanmean( newObj1{i,j}(k,ix_0:ix_30),2);




        end
    end
end


%%


col = linspecer(5,'heat');
close all
for j  =  1:size(u_obj,1)
    figure
for i = 1:5
plot(u_obj{j,i}(1),u_obj{j,i}(2),'o','color',col{i})
shg
xlim([-.25 1])
ylim([-.25 1])
hold on
plot(-.15:.01:1,-.15:.01:1)


plot(u_obj{j,i}(1),u_obj{j,i}(3),'x','color',col{i})
%plot(u_obj{j,i}(1),u_obj{j,i}(4),'d','color',col{i})


xlabel('GRABDA 1st encounter')
ylabel('GRABDA Nth encounter')
end

end

%%
%close all
figure
plot(ts_PETH,nanmean(loopOn1),'b')
hold on
plot(ts_PETH,nanmean(loopOff1), 'k')
xlim([-300 300])

xlabel ('Time (s)') 
ylabel ('Z-Score') 
legend ('Start','Finish')
set(gca,'box','off')
%%

%get index at time 0
[a,ix_0] = bestmatch(0,ts_PETH);

%get index at time 100
[a,ix_100] = bestmatch(60,ts_PETH);


u_response_On = sum(loopOn1(:,ix_0:ix_100),2);
u_response_Off = sum(loopOff1(:,ix_0:ix_100),2);


%%


% make scatter plot of homecage vs novel

figure
 plot(mean(loopOff1(:,ix_0:ix_100),2),mean(loopOn1(:,ix_0:ix_100),2),'.')
 hold on
  plot(-.3:.1:1,-.3:.1:1)
  xlim([-.3 1])
  ylim([-.3 1])

 
  

xlabel ('Start') 
ylabel ('Finish') 
set(gca,'box','off')



