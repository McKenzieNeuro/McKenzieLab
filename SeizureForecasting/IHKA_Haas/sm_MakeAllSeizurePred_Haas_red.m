% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
ClassifierFileOutput =  'G:\data\IHKA_Haas\classification_red.mat';
load(ClassifierFileOutput)

FeatureFileOutput = 'G:\data\IHKA_Haas\features_red.mat';
load(FeatureFileOutput,'sesID')

%%

warning off
% loop over files to predict
fils = getAllExtFiles('G:\data\IHKA_Haas\LFS','dat',1);
sessions  =unique(fileparts(fils));


fils = getAllExtFiles('R:\IHKA_Haas\LFS','csv',1);
sz_annots = fils(contains(fils,'Hz'));
und = regexp(sz_annots,'_');
sl = regexp(sz_annots,'\');

basename = cellfun(@(a,b,c) a(c(end)+1:b(4)-1),sz_annots,und,sl,'UniformOutput',false);
clear gd_ses
for i = 1:length(basename)
     gd_ses{i,1} =  sz_annots{i};
   tmp =  sessions(cellfun(@any,regexp(sessions,basename{i})));
   
   if ~isempty(tmp)
       
       [~,bb] = fileparts( tmp{1});
       gd_ses{i,2}  = [tmp{1} filesep bb];
   end
end

kp = all(~cellfun(@isempty,gd_ses),2);
gd_ses = gd_ses(kp,:);

%%
for i = 10:size(gd_ses,1)
    
    %get times used in training
    
    [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred_Haas_red(gd_ses{i,2},gd_ses{i,1},rusTree,ops);
    [dirOut,fileOut] = fileparts(gd_ses{i,2});
    save([dirOut filesep fileOut '_predict_red.mat'],'estimateLabel','trueLabel','time2seizure','inTrainingSet','seizure_start')
    disp([' saved: ' dirOut filesep fileOut '_predict_red.mat'])
end
fils = getAllExtFiles('G:\data\IHKA_Haas\LFS','dat',1);
incomplete_ses = cellfun(@(a) a(1:end-6),fils(contains(fils,'pre')),'uni',0);

for i = 2:length(incomplete_ses)
    
    [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred_Haas_red(incomplete_ses{i},[],rusTree,ops);
    [dirOut,fileOut] = fileparts(incomplete_ses{i});
    save([dirOut filesep fileOut '_predict_red.mat'],'estimateLabel','trueLabel','time2seizure','inTrainingSet','seizure_start')
    disp([' saved: ' dirOut filesep fileOut '_predict_red.mat'])
end
%%
fils = getAllExtFiles(fileparts(dirOut),'mat',1);
fils = fils(contains(fils,'predict_red')&~contains(fils,'pre-cycle'));
k  = gaussian2Dfilter([10000 1],10);

ts = 1:size(estimateLabel_all,2);
LFS = [1800+30 2400; 3001+30 3600;4201+30 4800;5401+30 6000;6601+30 10000];
stim = excludeEpochs([0 10000],LFS );
in_stim = InIntervals(ts,stim);

for j = 1:6
for i = 1:length(fils)
    v = load(fils{i});
    
    estimateLabel_all(i,:,j) = nanconvn(v.estimateLabel(1:10000,1)==j,k,'nanout',true);
end
end

estimateLabel_all(:,in_stim,:) = nan;
fils = getAllExtFiles(fileparts(dirOut),'mat',1);
fils = fils(contains(fils,'predict_red')&contains(fils,'pre-cycle'));
k  = gaussian2Dfilter([10000 1],1);
clear estimateLabel_pre
for i = 1:length(fils)
       v = load(fils{i});
       
for j = 1:6

 
    estimateLabel_pre(i,j) = nanmean(nanconvn(v.estimateLabel(:,1)==j,k));
end
end



clear bl_mean
for i = 1:size(estimateLabel_all,1)
    for j = 1:size(estimateLabel_all,3)
        [mat,vs] = Epoch2Mat(ts,LFS,estimateLabel_all(i,:,j));
        
        bl_mean(i,:,j) = cellfun(@mean,vs);
    end
end

%%
ix = 5;
close all
plotMeanSEM(1:10000,estimateLabel_all(:,:,ix),'r')
hold on

shg
plot([0 10000],[nanmean(estimateLabel_pre(:,ix)) nanmean(estimateLabel_pre(:,ix))],'--','color','r')
plot([LFS(:) LFS(:)],[0 1],'k','linewidth',4)
ylim([0 .4])

figure
ix = 1:2;

plotMeanSEM(1:10000,nanmean(estimateLabel_all(:,:,ix),3),'b')
hold on
shg
plot([0 10000],[nanmean(nanmean(estimateLabel_pre(:,ix),2)) nanmean(nanmean(estimateLabel_pre(:,ix),2))],'--','color','b')
plot([LFS(:) LFS(:)],[0 1],'k','linewidth',4)
ylim([.15 .4])
%%
close all
figure
k = gaussian2Dfilter([1 100],[ 1 50]);
%estimateLabel1 = estimateLabel;
%estimateLabel1(inTrainingSet) = nan;
for i = 1:6
    subplot(6,1,i)
   
    hold on
    plot([seizure_start seizure_start],[0 1],'--','color','r')
     plot(nanconvn(estimateLabel==i,k'),'k')
end
    %%
    
    alldat =  cell2mat(dat');
    group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');
    testix = ops.rix(ops.N+1:end);
    training = alldat(testix,:);
    
    [pred,score] = predict(rusTree,training);

actual_Y = group(testix);
predicted_Y = pred;
C = confusionmat(actual_Y,predicted_Y);

%%
figure
imagesc(C./sum(C,2),[0 1])
set(gca,'yticklabel',{'3hrs','1hr','10min','10s','sz','post'})
set(gca,'xticklabel',{'3hrs','1hr','10min','10s','sz','post'})
xlabel('Predicted class')
ylabel('Real class')

set(gca,'ydir','normal')

%%

Cr = nan(6,6,1000);
for i = 1:1000
    
    % get conf
   actual_Yr =  actual_Y(randsample(1:length(actual_Y),length(actual_Y)));
   tmp = confusionmat(actual_Yr,predicted_Y);
   
   tmp = tmp./sum(tmp,2);
    Cr(:,:,i) =tmp;
    i
end
%%
 

lab = [num2cell(ops.bins),{'sz'},{'post'}];
ok = C./sum(C,2);
ok1 = eye(6)';
okRu = nanmean(Cr,3);
close all
ax  =tight_subplot(6,1);
ixx = 1;
for i = [1:6]
    
    
    axes(ax(i))

semilogx([-ops.bins(1:4) -.1 -.01], 100* ((ok(i,1:6)- okRu(i,1:6))./okRu(i,1:6)),'o-')
    hold on
   % semilogx([-ops.bins(1:4) -.1], 100* ((ok(i,1:5)- ok1(i,1:5))./ok(i,1:5)),'x-')
    x = [-ops.bins(1:4) -.1 -.01]';
    z1 = -prctile(Cr(i,:,:),5,3);
z2 = prctile(Cr(i,:,:),95,3);
    xx = [x;flipud(x)];
zz = 100*[z1';z2'];

col = 'k';
plot([-ops.bins(1:4) -.1 -.01],[0 0 0 0 0 0],'k')
patch(xx,zz,col,'facealpha',.25,'edgecolor','none')

    % semilogx([-ops.bins(1:4) -.1], okRu(i,1:5),'--','color','r')

if i<=4
title(['Model predicts ' num2str(lab{i}) 's to seiz.'])
elseif i==5
title(['Model predicts seiz.'])    
else
    title(['Model predicts post ict.'])    
end
if i ==6
xlabel('Prediction')
end
%xlim([-700 1])
xlim([-1800 0])
if i <5
%ylim([0 .45])

set(gca,'xtick',[-ops.bins(1:4) -.1 -.01],'xticklabel',[])
else
   %ylim([0 .8])
%   ylim(100*[-1 2])
   set(gca,'xtick',[-ops.bins(1:4) -.1 -.01],'xticklabel',lab(1:6))
   xlabel('Real time to seizure')
end
ylim(100*[-1 3.5])
ixx = ixx+1;
end


%%

