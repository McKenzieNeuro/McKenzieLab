% this function will load a feature file to train the RUSboost model to
% predict which time bin the data originated from relative to the nearest
% seizure


% relies on the output of sm_PredictIHKA_getAllFeatures which calculates 
% and saves the feature file 
% features are loaded into the variable 'dat' which is is a cell array
% with each element as a category (time bin) to predict




FeatureFileOutput = 'G:\data\IHKA_Haas\features.mat';
load(FeatureFileOutput)
ops.ClassifierFileOutput =  'G:\data\IHKA_Haas\classification.mat';

%%

% extract the features (each group is an element of the cell array)
ops.nGroup  = length(dat);
training = cell2mat(dat');

%define the groups (1:length(dat))
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');


% randomly sort for cross validation
ops.rix  = randsample(1:length(group),length(group));
training = training(ops.rix,:);
group = group(ops.rix);


%set up classifer


ops.N = round(size(training,1)/2);         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';

%train model
rusTree = fitcensemble(training(1:ops.N,:),group(1:ops.N,:),'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate);




%save (does not resave features, does save filename)
save(ops.ClassifierFileOutput,'ops','rusTree','training','group','sessions','-v7.3')
%%
%C1 = [];
ses = cell2mat(sesID(1:6)');
subjID = ses(:,1);
b = nchoosek(1:19,5);
%C2 =[];
for i = 1:38


tic
ops.nGroup  = length(dat(1:6));
training = cell2mat(dat(1:6)');

%define the groups (1:length(dat))
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat(1:6),num2cell(1:length(dat(1:6))),'uni',0)');
kp = ~ismember(subjID,i);

training = training(kp,:);
group = group(kp);


% randomly sort for cross validation
%ops.rix  = randsample(1:length(group),length(group));
%training = training(ops.rix,:);
%group = group(ops.rix);


%set up classifer


ops.N = round(size(training,1));         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';

%train model
rusTree = fitcensemble(training,group,'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate);


kp = ~kp;


    alldat =  cell2mat(dat(1:6)');
    group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat(1:6),num2cell(1:length(dat(1:6))),'uni',0)');
    
    alldat = alldat(kp,:);
    group = group(kp);
    testix = ops.rix(ops.N+1:end);
    training = alldat(:,:);
    
    [pred,score] = predict(rusTree,training);

actual_Y = group(:);
predicted_Y = pred;
C = confusionmat(actual_Y,predicted_Y);
if size(C,1) ==5 && size(C,2) ==5
    C = [nan(1,5);C];
    C = [nan(6,1) C];
end

C2(:,:,i) = C./nansum(C,2);
toc 
i
end
%%
%close all
figure

imagesc(C1,[0 1])
set(gca,'yticklabel',{'1800','600','300','10s','sz','post'})
set(gca,'xticklabel',{'1800','600','300','10s','sz','post'})
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
ylim(100*[-1 4.5])
ixx = ixx+1;
end