% get all features

featureDir = 'G:\data\isip\oneTreeData\features1';
fils = getAllExtFiles(featureDir,'mat',1);


%%
dat = cell(1,6);
sesID = cell(1,6);

%initialize

for j = 1:6
   
    dat{j} = nan(1e5,6060);
      sesID{j} = nan(1e5,2);
end
%%


% load all features (v.dat) and session IDs (v.sesID)
N = zeros(1,6);
for i = 1:length(fils)
    try
    v=load(fils{i});
    for j = 1:6
        
        ix = (N(j)+1 ):  ( N(j)  + length(v.sesID{j}));
        dat{j}(ix,:) = v.dat{j};
        sesID{j}(ix,:) = v.sesID{j};
        
        N(j) =   N(j)  + length(v.sesID{j});
    end
    
    
    i
    end
    
end


%%

%delete all initialized rows that weren't used
for j  =1:6
    kp = ~all(isnan(dat{j}),2);
    
    dat{j} = dat{j}(kp,:);
      sesID{j} = sesID{j}(kp,:);
end




%%

%set up features and labels for classification (ses I exlcude the post
%ictal state (6th element of the cell array)

% get session IDs
sesID1 = cell2mat(sesID(1:5)');

%define how many classes
ops.nGroup  = length(dat(1:5));

%get features
training = cell2mat(dat(1:5)');

%sanitize features
training(isinf(training)) = nan;

training = training(:,~all(isnan(training)));

 u_dim = nanmean(training);
for i = 1:size(training,2)
    
    training(isnan(training(:,i)),i) = training(i);
end



%define the groups labels
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat(1:5),num2cell(1:length(dat(1:5))),'uni',0)');


% randomly sort for cross validation
ops.rix  = randsample(1:length(group),length(group));
training = training(ops.rix,:);
group = group(ops.rix);
sesID1 = sesID1(ops.rix,:);
subj = cellfun(@(a) str2num(a(end-24:end-17)),v.sessions(sesID1(:,1),1));

%define subset of data to train (odd sessions)
kp = ismember(sesID1(:,1),[1:2:max(sesID1(:,1))]);


%set up classifer
ops.N = round(size(training,1)/2);         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N/10);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';


 % train model
 rusTree = fitcensemble(training(kp,:),group(kp,:),'Method',ops.Method, ...
'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate,'ScoreTransform','logit');

%get prediction on held out data

[pred,score] = predict(rusTree,training(~kp,:));


%get confusion matrix
actual_Y = group(~kp);
predicted_Y = pred;
C = confusionmat(actual_Y,predicted_Y);
%%
Cr = nan(5,5,1000);
for i = 1:1000
    
    % get conf
   actual_Yr =  actual_Y(randsample(1:length(actual_Y),length(actual_Y)));
   tmp = confusionmat(actual_Yr,predicted_Y);
   
   tmp = tmp./sum(tmp,2);
    Cr(:,:,i) =tmp;
    i
end

%%
figure
for i = 1:5
    hold on
    ok = accumarray(sesIDkp(actual_Y==i),predicted_Y(actual_Y==i)==actual_Y(actual_Y==i),[],@nanmean,[]);
    ok = ok(ok>0);
    plot(histc(ok,0:.1:1)/sum(ok))
end
    
%%

% project in lowD space to explore feature PDF
dat1  =cell2mat(dat');
 a = tsne(dat1);
 ix = cell2mat(cellfun(@(a,b) a*ones(sum(~any(isnan(b),2)),1),num2cell(1:6),dat,'uni',0)');
 
 
 % Haojie stop here!
 %%%%%%%%%%%
 
%%

% if yHat are your predictions and yval are your y true then

for i = 1:5
tp = sum((predicted_Y == i) & (actual_Y == i));
fp = sum((predicted_Y == i) & (actual_Y ~= i));
fn = sum((predicted_Y ~= i) & (actual_Y == i));
tn = sum((predicted_Y ~= i) & (actual_Y ~= i));
precision(i) = tp / (tp + fp);
sensitivity(i) = tp/(tp + fn);  %TPR;
specificity(i) = tn/(tn + fp);  %TNR
FPR(i) = fp/(tn+fp);
Accuracy(i) = (tp+tn)./(tp+fp+tn+fn);
recall(i) = tp / (tp + fn);
F1(i) = (2 * precision(i) * recall(i)) / (precision(i) + recall(i));
end
%%
    ops.bins = [ 600 300 100 10];

lab = [num2cell(ops.bins),{'sz'},{'post'}];
ok = C./sum(C,2);
ok1 = eye(5)';
okRu = nanmean(Cr,3);
close all
ax  =tight_subplot(5,1);
ixx = 1;
for i = [1:5]
    
    
    axes(ax(i))

semilogx([-ops.bins(1:4) -.1], 100* ((ok(i,1:5)- okRu(i,1:5))./okRu(i,1:5)),'o-')
    hold on
   % semilogx([-ops.bins(1:4) -.1], 100* ((ok(i,1:5)- ok1(i,1:5))./ok(i,1:5)),'x-')
    x = [-ops.bins(1:4) -.1]';
    z1 = -prctile(Cr(i,:,:),5,3);
z2 = prctile(Cr(i,:,:),95,3);
    xx = [x;flipud(x)];
zz = 100*[z1';z2'];

col = 'k';
plot([-ops.bins(1:4) -.1],[0 0 0 0 0],'k')
patch(xx,zz,col,'facealpha',.25,'edgecolor','none')

    % semilogx([-ops.bins(1:4) -.1], okRu(i,1:5),'--','color','r')

if i<=4
title(['Model predicts ' num2str(lab{i}) 's to seiz.'])
else
title(['Model predicts seiz.'])    
end
if i ==6
xlabel('Prediction')
end
xlim([-700 1])

if i <5
%ylim([0 .45])
%ylim(100*[-.5 1.5])
set(gca,'xtick',[-ops.bins(1:4) -.1],'xticklabel',[])
else
   %ylim([0 .8])
%   ylim(100*[-1 2])
   set(gca,'xtick',[-ops.bins(1:4) -.1],'xticklabel',lab(1:5))
   xlabel('Real time to seizure')
end
ixx = ixx+1;
end
%%


figure
imagesc(C./sum(C,2),[0 .4])
set(gca,'xtick',1:5,'xticklabel',[num2cell(ops.bins),{'sz'}],'ytick',1:5,'yticklabel',[num2cell(ops.bins),{'sz'}])
xlabel('Estimate')
ylabel('True')

%%
dat1  =cell2mat(dat');
dat1 = dat1(:,~all(isnan(dat1)));
dat1(isinf(dat1)) = nan;
 u_dim = nanmean(dat1);
for i = 1:size(dat1,2)
    
    dat1(isnan(dat1(:,i)),i) = u_dim(i);
end

[a,b,c] = pca(dat1);
 a = tsne(dat1);
 ix = cell2mat(cellfun(@(a,b) a*ones(sum(~all(isnan(b),2)),1),num2cell(1:6),dat,'uni',0)');
%%



%for SAP

subj = cellfun(@(a) str2num(a(end-24:end-17)),v.sessions(sesID1(:,1),1));
%%

[N,~,~,nIX] = histcn(a(: ,1:2),-100:100,-100:100);


ID  =1;
thres=prctile(linearize(nanconvn(histcn(a(ix==ID,:),-100:100,-100:100),k)/sum(ix==ID)),99);
bnds = bwboundaries(nanconvn(histcn(a(ix==ID ,:),-100:100,-100:100),k)/sum(ix==ID)>thres);
in = inpolygon(nIX(:,1),nIX(:,2),bnds{6}(:,1),bnds{6}(:,2));

goodbnds1 = bnds{6};

ID  =4;
thres=prctile(linearize(nanconvn(histcn(a(ix==ID,:),-100:100,-100:100),k)/sum(ix==ID)),99);
bnds = bwboundaries(nanconvn(histcn(a(ix==ID ,:),-100:100,-100:100),k)/sum(ix==ID)>thres);
goodbnds4 = bnds{2};
in1 = inpolygon(nIX(:,1),nIX(:,2),bnds{2}(:,1),bnds{2}(:,2));



%%

sesID1 = cell2mat(sesID(1:6)');
% get boundaries

figure
k  = gaussian2Dfilter([1000 1000],[5 5]);
ax = tight_subplot(1,6); 
ID  =4;
thres=prctile(linearize(nanconvn(histcn(a(ix==ID,:),-100:100,-100:100),k)/sum(ix==ID)),99);
bnds = bwboundaries(nanconvn(histcn(a(ix==ID ,:),-100:100,-100:100),k)/sum(ix==ID)>thres);
for i = 1:6
    axes(ax(i))
    imagesc(nanconvn(histcn(a(ix==i ,1:2),-100:100,-100:100),k,'nanout',true)/sum(ix==i),[0 .00035])
hold on
% for j = 1:length(bnds)
% plot(bnds{j}(:,2),bnds{j}(:,1),'w')
% end
plot(goodbnds1(:,2),goodbnds1(:,1),'r')
plot(goodbnds4(:,2),goodbnds4(:,1),'w')
axis off
axis square
colormap('jet')
end
%%


%%


channel  = [ ...
    {'fp1'} ; ...
    {'fp2'} ; ...
    {'f7'} ; ...
    {'f8'} ; ...
    {'t3'} ; ...
    {'t4'} ; ...
    {'t5'} ; ...
    {'t6'} ; ...
    {'o1'} ; ...
    {'o2'} ; ...
    {'f3'} ; ...
    {'f4'} ; ...
    {'c3'} ; ...
    {'c4'} ; ...
    {'p3'} ; ...
    {'p4'} ; ...
    {'a1'} ; ...
    {'a2'} ; ...
    {'fpz'} ; ...
    {'fz'} ; ...
    {'cz'} ; ...
    {'pz'} ; ...
    {'oz'} ; ...
    {'ekg'} ; ...
    
    ];
notEEGCh = [ ...
    {'a1'} ; ...
    {'a2'} ; ...
    {'ekg'} ; ...
 
    ];
goodCh = (~ismember(channel,notEEGCh));
%%

ok = predictorImportance(rusTree);
%%
%ch 18 is missing

for j = 1:20
    for i =1:24
        if i~=18
            fr(i) = ok((i-1)*20+j);
        end
    end
    figure
    [z,map]=eegplot(fr(goodCh),[],[],[],'cubic' ,[],channel(goodCh),[-.5 1.5]);
    waitforbuttonpress
    close all
    title(num2str(j))
end
%%


% plot low freq
close all
ops.freqs = logspace(log10(.5),log10(200),20);
dat1(isinf(dat1)) = nan;
figure
ax = tight_subplot(1,2);
i=2
fr1a =[]; fr2a =[];
for i = 2 
    fr = dat1(:,(i-1)*20+1:20:(i-1)*20+20*24);
    fr = fr(:,goodCh);
    fr1 = fr(in,:);
    fr2 = fr(in1,:);
    fr1(isnan(fr1)) = 0;
    fr2(isnan(fr2)) = 0;
    fr1 = nanmean(fr1);
    fr2 = nanmean(fr2);
    fr1a = [fr1a;fr1];
    fr2a = [fr2a;fr2];
end
fr1a = nanmean(fr1a,1);
fr2a = nanmean(fr2a,1);
axes(ax(1));
[z,map]=eegplot(fr1a,[],[],[],'cubic' ,[],channel(goodCh),[-.5 1.5]);
axis off
axes(ax(2));
[z,map]=eegplot(fr2a,[],[],[],'cubic' ,[],channel(goodCh),[-.5 1.5]);
axis off


%%
usubjs = unique(subj);
pp = nan(188,size(dat1,2));

    
    for va = 1:size(dat1,2)
        datTable = table(dat1(:,va), sesID1(:,1),subj,double(in6),double(in4), ...
            'VariableNames',{'neural','subjSes','subj','in6','in4'} );
        for n = 1:188
        
        
        %neur
      
        try
            kpSubjs = [randsample(unique(subj(in6)),n);randsample(unique(subj(in6)),n)];
            kp = (in6 | in4) & ismember(subj,kpSubjs);
            glme1 = fitglme(datTable(kp,:),'in6 ~ 1   + neural +   (1|subjSes:subj) + (1|subj)','Distribution','Binomial');
            pp(n,va) = glme1.Coefficients.pValue(2);
        end
     
        end
    va
end
    
%%
k  = gaussian2Dfilter([1 100000],1);

figure
semilogx(0:5000,fliplr(avghist(time2seize,nanconvn(double(estimateLabel==1),k'),-5000:0)))

hold on
semilogx(0:5000,fliplr(avghist(time2seize,nanconvn(double(estimateLabel==3),k'),-5000:0)))
semilogx(0:5000,fliplr(avghist(time2seize,nanconvn(double(estimateLabel==4),k'),-5000:0)))
semilogx(0:5000,fliplr(avghist(time2seize,nanconvn(double(estimateLabel==5),k'),-5000:0)))