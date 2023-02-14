% get all predictions


topDirs  = [ ...
    {'E:\data\IHKA'} ; ...
    {'F:\data1\IHKA'} ; ...
    ];

fils = [];

for i = 1:length(topDirs)
    fils  = [fils;getAllExtFiles(topDirs{i},'mat',1)];
end

kp = contains(fils,'predict');

fils = fils(kp);
%%
estimateLabel = [];trueLabel =[];time2seize = [];
for i = 1:length(fils)
    
    v = load(fils{i});
    
    siz = min(length(v.trueLabel),length(v.estimateLabel));
    ix  = [1:siz]';
    estimateLabel = [estimateLabel;linearize(v.estimateLabel(ix))];
    if isfield(v,'time2seize')
      time2seize = [time2seize;linearize(v.time2seize(ix))];
    else
           time2seize = [time2seize;linearize(v.time2seizure(ix))];
    end
   
     trueLabel = [trueLabel;linearize(v.trueLabel(ix))];
    i
end

%%
clear ok
for i = 1:6
    
    ok{i} = CrossCorr(find(double(estimateLabel==i)),find(double(estimateLabel==i)),5,1001)/sum(estimateLabel==i)/5;
    ok{i}(501) = nan;
end
%%


clear ok1
for i = 1:6
    
    ok1{i} = CrossCorr(find(double(trueLabel==6)),find(double(estimateLabel==i)),5,1001)/sum(trueLabel==6)/5;
   
end


%%
kp = estimateLabel>1 & estimateLabel<5;
seizHour =  CrossCorr(find(double(kp )),find(double(kp)),5,1001)/sum(kp)/5;

kpt = trueLabel>2 & trueLabel<5;
ok1t =  CrossCorr(find(double(kpt )),find(double(kpt)),5,1001)/sum(kpt)/5;


ok2 =  CrossCorr(find(double(estimateLabel==4 )),find(double(kp)),5,1001)/sum(estimateLabel==4)/5;
%%

figure
col  = linspecer(2,'jet');
for i = [1 4 5 6]

            plot((-500:500)*5,ok1{i},'linewidth',4)
   
    hold on
end

%%
kpt = trueLabel>1 & trueLabel<6;
close all
figure
col  = linspecer(2,'jet');
for i = [1 4]
    switch i
        case 1
            semilogx((-500:500)*5,ok1{4}-nanmean(trueLabel==4),'color','k')
            
        case 4
            semilogx((-500:500)*5,ok1t-nanmean(kpt),'color','r')
    end
    hold on
end
%%
kp = estimateLabel>1 & estimateLabel<5;
ok2 =  CrossCorr(find(double(trueLabel==6 )),find(estimateLabel==4),5,1001)/sum(trueLabel==6 )/5;
ok3 =  CrossCorr(find(double(trueLabel==6 )),find(kp),5,1001)/sum(trueLabel==6 )/5;
figure
plot((-500:500)*5,ok2)
hold on
plot((-500:500)*5,ok3)

%%

figure
ax = tight_subplot(6,1);
for i  = 1:6
    o =  CrossCorr(find(double(trueLabel==6 )),find(estimateLabel==i),5,10001)/sum(trueLabel==6 )/5;
    axes(ax(i))
    plot(((-5000:5000)*5)/60,o)
    ylim([0 1])
end

%%

FP = CrossCorr(find(trueLabel~=5 & estimateLabel==4),find(trueLabel~=5 & estimateLabel==4),5,1001)/sum(trueLabel~=5 & estimateLabel==4)/5;
FP(501) = nan;
figure

semilogx((-500:500)*5,ok{4}-nanmean(estimateLabel==4),'color','k')

hold on
semilogx((-500:500)*5,FP - -nanmean(trueLabel~=5 & estimateLabel==4),'color','r')

%%
close all
thres = .75;
%t = nanconvn(kp,k');
k  = gaussian2Dfilter([1 100000],5);

t = nanconvn(estimateLabel==4,k');
t_on = find(diff(t>thres)>0);
t_off = find(diff(t>thres)<0);

dur = t_off-t_on;
bins =  sort([0 -logspace(log10(1),log10(1000),200)]);
semilogx(bins,(avghist(time2seize(diff(t>thres)>0),dur,bins,@nanmedian)),'o')

ylim([0 50])

%%
load('G:\data\isip\oneTreeData\features1.mat')
sesID  = cell2mat(sesID(1:5)');
%%

ops.nGroup  = length(dat(1:5));
training = cell2mat(dat(1:5)');
%define the groups (1:length(dat))
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat(1:5),num2cell(1:length(dat(1:5))),'uni',0)');
% randomly sort for cross validation
ops.rix  = randsample(1:length(group),length(group));
training = training(ops.rix,:);
group = group(ops.rix);
sesID = sesID(ops.rix,:);
kp = ismember(sesID(:,1),[1:2:max(sesID(:,1))]);
%set up classifer
ops.N = round(size(training,1)/2);         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';


 
 rusTree = fitcensemble(training(kp,:),group(kp,:),'Method',ops.Method, ...
'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate,'ScoreTransform','logit');
[pred,score] = predict(rusTree,training(~kp,:));

actual_Y = group(~kp);
predicted_Y = pred;
C = confusionmat(actual_Y,predicted_Y);


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

close all
lab = [num2cell(ops.bins),{'sz'},{'post'}];
ok = C./sum(C,2);
figure
ax  =tight_subplot(5,1);
for i = 1:5
    axes(ax(i))
plot([-ops.bins(1:4) -.1], ok(i,1:5),'o-')
set(gca,'xtick',[-ops.bins(1:4) -.1],'xticklabel',lab(1:5))
if i<=4
title(['Acutal: ' num2str(lab{i})])
else
title(['Acutal: ' lab{i}])    
end
if i ==6
xlabel('Prediction')
end
xlim([-700 1])

if i <5
ylim([0 .5])
else
    ylim([0 .8])
end
end
%%


figure
imagesc(C./sum(C,2),[0 .5])
set(gca,'xtick',1:5,'xticklabel',[num2cell(ops.bins),{'sz'}],'ytick',1:5,'yticklabel',[num2cell(ops.bins),{'sz'}])
xlabel('Estimate')
ylabel('True')

%%
dat1  =cell2mat(dat');
 a = tsne(dat1);
 ix = cell2mat(cellfun(@(a,b) a*ones(sum(~any(isnan(b),2)),1),num2cell(1:6),dat,'uni',0)');
%%

% get boundaries
close all
figure
k  = gaussian2Dfilter([1000 1000],[3 3]);
ax = tight_subplot(1,6); 

thres= .00022;% prctile(linearize(nanconvn(histcn(a(ix==4,:),-100:100,-100:100),k)/sum(ix==4)),99);
bnds = bwboundaries(nanconvn(histcn(a(ix==4,:),-100:100,-100:100),k)/sum(ix==4)>thres);
for i = 1:6
    axes(ax(i))
    imagesc(nanconvn(histcn(a(ix==i,:),-100:100,-100:100),k,'nanout',true)/sum(ix==i),[0 .0004])
hold on
for j = 1:length(bnds)
plot(bnds{j}(:,2),bnds{j}(:,1),'w')
end
axis off
axis square

end


%%
k  = gaussian2Dfilter([1 100000],1);

figure
semilogx(0:5000,fliplr(avghist(time2seize,nanconvn(double(estimateLabel==1),k'),-5000:0)))

hold on
semilogx(0:5000,fliplr(avghist(time2seize,nanconvn(double(estimateLabel==3),k'),-5000:0)))
semilogx(0:5000,fliplr(avghist(time2seize,nanconvn(double(estimateLabel==4),k'),-5000:0)))
semilogx(0:5000,fliplr(avghist(time2seize,nanconvn(double(estimateLabel==5),k'),-5000:0)))




%%
load('F:\data1\IHKA\classification.mat')

kp = ismember(1:length(group),1:ops.N);
[pred,score] = predict(rusTree,training(~kp,:));

actual_Y = group(~kp);
predicted_Y = pred;
C = confusionmat(actual_Y,predicted_Y);
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

semilogx([-ops.bins(1:4) -5 -1], 100* ((ok(i,1:6)- okRu(i,1:6))./okRu(i,1:6)),'o-')
    hold on
   % semilogx([-ops.bins(1:4) -.1], 100* ((ok(i,1:5)- ok1(i,1:5))./ok(i,1:5)),'x-')
    x = [-ops.bins(1:4) -5 -1]';
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
%xlim([-700 1])
xlim([-10800 -1])
if i <5
%ylim([0 .45])

set(gca,'xtick',[-ops.bins(1:4) -.1],'xticklabel',[])
else
   %ylim([0 .8])
%   ylim(100*[-1 2])
   set(gca,'xtick',[-ops.bins(1:4) -.1],'xticklabel',lab(1:5))
   xlabel('Real time to seizure')
end
ylim(100*[-.5 2.5])
ixx = ixx+1;
end
