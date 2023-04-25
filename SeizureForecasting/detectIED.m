% train IED detector

topDir = 'R:\DGregg\NeuralData\LTD 10.0';
fils = getAllExtFiles(topDir,'IED',1);
kp = contains(fils,'autoDetect');
fils = fils(kp);
[a,b ] =fileparts(fils);
%%
tim=0;
IED_ex = [];
for i = 1:length(a)
    cd(a{i})
    xml  = LoadXml('amplifier.xml');
    ev = LoadEvents('autoDetect.evt.IED');
    nCh = xml.nChannels;
    for j = 1:length(ev.time)
        d = LoadBinary('amplifier.lfp','nchannels',nCh,'channels',1:nCh,'start',ev.time(j)-.2,'duration',.5,'frequency',1250);
        [~,b] = max(max(abs(d-repmat(d(:,1),1,nCh))));
        
        IED_ex = [IED_ex;d(:,b)'];
    end
    i
end
%%
IED1 = double(IED_ex);

for i = 1:size(IED_ex,1)
    
    
    IED1(i,:) =  BandpassFilter(double(IED_ex(i,:)),1250,[10 300]);
    i
end
ix = kmeans(double(IED1),25);
%%

close all
figure
for i = 1:25
    
    subplot(5,5,i)
    plot(nanmean(IED1(ix==i,:)))
    title(num2str(i))
end

kp = ismember(ix,[2 3 6 7 10 11 13 17 18 ]);

%%
IED_good = IED1(kp,:);
IED_bad = IED1(~kp,:);


training = [IED_good;IED_bad];
group = [true(sum(kp),1);false(sum(~kp),1)];
kp1 = ismember(1:length(group),1:floor(length(group)/2)) ;
idx  = randsample(1:length(group),length(group));

training = training(idx,:);
group = group(idx);
%set up classifer
ops.N = round(size(training,1)/2);         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N/10);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';


% train model
rusTree = fitcensemble(training(kp1,:),group(kp1,:),'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate,'ScoreTransform','logit');

%get prediction on held out data

[pred,score] = predict(rusTree,training(~kp1,:));
actual_Y = group(~kp1 >0);
predicted_Y = pred;
C = confusionmat(actual_Y,predicted_Y);

%%
%%
%dat = nan(50,1250,14832);
for k1 = 11956:length(IED)
    tmp = LoadBinary('continuous.lfp','nchannels',43','channels',14,'frequency',1250,'start',IED(k1)-5,'duration',10);
    datT = awt_freqlist(double(tmp),1250,logspace(log10(1),log10(300),50));
    datT = abs(datT(1:10:end,:));
    dat(:,:,k1) = datT';
    k1
    
end
%%

ok  =dat(:,:,IED<3600*12);
ok1  = dat(:,:,IED>3600*12);


close all
figure
imagesc((1:1250)/125 -5,[],nanmean(ok,3),[0 400])
set(gca,'ytick',1:10:50','yticklabel',freqs(1:10:end),'ydir','normal')

figure
imagesc((1:1250)/125 -5,[],nanmean(ok1,3),[0 400])
set(gca,'ytick',1:10:50','yticklabel',freqs(1:10:end),'ydir','normal')


%%
pp = nan(size(ok,1),size(ok,2));
for i = 1:size(ok,1)
    for j = 1:size(ok,2)
        
        
     pp(i,j) =    ranksum(squeeze(ok(i,j,:)),squeeze(ok1(i,j,:)));
        
    end
    i
end