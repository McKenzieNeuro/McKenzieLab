% this function takes a trained model with accompanying feature definition
% and calculates the predicted time to seizure. pulls data from both the
% raw time series, and the pre-calculated feature space
%
%
%
% see: sm_MakeAll_getPowerPerChannel,sm_PredictIHKA_getAllFeatures , sm_PredictIHKA



%%
%load classifier, loads 'ops','rusTree','sessions'
ClassifierFileOutputDir =  'E:\Dropbox\Scn1a_EEG_for_SM\analysis\classification';


FeatureFileOutput = 'E:\Dropbox\Scn1a_EEG_for_SM\analysis\features.mat';
load(FeatureFileOutput)

%%

warning off
% loop over files to predict
for i = 1:10%size(sessions,1)
    ClassifierFileOutput = [ClassifierFileOutputDir filesep 'classification_' num2str(i) '.mat'];
    load(ClassifierFileOutput)
    featureFile =  sessions{i,2};
    seizureFile = sessions{i,1};
    ev = LoadEvents(seizureFile);
    sz_off = ev.time(contains(ev.description,'sz_off'));
    %get times used in training
   % trainingTime = sort(cell2mat(cellfun(@(a) a(a(:,1)==i,2),sesID,'UniformOutput',false)'));
    trainingTime = [];
    estimateLabel = sm_getSeizurePred_Mattis(featureFile,rusTree,sz_off,ops);
    outfil = [ClassifierFileOutputDir 'predict_' num2str(i) '.mat'];
    save(outfil,'estimateLabel')
    disp([' saved: ' outfil])
end

%%
%%

% load all predictions
fils = getAllExtFiles(ClassifierFileOutputDir,'mat',1);
kp = contains(fils,'predict') & ~contains(fils,'before');
fils  = fils(kp);

estimateLabel =[];
k = gaussian2Dfilter([1000 1],1);
for i = 1:length(fils)
    
    v= load(fils{i});

    tmp = v.estimateLabel(:,1);
    tmp = nanconvn(tmp==1,k);
    tmp = nanPad(tmp',601);
    estimateLabel = [estimateLabel;tmp];
    i
end

%%



estimateLabel1 = estimateLabel;
estimateLabel1(isnan(estimateLabel)) = -1;
figure
plotMeanSEM(0:600,estimateLabel,'k')
xlim([0 500])
xlabel('Time from seizure end (s)')
ylabel('Pr(inter-ictal | LFP)')
set(gca,'fontsize',16)
%%
idx = [];
for i = 1:size(estimateLabel,1)
    
    tmp = find(isnan(estimateLabel(i,:)),1,'first');
    
    if ~isempty(tmp)
        
        idx(i) = tmp;
    else
        idx(i) = nan;
    end
end

[~,b]  = sort(idx);
%%
close all
figure
imagesc(estimateLabel1(b,:))
caxis([-.01, 1]);
x = colormap;
x = [.5 .5 .5;x];
colormap(x)
xlim([0 500])
xlabel('Time from seizure end (s)')
ylabel('Subject')
set(gca,'fontsize',16)

%%
k = gaussian2Dfilter([1000 1],10);
estimateLabel2 =[];
for i = 1:length(fils)
    
    v= load(fils{i});

    tmp = v.estimateLabel(:,1);
    tmp = nanconvn(tmp==1,k);
    tmp = nanPad(tmp',601);
    estimateLabel2 = [estimateLabel2;tmp];
    i
end



s = fittype('1/(1+exp(-a*(x-b)))', 'Coefficients',{'a','b'}, 'Independent','x', 'Dependent','y')

%%
x = [1:500]';
for i = 1:size(estimateLabel2,1)
y = estimateLabel2(i,1:500)';
kp = ~isnan(y);
sfit =  fit(x(kp), y(kp), s,'lower',[0  0],'upper',[10  600]);
beta(i,:) = coeffvalues(sfit);
end

%%
mu = mean(beta(:,2));
sigma = std(beta(:,2));
pct = [1 -.2];
ratio_within_across = 4;
nsub =  10;

alpha=sm_power_repeated_anova(mu,sigma,nsub,pct,ratio_within_across)
 
 