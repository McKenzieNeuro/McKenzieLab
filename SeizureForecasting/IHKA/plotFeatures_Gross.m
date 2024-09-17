%%
%N = round(size(training,1)/2);         % Number of observations in the training sample
%label = predict(rusTree,training); % only N+1 CV

%zT  = zscore(training);

%Y = tsne(zT);
%%
% 8/21/2024
load('G:\data\IHKA_gross\features1.mat')

%%




warning off
%get all relevant timepoints around the seizure start and end

training = cell2mat(dat');
ses = cell2mat(sesID');
group = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),dat,num2cell(1:length(dat)),'uni',0)');
freqs = ops.freqs;
bins = ops.bins;

%%

ix=1;
subject_IDs = {'KA5','KA6','KA8','KA12','KA13','KA14','KA11'};
tot = 0;
for k1 =1:length(subject_IDs) %loops through all the animal IDs you provide at the top of the code
    
    
    [rows_to_extract, stim_metatable] = read_metatable_KAsz(subject_IDs{k1}, true);
    
    sz_o = cellfun(@str2num,stim_metatable.Seizures_sec_(rows_to_extract),'UniformOutput',false);
    %sz_off = cellfun(@str2num,stim_metatable.SeizuresEnd_sec_(rows_to_extract),'UniformOutput',false);
    ID = [stim_metatable(rows_to_extract,:).AnimalIdentification stim_metatable(rows_to_extract,:).TrialNumber];
    for k2 = 1:size(ID,1)
        
        str = [ID{k2,1} filesep ID{k2,1} '_' ID{k2,2}];
        idx =  find(contains(sessions,str));
        tmp = sz_o{k2} ;
        seizure_start{idx}(:,1) = tmp;
        seizure_start{idx}(:,2) = (1:length(tmp))+tot;
        tot = tot+length(tmp);
        
    end
    
    
end



%%
% get subject info

[~,subj] = fileparts(sessions);
d = regexp(subj,'_');
subj = cellfun(@(a,b) a(1:b(1)-1),subj,d,'uni',0);
[a,b,subjID] = unique(subj);
%%
ses = [ses nan(size(ses,1),2)];
for i = 1:length(sessions)
    idx = find(ses(:,1)==i);
    st = ses(idx,2);
    
    
    for j = 1:size(seizure_start{i},1)
        
        if j ==1
            kp = st>0 & st <seizure_start{i}(j,1)+700;
        else
            kp = st>seizure_start{i}(j-1,1)+700 & st <seizure_start{i}(j,1)+700;
            
        end
        
        ses(idx(kp),3) = seizure_start{i}(j,2);
        ses(idx(kp),4)  = subjID(i);
        
        
    end
    i
    
end
%%

% get mean per seizure
nsz = max(ses(:,3));
trainU = nan(nsz,size(training,2),6);
for i = 1:nsz
    
    for j = 1:6
        kp = ses(:,3) == i & group==j;
        
        trainU(i,:,j) = nanmean(training(kp,:));
        meanID(i,:) = nanmean(ses(kp,3:4));
    end
end
%%
close all
figure

c = [bins 0 1 2];
apr = ops.amp_idx;
idx1 = ops.ph_idx;

col = linspecer(4,'jet');
for j = 1:6
    subplot(3,2,j)
    for i = 1:3
        
        ix = (i-1)*20+(1:20);
        kp = group==j;
        %kp = kp2;
        semilogx(1,1)
        hold on
        plotMeanSEM(freqs(1:20), (trainU(:,ix,j)),col{i})
        hold on
        
        
        
        
        
    end
    
    
    if j <5
        title([num2str(round(c(j))) ' to ' num2str(round(c(j+1)))])
        ylim([-.4 .4])
    elseif j ==5
        title('seizure')
        ylim([-.4 4])
    else
        title('post')
        ylim([-.4 2])
    end
    if j ==5
        legend({'','C Cort','','I Cort','I H','','C H'},'Location','northwest')
    end
end

%%


nFeat = size(training,2);
Interval = (upSample(1:4,nsz));
colNames = {'SeizureID','Animal', 'Interval'};
featureColNames = arrayfun(@(x) sprintf('Feature_%d', x), 1:nFeat, 'UniformOutput', false);
allColNames = [colNames, featureColNames];
variableTypes = ['categorical' 'categorical' 'categorical' repmat({'double'}, 1, nFeat)];

nRows = size(trainU,1)*4;
tempTbl = table('Size', [nRows, length(allColNames)], 'VariableTypes', variableTypes, 'VariableNames', allColNames);



% Split features into separate columns
for i = 1:nFeat
    tempTbl.(featureColNames{i}) = linearize(squeeze(trainU(:,i,1:4)));
end

tempTbl.SeizureID = categorical(repmat(meanID(:,1),4,1));
tempTbl.Animal = categorical(repmat(meanID(:,2),4,1));

tempTbl.Interval = categorical(Interval);


for i = 1:nFeat
    featname = ['Feature_' num2str(i)];
    
    formula = ['Feature_' num2str(i) ' ~ Interval +  (1|Animal) + (1|Animal:SeizureID)'];
    lme = fitlme(tempTbl, formula);
    MSE(i) = lme.MSE;
    
    % Perform analysis of variance
    anovaTable = anova(lme, 'DFMethod','satterthwaite');
    
    % Extract the F ratio for the specific feature
    F_ratio_Feature(i) = anovaTable.FStat(   ismember(anovaTable.Term, 'Interval'   ));
    %F_ratio_Intercept(i) = anovaTable.FStat(   ismember(anovaTable.Term, '(Intercept)')   );
    i
end

%%
figure
for i = 1:3
    
    ix = (i-1)*20+(1:20);
    semilogx(freqs,(F_ratio_Feature(ix)))
    hold on
end
%%

idx = (20*3)+ ( 1: ( length(apr)*length(idx1))) ; % phase amp coupling

tt = nan(6,10,length(training));
for i = 1:length(training)
    
    tt(:,:,i) =  reshape(training(i,idx),length(apr),length(idx1));
    
end


%%
close all

figure
for j = 1:6
    subplot(3,2,j)
    kp = group ==j;
    % kp=kp1;
    imagesc(nanmedian(tt(:,:,kp),3),[0 .25])
    set(gca,'ytick',1:10,'yticklabel',round(freqs(apr)),'xtick',1:10,'xticklabel',round(freqs(idx1)*10)/10)
    colorbar
    xlabel('Phase freq.')
    ylabel('Amplitude freq.')
    
    if j <5
        title([num2str(round(c(j))) ' to ' num2str(round(c(j+1)))])
    elseif j ==5
        title('seizure')
    else
        title('post')
    end
    
end

%%
F_ratio_Feature_phase_amp = reshape(F_ratio_Feature(idx),length(apr),length(idx1));
%close all

figure
    kp = group ==j;
    % kp=kp1;
    imagesc(F_ratio_Feature_phase_amp)
    set(gca,'ytick',1:10,'yticklabel',round(freqs(apr)),'xtick',1:10,'xticklabel',round(freqs(idx1)*10)/10)
    colorbar
    xlabel('Phase freq.')
    ylabel('Amplitude freq.')
    
   

%%
figure
col = linspecer(6,'jet');
for j = 1:6
    subplot(3,2,j)
    for i = 1:2
        
        ix = (120 + (i-1)*20)  + (1:20);
        kp = group==j;
        %kp=kp1;
        semilogx(freqs(1:20), nanmean(training(kp ,ix)),'linewidth',3,'color',col{i})
        hold on
        
        
        ylim([0 1])
    end
    
    
    if j <5
        title([num2str(round(c(j))) ' to ' num2str(round(c(j+1)))])
    elseif j ==5
        title('seizure')
    else
        title('post')
    end
    if j ==1
        
        legend({'C Cort/I Cort','C Cort/I H','C Cort/C H','I Cort/I H','I Cort/C H','I H/C H'},'Location','northwest')
        
        
    end
end
%%
col = linspecer(6,'jet');

close all
figure
 for i = 1:2
        
        ix = (140 + (i-1)*20)  + (1:20);
        F_ratio_Feature_phase_coh = F_ratio_Feature(ix);

        %kp=kp1;
        semilogx(freqs(1:20), F_ratio_Feature_phase_coh,'linewidth',3,'color',col{i})
        hold on
        
        
       
    end
%get coherence
%%
figure
mat = histcn([ group((N+1):end) label(N+1:end)],1:6,1:6)./repmat(histc( group((N+1):end),1:6),1,6);
imagesc(mat)
colorbar
set(gca,'ydir','normal','xtick',1:6,'xticklabel',{'3hrs','1hr','10min','10s','Sz','PI'},'ytick',1:6,'yticklabel',{'3hrs','1hr','10min','10s','Sz','PI'})
title('Prob(Estimate| True)')
ylabel('True')
xlabel('Estimate')

%%
for i = 1:6
    TP = group ==i & label ==i;
    FP = group~=i & label ==i;
    FN = group ==i & label~=i;
    TN = group~=i & label~=i;
    
    Sp(i) = sum(TN)/(sum(TN) + sum(FP));
    
    Se(i) = sum(TP)/(sum(TP) + sum(FN));
end
figure
plot(Sp,'o')
hold on
plot(Se,'o')
legend('Specificity','Sensitivity')
set(gca,'ydir','normal','xtick',1:6,'xticklabel',{'3hrs','1hr','10min','10s','Sz','PI'})

ylabel('Sensitivity/Specificity')
xlabel('time to seizure')
xlim([0 7])

%%

k = gaussian2Dfilter([100 100],[ 2 2]);
figure
ix = -70:70;
for i = 1:6
    
    subplot(3,2,i)
    imagesc(ix,ix,nanconvn(histcn(Y(group==i,:),ix,ix),k)/sum(group==i))
end

%%
kp1 = Y(:,1)>-40 & Y(:,1) < -30 & Y(:,2)> -5 & Y(:,2)< 5;
kp2 = Y(:,1)>20 & Y(:,1) < 30 & Y(:,2)> -30 & Y(:,2)< -20;
%%


dirN = 'E:\Dropbox\UNM\Analysis\IHKA\data';
fils = getAllExtFiles(dirN,'mat',0);

fils = fils(1:end);



power_dirN = 'F:\data1\IHKA';
filP = dir(power_dirN);
filP = {filP.name}';

filenameps = 'E:\Dropbox\UNM\Analysis\IHKA\allPredict.ps';

O =[];
ts =[];
TMP =[];
sz_start =0;
for j = 1:length(fils)
    load(fils{j})
    [a,bn] = fileparts(fils{j});
    
    %get dat file
    % h= figure;
    % ax = tight_subplot(6,1);
    
    
    ok = cell2mat(cellfun(@(a) estimateLabel' ==a,num2cell(1:6),'uni',0)');
    
    
    sz_start = length(find(diff((trueLabel == 5) >0)))+sz_start;
    
    [~,b] = histc(ts1,-6000:0);
    k = gaussian2Dfilter([1 1000],[ 1 2]);
    ok = double(ok);
    
    for i = 1:6
        
        
        ok(i,:) = (nanconvn(ok(i,:) , k));
        
        
        
        
    end
    
    mix = min(size(ok,2),length(ts1));
    ts1 = ts1(1:mix);
    ok = ok(:,1:mix);
    
    O = [O ok];
    
    ts = [ts ts1(1:size(ok,2))];
    j
end


%%
dirN = 'E:\Dropbox\UNM\Analysis\IHKA\data';
fils = getAllExtFiles(dirN,'mat',0);

fils = fils(1:end);



power_dirN = 'F:\data1\IHKA';
filP = dir(power_dirN);
filP = {filP.name}';

filenameps = 'E:\Dropbox\UNM\Analysis\IHKA\allPredict.ps';

O =[];
ts =[];
TMP =[];

for j = 1:length(fils)
    load(fils{j})
    [a,bn] = fileparts(fils{j});
    
    %get dat file
    h= figure;
    ax = tight_subplot(6,1);
    
    
    ok = cell2mat(cellfun(@(a) estimateLabel' ==a,num2cell(1:6),'uni',0)');
    
    
    sz_start = find(diff((trueLabel == 5) >0));
    
    [~,b] = histc(ts1,-6000:0);
    k = gaussian2Dfilter([1 1000],[ 1 10]);
    ok = double(ok);
    
    for i = 1:6
        
        axes(ax(i))
        ok(i,:) = (nanconvn(ok(i,:) , k));
        
        
        hold on
        plot([sz_start' sz_start'],[0 1],'color',[.7 .7 .7])
        plot(ok(i,:),'k')
        
    end
    
    mix = min(size(ok,2),length(ts1));
    ts1 = ts1(1:mix);
    ok = ok(:,1:mix);
    
    
    %    print(h, '-dpsc2',filenameps ,'-append');
    close all
    onset = find(diff(ok(3,:)>.7)>0);
    offset = find(diff(ok(3,:)<.7)>0);
    
    %get seizure time
    sz_ix = find(ts1>-1 & ts1<=0);
    
    % load power
    basename =  filP{find(cellfun(@any,regexp(filP,bn(1:end-8))))};
    ff = [power_dirN '/' basename '/' basename '_3.dat'];
    
    % gd = onset(find((bestmatch(onset,sz_ix) - onset') < 1200 & (bestmatch(onset,sz_ix) - onset') > 600));
    %onset = sz_ix;
    onset = onset( (onset>600) & (onset<(size(ok,2)-600)));
    
    tmp1 = nan(240000,20,length(onset));
    for n = 1:length(onset)
        
        tp = LoadBinary(ff,'nchannels',41,'frequency',2000,'channels',2:2:40,'duration', 1200,'start',onset(n)-600);
        %      tmp(:,:,n) = abs(awt_freqlist(double(tp(1:10:end)),2000,1:200));
        tmp1(:,:,n) = tp(1:10:end,:);
    end
    
    TMP(:,:,j) = nanmean(tmp1,3);
    
    O = [O ok];
    
    ts = [ts ts1(1:size(ok,2))];
    j
end

%ps2pdf('psfile', filenameps ,'pdffile','E:\Dropbox\UNM\Analysis\IHKA\allPredict.pdf','gscommand','C:\Program Files (x86)\gs\gs9.54.0\bin\gswin32.exe')
%%

s1 = 1:10;
s2 = 11:13;
s3 = 14:23;
freqs = logspace(log10(.5),log10(200),20);
figure
imagesc(-600:(1/200):600,[],nanmean(TMP(:,:,s2),3)'/1000,[-.25 .25])
%xlim([-50 100])
set(gca,'ydir','normal','ytick',1:20','yticklabel',round(10*freqs)/10)
xlabel('time to 10-600s warning onset(s)')
ylabel('frequency')
colormap('jet')
%%
figure
imagesc(ok)


%%

figure
%    idx = find((diff(O(3,:)>.5)>0));
idx2 = find((diff(O(3,:)<.5)>0));

plot(0:60:3600,histc(idx2 - idx,0:60:3600)/length(idx),'k')
xlim([0 1000])

%semilogx(0:60:3600,histc(idx2 - idx(1:end-1),0:60:3600)/length(idx),'r')



%plot(0:60:3600,histc(idx2 - idx,0:60:3600)/length(idx),'k')

sz_ix = find(ts>-1 & ts<=0);
xx=[];
for i = 1:length(idx)
    xx = [xx;sz_ix(find( sz_ix>idx(i),1,'first'))-idx(i)];
    
    
    
end
hold on
plot(0:300:3600,histc(xx,0:300:3600)/sum(xx<3600))

%%

dt = 1;
bins = -4*3600:dt:0;
NN = length(bins);
[~,b] = histc(ts,bins);
figure
semilogx(bins,accumarray(b(b~=0)',O(1,b~=0)',[NN 1],@nanmean,nan),'b')
hold on
semilogx(bins,accumarray(b(b~=0)',O(3,b~=0)',[NN 1],@nanmean,nan),'g')

semilogx(bins,accumarray(b(b~=0)',O(4,b~=0)',[NN 1],@nanmean,nan),'r')
semilogx(bins,accumarray(b(b~=0)',O(5,b~=0)',[NN 1],@nanmean,nan),'k')

legend('>1hr','10-600s','0-10s','seizure')
%%

