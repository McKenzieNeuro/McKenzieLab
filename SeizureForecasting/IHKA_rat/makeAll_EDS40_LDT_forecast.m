LDT_stim = [...
    
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction1\11-3-2023(14.26)\RHS_231103_142651'} ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction2\11-8-2023(13.41)\RHS_231108_134150'} ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction3\11-10-2023(13.7)\RHS_231110_130745'} ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction3\11-12-2023(12.58)\RHS_231112_130000'} ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\Induction1\11-1-2023(13.22)\RHS_231101_132251'} ;...

];

%%

vv=load('R:\Analysis\SeizureForecasting\IHKA_rat_RF\features.mat');
stim_ses = find(ismember(fileparts(vv.sessions(:,1)),LDT_stim));
%%
ClassifierFileOutputDir = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\Classification\';

k  = gaussian2Dfilter([10000 1],10);
for i = 1:length(stim_ses)
    
    outfil = [ClassifierFileOutputDir 'predict_' num2str(stim_ses(i)) '.mat'];
    v=load(outfil);
    % find stim times
    
    lfpfil = vv.sessions{stim_ses(i),2};
    dat  = LoadBinary(lfpfil,'frequency',1250,'nchannels',32,'channels',7);
    stims = find(diff(dat>1.5e4)>0)/1250;
    estInStim = InIntervals(1:length(v.estimateLabel),[stims(1) stims(end)]);
    pr_far(i,1) = nanmean(v.estimateLabel(estInStim,1)==1);
    pr_far(i,2) = nanmean(v.estimateLabel(~estInStim,1)==1);
    i
end
%%
%get baseline

BL_ses= find(contains(sessions(:,2),'BL'));


for i = 1:length(BL_ses)
    
    outfil = [ClassifierFileOutputDir 'predict_' num2str(BL_ses(i)) '.mat'];
    v=load(outfil);
    BL_far(i) = nanmean(v.estimateLabel(:,1)==1);
    
end

%%
%example
close all
for i = 5
    outfil = [ClassifierFileOutputDir 'predict_' num2str(stim_ses(i)) '.mat'];
    v=load(outfil);
    % find stim times
    
    lfpfil = vv.sessions{stim_ses(i),2};
    
     
    dat  = LoadBinary(lfpfil,'frequency',1250,'nchannels',32,'channels',7);
    k  = gaussian2Dfilter([10000 1],10);
    stims = find(diff(dat>1.5e4)>0)/1250;
    figure
    plot(nanconvn(v.estimateLabel==1,k),'k')
    hold on
    plot([stims ],[ 1.1],'.','color','r')
    %plot([v.seizure_start,1.5,'x','b')
end

%%

%CA1 stim
CA1_stim = [...
    
{'R:\DGregg\NeuralData\EDS_Cohort2\LFHS1\11-14-2023(14.46)\RHS_231114_144651'} ;...
{'R:\DGregg\NeuralData\EDS_Cohort2\LFHS1\11-20-2023(12.58)\RHS_231120_130000'};...
{'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-22-2023(12.58)\RHS_231122_130000'};...
{'R:\DGregg\NeuralData\EDS_Cohort2\LFHSP1\11-24-2023(12.58)\RHS_231124_130000'};...

];

stim_ses_CA1 = find(ismember(fileparts(vv.sessions(:,1)),CA1_stim));
%%

close all
ix=1;
for i = 2:4
    outfil = [ClassifierFileOutputDir 'predict_' num2str(stim_ses_CA1(i)) '.mat'];
    v=load(outfil);
    % find stim times
    
    lfpfil = vv.sessions{stim_ses_CA1(i),2};
    
     
    dat  = LoadBinary(lfpfil,'frequency',1250,'nchannels',32,'channels',5);
    k  = gaussian2Dfilter([10000 1],10);
    stims = find(diff(dat>1e4)>0)/1250;
    ok(ix,:) = nanconvn(v.estimateLabel(1:1e4,1)==2,k);
   % figure
   % plot(nanconvn(v.estimateLabel(:,1)==1,k),'k')
   % hold on
   %plot([stims ],[ 1.1],'.','color','r')
    %plot([v.seizure_start,1.5,'x','b')
    ix=ix+1;
end

%%
plotMeanSEM(1:1e4,ok,'b')
hold on
plot([stims ],[ 1.1],'.','color','r')


%%

% get example seizure

cd('R:\DGregg\NeuralData\EDS_Cohort2\LFHS1\11-17-2023(23.35)_BL\RHS_231117_233532')

d = LoadBinary('amplifier.lfp','frequency',1250,'start',1117,'duration',40,'channels',1:8,'nchannels',32);
%%
close all
figure
for i = 1:7
    
    plot(ts,double(d(:,i))*.195 - i*2000,'k')
    hold on
end