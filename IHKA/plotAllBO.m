% get all files
load('R:\DGregg\NeuralData\PCP\Recordings\BayesOpt\Box1_gridData_TD.mat')
topDir = 'R:\DGregg\NeuralData\PCP\Recordings\BayesOpt\rec';
fils = getAllExtFiles(topDir,'dat',1);
kp = contains(fils,'amplifier.dat');
fils = fils(kp);
[~,b] = sort(cellfun(@(a) a(end-26:end-14),fils,'uni',0));
fils =  fils(b);

%%
clear ses
for i = 1:length(fils)
    subDir = fileparts(fils{i});
    xmlfil = [subDir filesep 'amplifier.xml'];
    if ~exist(xmlfil)
        copyfile([topDir filesep 'amplifier.xml'],xmlfil)
    end
    
    cd(subDir)
    
    % get stim
    s = LoadBinary('stim.dat','frequency',20000,'nchannels',8,'channels',7);
    ses(i).stim(1) = find(s>5000,1,'first')/20000;
    ses(i).stim(2) = find(s>5000,1,'last')/20000;
    %get ied
    
    if ~exist('autoDetect.evt.IED')
        st = sm_detectIED(pwd,'basename','amplifier','ch',1:8,'pass_band',[10 200]);
    end
    ev = LoadEvents('autoDetect.evt.IED');
    k = gaussian2Dfilter([ 1 10000],1250);
    tmp =  LoadBinary('amplifier.lfp','nchannels',8,'channels',6,'frequency',1250);
    
    theta = BandpassFilter(double(tmp),1250,[4 12]);
    delta = BandpassFilter(double(tmp),1250,[1 4]);
    deltaP = nanconvn(InstAmplitude(delta),k');
    thetaP = nanconvn(InstAmplitude(theta),k');
    
    
    
    TD = thetaP./deltaP;
    ts = (1:length(TD))/1250;
    
    
    ied = histc(ev.time,ts);
    ses(i).TD = avghist(ts,TD,ses(i).stim(2)-300:ses(i).stim(2)+300);
    ses(i).ied = avghist(ts,ied,ses(i).stim(2)-300:ses(i).stim(2)+300);
    ses(i).obs = obs(i,:);
    ses(i).search = search(i,:);
    i
end


%%
% get all files
load('R:\DGregg\NeuralData\PCP\Recordings\BayesOpt2\Box1_gridData_TD.mat')
topDir = 'R:\DGregg\NeuralData\PCP\Recordings\BayesOpt2\rec';
fils = getAllExtFiles(topDir,'dat',1);
kp = contains(fils,'amplifier.dat');
fils = fils(kp);
[~,b] = sort(cellfun(@(a) a(end-26:end-14),fils,'uni',0));
fils =  fils(b);

idx = length(ses)+1;
%%
idx=1
for i = 144:length(fils)
    subDir = fileparts(fils{i});
    xmlfil = [subDir filesep 'amplifier.xml'];
    if ~exist(xmlfil)
        copyfile([topDir filesep 'amplifier.xml'],xmlfil)
    end
    
    cd(subDir)
    sl = regexp(fils{i},'\');
    ses_name = fils{i}(sl(end-2)+1:sl(end-1)-1);
    [~,ix] = ismember(ses_name,folders);
    % get stim
    s = LoadBinary('stim.dat','frequency',20000,'nchannels',8,'channels',7);
    ses(idx).stim(1) = find(s>5000,1,'first')/20000;
    ses(idx).stim(2) = find(s>5000,1,'last')/20000;
    %get ied
    
    if ~exist('autoDetect.evt.IED')
        st = sm_detectIED(pwd,'basename','amplifier','ch',1:8,'pass_band',[10 200]);
    end
    ev = LoadEvents('autoDetect.evt.IED');
    k = gaussian2Dfilter([ 1 10000],1250);
    tmp =  LoadBinary('amplifier.lfp','nchannels',8,'channels',6,'frequency',1250);
    
    theta = BandpassFilter(double(tmp),1250,[4 12]);
    delta = BandpassFilter(double(tmp),1250,[1 4]);
    deltaP = nanconvn(InstAmplitude(delta),k');
    thetaP = nanconvn(InstAmplitude(theta),k');
    
    
    
    TD = thetaP./deltaP;
    ts = (1:length(TD))/1250;
    
    
    ied = histc(ev.time,ts);
    ses(idx).TD = avghist(ts,TD,ses(idx).stim(2)-300:ses(idx).stim(2)+300);
    ses(idx).ied = avghist(ts,ied,ses(idx).stim(2)-300:ses(idx).stim(2)+300);
    if ix~=0
        ses(idx).obs = obs(ix,:);
        ses(idx).search = search(ix,:);
    else
        ses(idx).obs = nan;
        ses(idx).search  = nan;
        
    end
    idx = idx+1;
    i
end

%%
kp1 = arrayfun(@(a) ~any(isnan(a.search)),ses);
kp2 = arrayfun(@(a) ~isempty(a.ied),ses);
kp = kp1&kp2;
TD = cell2mat({ses(kp).TD}');
ied = cell2mat({ses(kp).ied}')*1250;
stim  = cell2mat({ses.stim}');
stim = stim(kp,:);
search =  cell2mat({ses(kp).search}');
TD_diff = nanmean(TD(:,301:end),2)-nanmean(TD(:,1:290),2);


%search = search(kp,:);
%obs = obs(kp,:);
%%
clear on ied2
for j =1:size(TD,1)
    ix = find(TD(j,301:end)>1,1,'first');
    
    if any(ix)
        on(j,1) = ix+300;
        ix_off  = find(TD(j,ix+300:end)<1,1,'first');
        if any(ix_off)
            on(j,2) = ix_off+300+ix;
        else
            on(j,2) = 601;
        end
    else
        on(j,:) = nan(1,2);
    end
end


for i = 1:size(ied,1)
    k = gaussian2Dfilter([ 1 10000],3);
    ied2(i,:) = nanconvn(ied(i,:),k);
end

%%
close all
kp = obs(:,3)~=1000;
kp = ~isnan(diff(on,[],2));
[~,b] = sort(nanmean(TD(:,301:end),2)-nanmean(TD(:,1:290),2));
[~,b] = sort(nanmean(TD(:,301:600),2));

[a,b] = sort(diff(on(kp,:),[],2));
TD1 = TD(kp,:);
ied3 = ied2(kp,:);
%fils1 = fils(kp);

figure
imagesc(-300:300,[],TD1(b,:),[0 1.5])
xlabel('Time from stim off (s)')
ylabel('Trial #')
xlim([-100 100])
figure
imagesc(-300:300,[],ied3(b,:),[0 1e-3])
xlabel('Time from stim off (s)')
ylabel('Trial #')
xlim([-100 100])

%close all
%[~,b] = sort(nanmean(TD(:,301:600),2));


ok1=  ied3(b,:);
figure
plotMeanSEM(-300:300,(ok1(end-200:end,:)),'r')
hold on
plotMeanSEM(-300:300,(ok1(1:200,:)),'k')
xlabel('Time from stim off (s)')
ylabel('IED rate')
xlim([-100 100])

%%



%%

gd_idx = 379:size(TD,1);
gd_idx = 1:size(TD,1);
TD_diff = nanmean(TD(gd_idx,301:501),2)-nanmean(TD(gd_idx,1:290),2);
%TD_diff = nanmean(TD(:,301:end),2);
%TD_diff = nanmean(ied(gd_idx,301:end),2)-nanmean(ied(gd_idx,1:290),2);
close all
k = gaussian2Dfilter([ 100 100],2);
figure
[~,~,~,b] = histcn([search(gd_idx,2),search(gd_idx,3)],0:50:200,0:50:300);
imagesc(0:50:300,0:50:200,nanconvn(accumarray(b,TD_diff,[],@nanmean,nan),k,'nanout',true))
xlabel('duration')
ylabel('frequency')

%%
TD_diff = nanmean(TD(:,301:501),2)-nanmean(TD(:,1:290),2);
x_observed = search(gd_idx,2:3);
y_observed1 = -TD_diff(gd_idx);
model = fitrgp(x_observed,-y_observed1);
 [xx,yy] = meshgrid(.1:10:200,1:10:300);
pr = predict(model,[xx(:) yy(:)]);
figure
pr = reshape(pr,size(xx,1),size(xx,2));
imagesc(1:10:300,.1:10:200,pr')

xlabel('duration')
ylabel('frequency')
%%
x1 = optimizableVariable('x1',[.1,200]);

x2 = optimizableVariable('x2',[1,300]);
fun = @(x) myPred1(x_observed,y_observed1,table2array([x]));

results = bayesopt(fun,[x1,x2],'Verbose',0,...
    'AcquisitionFunctionName','expected-improvement-plus','ExplorationRatio',.75,'IsObjectiveDeterministic',false,...
    'UseParallel',true,'MaxObjectiveEvaluations',100);





%%

figure
 plot(-300:300,ses(562).TD)
hold on
plot(-300:300,ses(562).ied*1250)
xlim([-100 100])

%%
ev = LoadEvents('R:\DGregg\NeuralData\PCP\Recordings\BayesOpt2\rec\Stim0185_9-17-2024(9.16)\RHS_240917_091711\autoDetect.evt.IED');
d = LoadBinary('R:\DGregg\NeuralData\PCP\Recordings\BayesOpt2\rec\Stim0185_9-17-2024(9.16)\RHS_240917_091711\amplifier.dat',...
    'frequency',20000,'nchannels',8,'channels',1:8,'start',304,'duration',70);
%%
figure
ts = ((1:length(d(:,5)))/20000) - (ses(562).stim(2)-304);
plot(ts,d(:,5)*.195,'k','linewidth',2)
ylim([-2000 1500])
set(gca,'fontsize',16)

hold on
plot([5 5],[-1500 -1000],'k','linewidth',6)
axis off
plot(ev.time-ses(562).stim(2),1000,'*','color','r')
xlim([-14 60])