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
kp = arrayfun(@(a) ~isempty(a.ied),ses);
TD = cell2mat({ses(kp).TD}');
ied = cell2mat({ses(kp).ied}');
stim  = cell2mat({ses.stim}');
stim = stim(kp,:);
search = search(kp,:);
obs = obs(kp,:);
%%

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

kp = obs(:,3)~=1000;
[~,b] = sort(nanmean(TD(:,301:end),2)-nanmean(TD(:,1:290),2));
[~,b] = sort(diff(on(kp,:),[],2));
TD1 = TD(kp,:);
ied3 = ied2(kp,:);
fils1 = fils(kp);

figure
imagesc(-300:300,[],TD1(b,:),[0 1.5])
xlabel('Time from stim off (s)')
ylabel('Trial #')
figure
imagesc(-300:300,[],ied3(b,:),[0 1e-3])
xlabel('Time from stim off (s)')
ylabel('Trial #')
%%
close all
[~,b] = sort(nanmean(TD(:,301:600),2));

ok1=  ied2(b,:);
figure
 plotMeanSEM(-300:300,(ok1(end-50:end,:)),'r')
hold on
plotMeanSEM(-300:300,(ok1(1:50,:)),'k')
xlabel('Time from stim off (s)')
ylabel('IED rate')
