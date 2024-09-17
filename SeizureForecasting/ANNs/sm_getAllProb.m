ops.numfreq  = [4.9334 7.7869 12.2910 19.4002 30.6214 48.3330 76.2892 120.4155 190.0649 300.0000]; % frequencies for wavelet decomp.
ops.nchanFil= 8; % number of channel in dat file
ops.nChanSubj = 8; % # channels for subj
ops.channels = 1:8; % channels for subject
ops.Fs =  20000; % sampling rate
ops.nTempBin =  32; % ?
ops.data =  15; % ?
ops.Twin =  2; % window size
ops.feature_fun = @sm_GetDataFeature2; % function for featurizing the data




%for retrospective
fil = 'R:\DGregg\NeuralData\EDS\ClosedLoop\EDS2.3\9-17-2023(15.0)\RHS_230917_150102\amplifier.dat';
tmp = dir(fil);

filedur = tmp.bytes/2/ops.nchanFil/ops.Fs;
%%

stimfil = 'G:\RecordingData\EDS\closedLoop\EDS2.3\9-17-2023(15.0)\RHS_230917_150102\stim.dat';
data = LoadBinary(stimfil,'nchannels',8,'frequency',ops.Fs,'channels',7);
pulse = find(diff(data>4000)>0);
ts = (1:length(data))/20000;
stim_on = ts(pulse);
kp = diff([0 stim_on])>1;
stim_on = stim_on(kp);
%%
%LOAD MODEL
model_fname  = 'R:\DGregg\SeizureForecast\Seizuredetect_demo\Training\netCNN2.mat';
model_fname  = 'R:\DGregg\SeizureForecast\Seizuredetect_demo\Training\randomForest.mat';
v = load(model_fname);
%%


%get sz prob in 1s increments
clear tim

Pr = nan(ceil(filedur-ops.Twin/.1),1);
ix = 1;
Fstest = 1;
for i = 1:(1/Fstest):filedur-ops.Twin
    tic
data = LoadBinary(fil,'nchannels',ops.nchanFil,'frequency',ops.Fs,'channels',ops.channels,'duration',2 ,'start',i);
Prob = sm_getSeizProb(data,v.rusTree,ops);
Pr(ix) = Prob(1);
tim(ix) = toc;
ix = ix+1

end

%%

% save auto-detected seizures
sz_ep = [find(diff([0;Pr]==2)>0) find(diff([0;Pr]==2)<0)]/Fstest;
ev.time  = [sz_ep(:,1);sz_ep(:,2)];
ev.description = [repmat({'sz_on'},size(sz_ep,1),1);repmat({'sz_off'},size(sz_ep,1),1)];

[~,b] = sort(ev.time(:,1));
ev.time = ev.time(b,:);
ev.description = ev.description(b);
filename = 'R:\DGregg\NeuralData\EDS\ClosedLoop\EDS2.3\9-17-2023(15.0)\RHS_230917_150102\autodetect.evt.szr';
SaveEvents(filename,ev);

%%

% plot  real vs dummy CCG
sz_ep1 = sz_ep(diff(sz_ep,[],2)>3,:);

ts  = (1:length(Pr))/Fstest;

stimblock = mod(floor(ts/3600),2)==0;
stimblock(ts>ev.time(end)) = false;

kp = stimblock(:)==0 & Pr(:)==2;
szt = find(kp);


% get fake sz
stim_dummy = szt(1);
for i = 2:length(szt)
    if szt(i)>stim_dummy(end)+15
        stim_dummy = [stim_dummy;szt(i)];
    end
    
end
k  = gaussian2Dfilter([10000 1],.1);
kp = diff(ev.time)>1;
lag  = -500:500;
cc_real  = CrossCorr(ev.time(kp),sz_ep1(:,1),2,1001)/length(ev.time(kp));
cc_dummy  = CrossCorr(stim_dummy,sz_ep1(:,1),2,1001)/length(stim_dummy);
%cc_dummy(501) = nan;
close all
figure
 hold on
plot(lag,nanconvn(cc_real,k),'k')
plot(lag,nanconvn(cc_dummy,k,'nanout',true),'r')
%%
warning off
close all
figure
% check events
ev = LoadEvents(filename);
fil = 'R:\DGregg\NeuralData\EDS\Prophylactic3\4-13-2023(12.59)\RHS_230413_130000\amplifier.dat';
ts = (1:size(data,1))/20000 - 5;
d = [];
for i = 1:length(ev.time)






    data = LoadBinary(fil,'nchannels',ops.nchanFil,'frequency',ops.Fs,'channels',ops.channels,'duration',10 ,'start',ev.time(i)-5);

    for j = 1:size(data,2)
        plot(ts,double(data(:,j))-j*6000,'k')
        hold on
    end
    waitforbuttonpress
    close all

%d(:,:,i) = abs(awt_freqlist(double(data(:,1)),20000,logspace(log(1),log10(300),50))');
i
end