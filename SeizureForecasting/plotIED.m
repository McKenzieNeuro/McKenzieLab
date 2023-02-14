topDir = 'R:\DGregg\NeuralData\LTD 10.0';

fils = getAllExtFiles(topDir,'IED',1);
d= nan(984,40000);
ix=1;
for i = 1:length(fils)
    [a,b] = fileparts(fils{i});
    cd(a)
    
    
    xml  = LoadXml('amplifier.xml');
    ev = LoadEvents(fils{i});
    
    kp = contains(ev.description,'IED');
    
    ts = ev.time(kp);
    
    if xml.nChannels>8
        ch = 33:40;
    else
        ch = 1:8;
    end
    
    for j = 1:length(ts)
        
        tmp = LoadBinary('amplifier.dat','nchannels',xml.nChannels,'channels',ch,'frequency',20000,'start',ts(j)-.25,'duration',2);
        if size(tmp,1)==40000
            [~,ix1] = max(max(abs(tmp(4500:5500,:))));
            tmp = tmp(:,ix1);
            d(ix,:) = tmp;
            ix = ix+1;
        end
    end
    i
end


%%

topDir = 'R:\DGregg\NeuralData\LTD 10.0\10-6-2022(12.41)_LTD(L)_RHS';
fils = dir(topDir);

fils = {fils(cell2mat({fils.isdir})).name}';
fils = fils(~contains(fils,'.'));

ix=1;
cd(topDir)
for i = 1:length(fils)
    cd(topDir)
  
    cd(fils{i})
    xml  = LoadXml('amplifier.xml');
    
    
    if xml.nChannels>8
        sm_SplitDatChannel('amplifier.dat','amplifier1.dat','channelix',1:8);
        
        
        
    end
    i
end

%%

topDir = 'R:\DGregg\NeuralData\LTD 10.0\10-6-2022(12.41)_LTD(L)_RHS';
fils = dir(topDir);

fils = {fils(cell2mat({fils.isdir})).name}';
fils = fils(~contains(fils,'.'));

for i = 1:length(fils)
    cd(topDir)
    cd(fils{i})
    
    bz_LFPfromDat(pwd,'basename','amplifier')
    
    i
end

%%


fils = fils([3 1 2]);
%%
% relate IED to theta power

for i = 1:length(fils)
    cd(topDir)
    cd(fils{i})
    
    d = LoadBinary('amplifier.lfp','nchannels',8','channels',6,'frequency',1250);
    
    
    theta = BandpassFilter(double(d),1250,[5 12]);
    power_theta = InstAmplitude(theta);
    
    delta = BandpassFilter(double(d),1250,[1 4]);
    power_delta = InstAmplitude(delta);
    
    
    k = gaussian2Dfilter([100000 1],5000);
    power_theta = nanconvn(power_theta,k);
    power_delta = nanconvn(power_delta,k);
    
    
    theta_delta{i} = power_theta./power_delta;
    
    if i ==1
        ts{i} = (1:length(theta))/1250;
    else
        ts{i} = ((1:length(theta))/1250) + ts{i-1}(end);
    end
    
    %   ev = LoadEvents(fils{i});
    
    % kp = contains(ev.description,'IED');
    
    % ts_IED = ev.time(kp);
    
    %     bins = prctile(power_theta,0:10:100);
    %      IED = histc(ts_IED,ts);
    %      power1 = power_theta(1250:end);
    %      IED1 = IED(1:end-1249);
    %
    %      [~,b] = histc(power1,bins);
    %      ok = accumarray(b(b>0),IED1(b>0),[],@sum,nan);
    %      % bin power
    i
    
end

%%

figure

plot(cell2mat(ts),cell2mat(theta_delta'))
hold on
plot([ts{1}(end) ts{1}(end)],[0 7],'k')
plot([ts{2}(end) ts{2}(end)],[0 7],'k')
%%

cd('R:\DGregg\NeuralData\LTD 10.0\9-26-2022(12.53)_LTD(L)_RHS\initialBaseline__220926_130000')
d = LoadBinary('amplifier.lfp','nchannels',8,'channels',[3 5 6 7 8],'frequency',1250,'start',5775,'duration',60);

ts = (1:size(d,1))/1250;
%%
figure
for i = 1:5
    hold on
    
plot(ts,d(:,i)*.195 - (i*4000),'k')

plot([5 5],[-2200 -200],'k')
end

%%

% plot example
cd('R:\DGregg\NeuralData\LTD 10.0\10-24-2022(12.56)_LTD(L)_RHS\(s1-100)uA__221024_130000')

d = LoadBinary('amplifier.lfp','nchannels',8,'channels',8,'frequency',1250);
pulse = find(diff(d)>1e4);
kp = diff(pulse)>100;
clear d

d1 = LoadBinary('amplifier.lfp','nchannels',8,'channels',6,'frequency',1250);
theta = BandpassFilter(double(d1),1250,[5 12]);
power_theta = InstAmplitude(theta);
delta = BandpassFilter(double(d1),1250,[1 4]);
power_delta = InstAmplitude(delta);
k = gaussian2Dfilter([100000 1],5000);
power_theta = nanconvn(power_theta,k);
power_delta = nanconvn(power_delta,k);
TD =  power_theta./power_delta;
%%
close all
figure
hold on
plot((1:length(TD(1:1000:end)))/1.25,TD(1:1000:end))
plot([pulse(kp)/1250 pulse(kp)/1250],[7 7.2],'k')
