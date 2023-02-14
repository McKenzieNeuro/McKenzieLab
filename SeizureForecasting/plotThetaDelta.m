
topDir  ='R:\DGregg\NeuralData\LTD 10.0';
fils = getAllExtFiles(topDir,'mp4',1);
kp = (contains(fils,'s1') | contains(fils,'100uA')) & contains(fils,'LTD(L)') & contains(fils,'uA');

dirs = fileparts(fils(kp));
%theta_delta = nan(length(dirs),3);
%%
for i = 5:length(dirs)
    
    cd(topDir)
    fils = dir(dirs{i});
    
    fils = {fils(cell2mat({fils.isdir})).name}';
    fils = fils(~contains(fils,'.'));
    baselineF = contains(fils,'initial');
    stimF = contains(fils,'100uA') |  contains(fils,'s1');
    postF = contains(fils,'final');
    for j = 1:3
        cd(dirs{i})
        switch j
            
            case 1
                fil = fils{baselineF};
            case 2
                fil = fils{stimF};
                
                if ~any(postF)
                    
                    
                    noPost = true;
                else
                    noPost = false;
                end
                
                
            case 3
                
                if any(postF)
                    
                    
                    fil = fils{postF};
                else
                    
                    fil = [];
                end
                
        end
        
        
        if ~isempty(fil)
            cd(fil)
            
            if ~exist('amplifier.lfp')
                bz_LFPfromDat(pwd,'basename','amplifier')
            end
            
            d = LoadBinary('amplifier.lfp','nchannels',8','channels',6,'frequency',1250);
            
            
            theta = BandpassFilter(double(d),1250,[5 12]);
            power_theta = InstAmplitude(theta);
            
            delta = BandpassFilter(double(d),1250,[1 4]);
            power_delta = InstAmplitude(delta);
            
            
            k = gaussian2Dfilter([10000 1],500);
            power_theta = nanconvn(power_theta(1:10:end),k);
            power_delta = nanconvn(power_delta(1:10:end),k);
            
            
            
            
            if j ==1 || j ==3
                
                theta_delta(i,j) = nanmean(power_theta./power_delta);
                theta_delta_acg{i,j} = power_theta./power_delta;
            elseif j ==2
                
                if ~noPost
                    theta_delta(i,2) = nanmean(power_theta./power_delta);
                    theta_delta_acg{i,2} = power_theta./power_delta;
                    
                else
                    ev = LoadEvents('amplifier.evt.sti');
                    offTime= ev.time(contains(ev.description,'off'));
                    status = InIntervals(((1:10:length(theta))/1250),[0 offTime]);
                    theta_delta(i,2) = nanmean(power_theta(status)./power_delta(status));
                    theta_delta(i,3) = nanmean(power_theta(~status)./power_delta(~status));
                    
                    
                    theta_delta_acg{i,2} = power_theta(status)./power_delta(status);
                    theta_delta_acg{i,3} = power_theta(~status)./power_delta(~status);
                    
                end
                
                
                
            end
            
            
        end
    end
    i
end

%%

theta_delta_acg(:,1) = cellfun(@(a) nanPad(a',400000,1),theta_delta_acg(:,1),'uni',0);
theta_delta_acg(:,2) = cellfun(@(a) nanPad(a',400000),theta_delta_acg(:,2),'uni',0);
theta_delta_acg(:,3) = cellfun(@(a) nanPad(a',400000),theta_delta_acg(:,3),'uni',0);

%%
close all

figure
hold on
plotMeanSEM((1:400000)/125,(cell2mat(cellfun(@(a) a(end-400000:end),theta_delta_acg(:,1),'UniformOutput',false))),'r')
plotMeanSEM((400001:400000+400000)/125,(cell2mat(cellfun(@(a) a(1:400000),theta_delta_acg(:,2),'UniformOutput',false))),'k')


plotMeanSEM((400000+400001:400000+400000+400000)/125,(cell2mat(cellfun(@(a) a(1:400000),theta_delta_acg(:,3),'UniformOutput',false))),'r')

%%
nSTD = 5;
k  = gaussian2Dfilter([1 1000],10);
ch = 2;
for i = 1:length(dirs)
    
    cd(topDir)
    fils = dir(dirs{i});
    
    fils = {fils(cell2mat({fils.isdir})).name}';
    fils = fils(~contains(fils,'.'));
    baselineF = contains(fils,'initial');
    stimF = contains(fils,'100uA') |  contains(fils,'s1');
    postF = contains(fils,'final');
    for j = 1:3
        cd(dirs{i})
        switch j
            
            case 1
                fil = fils{baselineF};
            case 2
                fil = fils{stimF};
                
                if ~any(postF)
                    
                    
                    noPost = true;
                else
                    noPost = false;
                end
                
                
            case 3
                
                if any(postF)
                    
                    
                    fil = fils{postF};
                else
                    
                    fil = [];
                end
                
        end
        
        
        if ~isempty(fil)
            cd(fil)
            
            if ~exist('amplifier.lfp')
                bz_LFPfromDat(pwd,'basename','amplifier')
            end
            
            tmp = LoadBinary('amplifier.lfp','nchannels',8','channels',ch,'frequency',1250);
            
            
            
            
            
            ts = (1:length(tmp))/1250;
            st = std(abs(double(tmp)));
            high_tres = ts(abs(tmp)>nSTD*st);
            kp = diff(ts(abs(tmp)>nSTD*st))>.01;
            IED = high_tres(kp);
            
            
            %get width
            width = nan(length(IED),1);
            for k1 = 1:length(IED)
                tmp = LoadBinary('amplifier.lfp','nchannels',8','channels',ch,'frequency',1250,'start',IED(k1)-.25,'duration',.5);
                tmp = nanconvn(tmp,k');
                ii = round(length(tmp)/2);
                %determine if its a peak or valley
                if tmp(ii)>0
                    [height,loc,w] = findpeaks(double(tmp));
                    ix = bestmatch(ii,loc);
                    if any(loc)
                        width(k1) =  w(loc==ix)/fs;
                    end
                else
                    [height,loc,w] = findpeaks(-double(tmp));
                    ix = bestmatch(ii,loc);
                    if any(loc)
                        width(k1) =  w(loc==ix)/fs;
                    end
                end
                
            end
            
            IED = IED(width<.075);
            
            kp = diff(IED)>.2;
            IED = IED(kp);
            
            events.description = repmat({'IED'},length(IED),1);
            events.time = IED;
            SaveEvents('autoDetect.evt.IED',events)
            
            
            
            
            j
            
        end
        
        
    end
    i
end
%%


nSTD = 5;
k  = gaussian2Dfilter([1 1000],10);
ch = 2;
IED = nan(length(dirs),3);
for i = 1:length(dirs)
    
    cd(topDir)
    fils = dir(dirs{i});
    
    fils = {fils(cell2mat({fils.isdir})).name}';
    fils = fils(~contains(fils,'.'));
    baselineF = contains(fils,'initial');
    stimF = contains(fils,'100uA') |  contains(fils,'s1');
    postF = contains(fils,'final');
    for j = 1:3
        cd(dirs{i})
        switch j
            
            case 1
                fil = fils{baselineF};
            case 2
                fil = fils{stimF};
                
                if ~any(postF)
                    
                    
                    noPost = true;
                else
                    noPost = false;
                end
                
                
            case 3
                
                if any(postF)
                    
                    
                    fil = fils{postF};
                else
                    
                    fil = [];
                end
                
        end
        
        if ~isempty(fil)
            cd(fil)
            
            ev = LoadEvents('autoDetect.evt.IED');
            ok = dir('amplifier.lfp');
            len = ok.bytes/1250/2/8;
            if j ==1 || j ==3
                
                IEDRate = length(ev.time)/len;
                IED(i,j) = IEDRate;
            elseif j ==2
                
                if ~noPost
                    IEDRate = length(ev.time)/len;
                    IED(i,j) = IEDRate;
                    
                else
                    ev1 = LoadEvents('amplifier.evt.sti');
                    offTime= ev1.time(contains(ev1.description,'off'));
                    status = InIntervals(ev.time,[0 offTime]);
                    
                    IEDRate = length(ev.time(status))/offTime;
                    IED(i,2) = IEDRate;
                    
                    
                    len = len - offTime;
                    IEDRate = length(ev.time(~status))/len;
                    IED(i,3) = IEDRate;
                end
                
                
            end
        end
        
    end
    
    
end


%%

clear IED_theta
for i = 1:length(dirs)
    
    cd(topDir)
    fils = dir(dirs{i});
    
    fils = {fils(cell2mat({fils.isdir})).name}';
    fils = fils(~contains(fils,'.'));
    baselineF = contains(fils,'initial');
    stimF = contains(fils,'100uA') |  contains(fils,'s1');
    postF = contains(fils,'final');
    for j = 1:3
        cd(dirs{i})
        switch j
            
            case 1
                fil = fils{baselineF};
            case 2
                fil = fils{stimF};
                
                if ~any(postF)
                    
                    
                    noPost = true;
                else
                    noPost = false;
                end
                
                
            case 3
                
                if any(postF)
                    
                    
                    fil = fils{postF};
                else
                    
                    fil = [];
                end
                
        end
        
        
        if ~isempty(fil)
            cd(fil)
            
            if ~exist('amplifier.lfp')
                bz_LFPfromDat(pwd,'basename','amplifier')
            end
            
            d = LoadBinary('amplifier.lfp','nchannels',8','channels',6,'frequency',1250);
            
            
            theta = BandpassFilter(double(d),1250,[5 12]);
            power_theta = InstAmplitude(theta);
            ts = (1:length(theta))/1250;
            delta = BandpassFilter(double(d),1250,[1 4]);
            power_delta = InstAmplitude(delta);
            
            
            k = gaussian2Dfilter([100000 1],5000);
            power_theta = nanconvn(power_theta,k);
            power_delta = nanconvn(power_delta,k);
            
            
            
            
            
            ev = LoadEvents('autoDetect.evt.IED');
            
            kp = contains(ev.description,'IED');
            
            ts_IED = ev.time(kp);
            
            bins = prctile(power_theta./power_delta,0:10:100);
            bins = 0.5:.1:2;
            % bin power
            i
            power1 = power_theta./power_delta;
            
            if j ==1 || j ==3
                [n,b] = histc(power1,bins);
                IED = histc(ts_IED,ts);
                ok = accumarray(b(b>0),IED(b>0),[length(bins) 1],@sum,nan)./n * 1250;
                IED_theta{j}(i,:) =ok;
            elseif j ==2
                
                if ~noPost
                    [n,b] = histc(power1,bins);
                    IED = histc(ts_IED,ts);
                    ok = accumarray(b(b>0),IED(b>0),[length(bins) 1],@sum,nan)./n * 1250;
                    IED_theta{j}(i,:) =ok;
                    
                else
                    ev1 = LoadEvents('amplifier.evt.sti');
                    offTime= ev1.time(contains(ev1.description,'off'));
                    status = InIntervals(((1:length(theta))/1250),[0 offTime]);
                    
                    
                    [n,b] = histc(power1(status),bins);
                    IED = histc(ts_IED(ts_IED<offTime),ts(status));
                    ok = accumarray(b(b>0),IED(b>0),[length(bins) 1],@sum,nan)./n * 1250;
                    IED_theta{2}(i,:) =ok;
                    
                    [n,b] = histc(power1(~status),bins);
                    IED = histc(ts_IED(ts_IED>offTime),ts(~status));
                    ok = accumarray(b(b>0),IED(b>0),[length(bins) 1],@sum,nan)./n * 1250;
                    IED_theta{3}(i,:) =ok;
                end
                
                
                
            end
            
            
        end
    end
end





%%
figure
boxplot(theta_delta, 'datalim' ,[-inf inf], 'whisker' ,0,'symbol','')
hold on

plot(1:3,theta_delta')
set(gca,'xtick',1:3,'xticklabel',{'Pre','LTD stim','Post'})
ylabel('theta/delta ratio')

%%
figure
boxplot(IED, 'datalim' ,[-inf inf], 'whisker' ,0,'symbol','')
hold on

plot(1:3,IED')
set(gca,'xtick',1:3,'xticklabel',{'Pre','LTD stim','Post'})
ylabel('IED rate (Hz)')