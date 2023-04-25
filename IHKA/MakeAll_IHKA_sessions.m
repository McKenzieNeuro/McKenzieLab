%handle single session

chReg = [...
    {'PrL(L)' };...
    {'PrL(R)' };...
    {'AVT(L)' };...
    {'BLA(R)' };...
    {'CA1(L)' };...
    {'CA1(R)' };...
    {'LDT(L1)'};...
    {'LDT(L2)'}];

topdir = 'R:\DGregg\NeuralData\EDS\';
files = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files,'recInfo') & ~contains(files,'WASH');
files = files(keepFiles);
recInfoSubDir = fileparts(files);


lfpfiles = getAllExtFiles(topdir,'lfp',1);
keepFiles = ~contains(lfpfiles,'WASH');
lfpfiles = lfpfiles(keepFiles);

lfpsubdir = unique(fileparts(fileparts(lfpfiles)));


kp = ismember(lfpsubdir,recInfoSubDir);



%%
idx = 1;
for kk = 1:length(lfpsubdir)
    
    cd(lfpsubdir{kk});
    clear recInfo
    load('recInfo.mat')
    
    ISI = load('config.mat','ISI');
    ISI = mean(ISI.ISI);
    
    
    
    
    %loop through systems (RHD/RHS)
    for jj = 1:2
        
        switch jj
            case 1
                syst = find(contains(recInfo(:,2),'RHD'));
                kp = contains(lfpfiles,lfpsubdir{kk}) & contains(lfpfiles,'RHD');
                
            case 2
                syst = find(contains(recInfo(:,2),'RHS'));
                kp = contains(lfpfiles,lfpsubdir{kk}) & contains(lfpfiles,'RHS');
        end
        
        
        if sum(kp) ==1
            
            lfp = lfpfiles{kp};
            
            
            
            xmlfil = [lfp(1:end-3) 'xml'];
            xml = LoadXml(xmlfil);
            %loop through animals on each system (2per)
            totCh =0;
            
            % check if recInfo is made
            
            if size(recInfo,2)>3
                
                recInfo = sm_recfileAnalysis(lfpsubdir{kk});
                
            end
            
            [outDir] = fileparts(lfp);
            for ii = 1:2
                clear sessiondata
                
                subjectName = recInfo{syst,3}{1+ii,2}{1};
                
                
                outfil = fullfile(outDir,[subjectName '.sessiondata.mat']);
                
                
                if ~exist(outfil)
                    if isempty(recInfo{syst,3}{2,7})
                        
                        %make sure the session has 8 channels
                        if xml.nChannels ~=16
                            error('non default channel list')
                        end
                        
                        if jj==2
                            
                            error('need to make recInfofile')
                        end
                        
                        stimAmp = nan;
                        stimON = nan;
                        stimType = 'none';
                        % hard code 8 channel with standard labels
                        chNum = (0:7) + totCh;
                        chName = chReg;
                    else
                        chNum = totCh + cellfun(@(a,b) str2num(a(b)) ,recInfo{syst,3}{1+ii,end}(:,1),regexp(recInfo{syst,3}{1+ii,end}(:,1),'[0-9]'));
                        chName = recInfo{syst,3}{1+ii,end}(:,2);
                        
                        %get the stim information
                        if jj==2
                            stimChannels = recInfo{syst,3}{1+ii,4}; % load the active stimulation channels for comparison
                            if stimChannels(2) == 0
                                stimType = 'mono';
                            elseif stimChannels(2) == 8
                                stimType = 'bi';
                            else
                                stimType = 'none';
                            end
                            stimAmp = recInfo{syst,3}{1+ii,3}; % load the stimulation amplitude
                            stimON = recInfo{syst,3}{1+ii, 6}/1000;
                            
                        else
                            stimAmp = nan;
                            stimON = nan;
                            stimType = 'none';
                        end
                        % save the data
                        
                        
                        
                    end
                    
                    sessiondata.subject = subjectName;
                    sessiondata.StimType = stimType;
                    sessiondata.StimAmp = stimAmp;
                    sessiondata.stimON = stimON;% stim data encoded in ms. convert to seconds
                    
                    
                    
                    
                    sessiondata.channel = chName;
                    sessiondata.channelID = chNum+1;
                    sessiondata.lfpFile = lfp;
                    
                    save(outfil,'sessiondata')
                    totCh = length(chNum);
                    idx = idx+1;
                    
                end
            end
        end
    end
    
end
%%




files = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files,'sessiondata');

files = files(keepFiles);

for i = 1:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    
    if ~isfield(sessiondata,'theta_delta') || (exist('autoDetect.evt.IED') && ~isfield(sessiondata,'IED'))
        % get TD
        xmlfil = [sessiondata.lfpFile(1:end-3) 'xml'];
        xml = LoadXml(xmlfil);
        fs = xml.lfpSampleRate;
        k = gaussian2Dfilter([10000 1],1250/2);
        ch = sessiondata.channelID(contains(sessiondata.channel,'BLA'));
        
        if isempty(ch)
            
            ch = sessiondata.channelID(contains(sessiondata.channel,'CA1(R)'));
        end
        d = LoadBinary(sessiondata.lfpFile,'nchannels',xml.nChannels,'channels',ch,'frequency',fs);
        
        theta = BandpassFilter(double(d),fs,[5 12]);
        power_theta = InstAmplitude(theta);
        
        delta = BandpassFilter(double(d),fs,[1 4]);
        power_delta = InstAmplitude(delta);
        
        
        
        power_theta = nanconvn(power_theta(1:10:end),k);
        power_delta = nanconvn(power_delta(1:10:end),k);
        
        sessiondata.theta_delta  = power_theta./power_delta;
        sessiondata.ts = (1:length(power_delta))/(fs/10);
        % get IED
        
        if exist('autoDetect.evt.IED')
            evs = LoadEvents('autoDetect.evt.IED');
            
            
            for j =1:length(sessiondata.channelID)
                
                
                str = ['IED: Ch' num2str(sessiondata.channelID(j))];
                sessiondata.IED{j} = evs.time(contains(evs.description,str));
            end
        end
        save(files{i},'sessiondata')
    end
    i
end
%%


%handle merged sessions


chReg = [...
    {'PrL(L)' };...
    {'PrL(R)' };...
    {'AVT(L)' };...
    {'BLA(R)' };...
    {'CA1(L)' };...
    {'CA1(R)' };...
    {'gRSC (L)'};...
    {'LDT(L2)'}];

topdir = 'R:\DGregg\NeuralData\LTD 10.0';
fils = getAllExtFiles(topdir,'mp4',1);
kp = (contains(fils,'s1') | contains(fils,'100uA')) & contains(fils,'LTD(L)') & contains(fils,'uA');
files = fils(kp);
recInfoSubDir = fileparts(fils(kp));



lfpfiles = getAllExtFiles(topdir,'lfp',1);
keepFiles = ~contains(lfpfiles,'WASH');
lfpfiles = lfpfiles(keepFiles);

lfpsubdir = unique(fileparts(fileparts(lfpfiles)));


kp = ismember(lfpsubdir,recInfoSubDir);



%%
idx = 1;
for kk = 1:6%length(lfpsubdir)
    
    cd(lfpsubdir{kk});
    
    
    sl = regexp(lfpsubdir{kk},filesep);
    subjectName = lfpsubdir{kk}(sl(3)+1:sl(4)-1);
    
    
    
    kp = contains(lfpfiles,lfpsubdir{kk}) & contains(lfpfiles,'RHS');
    
    if exist('recInfo.mat')
        load('recInfo.mat')
        
    else
        
        recInfo{1,2} = 100;
    end
    clear sessiondata
    
    outfil = fullfile(lfpsubdir{kk},[subjectName '.sessiondata.mat']);
    if sum(kp) ==1
        
        lfp = lfpfiles{kp};
        
        
        
        xmlfil = [lfp(1:end-3) 'xml'];
        xml = LoadXml(xmlfil);
        
        
        
        
        
        
        
        
        
        
        
        
        
        %  if ~exist(outfil)
        chNum = 1:xml.nChannels;
        chName = chReg;
        stimType = 'mono';
        
        if size(recInfo,1) ==1
            stimAmp = recInfo{1,2};
        else
            stimAmp = recInfo{2,2};
        end
        pulse_dir = fileparts(lfp);
        if exist(fullfile(pulse_dir,'TTL_pulse.mat'))
            load(fullfile(pulse_dir,'TTL_pulse.mat'))
        else
            [ups,dwns]  = sm_getDigitalin(pulse_dir,'digitalin',xml.SampleRate);
        end
        stimON = ups{2};
        
        
        sessiondata.subject = subjectName;
        sessiondata.StimType = stimType;
        sessiondata.StimAmp = stimAmp;
        sessiondata.stimON = stimON ;% stim data from start of session
        
        
        
        
        sessiondata.channel = chName;
        sessiondata.channelID = chNum;
        sessiondata.lfpFile = lfp;
        
        
        
        save(outfil,'sessiondata')
        
        %  end
        
        
    elseif sum(kp)==3
        
        gd_fils = lfpfiles(kp);
        xmlfil = [gd_fils{1}(1:end-3) 'xml'];
        xml = LoadXml(xmlfil);
        f1 = gd_fils{contains(gd_fils,'initial')};
        f2 = gd_fils{contains(gd_fils,'uA')};
        f3 = gd_fils{contains(gd_fils,'final')};
        
        
        
        
        
        
        
        % if ~exist(outfil)
        
        chNum = 1:xml.nChannels;
        chName = chReg;
        stimType = 'mono';
        if size(recInfo,1) ==1
            stimAmp = recInfo{1,2};
        else
            stimAmp = recInfo{2,2};
        end
        
        % get stim times
        
        
        pulse_dir = fileparts(f2);
        % Digital in
        if exist(fullfile(pulse_dir,'TTL_pulse.mat'))
            load(fullfile(pulse_dir,'TTL_pulse.mat'))
        else
            [ups,dwns]  = sm_getDigitalin(fileparts(f2),'digitalin',xml.SampleRate);
        end
        
        stimON = ups{1};
        
        
        
        
        tmp = dir(f1);
        dur(1) = tmp.bytes/xml.nChannels/xml.lfpSampleRate/2;
        
        tmp = dir(f2);
        dur(2) = tmp.bytes/xml.nChannels/xml.lfpSampleRate/2;
        
        tmp = dir(f3);
        dur(3) = tmp.bytes/xml.nChannels/xml.lfpSampleRate/2;
        
        
        sessiondata.subject = subjectName;
        sessiondata.StimType = stimType;
        sessiondata.StimAmp = stimAmp;
        sessiondata.stimON = stimON + dur(1);% stim data from start of session
        
        
        
        
        sessiondata.channel = chName;
        sessiondata.channelID = chNum;
        sessiondata.lfpFile{1} = f1;
        sessiondata.lfpFile{2} = f2;
        sessiondata.lfpFile{3} = f3;
        
        sessiondata.fileDur = dur;
        
        save(outfil,'sessiondata')
        
        
        
        
        
        %  end
    else
        
        disp('here')
        
    end
end


%%




files = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files,'sessiondata');

files = files(keepFiles);

for i = 1:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    if ~isfield(sessiondata,'theta_delta')
    % get TD
    
    if iscell(sessiondata.lfpFile)
        
        xmlfil = [sessiondata.lfpFile{1}(1:end-3) 'xml'];
        xml = LoadXml(xmlfil);
        fs = xml.lfpSampleRate;
        k = gaussian2Dfilter([10000 1],1250/2);
        ch = sessiondata.channelID(contains(sessiondata.channel,'BLA'));
        
        if isempty(ch)
            
            ch = sessiondata.channelID(contains(sessiondata.channel,'CA1(R)'));
        end
        totDur =0;
        sessiondata.ts =[];
        sessiondata.theta_delta =[];
        for j = 1:3
            
            dName = fileparts(sessiondata.lfpFile{j});
            d = LoadBinary(sessiondata.lfpFile{j},'nchannels',xml.nChannels,'channels',ch,'frequency',fs);
            
            theta = BandpassFilter(double(d),fs,[5 12]);
            power_theta = InstAmplitude(theta);
            
            delta = BandpassFilter(double(d),fs,[1 4]);
            power_delta = InstAmplitude(delta);
            
            
            
            power_theta = nanconvn(power_theta(1:10:end),k);
            power_delta = nanconvn(power_delta(1:10:end),k);
            ts = ((1:length(power_delta))/(fs/10)) + totDur;
            sessiondata.theta_delta  = [sessiondata.theta_delta;power_theta./power_delta];
            sessiondata.ts =[sessiondata.ts ts];
            
            
            
            if exist(fullfile(dName,'autoDetect.evt.IED'))
                evs = LoadEvents(fullfile(dName,'autoDetect.evt.IED'));
                
                
                for j1 =1:length(sessiondata.channelID)
                    
                    
                    str = ['IED: Ch' num2str(sessiondata.channelID(j1))];
                    
                    if j==1
                        
                        sessiondata.IED{j1} = [evs.time(contains(evs.description,str))];
                    else
                        sessiondata.IED{j1} = [sessiondata.IED{j1};evs.time(contains(evs.description,str))+totDur];
                    end
                end
            end
            
            
            
            
            totDur =  sessiondata.ts(end);
        end
    else
        
        xmlfil = [sessiondata.lfpFile(1:end-3) 'xml'];
        xml = LoadXml(xmlfil);
        fs = xml.lfpSampleRate;
        k = gaussian2Dfilter([10000 1],1250/2);
        ch = sessiondata.channelID(contains(sessiondata.channel,'BLA'));
        
        if isempty(ch)
            
            ch = sessiondata.channelID(contains(sessiondata.channel,'CA1(R)'));
        end
        
        
        dName = fileparts(sessiondata.lfpFile);
        d = LoadBinary(sessiondata.lfpFile,'nchannels',xml.nChannels,'channels',ch,'frequency',fs);
        
        theta = BandpassFilter(double(d),fs,[5 12]);
        power_theta = InstAmplitude(theta);
        
        delta = BandpassFilter(double(d),fs,[1 4]);
        power_delta = InstAmplitude(delta);
        
        
        
        power_theta = nanconvn(power_theta(1:10:end),k);
        power_delta = nanconvn(power_delta(1:10:end),k);
        ts = ((1:length(power_delta))/(fs/10)) ;
        sessiondata.theta_delta  = power_theta./power_delta;
        sessiondata.ts =  ts;
        
        
        
        if exist(fullfile(dName,'autoDetect.evt.IED'))
            evs = LoadEvents(fullfile(dName,'autoDetect.evt.IED'));
            
            
            for j1 =1:length(sessiondata.channelID)
                
                
                str = ['IED: Ch' num2str(sessiondata.channelID(j1))];
                
                
                
                sessiondata.IED{j1} = evs.time(contains(evs.description,str));
                
            end
        end
        
        
    end
    save(files{i},'sessiondata')
    i
end

end
    
    
     

%%

files = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files,'sessiondata');

files = files(keepFiles);

for i = 6
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    sessiondata.stimON = sessiondata.stimON(:);
    if isfield(sessiondata,'IED')
        sessiondata.IED_rate =[];
        for j = 1:3
            
            if ~isnan(sessiondata.stimON)
                
                
                preEpoch = [0 sessiondata.stimON(1)];
                postEpoch = [sessiondata.stimON(end) sessiondata.ts(end)];
                stimEpoch =  [sessiondata.stimON(1:end-1)+2 sessiondata.stimON(2:end)-.25];
                
            else
                preEpoch = [0 3600];
                
                stimEpoch  = 3600 + [[0:30:3000]' [29:30:3029]'];
                postEpoch = [max(stimEpoch(:)) max(stimEpoch(:))+3600];
            end
            
            
            
            switch j
                case 1
                    ep = preEpoch;
                    
                case 2
                    ep = stimEpoch;
                    
                case 3
                    ep = postEpoch;
                    
            end
            ep= MergeEpochs2(ep);
            totDur = sum(diff(ep,[],2));
            status = InIntervals(sessiondata.ts,ep);
            sessiondata.TD(j) = nanmean(sessiondata.theta_delta(status));
            
            
            IED_sub =[];
            for k = 1:length(sessiondata.IED)
                kp = InIntervals(sessiondata.IED{k},ep);
                sessiondata.IED_rate(k,j) = sum(kp)/totDur;
                IED_sub{k} = sessiondata.IED{k}(kp);
            end
            
            mm = nan(8);
            for k1 = 1:length(sessiondata.IED)
                [~,idx1] = ismember(sessiondata.channel{k1},chReg);
                for k2 = 1:length(sessiondata.IED)
                    [~,idx2] = ismember(sessiondata.channel{k2},chReg);
                    if k1~=k2
                        mm(idx1,idx2) = mean(abs(bestmatch(IED_sub{k1},IED_sub{k2}) - IED_sub{k1})<.025);
                    end
                end
            end
            
            sessiondata.IED_syn(:,:,j) = mm;
        end
        spikes.times = sessiondata.IED;
        [binnedPopStart,bin_times]=populationMatrix(spikes,3600,1800,5400,stimEpoch(1,1));
        
        % eliminate stim times
        spikes1.times{1} = stimEpoch(:,1);
        [spikePop,bin_times]=populationMatrix(spikes1,3600,1800,5400,stimEpoch(1,1));
        
        binnedPopStart(:,spikePop>0) = nan;
        sessiondata.binnedPopStart = sum(binnedPopStart);
        
        [binnedPopEnd,bin_times]=populationMatrix(spikes,1800,3600,5400,stimEpoch(end,end));
        
        % eliminate stim times
        spikes1.times{1} =stimEpoch(:,1);
        [spikePop,bin_times]=populationMatrix(spikes1,1800,3600,5400,stimEpoch(end,end));
        
        binnedPopEnd(:,spikePop>0) = nan;
        sessiondata.binnedPopEnd = sum(binnedPopEnd);
        
        
        
    end
    save(files{i},'sessiondata')
    i
end


%%

files = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files,'sessiondata');

files = files(keepFiles);
%cellfun(@(a,b) (a(b+4:b+9)),files,ixx,'UniformOutput',false)
datInj(:,1) = [...
    {'EDS 1.0'} ; ...
    {'EDS 1.3'}; ...
    {'EDS 2.0'}; ...
    {'EDS 2.1'}; ...
    {'EDS 2.2'}; ...
    {'EDS 2.3'}; ...
    {'LDT 10.0'}; ...
    
    ];


datInj{1,2} = datenum(2022,6,17);
datInj{2,2} = datenum(2022,6,18);
datInj{3,2} = datenum(2022,7,7);
datInj{4,2} = datenum(2022,7,7);
datInj{5,2} = datenum(2022,7,8);
datInj{6,2} = datenum(2022,7,8);
datInj{7,2} = datenum(2022,3,25);

datInj{1,3} = '6/17/2022';
datInj{2,3} = '6/18/2022';
datInj{3,3} = '7/7/2022';
datInj{4,3} = '7/7/2022';
datInj{5,3} = '7/8/2022';
datInj{6,3} = '7/8/2022';
datInj{7,3} = '3/24/2022';


for i = 1:length(files)
    load(files{i})
    
    if ~isfield(sessiondata,'days_post_IHKA')
    ixx = regexp(files{i},'RH');
    d = files{i}(ixx+4:ixx+9);
  
    
    if iscell(sessiondata.lfpFile)
           ixx = regexp(sessiondata.lfpFile{1},'22');
       
    d=sessiondata.lfpFile{1}(ixx(2):ixx(2)+5);
    else
        ixx = regexp(sessiondata.lfpFile,'22');
       
    d=sessiondata.lfpFile(ixx(2):ixx(2)+5);
    end
    
      Y = str2num(['20' num2str(d(1:2))]);
    M = str2num(d(3:4));
    D = str2num(d(5:6));
    
    
    
    dat= datenum(Y,M,D);
    
    [~,b] = ismember(sessiondata.subject,datInj(:,1));
    if b>0
        sessiondata.date_injected =  datInj{b,3};
        sessiondata.date_recorded = [num2str(M) '/' num2str(D) '/'  num2str(Y)];
        sessiondata.days_post_IHKA = dat - datInj{b,2};
        save(files{i},'sessiondata')
    end
    i
    end
end
% save days post KA




%%
chReg = [...
    {'PrL(L)' };...
    {'PrL(R)' };...
    {'AVT(L)' };...
    {'BLA(R)' };...
    {'CA1(L)' };...
    {'CA1(R)' };...
    {'gRSC (L)'};...
    {'LDT(L1)'}];

% collect all data
topdir = 'R:\DGregg\NeuralData\LTD 10.0';
files1 = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files1,'sessiondata');

files1 = files1(keepFiles);

topdir = 'R:\DGregg\NeuralData\EDS';
files2 = getAllExtFiles(topdir,'mat',1);
keepFiles = contains(files2,'sessiondata');

files2 = files2(keepFiles);

files = [files1;files2];
files = files2;
IED_rate_control =[];
TD_control =[];
IED_syn_stim =[];
IED_syn_con =[];
IED_rate_stim =[];
TD_stim =[];
stim_sub =[];
con_sub =[];
binnedPopStart =[];
binnedPopEnd =[];
ISI_stim =[];
daysPost_stim = [];
daysPost_con =[];
f_stim= [];
IED_rate_pre =[];
binnedTDEnd =[];
binnedTDStart =[];

for i = 1:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    [~,chIDX] = ismember(chReg,sessiondata.channel);
    if isfield(sessiondata,'IED')
        if strmatch(sessiondata.StimType,'mono') & sessiondata.StimAmp >0
            
            tmp = nan(8,3);
            tmp(chIDX>0,:) = sessiondata.IED_rate(chIDX(chIDX>0),:);
            f_stim= [f_stim; {files{i}}];
            IED_syn_stim = cat(4,IED_syn_stim,sessiondata.IED_syn);
            IED_rate_stim = cat(3,IED_rate_stim,tmp);
            TD_stim = [TD_stim;sessiondata.TD];
            stim_sub = [stim_sub;{sessiondata.subject}];
            binnedPopEnd = [binnedPopEnd;sessiondata.binnedPopEnd];
            binnedPopStart = [binnedPopStart;sessiondata.binnedPopStart];
            
            
            tmp = avghist(sessiondata.ts-sessiondata.stimON(1),sessiondata.theta_delta',-3600:1800);
            binnedTDStart = [binnedTDStart; tmp];
            
            tmp = avghist(sessiondata.ts-sessiondata.stimON(end),sessiondata.theta_delta',-1800:3600);
            binnedTDEnd = [binnedTDEnd; tmp];
            
            
            daysPost_stim = [daysPost_stim;sessiondata.days_post_IHKA];
            ISI_stim = [ISI_stim; median(diff(sessiondata.stimON)) sessiondata.StimAmp];
            
            
            [~,b] = histc(sessiondata.theta_delta(sessiondata.ts<floor(sessiondata.stimON(1))),0:.05:3);
            ok = sum(cell2mat(cellfun(@(a) histoc(a,sessiondata.ts(sessiondata.ts<floor(sessiondata.stimON(1))))',sessiondata.IED,'uni',0)'));
         
            IED_rate_pret =  accumarray(b(b>0),ok(b>0)',[length(0:.05:3) 1],@nanmean,nan);
            
            
            
            IED_rate_pre = [IED_rate_pre; IED_rate_pret'];
            
            
            
        elseif strmatch(sessiondata.StimType,'none')
            tmp = nan(8,3);
            tmp(chIDX>0,:) = sessiondata.IED_rate(chIDX(chIDX>0),:);
            
            IED_syn_con = cat(4,IED_syn_con,sessiondata.IED_syn);
            IED_rate_control = cat(3,IED_rate_control,tmp);
            TD_control = [TD_control;sessiondata.TD];
            con_sub = [con_sub;{sessiondata.subject}];
            daysPost_con = [daysPost_con;sessiondata.days_post_IHKA];
        end
        
    end
    i
end
%%

close all
td1 = binnedTDStart(:,1:3600);
IED1 = binnedPopStart(:,1:3600);

pre =[];
for i = 1:size(td1,1)
  pre(i,:) =  avghist(td1(i,:),IED1(i,:),0:.1:3) ;
end



td2 = [binnedTDStart(:,3600:end-1) binnedTDEnd(:,1:1799)];
IED2 = [binnedPopStart(:,3600:end-1) binnedPopEnd(:,1:1800)];

stim =[];
for i = 1:size(td1,1)
  stim(i,:) =  avghist(td2(i,:),IED2(i,:),0:.1:3) ;
end

p= [];
for i = 1:size(stim,2)
    
    [p(i)] = signtest(pre(:,i),stim(:,i));
end
p(p>.05) = nan;

p(~isnan(p)) = .6;
figure
plotMeanSEM(0:.1:3,pre,'k','yAxisFunction',@nanmean)
plotMeanSEM(0:.1:3,stim,'r','yAxisFunction',@nanmean)

%plot(0:.1:3,p,'o')
xlim([.6 1.3])
%%
close all
k = gaussian2Dfilter([10000 1],1250/16);
% typical example
144
154
for i = 144:length(files)
    load(files{i})
    if isfield(sessiondata,'IED') & ~isempty( sessiondata.stimON)
%load(  'R:\DGregg\NeuralData\EDS\OL3\2-7-2023(12.59)\RHS_230207_130000\EDS 2.2.sessiondata.mat')
figure
ts1 = (0: sessiondata.ts(end)) - sessiondata.stimON(1);
set(gcf,'position',[697    57   560   479]);
plot(sessiondata.ts-sessiondata.stimON(1),sessiondata.theta_delta)
hold on
plot(ts1,nanconvn(histc( sessiondata.IED{2},0: sessiondata.ts(end)),k)*50)
plot([sessiondata.stimON(1) sessiondata.stimON(1)]- sessiondata.stimON(1),[0 7])
plot([sessiondata.stimON(end) sessiondata.stimON(end)]- sessiondata.stimON(1),[0 7])
xlim([-2100 7200])
waitforbuttonpress
close all

    end
end
%%

ISI_stim =[];CA1_AMYG =[];

for i = 1:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    
    if sessiondata.StimAmp>0 && isfield(sessiondata,'IED') && strcmp(sessiondata.StimType,'mono')
        
        
        preEpoch = [0 sessiondata.stimON(1)];
        postEpoch = [sessiondata.stimON(end) sessiondata.ts(end)];
        stimEpoch =  [sessiondata.stimON(1:end-1)+3 sessiondata.stimON(2:end)-.25];
        ok =[];
        for j = 1:3
            switch j
                case 1
                    ep = preEpoch;
                    
                case 2
                    ep = stimEpoch;
                    
                case 3
                    ep = postEpoch;
                    
            end
            
            totDur = sum(diff(ep,[],2));
            
            
            
            [chIDX2] = ismember(sessiondata.channel,'BLA(R)');
            [chIDX1] = ismember(sessiondata.channel,'CA1(R)');
            
            if any(chIDX2)
                
                kp1 = InIntervals(sessiondata.IED{chIDX1},ep);
                kp2 = InIntervals(sessiondata.IED{chIDX2},ep);
                ok(j,:) = CrossCorr(sessiondata.IED{chIDX1}(kp1),sessiondata.IED{chIDX2}(kp2),1,101)/sum(kp1);
                
            end
        end
        CA1_AMYG = cat(3,CA1_AMYG,ok);
        
    end
    i
end

%%

close all
ax  = tight_subplot(2,2)

for i = 1:4
    axes(ax(i))
    if i<4
        imagesc(nanmean(IED_syn_stim(:,:,i,:) ,4),[0 .5])
        %  colormap('bluewhitered')
    else
        imagesc(nanmean(IED_syn_con(:,:,2,:),4),[0 .5])
        %colormap('bluewhitered')
    end
    
end
%%
k  = gaussian2Dfilter([ 1 10000 ],400);
k1 =k;
k2 = k;
k1(1:5000) = 0;

k1 = k1*2;
k2(5001:end) = 0;

k2 = k2*2;

binnedPopStart1 =[];
binnedPopEnd1 =[];
close all
for i = 1:size(binnedPopStart,1)
    
    binnedPopStart1(i,:) = (nanconvn(binnedPopStart(i,:) - nanmean(binnedPopStart(i,1800:3600)),k1));%./ nanmedian(binnedPopStart(i,1:3600));
    binnedPopEnd1(i,:) = nanconvn(binnedPopEnd(i,:) - nanmean(binnedPopStart(i,1800:3600)),k1);%./ nanmedian(binnedPopStart(i,1:3600));
binnedTDStart1(i,:) = nanconvn(binnedTDStart(i,:) - nanmean(binnedTDStart(i,1800:3600)),k1);
binnedTDEnd1(i,:) = nanconvn(binnedTDEnd(i,:) - nanmean(binnedTDEnd(i,1800:3600)),k1);

end
%%
binnedPopStart1(isinf(binnedPopStart1)) = nan;
binnedPopEnd1(isinf(binnedPopEnd1)) = nan;
close all
figure
plotMeanSEM(-3600:1799-60,binnedPopStart1(:,1:end-60),'k')
hold on
plotMeanSEM(1801+100:7200,(binnedPopEnd1(:,101:end)),'k')


plotMeanSEM(-3600:1799-60,binnedTDStart1(:,1:end-61),'r')
hold on
plotMeanSEM(1801+100:7200,(binnedTDEnd1(:,102:end)),'r')


set(gca,'xtick',[-3600:600:1730 1800:600:7200],'xticklabel',[-3600:600:1730 -1800:600:3600])
xlim([-1800 5400])
pp=[];
 for i = 1:5400
pp(i) = signtest(binnedPopStart1(:,i));
 end

 pp(pp>.05) = nan;
 pp(~isnan(pp)) = .27;
 
 plot(-3600:1799-60,pp(1:end-60),'k','linewidth',6)
 
 pp =[];
 for i = 1:5400
pp(i) = signtest(binnedPopEnd1(:,i));
 end

 pp(pp>.05) = nan;
 pp(~isnan(pp)) = .27;
 
 plot(1801+100:7200,pp(101:end),'k','linewidth',6)
 
 
 
 pp=[];
 for i = 1:5400
pp(i) = signtest(binnedTDStart1(:,i));
 end

 pp(pp>.05) = nan;
 pp(~isnan(pp)) = .25;
 
 plot(-3600:1799-60,pp(1:end-60),'r','linewidth',6)
 
 pp =[];
 for i = 1:5400
pp(i) = signtest(binnedTDEnd1(:,i));
 end

 pp(pp>.05) = nan;
 pp(~isnan(pp)) = .25;
 
 plot(1801+100:7200,pp(101:end),'r','linewidth',6)
 
 
 plot([0 0],[-.1 .05],'k')
 plot([3600 3600],[-.1 .05],'k')
%%
IED_mean_stim_12 =[];
for i = 1:7
    
    if i<7
        IED_mean_stim_12(i,:) = (IED_rate_stim(i,2,:) - IED_rate_stim(i,1,:))./(IED_rate_stim(i,1,:));
    else
        IED_mean_stim_12(i,:) = (sum(IED_rate_stim(7:8,2,:),1) - sum(IED_rate_stim(7:8,1,:),1))./( sum(IED_rate_stim(7:8,1,:),1));
    end
end


IED_mean_stim_23 =[];
for i = 1:7
    
    
    IED_mean_stim_23(i,:) = (IED_rate_stim(i,2,:) - IED_rate_stim(i,3,:))./(IED_rate_stim(i,2,:) + IED_rate_stim(i,3,:));
end

IED_mean_stim_13 =[];
for i = 1:7
    
    if i<7
        IED_mean_stim_13(i,:) = (IED_rate_stim(i,3,:) - IED_rate_stim(i,1,:))./(IED_rate_stim(i,3,:) + IED_rate_stim(i,1,:));
    else
        IED_mean_stim_13(i,:) = (sum(IED_rate_stim(7:8,3,:),1) - sum(IED_rate_stim(7:8,1,:),1))./(sum(IED_rate_stim(7:8,3,:),1) + sum(IED_rate_stim(7:8,1,:),1));
    end
end

TD_stim_mean = (TD_stim(:,2) - TD_stim(:,1))./( TD_stim(:,1));
IED_mean_stim_12(isinf(IED_mean_stim_12)) = nan;


%%
usub = unique(stim_sub);

[~,b] = ismember(stim_sub,usub);

close all
for i = 1:6
    figure
    
    violin(num2cell(IED_mean_stim_12(:,b==i),2)')
    subj_mean(i,:) = mean(IED_mean_stim_12(:,b==i),2);
    title(usub{i})
end

%%
figure


x = [IED_mean_stim_12(:)];
g = [linearize(repmat([1:7]',1,size(IED_mean_stim_12,2)))];


boxplot(x,g,'notch','on')
hold on
plot([0 8],[0 0],'k')



%%
ok = squeeze((nansum(IED_rate_stim(:,:,:))- repmat(nansum(IED_rate_stim(:,1,:)),1,3,1))./repmat(nansum(IED_rate_stim(:,1,:)),1,3,1));
ok1 = squeeze((nansum(IED_rate_control(:,:,:))- repmat(nansum(IED_rate_control(:,1,:)),1,3,1))./repmat(nansum(IED_rate_control(:,1,:)),1,3,1));


ok = squeeze((nansum(IED_rate_stim(:,:,:),1)));
%ok1 = squeeze((nansum(IED_rate_control(:,:,:),1)));


%ok = ok(:,daysPost_stim>190);
%ok1 = ok1(:,daysPost_con>190);
x = [ok(:);ok1(:)];
g = [linearize(repmat([1:3]',1,size(ok,2))); linearize(repmat([4:6]',1,size(ok1,2)))];
figure
boxplot(x,g,'notch','on')
%%
close all
figure
ok = squeeze((nansum(IED_rate_stim(:,:,:),1)));
plot(0:.01:1,cumsum(histc(ok(1,:),0:.01:1))/57)
hold on
plot(0:.01:1,cumsum(histc(ok(2,:),0:.01:1))/57)

figure
plot(.6:.01:1.6,cumsum(histc(TD_stim(:,1),.6:.01:1.6))/57)
hold on
plot(.6:.01:1.6,cumsum(histc(TD_stim(:,2),.6:.01:1.6))/57)

%%
usub = unique(stim_sub);

[~,b] = ismember(stim_sub,usub);

IED_mean_stim_12_all = squeeze(nansum(IED_rate_stim));
IED_mean_stim_12_all = (IED_mean_stim_12_all(2,:)-IED_mean_stim_12_all(1,:))./(IED_mean_stim_12_all(1,:));
close all
figure
col = linspecer(7,'jet');
for i = 1:7
    
    
    plot(TD_stim_mean(b==i),IED_mean_stim_12_all(b==i),'.','color',col{i},'markersize',30)
    hold on
    if i ==1
        plot([-.15 .25],[0 0],'k')
        plot([0 0],[-1 1],'k')
    end
end

%%
N=0;
for i = 1:length(files)
    
    load(files{i})
    if strcmp(sessiondata.subject,'EDS 1.1')
        cd(fileparts(files{i}))
        disp('here')
    end
end

%%
[~,b] = ismember(stim_sub,usub);
[~,b1] = ismember(con_sub,usub);
close all
col = flipud(linspecer(6,'jet'))
for i = 1:7
    
    plot(daysPost_stim(b==i),squeeze(sum(IED_rate_stim(:,2,b==i),1)) ,'x','color',col{i})
    hold on
    plot(daysPost_con(b1==i),squeeze(sum(IED_rate_control(:,2,b1==i),1)) ,'o','color',col{i})
end
plot([160 240],[0 0],'k')

figure

col = flipud(linspecer(6,'jet'))
for i = 1:6
    
    plot(daysPost_stim(b==i),squeeze(sum(IED_rate_stim(:,1,b==i),1)) ,'x','color',col{i})
    hold on
    plot(daysPost_con(b1==i),squeeze(sum(IED_rate_control(:,1,b1==i),1)) ,'o','color',col{i})
end
plot([160 240],[0 0],'k')

%%
[~,b_stim] = ismember(stim_sub,usub);
[~,b_con] = ismember(con_sub,usub);

N_stim = size(TD_stim,1);
N_con = size(TD_control,1);

IEDtot_stim = squeeze((nansum(IED_rate_stim(:,:,:),1)))';
IEDtot_con = squeeze((nansum(IED_rate_control(:,:,:),1)))';
 
%b_stim = repmat(b_stim,1,3);
condition_stim = repmat(1:3,N_stim,1);
session_stim = repmat([1:N_stim]',1,3);
 


%b_con = repmat(b_con,1,3);
condition_con = repmat(1:3,N_con,1);
session_con = repmat([1:N_con]',1,3);
 
IED = [IEDtot_stim(:,2) - IEDtot_stim(:,1);IEDtot_con(:,2)-IEDtot_con(:,1)];
TD = [TD_stim(:,2)-TD_stim(:,1);TD_control(:,2)-TD_control(:,1)];
%condition = categorical([condition_stim(:);condition_con(:)]);
subject = usub([b_stim(:);b_con(:)]);
%session = [session_stim(:);session_con(:)];
stim_on = (categorical([ones(numel(b_stim(:)),1);2*ones(numel(b_con(:)),1)]));


tbl = table(IED,TD,stim_on,subject);
lme = fitlme(tbl,'IED ~  -1+ stim_on*TD+(1/subject) ')


