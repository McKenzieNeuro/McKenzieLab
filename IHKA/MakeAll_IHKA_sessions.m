



% make all LFP

dirs = readtable('R:\DGregg\NeuralData\thetaInduction_Stim2.xlsx');

dirsTop = strrep(dirs.RootFolder,'McKenzieLab\','');

kp = ~cellfun(@isempty,dirsTop) | ~contains(dirsTop,'hr');
dirsTop = dirsTop(kp);

subfolder = dirs.DataFolder(kp);
dirs1 = cellfun(@(a,b) [a filesep b],dirsTop,subfolder,'uni',0);

%%
for i = 135:length(dirs1)
    cd(dirs1{i})
    
    if ~exist('amplifier.lfp')
        %error('here')
        bz_LFPfromDat(pwd,'basename','amplifier')
    end
    i
end

lfpsubdir = unique(fileparts(dirs));

%%

%check recInfo
makeRec = false(length(makeRec),1);
for kk = 1:length(lfpsubdir)
     cd(lfpsubdir{kk});
     
     if ~exist('recInfo.mat')
         makeRec(kk) = true;
     end
end
    
%%


% make all LFP

dirs = readtable('R:\DGregg\NeuralData\thetaInduction.xlsx');

dirs = strrep(dirs.Folder,'McKenzieLab\','');
dirs = dirs(~cellfun(@isempty,dirs));
dirs = dirs(~contains(dirs,'hr'));
%%
for i = 1:length(dirs)
    cd(dirs{i})
    
    if ~exist('amplifier.lfp')
        %error('here')
        bz_LFPfromDat(pwd,'basename','amplifier')
    end
    i
end

lfpsubdir = unique(fileparts(dirs));
%%
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
%%
topdir = 'R:\DGregg\NeuralData\PCP\Recordings\OpenLoop\';
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

%make sessiondata struct base on lfp files
% store meta data about recording sessions

for kk = 1:length(lfpsubdir)
    
    cd(lfpsubdir{kk});
    clear recInfo
    load('recInfo.mat')
    
    %figure out which recInfo format
    
    if size(recInfo,2)==3
        oldInfo = true;
    else
        oldInfo = false;
    end
    
    
    %ISI = load('config.mat','ISI');
    %ISI = mean(ISI.ISI);
    
    lfpfiles = getAllExtFiles(lfpsubdir{kk},'lfp',1);
    
    
    %loop through systems (RHD/RHS)
    for jj = 1:2
        
        switch jj
            case 1
                if oldInfo
                    syst = find(contains(recInfo(:,2),'RHD'));
                    
                else
                    syst = find(contains({recInfo.system},'RHD'));
                end
                
                kp = contains(lfpfiles,lfpsubdir{kk}) & contains(lfpfiles,'RHD');
                
            case 2
                if oldInfo
                    syst = find(contains(recInfo(:,2),'RHS'));
                else
                    syst = find(contains({recInfo.system},'RHS'));
                end
                kp = contains(lfpfiles,lfpsubdir{kk}) & contains(lfpfiles,'RHS');
                
        end
        
        
        if sum(kp) ==1
            
            lfp = lfpfiles{kp};
            
            
            
            xmlfil = [lfp(1:end-3) 'xml'];
            xml = LoadXml(xmlfil);
            %loop through animals on each system (2per)
            totCh =0;
            
            % check if recInfo is made
            
            if size(recInfo,2)>3 && oldInfo
                
                recInfo = sm_recfileAnalysis(lfpsubdir{kk});
                
            end
            
            [outDir] = fileparts(lfp);
            
            
            if oldInfo
                nSubj = 2;
            else
                nSubj = length(recInfo.fileData);
            end
            
            for ii = 1:nSubj
                clear sessiondata
                skip = false;
                if oldInfo
                    subjectName = recInfo{syst,3}{1+ii,2}{1};
                elseif ~isempty(recInfo.fileData(ii).subject)
                    
                    if iscell(recInfo.fileData(ii).subject)
                        
                    subjectName = recInfo.fileData(ii).subject{1};
                    
                    else
                           subjectName = recInfo.fileData(ii).subject;
                    end
                else 
                    skip = true;
                end
                
                outfil = fullfile(outDir,[subjectName '.sessiondata.mat']);
                
                
                if ~exist(outfil) && ~skip
                    if oldInfo && isempty(recInfo{syst,3}{2,7})
                        
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
                    elseif oldInfo
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
                        
                        
                    else
                        
                            % check in the stim.dat
                            f = getAllExtFiles(lfpsubdir{kk},'dat',1);
                            
                            kp = (contains(f,'stim'));
                            
                            if sum(kp)==1
                                f = f{kp};
                                ch = recInfo.fileData(ii).stim1_chan(1)+1; % base 1
                                
                                if contains(f,'1KHz')
                                    st = LoadBinary(f,'nchannels',xml.nChannels,'frequency',1000,'channels',ch);
                                    fs = 1000;
                                else
                                    st = LoadBinary(f,'nchannels',xml.nChannels,'frequency',xml.SampleRate,'channels',ch);
                                  %st=[];
                                    fs = xml.SampleRate;
                                end
                                stimON = find(diff([0;(st>2000)])>0)/fs;
                                stimON = stimON(diff([0;stimON])>1);
                            else
                                disp('here')
                                
                            end
                            stimType = recInfo.fileData(ii).stim1_type;
                            stimAmp = recInfo.fileData(ii).stim1_uA;
                      
                            if isfield(recInfo.fileData(ii),'data')
                                chName = {recInfo.fileData(ii).data.site}';
                                tmp = {recInfo.fileData(ii).data.chan}';
                                chNum = cellfun(@(a) str2num(a(3:end)),tmp);
                                
                            else
                                chName =  [...
                                    {'PrL(L)' };...
                                    {'PrL(R)' };...
                                    {'AVT(L)' };...
                                    {'BLA(R)' };...
                                    {'CA1(L)' };...
                                    {'CA1(R)' };...
                                    {'LDT(L1)'};...
                                    {'LDT(L2)'}];
                                chNum = 0:7;
                            end
                        
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

% select session data files


files = getAllExtFiles('R:\DGregg\NeuralData','mat',1);
%files = getAllExtFiles('R:\DGregg\NeuralData\PCP\Recordings\OpenLoop\','mat',1);
keepFiles = contains(files,'sessiondata');

files = files(keepFiles);

%%
% get all good sessions
files = [];
for i = 1:length(lfpsubdir)
    
    tmp = getAllExtFiles(lfpsubdir{i},'mat',1);
    keepFiles = contains(tmp,'sessiondata');
    
    tmp = tmp(keepFiles);
    files = [files;tmp];
end

%%

% save days post KA

clear datInj
datInj(:,1) = [...
    {'EDS 1.0'} ; ...
    {'EDS 1.3'}; ...
    {'EDS 2.0'}; ...
    {'EDS 2.1'}; ...
    {'EDS 2.2'}; ...
    {'EDS 2.3'}; ...
    {'LDT 10.0'}; ...
    {'EDS 3.0'}; ...
    {'EDS 3.2'}; ...
    {'EDS 4.0'}; ...
    {'EDS 4.1'}; ...
    {'EDS 4.2'}; ...
    {'EDS 5.1'}; ...
    {'EDS 1.1'} ; ...
    {'PCP 4.0'} ; ...
    ];


datInj{1,2} = datenum(2022,6,17);
datInj{2,2} = datenum(2022,6,18);
datInj{3,2} = datenum(2022,7,7);
datInj{4,2} = datenum(2022,7,7);
datInj{5,2} = datenum(2022,7,8);
datInj{6,2} = datenum(2022,7,8);
datInj{7,2} = datenum(2022,3,25);
datInj{8,2} = datenum(2023,3,21);
datInj{9,2} = datenum(2023,3,23);
datInj{10,2} = datenum(2023,3,27);
datInj{11,2} = datenum(2023,3,24);
datInj{12,2} = datenum(2023,3,27);
datInj{13,2} = datenum(2023,4,20);
datInj{14,2} = datenum(2022,6,17);
datInj{15,2} = datenum(2024,3,18);

datInj{1,3} = '6/17/2022';
datInj{2,3} = '6/18/2022';
datInj{3,3} = '7/7/2022';
datInj{4,3} = '7/7/2022';
datInj{5,3} = '7/8/2022';
datInj{6,3} = '7/8/2022';
datInj{7,3} = '3/24/2022';
datInj{8,3} = '3/21/2023';
datInj{9,3} = '3/23/2023';
datInj{10,3} = '3/27/2023';
datInj{11,3} = '3/24/2023';
datInj{12,3} = '3/27/2023';
datInj{13,3} = '4/20/2023';
datInj{14,3} = '6/17/2022';
datInj{15,3} = '3/18/2024';

%%
for i = 1:length(files)
    load(files{i})
    
    if ~isfield(sessiondata,'days_post_pilo')
        ixx = regexp(files{i},'RH');
        d = files{i}(ixx+4:ixx+9);
        
        
        if iscell(sessiondata.lfpFile)
            ixx = regexp(sessiondata.lfpFile{1},'_23');
            
            if isempty(ixx)
                ixx = regexp(sessiondata.lfpFile{1},'_22');
            end
            
            d=sessiondata.lfpFile{1}(ixx(1)+1:ixx(1)+6);
        else
            ixx = regexp(sessiondata.lfpFile,'_24');
            if isempty(ixx)
                ixx = regexp(sessiondata.lfpFile,'_25');
            end
            d=sessiondata.lfpFile(ixx(1)+1:ixx(1)+6);
        end
        
        Y = str2num(['20' num2str(d(1:2))]);
        M = str2num(d(3:4));
        D = str2num(d(5:6));
        
        
        
        dat= datenum(Y,M,D);
        
        [~,b] = ismember(sessiondata.subject,datInj(:,1));
        if b>0
            sessiondata.date_injected =  datInj{b,3};
            sessiondata.date_recorded = [num2str(M) '/' num2str(D) '/'  num2str(Y)];
            sessiondata.days_post_pilo = dat - datInj{b,2};
            save(files{i},'sessiondata')
        else
            disp('here')
        end
        
    end
    i
end




%%

for i = 182:length(files)
  cd(fileparts(files{i}))  
  if exist('events.evt.szr')
  movefile('events.evt.szr','manualDetect.evt.szr')
  end
  i
end

%%
%annotate seizures

for i = 1:length(files)
    cd(fileparts(files{i}))
    fils = getAllExtFiles(pwd,'szr',1);
    if isempty(fils)
        load(files{i})
        if iscell(sessiondata.lfpFile)
            for j = 1:length(sessiondata.lfpFile)
                
             %   sessiondata.lfpFile{j} = strrep(sessiondata.lfpFile{j},'NeuralData\','NeuralData\LTD\');
                cd(fileparts(sessiondata.lfpFile{j} ))
                if ~exist('manualDetect.evt.szr') && exist(files{i})
                    cmd = ['neuroscope "' sessiondata.lfpFile{j} '"'];
                    
                    system(cmd)
                end
            end
        else
            
            cmd = ['neuroscope "' sessiondata.lfpFile '"'];
            
            system(cmd)
        end
    end
    i
end


%%


% link theta delta and IED info with sessiondata


% get all IEDs

for i = 1:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    if ~isfield(sessiondata,'theta_delta') ||  ~isfield(sessiondata,'IED')
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
            
            
            
            if ~exist(fullfile(dName,'autoDetect.evt.IED'))
                
                st = sm_detectIED(dName);
                
            end
            
            
            
            evs = LoadEvents(fullfile(dName,'autoDetect.evt.IED'));
            
            
            for j1 =1:length(sessiondata.channelID)
                
                
                str = ['IED: Ch' num2str(sessiondata.channelID(j1))];
                
                
                
                sessiondata.IED{j1} = evs.time(contains(evs.description,str));
                
            end
        end
        save(files{i},'sessiondata')
        
    end
    
    i
end



%%

% link seizures with TD/IED
for i = 1:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
   
        
    if isfield(sessiondata,'stimON') %&& ~isfield(sessiondata,'seizure_rate')
        
        
        if iscell(sessiondata.lfpFile)
            totDur = 0;
            sessiondata.seizure = [];
            for j = 1:length(sessiondata.lfpFile)
                %sessiondata.lfpFile{j} =  strrep(sessiondata.lfpFile{j},'NeuralData\','NeuralData\LTD\');
                cd(fileparts(sessiondata.lfpFile{j}))
                if  exist('manualDetect.evt.szr')
                    ev = LoadEvents('manualDetect.evt.szr');
                    sub = strrep(sessiondata.subject,' ','');
                    
                    if ~isempty(ev.description)
                        kp_on = contains(ev.description,['sz_on_' sub]);
                        kp_off = contains(ev.description,['sz_off_' sub]);
                        
                        ep = [ev.time(kp_on) ev.time(kp_off)]+totDur;
                        
                        if any(ep(:,1)>ep(:,2))
                            
                            error('check events')
                            
                        else
                            sessiondata.seizure = [sessiondata.seizure ;ep];
                        end
                        
                        
                        
                        
                    end
                    
                    
                end
                totDur = totDur + sessiondata.fileDur(j);
            end
            
            
        else %if  ~isfield(sessiondata,'seizure_rate')
            
         
            
            
            
            fil = getAllExtFiles(pwd,'szr',1);
            if  ~isempty(fil) && length(fil) ==1
                ev = LoadEvents(fil{1});
                sub = strrep(sessiondata.subject,' ','');
                
                if ~isempty(ev.description)
                    
                    if ~contains(sessiondata.subject,'PCP')
                        kp_on = contains(ev.description,['sz_on_' sub]);
                        kp_off = contains(ev.description,['sz_off_' sub]);
                        
                        
                    else
                        
                        kp_on = contains(ev.description,'SeizureStart');
                        kp_off = contains(ev.description,'SeizureStop');
                    end
                    ep = [ev.time(kp_on) ev.time(kp_off)];
                    
                    if any(ep(:,1)>ep(:,2))
                        
                        error('check events')
                        
                    else
                        sessiondata.seizure = ep;
                    end
                    
                    
                else
                    sessiondata.seizure = [];
                end
                
            else
                
                error('annotate seizures')
            end
            
            
        end
        
        
          if ~isempty(sessiondata.stimON) && ~ all(isnan(sessiondata.stimON))
                eps = [0 sessiondata.stimON(1) sessiondata.stimON(end) sessiondata.ts(end)];
            else
                
                eps = [0 3600 2*3600 sessiondata.ts(end)];
          end
             dur = diff(eps);
          
        if ~isempty(sessiondata.seizure)
          
            
            [n,b] = histc(sessiondata.seizure(:,1),eps);
            n = n(1:3);
           sz_dur = accumarray(b,diff(sessiondata.seizure,[],2),[3 1],@nanmean,nan);
            sessiondata.seizure_rate = [n(:)./dur(:)]';
            sessiondata.seizure_num = n;
            sessiondata.seizure_duration = sz_dur';
        else
            
            sessiondata.seizure_rate = [0 0 0];
            sessiondata.seizure_num = [ 0 0 0];
              sessiondata.seizure_duration = [nan nan nan];
        end
       
        sessiondata.seizure_obs = dur;
    
        save(files{i},'sessiondata')
        i
    end
    
end



%%

% do some IED/TD analysis
for i = 182:length(files)
    dirN = fileparts(files{i});
    cd(dirN)
    
    load(files{i})
    sessiondata.stimON = sessiondata.stimON(:);
    
    % remove IEDs during seizure
    if ~isempty(sessiondata.seizure)
        
        for j = 1:length(sessiondata.IED)
            
            kp = ~InIntervals(sessiondata.IED{j},sessiondata.seizure);
            sessiondata.IED{j} = sessiondata.IED{j}(kp,:);
        end
        
    end
    
   % if isfield(sessiondata,'IED') && ~isfield(sessiondata,'binnedPopEnd')
        sessiondata.IED_rate =[];
        for j = 1:3
            
            if ~isnan(sessiondata.stimON)
                
                
                preEpoch = [0 sessiondata.stimON(1)];
                postEpoch = [sessiondata.stimON(end) sessiondata.ts(end)];
                stimEpoch =  [sessiondata.stimON(1:end-1)+2 sessiondata.stimON(2:end)-.25];
                
            else
                preEpoch = [0 3600];
                
                stimEpoch  = 3600 + [[0:30:3600-30]' [29:30:3629-30]'];
                postEpoch = [max(stimEpoch(:)) sessiondata.ts(end)];
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
            if ~isempty(sessiondata.seizure)
                ep = excludeEpochs(ep,sessiondata.seizure);
            end
            totDur = sum(diff(ep,[],2));
            status = InIntervals(sessiondata.ts,ep);
            sessiondata.TD(j) = nanmean(sessiondata.theta_delta(status));
            
            
            IED_sub =[];
            sessiondata.IED_dur(j) = totDur;
            
            for k = 1:length(sessiondata.IED)
                kp = InIntervals(sessiondata.IED{k},ep);
                sessiondata.IED_num(k,j) =  sum(kp);
                sessiondata.IED_rate(k,j) = sum(kp)/totDur;
                IED_sub{k} = sessiondata.IED{k}(kp);
            end
            
            mm = nan(8);
            for k1 = 1:length(sessiondata.IED)
                [~,idx1] = ismember(sessiondata.channel{k1},chReg);
                for k2 = 1:length(sessiondata.IED)
                    [~,idx2] = ismember(sessiondata.channel{k2},chReg);
                    if k1~=k2 && idx1>0 && idx2>0
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
         
         
         %eliminate seizures
         if ~isempty(sessiondata.seizure)
             tts = stimEpoch(1,1)+bin_times;
             kp = ~InIntervals(tts,sessiondata.seizure);
             binnedPopStart(:,kp) = nan;
         end
        
        
       
        sessiondata.binnedPopStart = nansum1(binnedPopStart);
        
        [binnedPopEnd,bin_times]=populationMatrix(spikes,1800,3600,5400,stimEpoch(end,end));
        
        % eliminate stim times
        spikes1.times{1} =stimEpoch(:,1);
        [spikePop,bin_times]=populationMatrix(spikes1,1800,3600,5400,stimEpoch(end,end));
        
        binnedPopEnd(:,spikePop>0) = nan;
      
         %eliminate seizures
         if ~isempty(sessiondata.seizure)
             tts = stimEpoch(end,end)+bin_times;
             kp = ~InIntervals(tts,sessiondata.seizure);
             binnedPopEnd(:,kp) = nan;
         end
          sessiondata.binnedPopEnd = nansum1(binnedPopEnd);
        save(files{i},'sessiondata')
   % end
    
    i
end





%%

%now link state and pose

for i = 1:length(files)
    cd(fileparts(files{i}))
    v = load(files{i});
    sessiondata = v.sessiondata;
    
    %if ~isfield(sessiondata,'sleepData')
    dirN = fileparts(files{i});
    
    %find associated video
    f= getAllExtFiles(dirN,'mp4',0);
    
    if isempty(f)
        f= getAllExtFiles(fileparts(dirN),'mp4',0);
    end
    
    % f = f(contains(f,sessiondata.subject) & ~contains(f,'Pose')  & ~contains(f,'Sync'));
    f = f(~contains(f,'Pose')  & ~contains(f,'Sync') & ~contains(f,'_O')); % for LDT
    if ~isempty(f)
        
        if (length(f)==1 || ~contains(files{i},'LTD'))
            
            if ~all(contains(files{i},'LTD'))  && any(contains(f,sessiondata.subject))
                
                
                f = f(contains(f,sessiondata.subject));
            end
            if length(f)>1
                kp = contains(f,'(2)');
                f= f(kp);
            end
            if length(f)~=1
                error('too many')
            end
            
            %
            
            
            
            sessiondata.videoFile = f{1};
            poseFile = strrep(sessiondata.videoFile,'.mp4','_poseData.mat');
            sleepFile = strrep(sessiondata.videoFile,'.mp4','_sleepData2.mat');
            if exist(poseFile)
                sessiondata.poseData = load(poseFile);
                fprintf([num2str(i) ' poseFound '])
            end
            if exist(sleepFile)
                sessiondata.sleepData = load(sleepFile);
                fprintf(['sleepFound \n'])
            else
                fprintf(['\n'])
            end
            
            
            
            
            
            
            save(files{i},'sessiondata')
            
        else
            %deal with LDT
            clear idx
            idx(1) = find(contains(f,'Initial'));
            idx(2) = find(contains(f,'uA'));
            idx(3) = find(contains(f,'Final'));
            
            
            if any(idx==0)
                error('no video for LDT')
            else
                
                for ff = 1:3
                    sessiondata.videoFile{ff} = f{idx(ff)};
                    poseFile = strrep(sessiondata.videoFile{ff},'.mp4','_poseData.mat');
                    sleepFile = strrep(sessiondata.videoFile{ff},'.mp4','_sleepData2.mat');
                    if exist(poseFile)
                        
                        if ff ==1
                            sessiondata.poseData = load(poseFile);
                        else
                            tmp = load(poseFile);
                            offset = sum(sessiondata.fileDur(1:ff-1));
                            sessiondata.poseData.XY =  [sessiondata.poseData.XY; tmp.XY];
                            sessiondata.poseData.neuralTS =  [sessiondata.poseData.neuralTS; tmp.neuralTS+offset];
                            sessiondata.poseData.velocity =  [sessiondata.poseData.velocity; tmp.velocity];
                        end
                        
                        
                    end
                    
                    if exist(sleepFile)
                        
                        if ff ==1
                            sessiondata.sleepData = load(sleepFile);
                        else
                            tmp = load(sleepFile);
                            offset = sum(sessiondata.fileDur(1:ff-1));
                            sessiondata.sleepData.SleepState.ints.WAKEstate =  [sessiondata.sleepData.SleepState.ints.WAKEstate; tmp.SleepState.ints.WAKEstate+offset];
                            sessiondata.sleepData.SleepState.ints.NREMstate =  [sessiondata.sleepData.SleepState.ints.NREMstate; tmp.SleepState.ints.NREMstate+offset];
                            sessiondata.sleepData.SleepState.ints.REMstate =  [sessiondata.sleepData.SleepState.ints.REMstate; tmp.SleepState.ints.REMstate+offset];
                            sessiondata.sleepData.SleepState.idx.states =  [ sessiondata.sleepData.SleepState.idx.states ;  tmp.SleepState.idx.states ];
                            sessiondata.sleepData.SleepState.idx.timestamps =  [ sessiondata.sleepData.SleepState.idx.timestamps ;  tmp.SleepState.idx.timestamps+offset ];
                            
                        end
                        
                        
                    end
                    
                    
                end
                
                 save(files{i},'sessiondata')
            end
            
        end
        
        
    else
        warning('missing videos')
        
    end
    
    
end


    % end
%%

for i = 182:length(files)
    
    load(files{i})
    
    if contains(sessiondata.StimType,'mono') & (isempty(sessiondata.stimON) & sessiondata.StimAmp>0)
        
        sessiondata.StimType = 'none';
        save(files{i},'sessiondata')
    end
    i
end


%%