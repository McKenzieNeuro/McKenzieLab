clear all %#ok<CLALL> 

loopType = 'open'; % Experiment type, 'open' or 'closed'.

StartTime = [13,0]; % Time to begin recording in 24hr time, format = [hr,min]
nextDay = 0; % Start recording tomorrow?

% Block schedule. 
% BL = Baseline, S1 = Stim Loop, S2 = Stim Train, SX = Stim sham, OD = Optional delay
% [{Day1 block},{Day2 block},....];

% % Port A
% blockSched{1,1} = {'BL'};
% 
% % Port B
% blockSched{1,2} = {'BL'};
% 
% % Port C
% blockSched{1,3} = 'N/A';
% 
% % Port D
% blockSched{1,4} = 'N/A';
% 
% % Block schedule in seconds. 
% blockTime(1,1) = 3600*12;

% Port A
blockSched{1,1} = [{'BL'};{'BL'};{'BL'};{'BL'}];
blockSched{2,1} = [{'SX'};{'S1'};{'SX'};{'S1'}];
blockSched{3,1} = [{'S2'};{'S2'};{'S2'};{'S2'}];
blockSched{4,1} = [{'BL'};{'BL'};{'BL'};{'BL'}];

% Port B
blockSched{1,2} = [{'BL'};{'BL'};{'BL'};{'BL'}];
blockSched{2,2} = [{'S1'};{'SX'};{'S1'};{'SX'}];
blockSched{3,2} = [{'S2'};{'S2'};{'S2'};{'S2'}];
blockSched{4,2} = [{'BL'};{'BL'};{'BL'};{'BL'}];

% Port C
blockSched{1,3} = 'N/A';
blockSched{2,3} = 'N/A';
blockSched{3,3} = 'N/A';
blockSched{4,3} = 'N/A';

% Port D
blockSched{1,4} = 'N/A';
blockSched{2,4} = 'N/A';
blockSched{3,4} = 'N/A';
blockSched{4,4} = 'N/A';

% Block schedule in seconds. 
blockTime(1,1) = 3600;
blockTime(2,1) = 600;
blockTime(3,1) = 10;
blockTime(4,1) = 13790;

% Video Settings:
fps = 25;
RecVid = 1; % Record from webcams, 1 = yes, 0 = no
ShowVideo = 0; % Show video feed from cameras, 1 = yes, 0 = no

% Arduino Settings
ArdPin{1,1} = 'D13'; % Arduino pin for LED.
ArdPin{2,1} = 'D8'; % Arduino pin for Timestamp.
LEDtime = 10; % Time between LED blinks is seconds

% Intan TCP
useRHS = 1;
intanIP_RHS = '127.0.0.1'; % IP for RHS.
intanPort1_RHS = 5000; % Port for RHS commands.
intanPort2_RHS = 5001; % Port for RHS data output.

useRHD = 0;
intanIP_RHD = '127.0.0.2'; % IP for RHD.
intanPort1_RHD = 5000; % Port for RHD commands.
intanPort2_RHD = 5001; % Port for RHD data output.

% Intan channel ID
chanID(1,1:2) = [{'000'},{'PrL(L)'}];
chanID(2,1:2) = [{'001'},{'PrL(R)'}];
chanID(3,1:2) = [{'002'},{'AVT(L)'}];
chanID(4,1:2) = [{'003'},{'BLA(R)'}];
chanID(5,1:2) = [{'004'},{'CA1(L)'}];
chanID(6,1:2) = [{'005'},{'CA1(R)'}];
chanID(7,1:2) = [{'006'},{'LDT(L1)'}];
chanID(8,1:2) = [{'007'},{'LDT(L2)'}];

boxSubject{1,1} = 'EDS 2.2'; % Subject in box 1
boxSubject{2,1} = 'EDS 2.1'; % Subject in box 2
boxSubject{3,1} = 'EDS 1.1'; % Subject in box 3
boxSubject{4,1} = 'EDS 2.3'; % Subject in box 4

boxRHS(1,1) = 3; % Box for RHS port A
boxRHS(2,1) = 4; % Box for RHS port B
boxRHS(3,1) = 0; % Box for RHS port C
boxRHS(4,1) = 0; % Box for RHS port D

boxRHD(1,1) = 1; % Box for RHD port A
boxRHD(2,1) = 2; % Box for RHD port B
boxRHD(3,1) = 0; % Box for RHD port C
boxRHD(4,1) = 0; % Box for RHD port D

startIter = 2;

% Stim 1 parameters
ISI = {20:30}; % Inter stim interval in seconds, randomly selected.
thresholdTD = 1; % Theda/delta threshold for stim

% Amplitudes. [day1,day2...]
stimI(1,:,1) = [125,125,125,125]; % Port A stim amp in uA.
stimI(2,:,1) = [125,125,125,125]; % Port B stim amp in uA.
stimI(3,:,1) = [0,0,0,0]; % Port C stim amp in uA.
stimI(4,:,1) = [0,0,0,0]; % Port D stim amp in uA.

% Channel, [0,0] = unused, [c1,0] = monopolar, [c1,c2] = bipolar.
stimChan(1,:,1) = [7,0]; % Port A
stimChan(2,:,1) = [7,0]; % Port B
stimChan(3,:,1) = [0,0]; % Port C
stimChan(4,:,1) = [0,0]; % Port D

% Train period
stimP(1,1,1) = 10; % Port A stim train period in ms.
stimP(2,1,1) = 10; % Port B stim train period in ms.
stimP(3,1,1) = 10; % Port C stim train period in ms.
stimP(4,1,1) = 10; % Port D stim train period in ms.

% Train number of pulses
stimNum(1,1,1) = 4; % Port A stim train pulse number.
stimNum(2,1,1) = 4; % Port B stim train pulse number.
stimNum(3,1,1) = 4; % Port C stim train pulse number.
stimNum(4,1,1) = 4; % Port D stim train pulse number.

% Stim 2 parameters
ISI{2,1} = 10; % Inter stim interval in seconds, randomly selected.

% Amplitudes. [day1,day2...]
stimI(1,:,2) = [250,250,250,250]; % Port A stim amp in uA.
stimI(2,:,2) = [400,400,400,400]; % Port B stim amp in uA.
stimI(3,:,2) = [0,0,0,0]; % Port C stim amp in uA.
stimI(4,:,2) = [0,0,0,0]; % Port D stim amp in uA.

% Channel, [0,0] = unused, [c1,0] = monopolar, [c1,c2] = bipolar.
stimChan(1,:,2) = [3,0]; % Port A
stimChan(2,:,2) = [5,0]; % Port B
stimChan(3,:,2) = [0,0]; % Port C
stimChan(4,:,2) = [0,0]; % Port D

% Train period
stimP(1,1,2) = 100; % Port A stim train period in ms.
stimP(2,1,2) = 100; % Port B stim train period in ms.
stimP(3,1,2) = 100; % Port C stim train period in ms.
stimP(4,1,2) = 100; % Port D stim train period in ms.

% Train number of pulses
stimNum(1,1,2) = 100; % Port A stim train pulse number.
stimNum(2,1,2) = 100; % Port B stim train pulse number.
stimNum(3,1,2) = 100; % Port C stim train pulse number.
stimNum(4,1,2) = 100; % Port D stim train pulse number.

% Input channel for closed-loop, 0 = unused.
dataChan(1,1) = 0; % Port A
dataChan(2,1) = 0; % Port B
dataChan(3,1) = 0; % Port C
dataChan(4,1) = 0; % Port D

framesPerBlock = 128;
blocksPerRead = 125;
numBandsPerChan = 1; % Amplifier channels can be displayed as LOW, WIDE, or HIGH

% Analysis Settings
runAnalysis = 1;
sF = 20000; % Recording sample frequency
dSF = 1000; % Downsample frequency
stimTrig = 1;
use60HzFilt = 0;
ACfiltPath = 'G:\RecordingData\EDS\CNN60Hz.mat';
replaceFile = 1;

% Active channels for RHS
activeChan(1,1:3,1) = [{'a-000'},{'PrL(L)'},{1}];
activeChan(2,1:3,1) = [{'a-001'},{'PrL(R)'},{1}];
activeChan(3,1:3,1) = [{'a-002'},{'AVT(L)'},{1}];
activeChan(4,1:3,1) = [{'a-003'},{'BLA(R)'},{0}];
activeChan(5,1:3,1) = [{'a-004'},{'CA1(L)'},{1}];
activeChan(6,1:3,1) = [{'a-005'},{'CA1(R)'},{1}];
activeChan(7,1:3,1) = [{'a-006'},{'LDT(L1)'},{1}];
activeChan(8,1:3,1) = [{'a-007'},{'LDT(L2)'},{1}];

activeChan(1,1:3,2) = [{'b-000'},{'PrL(L)'},{1}];
activeChan(2,1:3,2) = [{'b-001'},{'PrL(R)'},{1}];
activeChan(3,1:3,2) = [{'b-002'},{'AVT(L)'},{1}];
activeChan(4,1:3,2) = [{'b-003'},{'BLA(R)'},{1}];
activeChan(5,1:3,2) = [{'b-004'},{'CA1(L)'},{1}];
activeChan(6,1:3,2) = [{'b-005'},{'CA1(R)'},{1}];
activeChan(7,1:3,2) = [{'b-006'},{'LDT(L1)'},{1}];
activeChan(8,1:3,2) = [{'b-007'},{'LDT(L2)'},{1}];

activeChan(1,1:3,3) = [{'c-000'},{'PrL(L)'},{0}];
activeChan(2,1:3,3) = [{'c-001'},{'PrL(R)'},{0}];
activeChan(3,1:3,3) = [{'c-002'},{'AVT(L)'},{0}];
activeChan(4,1:3,3) = [{'c-003'},{'BLA(R)'},{0}];
activeChan(5,1:3,3) = [{'c-004'},{'CA1(L)'},{0}];
activeChan(6,1:3,3) = [{'c-005'},{'CA1(R)'},{0}];
activeChan(7,1:3,3) = [{'c-006'},{'LDT(L1)'},{0}];
activeChan(8,1:3,3) = [{'c-007'},{'LDT(L2)'},{0}];

activeChan(1,1:3,4) = [{'d-000'},{'PrL(L)'},{0}];
activeChan(2,1:3,4) = [{'d-001'},{'PrL(R)'},{0}];
activeChan(3,1:3,4) = [{'d-002'},{'AVT(L)'},{0}];
activeChan(4,1:3,4) = [{'d-003'},{'BLA(R)'},{0}];
activeChan(5,1:3,4) = [{'d-004'},{'CA1(L)'},{0}];
activeChan(6,1:3,4) = [{'d-005'},{'CA1(R)'},{0}];
activeChan(7,1:3,4) = [{'d-006'},{'LDT(L1)'},{0}];
activeChan(8,1:3,4) = [{'d-007'},{'LDT(L2)'},{0}];

WL = webcamlist;
for i = 1:length(WL)
    
    cam = webcam(i);
    sp = rgb2gray(snapshot(cam));
    sp(:,:,2:5) = 0;
    for ii = 2:5
        sp(:,:,ii) = rgb2gray(snapshot(cam));
    end
    if mean(sp,'all') <= 50
        WL(i,2) = {'N/A'};
        continue
    end
    preview(cam)
    prompt = 'Label WebCam';
    dlgtitle = 'Input';
    WL{i,2} = str2double(inputdlg(prompt,dlgtitle,[1,40],WL(i,1)));
end

for i = 1:size(WL,1)

    if strcmp(WL{i,2},'N/A') == 0
        WL{i,3} = char(boxSubject{WL{i,2},1}+" ("+WL{i,2}+")");
    else
        WL{i,3} = 'N/A';
    end
end

[camChoice,~] = listdlg('PromptString','Select cameras.','ListString',WL(:,3),'SelectionMode','multiple');
numCam = length(camChoice);

Res = cam.AvailableResolutions;
[resChoice,~] = listdlg('PromptString','Select resolution.','ListString',Res,'SelectionMode','single');
Dim = symsepchar(Res{1,resChoice},'x');
Dim = [str2double(Dim{1,2}),str2double(Dim{1,1})];
clear cam

PathName = uigetdir(cd,'Select folder for Datastore');
cd(PathName);

promptMessage = sprintf('Would you like to copy to another directory?');
titleBarCaption = 'settings';
copyAction = questdlg(promptMessage, titleBarCaption, 'Yes','No','Yes');

if strcmp(copyAction,'Yes') == 1
    copyPathName = uigetdir(cd,'Select directory');
end

save([PathName,'\config.mat'])
%%

global dataFig videoFig %#ok<NUSED,GVMIS>

vStartTime = (StartTime(1,1)*3600)+(StartTime(1,2)*60);

port = {'a','b','c','d'};

hardwareTimer = 1;
intanWriter = 2;
intanReader = 3;
dataHandler = 4;
camReader = 5;
camSaver = 6;

% Calculations for accurate parsing
numChan = length(find(dataChan > 0));
chanStim = 1:4;
chanStim(1,:,2) = 1:4;
chanStim(2,sum(stimChan,2) > 0) = 1;
numChanStim = length(chanStim);
numDsFrames = (framesPerBlock*blocksPerRead)/(sF/dSF);
readTime = numDsFrames/dSF;
numAmplifierBands = numBandsPerChan*numChan;
waveformBytesPerFrame = 4+2*numAmplifierBands;
waveformBytesPerBlock = framesPerBlock*waveformBytesPerFrame+4;
framesPerRead = framesPerBlock*blocksPerRead;
waveformBytesBlocks = blocksPerRead*waveformBytesPerBlock;

textOut = parallel.pool.DataQueue;
afterEach(textOut,@disp);

videoOut = parallel.pool.DataQueue;
afterEach(videoOut,@showVideo);

dataOut = parallel.pool.DataQueue;
afterEach(dataOut,@plotTD);

for Iter = startIter:size(stimI,2)

    disp(['Will start at ',num2str(StartTime(1,1)),':',num2str(StartTime(1,2))]); 
    while 1
        if nextDay == 0
            if getTime >= vStartTime-90
                nextDay = 1;
                break
            end
        elseif nextDay == 1
            if getTime <= vStartTime
                nextDay = 0;
            end
        end
        pause(0.01);
    end

    Now = clock;
    Time = [num2str(Now(1,2)),'-',num2str(Now(1,3)),'-',num2str(Now(1,1)),'(',num2str(Now(1,4)),'.',num2str(Now(1,5)),')'];
    
    seshPathName = [PathName,'\',Time];
    eval(['mkdir ',seshPathName]);
    
    if isempty(gcp('nocreate')) == 1
        parpool('local',6);
        disp('Connecting to hardware...');        
    end 
    
    spmd (6)
        if labindex == hardwareTimer

            noCon_Ard = 0;
            if Iter == startIter

                t = 0;                
                while 1
                    try
                        Ard = arduino('com3','uno');
                        send(textOut,'Connected to Arduino')
                        break
                    catch
                        send(textOut,'Issue connecting to Arduino')
                        pause(1)
                        t = t+1;
                        if t == 60
                            send(textOut,'Connection to Arduino timed out')
                            noCon_Ard = 1;
                            break
                        end
                    end                    
                end
            end   

            labSend(noCon_Ard,[intanWriter,intanReader,dataHandler,camReader,camSaver],0);
            noCon_TCPwrite = labReceive(intanWriter,0);
            noCon_TCPread = labReceive(intanReader,0);

            if noCon_Ard == 0 && noCon_TCPwrite == 0 && noCon_TCPread == 0
    
                % sync with system clock
                send(textOut,'Syncing hardware');
                while 1
                    if getTime >= vStartTime
                        break
                    end
                end            
                
                labBarrier

                seshStart = getTime;
                labSend(seshStart,[intanWriter,intanReader,dataHandler,camReader],1);
    
                Count = 1;
                send(textOut,'----------')
                for i = 1:size(blockSched,1)
                    for ii = 1:size(blockSched,2)
                        if strcmp(blockSched{i,ii},'N/A') == 0
                            if strcmp(blockSched{i,ii}{Iter,1},'BL') == 1
                                send(textOut,['Box ',num2str(boxRHS(ii,1)),' Baseline']) 
                            elseif strcmp(blockSched{i,ii}{Iter,1},'S1') == 1
                                send(textOut,['Box ',num2str(boxRHS(ii,1)),' Stim 1']) 
                            elseif strcmp(blockSched{i,ii}{Iter,1},'S2') == 1
                                send(textOut,['Box ',num2str(boxRHS(ii,1)),' Stim 2']) 
                            elseif strcmp(blockSched{i,ii}{Iter,1},'SX') == 1
                                send(textOut,['Box ',num2str(boxRHS(ii,1)),' Sham Stim']) 
                            end
                        end                        
                    end
                    send(textOut,'----------')

                    % Turn on block marker
                    writeDigitalPin(Ard,ArdPin{2,1},1);
                    Ck = 1;
                    while getTime < seshStart+sum(blockTime(1:i,1))       
                        
                        % Turn on LED
                        if Count == 1                  
                            writeDigitalPin(Ard,ArdPin{1,1},1);                        
                        end 
    
                        pause(0.2)
    
                        % Turn off LED
                        if Count == 1    
                            writeDigitalPin(Ard,ArdPin{1,1},0);
                        end
                        
                        % Turn off block marker
                        if Ck == 1
                             writeDigitalPin(Ard,ArdPin{2,1},0);
                             Ck = 0;
                        end                  
                        
                        Count = Count+1;
    
                        if Count == LEDtime-1
                            Count = 1;
                        end
    
                        t = getTime;
                        % Wait for next second
                        while getTime-t <= 0.999
                            pause(0.001)
                        end                    
                    end
                end
            end
        elseif labindex == intanWriter
            
            noCon_TCPwrite = 0;
            if useRHS+useRHD > 0
                if Iter == startIter

                    intanSend = cell(0,1);
                    check = [0,0];
                    t = 0;                    
                    while true
                        try
                            if useRHS == 1 && check(1,1) == 0
                                intanSend{1,1} = tcpclient(intanIP_RHS,intanPort1_RHS);
                                check(1,1) = 1;
                                send(textOut,'Connected to Intan RHS commands')
                            end
                            if useRHD == 1 && check(1,2) == 0
                                intanSend{2,1} = tcpclient(intanIP_RHD,intanPort1_RHD);
                                check(1,2) = 1;
                                send(textOut,'Connected to Intan RHD commands')
                            end
                        catch
                            if useRHS == 1 && check(1,1) == 0
                                send(textOut,'Issue connecting to Intan RHS commands')
                            end
                            if useRHD == 1 && check(1,2) == 0
                                send(textOut,'Issue connecting to Intan RHD commands')
                            end

                            pause(1)
                            t = t+1;
                            if t == 60
                                send(textOut,'Connection to Intan timed out')
                                noCon_TCPwrite = 1;
                                break
                            end
                        end
                        if ((useRHS == 1) ~= (useRHD == 1) && sum(check) == 1) || sum(check) == 2
                            break
                        end
                    end
                end

                if check(1,1) == 1                    
                    write(intanSend{1,1},uint8(['set Filename.Path ',seshPathName]));
                    write(intanSend{1,1},uint8('set FileFormat OneFilePerSignalType'));
                    write(intanSend{1,1},uint8('set Filename.BaseFilename RHS'));
                end

                if check(1,2) == 1                    
                    write(intanSend{2,1},uint8(['set Filename.Path ',seshPathName]));
                    write(intanSend{2,1},uint8('set FileFormat OneFilePerSignalType'));
                    write(intanSend{2,1},uint8('set Filename.BaseFilename RHD'));
                end
                
                % Set up stim parameters
                if check(1,1) == 1    

                    Count = 1;
                    for i = 1:2
                        for ii = 1:size(stimChan,1)
                            if stimChan(ii,1,i) > 0
            
                                if dataChan(ii,1) > 0 && strcmp(loopType,'closed') == 1
                                    write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{dataChan(ii,1),1},'.tcpdataoutputenabled true;']));
                                end
                                
                                for i3 = 1:size(stimChan(ii,:,i),2)
                                    if stimChan(ii,i3,i) > 0
                        
                                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.stimenabled true']));
                                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.source KeyPressF',num2str(Count)]));
                                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.FirstPhaseAmplitudeMicroAmps ',num2str(stimI(ii,Iter,i))]));
                                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.SecondPhaseAmplitudeMicroAmps ',num2str(stimI(ii,Iter,i))]));
                                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.NumberOfStimPulses ',num2str(stimNum(ii,1,i))]));
                                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.PulseTrainPeriodMicroseconds ',num2str(stimP(ii,1,i)*1000)]));
                                        
                                        if i3 == 1
                                            write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.Polarity NegativeFirst']));
                                        elseif i3 == 2
                                            write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.Polarity PositiveFirst']));
                                        end
                                    end
                                end
                                send(textOut,['Stim ',num2str(i),' port',num2str(ii),' amp set to ',num2str(stimI(ii,Iter,i)),'uA.'])
                            end
                            Count = Count+1;
                        end
                    end                 
                    
                    if useRHS == 1 && sum(stimI(:,Iter)) > 0
                        write(intanSend{1,1},uint8('execute uploadstimparameters'));
                        pause(10)
                    end
                end                
            end            
            
            labSend(noCon_TCPwrite,[hardwareTimer,intanReader,dataHandler,camReader,camSaver],0);
            noCon_Ard = labReceive(hardwareTimer,0);
            noCon_TCPread = labReceive(intanReader,0);

            labBarrier

            seshStart = labReceive(hardwareTimer,1);
            if noCon_Ard == 0 && noCon_TCPwrite == 0 && noCon_TCPread == 0      
                
                if useRHS == 1
                    write(intanSend{1,1},uint8('set runMode record'));
                end
                if useRHD == 1
                    write(intanSend{2,1},uint8('set runMode record'));
                end
                
                while getTime < seshStart+sum(blockTime(1:end,1))
                    
                    stim = labReceive(dataHandler);
                    if isempty(stim) == 1
                        break
                    end
                    
                    Count = 1;
                    for ii = 1:size(chanStim,3)
                        for i3 = 1:size(chanStim,2)
                            if stim(ii,i3) == 1                            
                                write(intanSend{1,1},uint8(['execute ManualStimTriggerPulse F',num2str(Count)]));
                            end
                            Count = Count+1;
                        end
                    end

                    if labProbe(camSaver) == 1
                        labReceive(camSaver);

                        if useRHS == 1
                            write(intanSend{1,1},uint8('set runMode stop'));
                        end
                        if useRHD == 1
                            write(intanSend{2,1},uint8('set runMode stop'));
                        end

                        if check(1,1) == 1                    
                            write(intanSend{1,1},uint8(['set Filename.Path ',seshPathName]));
                            write(intanSend{1,1},uint8('set FileFormat OneFilePerSignalType'));
                            write(intanSend{1,1},uint8('set Filename.BaseFilename RHS_Dark'));
                        end
                        if check(1,2) == 1                    
                            write(intanSend{2,1},uint8(['set Filename.Path ',seshPathName]));
                            write(intanSend{2,1},uint8('set FileFormat OneFilePerSignalType'));
                            write(intanSend{2,1},uint8('set Filename.BaseFilename RHD_Dark'));
                        end

                        if useRHS == 1
                            write(intanSend{1,1},uint8('set runMode record'));
                        end
                        if useRHD == 1
                            write(intanSend{2,1},uint8('set runMode record'));
                        end
                    end
                end
                
                if useRHS == 1
                    write(intanSend{1,1},uint8('set runMode stop'));
                end
                if useRHD == 1
                    write(intanSend{2,1},uint8('set runMode stop'));
                end
            end
            
            if useRHS > 0
                if check(1,1) == 1
                    for i = 1:2
                        for ii = 1:size(stimChan,1)
                            if stimChan(ii,1,i) > 0
            
                                if dataChan(ii,1) > 0 && strcmp(loopType,'closed') == 1
                                    write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{dataChan(ii,1),1},'.tcpdataoutputenabled false;']));
                                end
                                
                                for i3 = 1:size(stimChan(ii,:,i),2)
                                    if stimChan(ii,i3,i) > 0                        
                                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3,i),1},'.stimenabled false']));                                        
                                    end
                                end
                            end
                        end
                    end                 
                    
                    if useRHS == 1 && sum(stimI(:,Iter)) > 0
                        write(intanSend{1,1},uint8('execute uploadstimparameters'));
                    end
                end
            end
            
        elseif labindex == intanReader
    
            noCon_TCPread = 0;
            if strcmp(loopType,'closed') == 1

                t = 0;                
                while 1
                    try
                        intanRead = tcpclient(intanIP_RHS,intanPort2_RHS);
                        send(textOut,'Connected to Intan RHS data output')
                        break
                    catch
                        disp('Issue connecting to Intan data output')
                        pause(1)
                        t = t+1;
                        if t == 60
                            send(textOut,'Connection to Intan data output timed out')
                            noCon_TCPread = 1;
                            break
                        end
                    end
                end
        
                flush(intanRead);
        
                % Pre-allocate memory for blocks of waveform data (the amount that's
                % plotted at once)
                amplifierData = 32768*ones(numChan,framesPerRead);
                ampTimewindow = zeros(numChan,numDsFrames*5);
            
                % Initialize counters
                chunkCount = 0;
                blockCount = 0;
                amplifierTimestampsIndex = 1;
                Stopper = 0;
            else
                intanRead = [];
            end

            labSend(noCon_TCPread,[hardwareTimer,intanWriter,dataHandler,camReader,camSaver],0);
            noCon_Ard = labReceive(hardwareTimer,0);
            noCon_TCPwrite = labReceive(intanWriter,0);   

            labBarrier

            seshStart = labReceive(hardwareTimer,1); 
            if noCon_Ard == 0 && noCon_TCPwrite == 0 && noCon_TCPread == 0
                if strcmp(loopType,'closed') == 1
                    for i = 1:size(blockSched,1)
                        
                        % Check if any blocks are S1
                        seshStim = zeros(1,size(blockSched,2));
                        for ii = 1:size(blockSched,2)
                            if ~all(strcmp(blockSched{i,ii},'N/A'))
                                if strcmp(blockSched{i,ii}{Iter,1},'S1') == 1
                                    seshStim(1,ii) = 1;
                                end
                            end
                        end
    
                        while getTime < seshStart+sum(blockTime(1:i,1))
                            if all(seshStim == 0)
                                pause(0.001)
                                flush(intanRead);
                            else    
    
                                HaveByte = intanRead.BytesAvailable;
                            
                                %Read waveform data in block chunks
                                if any(HaveByte >= waveformBytesBlocks)
                                    
                                    waveformArray = read(intanRead,waveformBytesBlocks);
                                    rawIndex = 1;
                    
                                    %Read all incoming blocks                   
                                    if HaveByte >= waveformBytesBlocks
                                        for block = 1:blocksPerRead
                    
                                            blockCount = blockCount+1;
                                            
                                            %Start here, blocksPerRead not being and integer is probly the issue
                                            %Expect 4 bytes to be TCP Magic Number as uint32. If not what's expected, print that there was an error.
                                            Bytes = waveformArray(1,rawIndex:rawIndex+3);
                                            magicNumber = typecast(uint8(Bytes),'uint32');
                                            rawIndex = rawIndex+4;  
                    
                                            if magicNumber ~= 0x2ef07a08
                                                fprintf(1,'Error... block %d magic number incorrect.\n',block);
                                            end
                    
                                            %Each block should contain 128 frames of data - process each of these one-by-one
                                            for frame = 1:framesPerBlock
                    
                                                %Expect 4 bytes to be timestamp as int32
                                                Bytes = waveformArray(1,rawIndex:rawIndex+3);
                                                rawIndex = rawIndex+4;
                    
                                                %Parse all bands of amplifier channels
                                                for channel = 1:numChan
                    
                                                    %2 bytes of wide            
                                                    Bytes = waveformArray(rawIndex:rawIndex+1);
                                                    amplifierData(channel,amplifierTimestampsIndex) = typecast(uint8(Bytes),'uint16');
                                                    rawIndex = rawIndex+2;
                                                end
                                                amplifierTimestampsIndex = amplifierTimestampsIndex+1;
                                            end
                                        end
                                    end               
                    
                                    % Scale these data blocks and downsample
                                    amplifierData = downsample(0.195*(amplifierData'-32768),sF/dSF)';
                    
                                    % Shift time window to next step
                                    ampTimewindow(:,1:numDsFrames) = [];
                                    ampTimewindow = [ampTimewindow,amplifierData]; %#ok<AGROW>
                    
                                    labSend(ampTimewindow,dataHandler);          
                                    
                                    % Reset index
                                    blockCount = 0;
                                    amplifierTimestampsIndex = 1;
                                    Stopper = 0;
                                else
                                    Stopper = Stopper+1;
                                    pause(0.001)
                                end
                                if Stopper > 5000
                                    break
                                end
                            end
                        end
                        if any(seshStim == 1)
                            labSend([],dataHandler);
                        end
                    end
                end
        
                while getTime < seshStart+sum(blockTime(1:end,1))
                    if strcmp(loopType,'closed') == 1
                        flush(intanRead);
                    end
                    pause(0.001)
                end   
            end
        elseif labindex == dataHandler     
            
            noCon_Ard = labReceive(hardwareTimer,0);
            noCon_TCPwrite = labReceive(intanWriter,0);
            noCon_TCPread = labReceive(intanReader,0);

            labBarrier

            seshStart = labReceive(hardwareTimer,1);
            if noCon_Ard == 0 && noCon_TCPwrite == 0 && noCon_TCPread == 0
                for i = 1:size(blockSched,1)  
    
                    % Check if any blocks are S1 or S2
                    seshStim = zeros(1,size(blockSched,2));
                    for ii = 1:size(blockSched,2)
                        if ~all(strcmp(blockSched{i,ii},'N/A'))
                            if strcmp(blockSched{i,ii}{Iter,1},'S1') == 1
                                seshStim(1,ii) = 1;
                            elseif strcmp(blockSched{i,ii}{Iter,1},'S2') == 1
                                seshStim(1,ii) = 2;
                            end
                        end
                    end
    
                    if strcmp(loopType,'closed') == 1
        
                        bankISI = (round(ISI{1,1}(1,1)/readTime)):round((ISI{1,1}(1,end)/readTime));
                        nextISI(1:numChan,1) = bankISI(1,1);
                        readCount(1:numChan,1) = 0;
                        theta_delta = zeros(numChan,1);
                    else
                        stimType = {'S1'};
                        stimType(2,1) = {'S2'};
                        stim = zeros(2,4);
                        Count = nan(2,4);
                        Count(1,chanStim(2,:,1) == 1) = 0;
                        Count(2,chanStim(2,:,2) == 1) = 0;
                        nextISI = nan(2,4);
                        nextISI(1,chanStim(2,:,1) == 1) = 0;
                        nextISI(2,chanStim(2,:,2) == 1) = 0;
                        numStim = 0;
                    end
    
                    while getTime < seshStart+sum(blockTime(1:i,1))
                        if sum(seshStim) == 0
                            pause(0.001)
                        elseif sum(seshStim) > 0
                            if any(seshStim == 1) && strcmp(loopType,'closed') == 1
                
                                ampTimewindow = labReceive(intanReader);
                                if isempty(ampTimewindow) == 1
                                    break
                                end
                                %send(dataOut,[{ampTimewindow};{dSF}]); 
    
                                readCount = readCount+1;                
                                if any(readCount >= nextISI)
                    
                                    % Calculate theta/delta ratio
                                    theta = bandpass(ampTimewindow',[5,12],dSF)';
                                    delta = bandpass(ampTimewindow',[1,4],dSF)';
                                    thetaAmp = abs(hilbert(theta')');
                                    deltaAmp = abs(hilbert(delta')');
                                    for ii = 1:numChan
                                        theta_delta(ii,1) = mean(thetaAmp(ii,end-dSF+1:end))/mean(deltaAmp(ii,end-dSF+1:end));
                                    end
                                    stim = theta_delta < thresholdTD;
                                    stim = stim.*(readCount >= nextISI);
                
                                    labSend(stim,intanWriter);
                                    for ii = 1:numChan
                                        if stim(ii,1) == 1 
                                            
                                            send(textOut,['Stim Ch',num2str(ii),' @ ',num2str(theta_delta(ii,1)),' ISI: ',num2str(readCount(ii,1))]);
                                            nextISI(ii,1) = datasample(bankISI,1);
                                            readCount(ii,1) = 0;
                                        end
                                    end
                                end
                            else
            
                                t = getTime;
                                while getTime-t <= 0.999
                                    % Wait for next second
                                end
    
                                for ii = 1:2
                                    for i3 = chanStim(1,chanStim(2,:,ii) == 1,ii)
                                        if all(strcmp(blockSched{i,i3},'N/A') == 0)
                                            if strcmp(blockSched{i,i3}{Iter,1},stimType{ii,1}) == 1
                                                if Count(ii,i3) == 0
                                                    stim(ii,i3) = 1;
                                                    nextISI(ii,i3) = datasample(ISI{ii,1},1);
                                                    Count(ii,i3) = Count(ii,i3)+1;
                                                elseif Count(ii,i3) == nextISI(ii,i3)
                                                    Count(ii,i3) = 0;
                                                else
                                                    stim(ii,i3) = 0;
                                                    Count(ii,i3) = Count(ii,i3)+1;
                                                end
                                            end
                                        end
                                    end
                                end                            
    
                                if any(stim == 1,'all')
                                    labSend(stim,intanWriter);
                                end
                            end
                        end
                    end
                end
                labSend([],intanWriter);
            end
        elseif labindex == camReader 
            
            % Create all webcam objects           
            cam = cell(numCam,1); 
            for ii = 1:numCam                
                cam{ii,1} = webcam(camChoice(1,ii));
            end           
            frames = uint8(zeros(Dim(1,1),Dim(1,2),3,40,numCam));     

            errorCode = uint8([0,255,0;255,0,255;0,255,0]);
            
            noCon_Ard = labReceive(hardwareTimer,0);
            noCon_TCPwrite = labReceive(intanWriter,0);
            noCon_TCPread = labReceive(intanReader,0);

            labBarrier

            seshStart = labReceive(hardwareTimer,1);
            if noCon_Ard == 0 && noCon_TCPwrite == 0 && noCon_TCPread == 0  

                t = getTime;
                ck = zeros(1,numCam);
                while t < seshStart+sum(blockTime(1:end,1))

                    Count = 1;
                    t = getTime;
                    while getTime-t <= 0.95
                        for ii = 1:numCam
                            try
                                frames(:,:,:,Count,ii) = snapshot(cam{ii,1});
                            catch
                                if ck(1,ii) == 0
                                    send(textOut,['Could not get frame from cam ',num2str(ii)])
                                    send(textOut,'Atempting to reconnect...')
                                    cam{ii,1} = webcam(camChoice(1,ii));
                                    ck(1,ii) = 1;
                                elseif ck(1,ii) > 0 && ck(1,ii) <= 5
                                    send(textOut,'Atempting to reconnect...')
                                    cam{ii,1} = webcam(camChoice(1,ii));
                                    ck(1,ii) = ck(1,ii)+1;
                                elseif ck(1,ii) == 6
                                    send(textOut,'Connection to cam ',num2str(ii),' lost.')
                                    ck(1,ii) = ck(1,ii)+1;
                                end
                                frames(1:3,1:3,1,Count,ii) = errorCode;                                
                            end
                        end
                        Count = Count+1;
                    end
        
                    labSend({frames,Count},camSaver);
                end 
                labSend([],camSaver);
            end
        elseif labindex == camSaver       
            
            VidFileName = cell(numCam,1);
            writerObj = cell(numCam,1);
    
            % Create all video writer objects
            for ii = 1:numCam   
                       
                VidFileName{ii,1} = [WL{camChoice(1,ii),3}(1,1:end-4),'.mp4'];
                writerObj{ii,1} = VideoWriter([seshPathName,'\',VidFileName{ii,1}],'MPEG-4'); %#ok<TNMLP>
                writerObj{ii,1}.FrameRate = fps;
                open(writerObj{ii,1});
            end   

            errorCode = uint8([0,255,0;255,0,255;0,255,0]);
            
            noCon_Ard = labReceive(hardwareTimer,0);
            noCon_TCPwrite = labReceive(intanWriter,0);
            noCon_TCPread = labReceive(intanReader,0);
            
            labBarrier

            if noCon_Ard == 0 && noCon_TCPwrite == 0 && noCon_TCPread == 0

                ck = zeros(1,numCam);
                while 1
        
                    frames = labReceive(camReader);
                    
                    if isempty(frames)
                        break
                    end 
    
                    Count = frames{1,2}-1;
                    frames = frames{1,1};
    
                    idx = round(linspace(1,Count,fps));
                    frames = frames(:,:,:,idx,:);
    
                    if ShowVideo == 1
                        send(videoOut,frames);
                    end
    
                    for ii = 1:numCam
                        for iii = 1:fps
                            if ~all(frames(1:3,1:3,1,iii,ii) == errorCode,'all')
                                if mean(sum(frames(:,:,:,iii,ii),3),'all') > 250
                                    writeVideo(writerObj{ii,1},frames(:,:,:,iii,ii));
                                    ck = zeros(1,numCam);
                                else                                    
                                    if all(ck < 5)
                                        ck(1,ii) = ck(1,ii)+1;
                                    else
                                        labSend(0,intanWriter);
                                    end
                                end
                            end
                        end
                    end
                end
                for ii = 1:numCam
                    close(writerObj{ii,1});
                end
            end
        end
    end
    close(gcf)

    if runAnalysis == 1

        disp('Starting analysis');       

        listing = returnSubfolder(seshPathName);
        numFolders = size(listing,1);

        recInfo = struct;
        for i = 1:numFolders

            recInfo(i).file = [listing(i).folder,'\',listing(i).name];

            if contains(listing(i).name,'RHS') == 1
                recInfo(i).system = 'RHS';
            elseif contains(listing(i).name,'RHD') == 1
                recInfo(i).system = 'RHD';
            end

            recInfo(i).fileData(1).port = 'A';
            recInfo(i).fileData(2).port = 'B';
            recInfo(i).fileData(3).port = 'C';
            recInfo(i).fileData(4).port = 'D';
        
            if contains(listing(i).name,'RHS') == 1
        
                for ii = 1:4
                    if boxRHS(ii,1) > 0
                        
                        timestamps = cell(length(blockTime),2);
                        x = 1;
                        y = 0;
                        for i3 = 1:length(blockTime)

                            timestamps{i3,1} = blockSched{i3,1}{Iter,1};
                            y = y+blockTime(i3,1);
                            timestamps{i3,2} = [x,y];
                            x = x+blockTime(i3,1);
                        end
                        recInfo(i).fileData(ii).blockDuration = timestamps;

                        recInfo(i).fileData(ii).subject = boxSubject(boxRHS(ii,1),1);

                        recInfo(i).fileData(ii).stim1_uA = stimI(ii,Iter,1);
                        recInfo(i).fileData(ii).stim1_chan = stimChan(ii,:,1);
                        if any(stimChan(ii,:,1) == 0)
                            recInfo(i).fileData(ii).stim1_type = 'monopole';
                        else
                            recInfo(i).fileData(ii).stim1_type = 'bipole';
                        end

                        recInfo(i).fileData(ii).stim2_uA = stimI(ii,Iter,2);
                        recInfo(i).fileData(ii).stim2_chan = stimChan(ii,:,2);
                        if any(stimChan(ii,:,2) == 0)
                            recInfo(i).fileData(ii).stim2_type = 'monopole';
                        else
                            recInfo(i).fileData(ii).stim2_type = 'bipole';
                        end                   
                    end
                end
            elseif contains(listing(i).name,'RHD') == 1
                for ii = 1:4
                    if boxRHD(ii,1) > 0
                        recInfo(i).fileData(ii).subject = boxSubject(boxRHD(ii,1),1);
                    end
                end
            end    
        end
        
        for i = 1:numFolders
        
            if strcmp(recInfo(i).system,'RHD') == 1
                continue
            end
        
            cd(recInfo(i).file)
        
            % Digital in
            disp('Working on digitalin.dat')
            s = dir([recInfo(i).file,'\digitalin.dat']);
            numSamples = s.bytes/2;
            
            fileID = fopen([recInfo(i).file,'\digitalin.dat']);
    
            fseek(fileID,0,'bof');
            L = floor(numSamples/(sF*600));
        
            if replaceFile == 1
                fID = fopen(recInfo(i).file+"\digitalin_1Ch_"+dSF/1000+"KHz_.dat",'w');
            end
        
            for ii = 1:L
                
                temp = fread(fileID,[1,sF*600],'int16');
                temp = downsample(temp,sF/dSF);
                if replaceFile == 1
                    fwrite(fID,temp,'int16');
                end
            end
        
            if replaceFile == 1
                fclose(fileID);
                fclose(fID);
                %eval(['delete ',recInfo{i,1},'\digitalin.dat']);
            end
        
            % Find number of active channels
            numChan = length(find(cell2mat(activeChan(:,3,:)) == 1));
            chan = cell(numChan,2);
            count = 1;
            for ii = 1:size(activeChan,3)
                for i3 = 1:size(activeChan,1)
                    if activeChan{i3,3,ii} == 1
                        chan(count,:) = activeChan(i3,1:2,ii);
                        count = count+1;
                    end
                end
            end
        
            % Stim TS
            disp('Working on stim.dat')
            s = dir([recInfo(i).file,'\stim.dat']);
            numSamples = s.bytes/(2*numChan);   
            
            fileID = fopen([recInfo(i).file,'\stim.dat']);
            fseek(fileID,0,'bof');
            L = floor(numSamples/(sF*600));
            stimData = nan(numChan,L*600*dSF);
        
            if replaceFile == 1
                fID = fopen([recInfo(i).file,'\stim_',num2str(numChan),'Ch_1KHz_.dat'],'w');
            end
            
            wb = waitbar(0,'Working on stim.dat');
            Count = 1;
            for ii = 1:L
                
                temp = fread(fileID,[numChan,sF*600],'int16');
                temp = downsample(temp',sF/dSF)';
                if replaceFile == 1
                    fwrite(fID,temp,'int16');
                end
                temp(temp > 0) = 1;
                temp(temp < 0) = 0;
                stimData(:,Count:Count+size(temp,2)-1) = temp;
                Count = Count+size(temp,2);
                waitbar(ii/L,wb,'Working on stim.dat');
            end   
            close(wb);
            
            fclose(fileID);
            if replaceFile == 1
                fclose(fID);
                %eval(['delete ',recInfo{i,1},'\stim.dat']);
            end
            
            stimTS = cell(size(activeChan,3),2);
            count = 1;
            for ii = 1:size(activeChan,3)
                for i3 = 1:size(activeChan,1)
                    if ~isempty(recInfo(i).fileData(ii).subject) && activeChan{i3,3,ii} == 1
                        for i4 = 1:2
                            if any(stimChan(ii,:,i4) == i3)
                                temp = stimData(count,:);
                                temp = find(diff(temp) > 0);
                                stimTS{ii,i4} = temp';
                            end
                            if i4 == 1
                                recInfo(i).fileData(ii).stim1_timeStamps = stimTS{ii,i4};
                            elseif i4 == 2
                                recInfo(i).fileData(ii).stim2_timeStamps = stimTS{ii,i4};
                            end
                        end                        
                        count = count+1;
                    end
                end
            end
            clear stimData
            
            % Amplifier in
            disp('Working on amplifier.dat')
            s = dir([recInfo(i).file,'\amplifier.dat']);
            numSamples = s.bytes/(2*numChan);   
            
            fileID = fopen([recInfo(i).file,'\amplifier.dat']);
            fseek(fileID,0,'bof');
            L = floor(numSamples/(sF*600));
            data = nan(numChan,L*600*dSF);
            thetaAmp = nan(numChan,L*600*dSF);
            deltaAmp = nan(numChan,L*600*dSF);
        
            if replaceFile == 1
                fID = fopen([recInfo(i).file,'\amplifier_',num2str(numChan),'Ch_1KHz.dat'],'w');
                if use60HzFilt == 1
                    load(ACfiltPath);
                    fID2 = fopen([recInfo(i).file,'\amplifier_',num2str(numChan),'Ch_1KHz_Filt.dat'],'w');
                end
            end
            
            wb = waitbar(0,'Working on amplifier.dat');
            Count = 1;
            for ii = 1:L
                
                temp = fread(fileID,[numChan,sF*600],'int16');
                temp = convLPFilt(temp',dSF/2,sF,[],sF/dSF)';

                %tempDS = zeros(size(temp,1),size(temp,2)/(sF/dSF));
                %parfor i3 = 1:size(temp,1)
                %    tempDS(i3,:) = decimate(temp(i3,:),sF/dSF);
                %end
                %temp = tempDS;

                %temp = lowpass(temp',dSF/2,sF)';
                %temp = downsample(temp',sF/dSF)';                

                if replaceFile == 1
                    fwrite(fID,temp,'int16');
                    if use60HzFilt == 1
                        for i3 = 1:numChan
                            temp(i3,:) = predict(CNN60Hz,temp(i3,:)/50000,'ExecutionEnvironment','gpu')*50000;
                        end
                        fwrite(fID2,temp,'int16');
                    end
                end
                data(:,Count:Count+size(temp,2)-1) = temp;
                thetaAmp(:,Count:Count+size(temp,2)-1) = abs(hilbert(bandpass(temp',[5,12],dSF)'));
                deltaAmp(:,Count:Count+size(temp,2)-1) = abs(hilbert(bandpass(temp',[1,4],dSF)'));
                Count = Count+size(temp,2);
                waitbar(ii/L,wb,'Working on amplifier.dat');
            end
            close(wb);
        
            if replaceFile == 1
                fclose(fileID);
                fclose(fID);
                if use60HzFilt == 1
                     fclose(fID2);
                end
                %eval(['delete ',recInfo{i,1},'\amplifier.dat']);
            end
            
            disp('Calculating theta/delta ratio')
            count = zeros(4,1);
            for ii = 1:size(data,1)
        
                if contains(chan{ii,1},'a')
                    idx = 1;
                elseif contains(chan{ii,1},'b')
                    idx = 2;
                elseif contains(chan{ii,1},'c')
                    idx = 3;
                elseif contains(chan{ii,1},'d')
                    idx = 4;
                end
                count(idx,1) = count(idx,1)+1;
                recInfo(i).fileData(idx).data(count(idx,1)).chan = chan{ii,1};
                recInfo(i).fileData(idx).data(count(idx,1)).site = chan{ii,2};
        
                theta_deltaMean = zeros(1,floor(length(thetaAmp)/dSF));
                for i3 = 1:length(theta_deltaMean)
                    try
                        tWin_thetaAmp = thetaAmp(ii,(dSF*(i3-1))+1:dSF*(i3));
                        tWin_deltaAmp = deltaAmp(ii,(dSF*(i3-1))+1:dSF*(i3));
                        theta_deltaMean(1,i3) = mean(tWin_thetaAmp)/mean(tWin_deltaAmp);
                    catch
                        tWin_thetaAmp = thetaAmp(ii,(dSF*(i3-1))+1:end);
                        tWin_deltaAmp = deltaAmp(ii,(dSF*(i3-1))+1:end);
                        theta_deltaMean(1,i3) = mean(tWin_thetaAmp)/mean(tWin_deltaAmp);
                    end
                end
                recInfo(i).fileData(idx).data(count(idx,1)).allTD = theta_deltaMean;
            end
        end
        save([seshPathName,'\recInfo.mat'],'recInfo')
        disp('Analysis complete');
    end

    if strcmp(copyAction,'Yes') == 1
        copyDir = [copyPathName,'\',Time];
        disp(['Copying data to ',copyPathName]);
        eval(['mkdir ',copyDir]);
        copyfile(seshPathName,copyDir);
    end
end

delete(gcp('nocreate'))
disp('Finsished recording');

%%

function timeVec = getTime()
    Now = clock;
    timeVec = (Now(1,4)*3600)+(Now(1,5)*60)+(Now(1,6));
end

function showVideo(frames)

    global videoFig
    matObj = matfile([tempdir,'Go.mat'],'Writable',true);
    numCam = size(frames,5);
    fps = size(frames,4);

    if isempty(videoFig) == 1 || any(ishandle(videoFig) == 0)
        videoFig = figure;
        videoFig = axes(videoFig);
    else
        for i = 1:numCam
            cla(videoFig(1,i))
        end
    end

    if matObj.Go == 1
        for i = 1:2:fps
            for ii = 1:numCam
                videoFig(ii) = subplot(1,numCam,ii);
                subplot(videoFig(ii))
                imshow(frames(:,:,:,i,ii));
            end
            drawnow
        end
    end
end

function plotTD(data)

    global dataFig
    
    if isempty(dataFig) == 1 || ishandle(dataFig) == 0
        dataFig = figure;
        dataFig = axes(dataFig);
    end

    theta = bandpass(data{1,1}',[5,12],data{2,1})';
    delta = bandpass(data{1,1}',[1,4],data{2,1})';
    thetaAmp = abs(hilbert(theta));
    deltaAmp = abs(hilbert(delta));
    theta_delta = thetaAmp./deltaAmp;
    plot(dataFig,theta_delta)
    title(num2str(mean(theta_delta(end-data{2,1}+1:end))))
    drawnow
end

function [Out] = symsepchar(StrIn,Sym)

    % Takes a string of characters (StrIn) seperated by any special character (Sym)
    % and outputs a double array. ex 'C:\Users\Data.txt' -> [C:\,Users\Data.txt] or '11/30/2016' -> [11,30,2016]

    symcount = 1;
    charcount = 1;
    for ci = 1:length(StrIn)
        if strcmp(StrIn(1,ci),Sym) == 1
            symcount = symcount+1;
            charcount = 1;
            continue
        elseif strcmp(StrIn(1,ci),Sym) == 0 
            Out{1,symcount}(1,charcount) = StrIn(1,ci); %#ok<AGROW>
            charcount = charcount+1;
        end
    end
end

function [Out] = returnSubfolder(dirIn)

    Out = dir(dirIn);
    dirRemove = false(length({Out.name}),1);
    for i = 1:length({Out.name})
        if strcmp(Out(i).name,'.') == 1 || strcmp(Out(i).name,'..') == 1 || Out(i).isdir == 0
            dirRemove(i,1) = true;
        end
    end
    Out(dirRemove) = [];
end

function out = convLPFilt(data,cutOff,sampleFreq,filtLength,downSampleRatio)

    % data = input data, processed column-wise
    % cutOff = lowpass filter cutoff frequency, must not be > sampleFreq/2
    % sampleFreq = sampling frequency of input data
    % filtLength = filter length, [] uses default
    % downSampleRatio = ratio to downsample filtered data, [] to not use

    cutOff = cutOff/(sampleFreq/2);

    if isempty(filtLength)
        filtLength = 1025;
    end

    % create filter kernel    
    n = -floor(filtLength/2):floor(filtLength/2); % create time base
    B = sinc(cutOff*n).*cutOff; % make kernel     

    for i = 1:size(data,2)

        P = numel(data(:,i));
        Q = numel(B);
        filtLength = P+Q-1;
        K = 2^nextpow2(filtLength);
    
        % do the convolution
        afft = fft(data(:,i),K);
        bfft = fft(B,K);
        c = ifft(afft(:).*bfft(:));
    
        range = [floor(length(B)/2)+1,filtLength-ceil(length(B)/2)+1];
        data(:,i) = c(range(1):range(2));
    end

    if ~isempty(downSampleRatio)
        data = downsample(data,downSampleRatio);
    end

    out = data;
end