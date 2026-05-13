clear all %#ok<CLALL> 

loopType = 'open'; % Experiment type, 'open' or 'closed'.

StartTime = [5,40]; % Time to begin recording in 24hr time, format = [hr,min]
nextDay = 1; % Start recording tomorrow?

% Block schedule 
firstBL = 60*2; % Pre-stim baseline duration in minutes
stimDur = 60*8; % Duration of stim period in minutes
finalBL = 60*2; % Post-stim baseline duration in minutes

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

boxSubject{1,1} = 'EDS 1.3'; % Subject in box 1
boxSubject{2,1} = 'EDS 2.3'; % Subject in box 2
boxSubject{3,1} = 'EDS 2.2'; % Subject in box 3
boxSubject{4,1} = 'EDS 1.1'; % Subject in box 4

boxRHS(1,1) = 3; % Box for RHS port A
boxRHS(2,1) = 4; % Box for RHS port B
boxRHS(3,1) = 0; % Box for RHS port C
boxRHS(4,1) = 0; % Box for RHS port D

boxRHD(1,1) = 1; % Box for RHD port A
boxRHD(2,1) = 2; % Box for RHD port B
boxRHD(3,1) = 0; % Box for RHD port C
boxRHD(4,1) = 0; % Box for RHD port D

% Stim parameters
ISI = 20:30; % Inter stim interval in seconds, randomly selected.
thresholdTD = 1; % Theda/delta threshold for stim

% Stim amplitudes
stimI(1,1) = 125; % Port A stim amp in uA.
stimI(2,1) = 0; % Port B stim amp in uA.
stimI(3,1) = 0; % Port C stim amp in uA.
stimI(4,1) = 0; % Port D stim amp in uA.

% Stim train period
stimP(1,1) = 10; % Port A stim train period in ms.
stimP(2,1) = 10; % Port B stim train period in ms.
stimP(3,1) = 10; % Port C stim train period in ms.
stimP(4,1) = 10; % Port D stim train period in ms.

% Stim train number of pulses
stimNum(1,1) = 4; % Port A stim train pulse number.
stimNum(2,1) = 4; % Port B stim train pulse number.
stimNum(3,1) = 4; % Port C stim train pulse number.
stimNum(4,1) = 4; % Port D stim train pulse number.

% Stim channel, [0,0] = unused, [c1,0] = monopolar, [c1,c2] = bipolar.
stimChan(1,:) = [7,0]; % Port A
stimChan(2,:) = [7,0]; % Port B
stimChan(3,:) = [0,0]; % Port C
stimChan(4,:) = [0,0]; % Port D

% Input channel for closed-loop, 0 = unused.
dataChan(1,1) = 0; % Port A
dataChan(2,1) = 0; % Port B
dataChan(3,1) = 0; % Port C
dataChan(4,1) = 0; % Port D

framesPerBlock = 128;
blocksPerRead = 125;
numBandsPerChan = 1; % Amplifier channels can be displayed as LOW, WIDE, or HIGH

% Analysis Settings
sF = 20000; % Recording sample frequency
dSF = 1000; % Downsample frequency
use60HzFilt = 1;
ACfiltPath = 'G:\RecordingData\EDS\CNN60Hz.mat';
stimTrig = 1;
Pre = 0;
Post = ISI(1,1);

% Active channels for RHS
activeChan(1,1:3,1) = [{'a-000'},{'PrL(L)'},{1}];
activeChan(2,1:3,1) = [{'a-001'},{'PrL(R)'},{1}];
activeChan(3,1:3,1) = [{'a-002'},{'AVT(L)'},{1}];
activeChan(4,1:3,1) = [{'a-003'},{'BLA(R)'},{1}];
activeChan(5,1:3,1) = [{'a-004'},{'CA1(L)'},{1}];
activeChan(6,1:3,1) = [{'a-005'},{'CA1(R)'},{1}];
activeChan(7,1:3,1) = [{'a-006'},{'LDT(L1)'},{1}];
activeChan(8,1:3,1) = [{'a-007'},{'LDT(L2)'},{1}];

activeChan(1,1:3,2) = [{'b-000'},{'PrL(L)'},{0}];
activeChan(2,1:3,2) = [{'b-001'},{'PrL(R)'},{0}];
activeChan(3,1:3,2) = [{'b-002'},{'AVT(L)'},{0}];
activeChan(4,1:3,2) = [{'b-003'},{'BLA(R)'},{0}];
activeChan(5,1:3,2) = [{'b-004'},{'CA1(L)'},{0}];
activeChan(6,1:3,2) = [{'b-005'},{'CA1(R)'},{0}];
activeChan(7,1:3,2) = [{'b-006'},{'LDT(L1)'},{0}];
activeChan(8,1:3,2) = [{'b-007'},{'LDT(L2)'},{0}];

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
        WL{i,3} = boxSubject{WL{i,2},1};
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
%%

PathName = uigetdir(cd,'Select folder for Datastore');
cd(PathName);

global dataFig videoFig %#ok<NUSED,GVMIS>

vStartTime = (StartTime(1,1)*3600)+(StartTime(1,2)*60);
matObj = matfile([tempdir,'Go.mat'],'Writable',true);

port = {'a','b','c','d'};

hardwareTimer = 1;
intanWriter = 2;
intanReader = 3;
dataHandler = 4;
camReader = 5;
camSaver = 6;

disp(['Will start at ',num2str(StartTime(1,1)),':',num2str(StartTime(1,2))]); 
while 1
    if nextDay == 0
        if getTime >= vStartTime-60
           break
        end
    elseif nextDay == 1
        if getTime <= vStartTime
            nextDay = 0;
        end
    end
    pause(0.001);
end

disp('Connecting to hardware...');

Now = clock;
Time = [num2str(Now(1,2)),'-',num2str(Now(1,3)),'-',num2str(Now(1,1)),'(',num2str(Now(1,4)),'.',num2str(Now(1,5)),')'];

PathName = [cd,'\',Time];
eval(['mkdir ',PathName]);
save([PathName,'\config.mat'])

% Calculations for accurate parsing
firstBL = firstBL*60;
stimDur = stimDur*60;
finalBL = finalBL*60; 

numChan = length(find(dataChan > 0));
chanStim = find(sum(stimChan,2) > 0);
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

if isempty(gcp('nocreate')) == 1
    parpool('local',6);
end 

spmd (6)
    if labindex == hardwareTimer        
        while 1
            try
                Ard = arduino('com3','uno');                          
                break
            catch
                disp('Issue connecting to Arduino')
                pause(3)
            end   
        end                  

        matObj.Go = 1; 

        % sync with system clock
        send(textOut,'Syncing hardware');
        while 1
            if getTime >= vStartTime
                break
            end
        end   

        labBarrier
        if firstBL > 0
            send(textOut,'Pre-stime baseline');
        end
        
        Count = 0;
        ck = [0,0];
        tS = 0;
        seshStart = getTime;
        labSend(seshStart,intanWriter);
        labSend(seshStart,intanReader);
        labSend(seshStart,dataHandler);
        while 1 
            t = getTime;
            while getTime-t <= 0.999
                % Wait for next second
            end
            Count = Count+1;             

            if Count == 1                  
                writeDigitalPin(Ard,ArdPin{1,1},1);
            elseif Count == LEDtime
                Count = 0;
            end       
            
            if t >= seshStart+firstBL && ck(1,1) == 0
                writeDigitalPin(Ard,ArdPin{2,1},1);
                send(textOut,'Begin stim');
                ck(1,1) = 1;
                tS = 1;
            elseif t >= seshStart+firstBL+stimDur && ck(1,2) == 0
                writeDigitalPin(Ard,ArdPin{2,1},1);                
                send(textOut,'Post-stime baseline');
                ck(1,2) = 1;
                tS = 1;
            elseif t >= seshStart+firstBL+stimDur+finalBL
                send(textOut,'Finished');
                break
            else
                tS = 0;
            end

            pause(0.2)
            if Count == 1    
                writeDigitalPin(Ard,ArdPin{1,1},0);
            end
            if tS == 1
                writeDigitalPin(Ard,ArdPin{2,1},0);
            end
        end
        
        matObj.Go = 0;

    elseif labindex == intanWriter
        if useRHS+useRHD > 0
            check = [0,0];
            while true
                try
                    if useRHS == 1 && check(1,1) == 0
                        intanSend{1,1} = tcpclient(intanIP_RHS,intanPort1_RHS);
                        write(intanSend{1,1},uint8(['set Filename.Path ',PathName]));
                        write(intanSend{1,1},uint8('set FileFormat OneFilePerSignalType'));
                        write(intanSend{1,1},uint8('set Filename.BaseFilename RHS'));
                        check(1,1) = 1;
                    end
                    if useRHD == 1 && check(1,2) == 0
                        intanSend{2,1} = tcpclient(intanIP_RHD,intanPort1_RHD);
                        write(intanSend{2,1},uint8(['set Filename.Path ',PathName]));
                        write(intanSend{2,1},uint8('set FileFormat OneFilePerSignalType'));
                        write(intanSend{2,1},uint8('set Filename.BaseFilename RHD'));
                        check(1,2) = 1;
                    end
                catch
                    if useRHD == 1 && check(1,2) == 0
                        disp('Issue connecting to Intan RHS commands')
                    end
                    if useRHD == 1 && check(1,2) == 0
                        disp('Issue connecting to Intan RHD commands')
                    end
                end
                if ((useRHS == 1) ~= (useRHD == 1) && sum(check) == 1) || sum(check) == 2
                    break
                end
                pause(1)
            end
        end
       
        for ii = 1:size(stimChan,1)
            if stimChan(ii,1) > 0 && stimI(ii,1) > 0

                if dataChan(ii,1) > 0 && strcmp(loopType,'closed') == 1
                    write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{dataChan(ii,1),1},'.tcpdataoutputenabled true;']));
                end
                
                for i3 = 1:size(stimChan(ii,:),2)
                    if stimChan(ii,i3) > 0
        
                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3),1},'.stimenabled true']));
                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3),1},'.source KeyPressF',num2str(ii)]));
                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3),1},'.FirstPhaseAmplitudeMicroAmps ',num2str(stimI(ii,1))]));
                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3),1},'.SecondPhaseAmplitudeMicroAmps ',num2str(stimI(ii,1))]));
                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3),1},'.NumberOfStimPulses ',num2str(stimNum(ii,1))]));
                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3),1},'.PulseTrainPeriodMicroseconds ',num2str(stimP(ii,1)*1000)]));
                        
                        if i3 == 1
                            write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3),1},'.Polarity NegativeFirst']));
                        elseif i3 == 2
                            write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{stimChan(ii,i3),1},'.Polarity PositiveFirst']));
                        end
                    end
                end
                send(textOut,['Stim port',num2str(ii),' amp set to ',num2str(stimI(ii,1)),'uA.'])
            end
        end
        if useRHS == 1 && sum(stimI) > 0
            write(intanSend{1,1},uint8('execute uploadstimparameters'));
            pause(10)
        end

        labBarrier
        seshStart = labReceive(hardwareTimer);
        
        if useRHS == 1
            write(intanSend{1,1},uint8('set runMode record'));
        end
        if useRHD == 1
            write(intanSend{2,1},uint8('set runMode record'));
        end
        
        % Pre baseline
        while getTime <= seshStart+firstBL
            pause(0.001)
        end
        
        % Stim
        while 1
            
            stim = labReceive(dataHandler);
            if isempty(stim) == 1
                break
            end
            for ii = chanStim'
                if stim(ii,1) == 1
                    write(intanSend{1,1},uint8(['execute ManualStimTriggerPulse F',num2str(ii)]));
                end
            end
        end
        
        % Post baseline
        while getTime <= seshStart+firstBL+stimDur+finalBL
            pause(0.001)
        end        
        
        if useRHS == 1
            write(intanSend{1,1},uint8('set runMode stop'));
        end
        if useRHD == 1
            write(intanSend{2,1},uint8('set runMode stop'));
        end
    elseif labindex == intanReader
        if strcmp(loopType,'closed') == 1
            while 1
                try
                    intanRead = tcpclient(intanIP_RHS,intanPort2_RHS);
                    break
                catch
                    disp('Issue connecting to Intan data output')
                    pause(3)
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
    
        labBarrier
        seshStart = labReceive(hardwareTimer);

        if strcmp(loopType,'closed') == 1
            
            % Pre baseline
            while getTime <= seshStart+firstBL
                flush(intanRead);
                pause(0.001)
            end
            
            flush(intanRead);
            while getTime <= seshStart+firstBL+stimDur
    
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
            labSend([],dataHandler);
        end

        while getTime <= seshStart+firstBL+stimDur+finalBL
            if strcmp(loopType,'closed') == 1
                flush(intanRead);
            end
            pause(0.001)
        end   
    elseif labindex == dataHandler
        
        if strcmp(loopType,'closed') == 1

            bankISI = (round(ISI(1,1)/readTime)):round((ISI(1,end)/readTime));
            nextISI(1:numChan,1) = bankISI(1,1);
            readCount(1:numChan,1) = 0;
            theta_delta = zeros(numChan,1);
        else

            Count = nan(4,1);
            Count(chanStim,1) = 0;
            nextISI = nan(4,1);
            nextISI(chanStim,1) = 0;
            StimIter = round((stimDur*60)/mean(ISI));
            numStim = 0;
        end

        labBarrier
        seshStart = labReceive(hardwareTimer);

        while getTime <= seshStart+firstBL
            pause(0.001)
        end

        while 1
            if strcmp(loopType,'closed') == 1

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
                
                stim = Count == 0;
                if any(stim == 1)
                    labSend(stim,intanWriter);
                    numStim = numStim+1;
                end

                if numStim >= StimIter || getTime >= seshStart+firstBL+stimDur
                    break
                end

                Count = Count+1;
                
                for ii = chanStim'
                    if stim(ii,1) == 1
                        nextISI(ii,1) = datasample(ISI,1);
                    end
                    if Count(ii,1) == nextISI(ii,1)
                        Count(ii,1) = 0;
                    end
                end
            end
        end
        labSend([],intanWriter);

        while getTime <= seshStart+firstBL+stimDur+finalBL
            pause(0.001)
        end     

    elseif labindex == camReader 
        
        % Create all webcam objects           
        cam = cell(numCam,1); 
        for ii = 1:numCam                
            cam{ii,1} = webcam(camChoice(1,ii));
        end           
        frames = uint8(zeros(Dim(1,1),Dim(1,2),3,fps,numCam));

        labBarrier
        
        while matObj.Go == 1
            
            Count = 1;
            t = getTime;
            while getTime-t <= 0.999 && Count <= fps
                for ii = 1:numCam
                    frames(:,:,:,Count,ii) = snapshot(cam{ii,1});
                end
                Count = Count+1;
            end

            labSend(frames,camSaver);
            if ShowVideo == 1
                send(videoOut,frames);
            end
        end 
        labSend([],camSaver);
    elseif labindex == camSaver       
        
        VidFileName = cell(numCam,1);
        writerObj = cell(numCam,1);

        % Create all video writer objects
        for ii = 1:numCam   
                   
            VidFileName{ii,1} = [WL{camChoice(1,ii),3},'.mp4'];
            writerObj{ii,1} = VideoWriter([PathName,'\',VidFileName{ii,1}],'MPEG-4'); %#ok<TNMLP>
            writerObj{ii,1}.FrameRate = fps;
            open(writerObj{ii,1});
        end

        videos = uint8(zeros(Dim(1,1),Dim(1,2),3,fps,numCam));

        labBarrier
        while 1

            frames = labReceive(camReader);
            if isempty(frames)
                break
            end 
            for ii = 1:numCam
                for iii = 1:fps
                    writeVideo(writerObj{ii,1},frames(:,:,:,iii,ii));
                end
            end
        end
        for ii = 1:numCam
             close(writerObj{ii,1});
        end
    end
end
close(gcf)

delete(gcp('nocreate'))
disp('Finsished recording');
%%
disp('Starting analysis');

recInfo = [{'File'},{'System'},{'Data'}];

listing = dir(PathName);
dirRemove = false(length({listing.name}),1);
for i = 1:length({listing.name})
    if strcmp(listing(i).name,'.') == 1 || strcmp(listing(i).name,'..') == 1 || listing(i).isdir == 0
        dirRemove(i,1) = true;
    end
end
listing(dirRemove) = [];

numFolders = size(listing,1);
for i = 1:size(listing,1)

    recInfo{i+1,3} = cell(5,7);
    recInfo{i+1,3}(1,:) = [{'Port'},{'Subject'},{'Stim amp (uI)'},{'Stim Chan'},{'Stim start/stop'},{'Stim TS'},{'Theta/delta'}];
    recInfo{i+1,3}{2,1} = 'A';
    recInfo{i+1,3}{3,1} = 'B';
    recInfo{i+1,3}{4,1} = 'C';
    recInfo{i+1,3}{5,1} = 'D';

    if contains(listing(i).name,'RHS') == 1

        recInfo{i+1,1} = [listing(i).folder,'\',listing(i).name];
        recInfo{i+1,2} = 'RHS';

        for ii = 1:4
            if boxRHS(ii,1) > 0
                recInfo{i+1,3}{ii+1,2} = boxSubject(boxRHS(ii,1),1);
                recInfo{i+1,3}{ii+1,3} = stimI(ii,1);
                recInfo{i+1,3}{ii+1,4} = stimChan(ii,:);
            end
        end
    elseif contains(listing(i).name,'RHD') == 1

        recInfo{i+1,1} = [listing(i).folder,'\',listing(i).name];
        recInfo{i+1,2} = 'RHD';
        
        for ii = 1:4
            if boxRHD(ii,1) > 0
                recInfo{i+1,3}{ii+1,2} = boxSubject(boxRHD(ii,1),1);
            end
        end
    end    
end

for i = 2:numFolders+1

    if strcmp(recInfo{i,2},'RHD') == 1
        continue
    end

    % Digital in
    disp('Working on digitalin.dat')
    s = dir([recInfo{i,1},'\digitalin.dat']);
    numSamples = s.bytes/2;

    fileID = fopen([recInfo{i,1},'\digitalin.dat']);
    fseek(fileID,0,'bof');
    L = floor(numSamples/(sF*600));
    digIN = nan(1,L*600*DsF);
    
    Count = 1;
    for ii = 1:L
        
        temp = fread(fileID,[1,sF*600],'int16');
        temp = downsample(temp,sF/dSF);
        temp(temp ~= stimTrig) = 0;
        digIN(:,Count:Count+size(temp,2)-1) = temp;
        Count = Count+size(temp,2);
    end

    timeStamps = 0;
    ck = 0;
    count = 1;
    for ii = 1:length(digIN)
        if digIN(1,ii) == stimTrig && ck == 0            

            timeStamps(1,count) = ii; %#ok<SAGROW> 
            count = count+1;
            ck = 1;
        elseif digIN(1,ii) == 0 && ck == 1
            ck = 0;
        end
    end
    recInfo{i,3}([false;stimI > 0],5) = {round(timeStamps/dSF)};
    clear digIN

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
    s = dir([recInfo{i,1},'\stim.dat']);
    numSamples = s.bytes/(2*numChan);   
    
    fileID = fopen([recInfo{i,1},'\stim.dat']);
    fseek(fileID,0,'bof');
    L = floor(numSamples/(sF*600));
    stimTS = nan(numChan,L*600*DsF);
    
    Count = 1;
    for ii = 1:L
        
        temp = fread(fileID,[numChan,sF*600],'int16');
        temp = downsample(temp',sF/dSF)';
        temp(temp > 0) = 1;
        stimTS(:,Count:Count+size(temp,2)-1) = temp;
        Count = Count+size(temp,2);
    end   

    count = 1;
    for ii = 1:size(activeChan,3)
        for i3 = 1:size(activeChan,1)
            if ~isempty(recInfo{i,3}{ii+1,3}) && activeChan{i3,3,ii} == 1
                if any(recInfo{i,3}{ii+1,4} == i3)

                    temp = stimTS(count,:);
                    temp = find(diff(temp) > 0);
                    recInfo{i,3}{ii+1,6} = temp';
                end
                count = count+1;
            end
        end
    end
    clear stimTS
    
    % Amplifier in
    disp('Working on amplifier.dat')
    s = dir([recInfo{i,1},'\amplifier.dat']);
    numSamples = s.bytes/(2*numChan);   
    
    fileID = fopen([recInfo{i,1},'\amplifier.dat']);
    fseek(fileID,0,'bof');
    L = floor(numSamples/(sF*600));
    data = nan(numChan,L*600*DsF);
    thetaAmp = nan(numChan,L*600*DsF);
    deltaAmp = nan(numChan,L*600*DsF);

    Count = 1;
    for ii = 1:L
        
        temp = fread(fileID,[numChan,sF*600],'int16');
        temp = downsample(temp',sF/dSF)';
        data(:,Count:Count+size(temp,2)-1) = temp;
        thetaAmp(:,Count:Count+size(temp,2)-1) = abs(hilbert(bandpass(temp',[5,12],dSF)'));
        deltaAmp(:,Count:Count+size(temp,2)-1) = abs(hilbert(bandpass(temp',[1,4],dSF)'));
        Count = Count+size(temp,2);
    end

    if use60HzFilt == 1
        
        disp('Applying ANN 60Hz filter')
        load(ACfiltPath);
        data2 = data;
        for ii = 1:dSF:size(data2,2)
            for i3 = 1:size(data2,1)
                data2(i3,ii:ii+dSF-1) = predict(CNN60Hz,data2(i3,ii:ii+dSF-1)/50000,'ExecutionEnvironment','gpu')*50000;
            end
            if rem((ii/size(data2,2))*100,1) == 0
                disp([num2str(Count),'%'])
                Count = Count+1;
            end
        end
        
        fID = fopen([recInfo{i,1},'\amplifier_60HzFilt.dat'],'w');
        Count = 1;
        for ii = 1:size(data2,2)
            for i3 = 1:size(data2,1)
                fwrite(fID,data2(i3,ii),'int16');
            end

            if rem((ii/size(data2,2))*100,1) == 0
                disp([num2str(Count),'%'])
                Count = Count+1;
            end
        end
        fclose(fID);
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
        recInfo{i,3}{idx+1,7}(count(idx,1),1:2) = chan(ii,:);

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
        recInfo{i,3}{idx+1,7}{count(idx,1),3} = theta_deltaMean;
    end
end
save([PathName,'\recInfo.mat'],'recInfo')

beep
disp('Finsished');
%%



for i = 1:size(recInfo{3,3}{3,7},1)

    subplot(size(recInfo{3,3}{3,7},1),1,i)
    hold on
    plot(recInfo{3,3}{3,7}{i,3})
    axis tight
    xline(recInfo{3,3}{3,5}(1,1),'LineWidth',3)
    xline(recInfo{3,3}{3,5}(1,2),'LineWidth',3)
    ylim([0,15])
    title(recInfo{3,3}{3,7}{i,2})
    hold off
end

%%
% d = recInfo{2,4}{6,2}(1:round(recInfo{2,3}(1,1)));
% numD9 = round(length(d)*0.9);
% h1 = histogram(d,50);
% bins = h1.BinEdges;
% h1 = h1.BinCounts;
% 
% d = recInfo{2,4}{6,2}(1,round(recInfo{2,3}(1,1)):round(recInfo{2,3}(2,1)));
% numD9 = round(length(d)*0.9);
% h2 = histogram(d,'BinEdges',bins);
% h2 = h2.BinCounts;
% 
% d = recInfo{2,4}{6,2}(1,round(recInfo{2,3}(2,1)):end);
% numD9 = round(length(d)*0.9);
% h3 = histogram(d,'BinEdges',bins);
% h3 = h3.BinCounts;
% 
% medianBin = zeros(1,50);
% count = 1;
% for i = 1:length(bins)-1
% 
%     medianBin(1,count) = mean(bins(1,i:i+1));
%     count = count+1;
% end
% 
% figure
% subplot(2,1,1)
% bar(medianBin,h2-h1)
% title('stim - pre-baseline hist')
% 
% subplot(2,1,2)
% bar(medianBin,h3-h1)
% title('post-baseline - pre-baseline hist')


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