clear all %#ok<CLALL>

loopType = 'open'; % Experiment type, 'open' or 'closed'.

StartTime = [12,0]; % Time to begin recording in 24hr time, format = [hr,min]
nextDay = 0; % Start recording tomorrow? 0 = no, 1 = yes


% define who is in which box

boxSubject{1,1} = 'EDS 5.1'; % Subject in box 1
boxSubject{2,1} = 'EDS 4.2'; % Subject in box 2
boxSubject{3,1} = 'EDS 4.0'; % Subject in box 3
boxSubject{4,1} = 'EDS 3.0'; % Subject in box 4

% Block schedule.
% BL = Baseline = record for observation no stim (RHS only)
% S1 = Stim Loop =  stim (many short bursts)
% S2 = Stim Train = induction stim (long duration stim)
% SX = Stim sham = like BL but with control TTLs for sham stims
% OD = Optional delay = block of time with no stim
%

% define experimental sequence per port
% e.g. % [{Day1 block},{Day2 block},....];

% Baseline only

% Port A
%blockSched{1,1} =  [{'BL'}];%,{'BL'},{'BL'},{'BL'},{'BL'},{'BL'},{'BL'}];

% Port B
%blockSched{1,2} =  [{'BL'}];%,{'BL'},{'BL'},{'BL'},{'BL'},{'BL'},{'BL'}];

% Port C
%blockSched{1,3} =   [{'BL'}];%,{'BL'},{'BL'},{'BL'},{'BL'},{'BL'},{'BL'}];

% Port D
%blockSched{1,4} =   [{'BL'}];%,{'BL'},{'BL'},{'BL'},{'BL'},{'BL'},{'BL'}];

%blockTime = 18*3600*ones(1,1);


% Main experimet

% Port A
blockSched{1,1} = [{'BL'};{'BL'};{'BL'};{'BL'}];
blockSched{2,1} = [{'S1'};{'SX'};{'S1'};{'SX'}];
blockSched{3,1} = [{'S2'};{'S2'};{'S2'};{'S2'}];
blockSched{4,1} = [{'BL'};{'BL'};{'BL'};{'BL'}];

% Port B
blockSched{1,2} = [{'BL'};{'BL'};{'BL'};{'BL'}];
blockSched{2,2} = [{'S1'};{'SX'};{'S1'};{'SX'}];
blockSched{3,2} = [{'S2'};{'S2'};{'S2'};{'S2'}];
blockSched{4,2} = [{'BL'};{'BL'};{'BL'};{'BL'}];

% Port C
blockSched{1,3} = [{'BL'};{'BL'};{'BL'};{'BL'}];
blockSched{2,3} = [{'S1'};{'SX'};{'S1'};{'SX'}];
blockSched{3,3} = [{'S2'};{'S2'};{'S2'};{'S2'}];
blockSched{4,3} = [{'BL'};{'BL'};{'BL'};{'BL'}];

% Port D
blockSched{1,4} = [{'BL'};{'BL'};{'BL'};{'BL'}];
blockSched{2,4} = [{'S1'};{'SX'};{'S1'};{'SX'}];
blockSched{3,4} = [{'S2'};{'S2'};{'S2'};{'S2'}];
blockSched{4,4} = [{'BL'};{'BL'};{'BL'};{'BL'}];

% Block schedule in seconds.
blockTime(1,1) = 3600;
blockTime(2,1) = 600;
blockTime(3,1) = 10;
blockTime(4,1) = 10190;


% Arduino Settings
ArdPin{1,1} = 'D13'; % Arduino pin for LED.
ArdPin{2,1} = 'D8'; % Arduino pin for Timestamp.
LEDtime = 10; % Time between LED blinks is seconds

pyrun('import pyautogui')

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
chanID(1,1:2) = [{'000'},{'M2(L)'}];
chanID(2,1:2) = [{'001'},{'M2(R)'}];
chanID(3,1:2) = [{'002'},{'AVT(L)'}];
chanID(4,1:2) = [{'003'},{'BLA(R)'}];
chanID(5,1:2) = [{'004'},{'CA1(L)'}];
chanID(6,1:2) = [{'005'},{'CA1(R)'}];
chanID(7,1:2) = [{'006'},{'LDT(L1)'}];
chanID(8,1:2) = [{'007'},{'LDT(L2)'}];


boxRHS(1,1) = 1; % Box for RHS port A
boxRHS(2,1) = 2; % Box for RHS port B
boxRHS(3,1) = 3; % Box for RHS port C
boxRHS(4,1) = 4; % Box for RHS port D

boxRHD(1,1) = 5; % Box for RHD port A
boxRHD(2,1) = 6; % Box for RHD port B
boxRHD(3,1) = 7; % Box for RHD port C
boxRHD(4,1) = 8; % Box for RHD port D

modDir = 'R:\McKenzieLab\DGregg\SeizureForecast\Seizuredetect_demo\Training\randomForest.mat';

% S1 (Stim 1) parameters

ISI = {20:30}; % Inter stim interval in seconds, randomly selected.
thresholdTD = 1; % Theda/delta threshold for stim

% Amplitudes =  (port, day, S1/S2)
stimI(1,:,1) = [0]; % Port A stim amp in uA.
stimI(2,:,1) = [0]; % Port B stim amp in uA.
stimI(3,:,1) = [0]; % Port C stim amp in uA.
stimI(4,:,1) = [0]; % Port D stim amp in uA.

% Channel, [0,0] = unused, [c1,0] = monopolar, [c1,c2] = bipolar.
stimChan(1,:,1) = [0,0]; % Port A
stimChan(2,:,1) = [0,0]; % Port B
stimChan(3,:,1) = [0,0]; % Port C
stimChan(4,:,1) = [0,0]; % Port D

% Train period (port, 1, S1/S2)
stimP(1,1,1) = 10; % Port A stim train period in ms.
stimP(2,1,1) = 10; % Port B stim train period in ms.
stimP(3,1,1) = 10; % Port C stim train period in ms.
stimP(4,1,1) = 10; % Port D stim train period in ms.

% Train number of pulses (port, 1, S1/S2)
stimNum(1,1,1) = 4; % Port A stim train pulse number.
stimNum(2,1,1) = 4; % Port B stim train pulse number.
stimNum(3,1,1) = 4; % Port C stim train pulse number.
stimNum(4,1,1) = 4; % Port D stim train pulse number.

% Stim 2 parameters
ISI{2,1} = 10; % Inter stim interval in seconds, randomly selected.

% Amplitudes.  (port, day, S1/S2)
stimI(1,:,2) = [0]; % Port A stim amp in uA.
stimI(2,:,2) = [0]; % Port B stim amp in uA.
stimI(3,:,2) = [0]; % Port C stim amp in uA.
stimI(4,:,2) = [0]; % Port D stim amp in uA.

% Channel, [0,0] = unused, [c1,0] = monopolar, [c1,c2] = bipolar.
stimChan(1,:,2) = [0,0]; % Port A
stimChan(2,:,2) = [0,0]; % Port B
stimChan(3,:,2) = [0,0]; % Port C
stimChan(4,:,2) = [0,0]; % Port D

% Train period (Train period (port, 1, S1/S2))
stimP(1,1,2) = 100; % Port A stim train period in ms.
stimP(2,1,2) = 100; % Port B stim train period in ms.
stimP(3,1,2) = 100; % Port C stim train period in ms.
stimP(4,1,2) = 100; % Port D stim train period in ms.

% Train number of pulses (Train number of pulses (port, 1, S1/S2)
stimNum(1,1,2) = 100; % Port A stim train pulse number.
stimNum(2,1,2) = 100; % Port B stim train pulse number.
stimNum(3,1,2) = 100; % Port C stim train pulse number.
stimNum(4,1,2) = 100; % Port D stim train pulse number.

% Input channel (base 1) for closed-loop, 0 = unused.
dataChan(1,:) = zeros(1,8); % Port A
dataChan(2,:) = 1:8; % Port B
dataChan(3,:) = zeros(1,8); % Port C
dataChan(4,:) = zeros(1,8); % Port D

winSize = 2;
framesPerBlock = 128;
blocksPerRead = 200;
numBandsPerChan = 1; % Amplifier channels can be displayed as LOW, WIDE, or HIGH

% Analysis Settings
runAnalysis = 1;
sF = 20000; % Recording sample frequency
dSF = 2000; % Downsample frequency
stimTrig = 1;
use60HzFilt = 0;
ACfiltPath = 'G:\RecordingData\EDS\CNN60Hz.mat'; % for adaptive line noise filter
replaceFile = 1;

% Active channels for RHS
activeChan(1,1:3,1) = [{'a-000'},{'M2(L)'},{1}];
activeChan(2,1:3,1) = [{'a-001'},{'M2(R)'},{1}];
activeChan(3,1:3,1) = [{'a-002'},{'AVT(L)'},{1}];
activeChan(4,1:3,1) = [{'a-003'},{'BLA(R)'},{1}];
activeChan(5,1:3,1) = [{'a-004'},{'CA1(L)'},{1}];
activeChan(6,1:3,1) = [{'a-005'},{'CA1(R)'},{1}];
activeChan(7,1:3,1) = [{'a-006'},{'LDT(L1)'},{1}];
activeChan(8,1:3,1) = [{'a-007'},{'LDT(L2)'},{1}];

activeChan(1,1:3,2) = [{'b-000'},{'M2(L)'},{1}];
activeChan(2,1:3,2) = [{'b-001'},{'M2(R)'},{1}];
activeChan(3,1:3,2) = [{'b-002'},{'AVT(L)'},{1}];
activeChan(4,1:3,2) = [{'b-003'},{'BLA(R)'},{1}];
activeChan(5,1:3,2) = [{'b-004'},{'CA1(L)'},{1}];
activeChan(6,1:3,2) = [{'b-005'},{'CA1(R)'},{1}];
activeChan(7,1:3,2) = [{'b-006'},{'LDT(L1)'},{1}];
activeChan(8,1:3,2) = [{'b-007'},{'LDT(L2)'},{1}];

activeChan(1,1:3,3) = [{'c-000'},{'M2(L)'},{0}];
activeChan(2,1:3,3) = [{'c-001'},{'M2(R)'},{0}];
activeChan(3,1:3,3) = [{'c-002'},{'AVT(L)'},{0}];
activeChan(4,1:3,3) = [{'c-003'},{'BLA(R)'},{0}];
activeChan(5,1:3,3) = [{'c-004'},{'CA1(L)'},{0}];
activeChan(6,1:3,3) = [{'c-005'},{'CA1(R)'},{0}];
activeChan(7,1:3,3) = [{'c-006'},{'LDT(L1)'},{0}];
activeChan(8,1:3,3) = [{'c-007'},{'LDT(L2)'},{0}];

activeChan(1,1:3,4) = [{'d-000'},{'M2(L)'},{0}];
activeChan(2,1:3,4) = [{'d-001'},{'M2(R)'},{0}];
activeChan(3,1:3,4) = [{'d-002'},{'AVT(L)'},{0}];
activeChan(4,1:3,4) = [{'d-003'},{'BLA(R)'},{0}];
activeChan(5,1:3,4) = [{'d-004'},{'CA1(L)'},{0}];
activeChan(6,1:3,4) = [{'d-005'},{'CA1(R)'},{0}];
activeChan(7,1:3,4) = [{'d-006'},{'LDT(L1)'},{0}];
activeChan(8,1:3,4) = [{'d-007'},{'LDT(L2)'},{0}];


%%
%model info
model_fname  = 'R:\McKenzieLab\DGregg\SeizureForecast\Seizuredetect_demo\Training\randomForest.mat';
v = load(model_fname);
ops.numfreq  = [4.9334 7.7869 12.2910 19.4002 30.6214 48.3330 76.2892 120.4155 190.0649 300.0000];
ops.nchanFil= 16;
ops.nChanSubj = 8;
ops.channels = 8:15;
ops.Fs =  2000;
ops.nTempBin =  32;
ops.data =  15;
ops.Twin =  2;
ops.feature_fun = @sm_GetDataFeature2;



PathName = uigetdir(cd,'Select folder for Datastore');
cd(PathName);

promptMessage = sprintf('Would you like to copy to another directory?');
titleBarCaption = 'settings';
copyAction = questdlg(promptMessage, titleBarCaption, 'Yes','No','Yes');

if strcmp(copyAction,'Yes') == 1
    copyPathName = uigetdir(cd,'Select directory');
    save([copyPathName,'\config.mat'])
end

save([PathName,'\config.mat'])
%%
%RUN MAIN EXPERIMENT - DO NOT EDIT BELOW
clc
global  dataFig %#ok<NUSED,GVMIS>

vStartTime = (StartTime(1,1)*3600)+(StartTime(1,2)*60);

port = {'a','b','c','d'};

hardwareTimer = 1;
intanWriter = 2;
intanReader = 3;
dataHandler = 4;


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
        parpool('local',4);
        disp('Connecting to hardware...');
    end
   

    pyrun('pyautogui.press(''F11'')') % start multicam
    send(textOut,'Recording video')
    pause(4)

    spmd (4)
        if spmdIndex == hardwareTimer

            noCon_Ard = 0;
            if Iter == startIter

                t = 0;
                while 1
                    try
                        Ard = arduino('com3','Uno');
                        send(textOut,'Connected to Arduino')

                      
                        break
                    catch
                        send(textOut,'Issue connecting to Arduino or multicam')
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

            spmdSend(noCon_Ard,[intanWriter,intanReader,dataHandler],0);
            noCon_TCPwrite = spmdReceive(intanWriter,0);
            noCon_TCPread = spmdReceive(intanReader,0);

            if noCon_Ard == 0 && noCon_TCPwrite == 0 && noCon_TCPread == 0

                % sync with system clock
                send(textOut,'Syncing hardware');
                while 1
                    if getTime >= vStartTime
                        break
                    end
                end

                spmdBarrier

                seshStart = getTime;
                spmdSend(seshStart,[intanWriter,intanReader,dataHandler],1);

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
        elseif spmdIndex == intanWriter

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
                                    for i3 = 1:size(dataChan,2)
                                        write(intanSend{1,1},uint8(['set ',port{1,ii},'-',chanID{dataChan(ii,i3),1},'.tcpdataoutputenabled true;']));
                                    end
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

            spmdSend(noCon_TCPwrite,[hardwareTimer,intanReader,dataHandler],0);
            noCon_Ard = spmdReceive(hardwareTimer,0);
            noCon_TCPread = spmdReceive(intanReader,0);

            spmdBarrier

            seshStart = spmdReceive(hardwareTimer,1);
            if noCon_Ard == 0 && noCon_TCPwrite == 0 && noCon_TCPread == 0

                if useRHS == 1
                    write(intanSend{1,1},uint8('set runMode record'));
                end
                if useRHD == 1
                    write(intanSend{2,1},uint8('set runMode record'));
                end

                while getTime < seshStart+sum(blockTime(1:end,1))

                    stim = spmdReceive(dataHandler);
                    if isempty(stim) == 1
                        break
                    end

                    if stim
                        write(intanSend{1,1},uint8('execute ManualStimTriggerPulse F2'));
                        send(textOut,'stim')
                    end


                    % fix this for stim on multiple boxes
                    % Count = 1;
                    % for ii = 1:size(chanStim,3)
                    %     for i3 = 1:size(chanStim,2)
                    %         if stim(ii,i3) == 1
                    %             write(intanSend{1,1},uint8(['execute ManualStimTriggerPulse F',num2str(Count)]));
                    %         end
                    %         Count = Count+1;
                    %     end
                    % end


                    %     if useRHS == 1
                    %         write(intanSend{1,1},uint8('set runMode stop'));
                    %     end
                    %     if useRHD == 1
                    %         write(intanSend{2,1},uint8('set runMode stop'));
                    %     end
                    %
                    %     if check(1,1) == 1
                    %         write(intanSend{1,1},uint8(['set Filename.Path ',seshPathName]));
                    %         write(intanSend{1,1},uint8('set FileFormat OneFilePerSignalType'));
                    %         write(intanSend{1,1},uint8('set Filename.BaseFilename RHS_Dark'));
                    %     end
                    %     if check(1,2) == 1
                    %         write(intanSend{2,1},uint8(['set Filename.Path ',seshPathName]));
                    %         write(intanSend{2,1},uint8('set FileFormat OneFilePerSignalType'));
                    %         write(intanSend{2,1},uint8('set Filename.BaseFilename RHD_Dark'));
                    %     end
                    %
                    %     if useRHS == 1
                    %         write(intanSend{1,1},uint8('set runMode record'));
                    %     end
                    %     if useRHD == 1
                    %         write(intanSend{2,1},uint8('set runMode record'));
                    %     end
                    % end
                end

                if useRHS == 1
                    write(intanSend{1,1},uint8('set runMode stop'));
                end
                if useRHD == 1
                    write(intanSend{2,1},uint8('set runMode stop'));
                end
                 pyrun('import pyautogui')
                 pyrun('pyautogui.press(''F11'')') % start multicam
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
                        pause(10)
                    end
                end
            end

        elseif spmdIndex == intanReader

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
                ampTimewindow = zeros(numChan,numDsFrames*winSize);

                % Initialize counters
                chunkCount = 0;
                blockCount = 0;
                amplifierTimestampsIndex = 1;
                Stopper = 0;
            else
                intanRead = [];
            end

            spmdSend(noCon_TCPread,[hardwareTimer,intanWriter,dataHandler],0);
            noCon_Ard = spmdReceive(hardwareTimer,0);
            noCon_TCPwrite = spmdReceive(intanWriter,0);

            spmdBarrier

            seshStart = spmdReceive(hardwareTimer,1);
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

                                    spmdSend(ampTimewindow,dataHandler);

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
                            spmdSend([],dataHandler);
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
        elseif spmdIndex == dataHandler

            noCon_Ard = spmdReceive(hardwareTimer,0);
            noCon_TCPwrite = spmdReceive(intanWriter,0);
            noCon_TCPread = spmdReceive(intanReader,0);

            spmdBarrier

            seshStart = spmdReceive(hardwareTimer,1);
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
                        nextISI(1,1) = bankISI(1,1); % for more boxes make this an array
                        readCount(1,1) = 0; % for more boxes make this an array
                        theta_delta = zeros(numChan,1);

                        %model info

                        v = load(modDir);
                        ops.numfreq  = [4.9334 7.7869 12.2910 19.4002 30.6214 48.3330 76.2892 120.4155 190.0649 300.0000];
                        ops.nchanFil= 16;
                        ops.nChanSubj = 8;
                        ops.channels = 8:15;
                        ops.Fs =  2000;
                        ops.nTempBin =  32;
                        ops.data =  15;
                        ops.Twin =  2;
                        ops.feature_fun = @sm_GetDataFeature2;

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

                                ampTimewindow = spmdReceive(intanReader);
                                if isempty(ampTimewindow) == 1
                                    break
                                end
                                %send(dataOut,[{ampTimewindow};{dSF}]);

                                readCount = readCount+1;
                                % send(textOut,readCount)
                                if any(readCount >= nextISI)

                                    % Calculate theta/delta ratio
                                    % theta = bandpass(ampTimewindow',[5,12],dSF)';
                                    % delta = bandpass(ampTimewindow',[1,4],dSF)';
                                    % thetaAmp = abs(hilbert(theta')');
                                    % deltaAmp = abs(hilbert(delta')');
                                    % for ii = 1:numChan
                                    %     theta_delta(ii,1) = mean(thetaAmp(ii,end-dSF+1:end))/mean(deltaAmp(ii,end-dSF+1:end));
                                    % end


                                    % stim = theta_delta < thresholdTD;

                                    Prob = sm_getSeizProb(ampTimewindow',v.rusTree,ops); % for multiple boxes need to loop and index into ampTimeWindow
                                    stim = ~isempty(strmatch(char(categorical(Prob)),'true'));
                                    %send(textOut,stim)
                                    stim = stim & (readCount >= nextISI);

                                    spmdSend(stim,intanWriter);
                                    %for ii = 1:numChan
                                    if stim

                                        %   send(textOut,['Stim Ch',num2str(7),' @ ',num2str(theta_delta(ii,1)),' ISI: ',num2str(readCount(ii,1))]);
                                        nextISI(1) = datasample(bankISI,1);
                                        readCount(1) = 0;
                                    end
                                    %end
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
                                    spmdSend(stim,intanWriter);
                                end
                            end
                        end
                    end
                end
                spmdSend([],intanWriter);
            end

        else

            noCon_Ard = spmdReceive(hardwareTimer,0);
            noCon_TCPwrite = spmdReceive(intanWriter,0);
            noCon_TCPread = spmdReceive(intanReader,0);

            spmdBarrier
        end
    end
end
close(gcf)

%%
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

  %%
    save([seshPathName,'\recInfo.mat'],'recInfo')

    disp('Analysis complete');
end

if strcmp(copyAction,'Yes') == 1
    copyDir = [copyPathName,'\',Time];
    disp(['Copying data to ',copyPathName]);
    eval(['mkdir ',copyDir]);
    copyfile(seshPathName,copyDir);
end


delete(gcp('nocreate'))
disp('Finsished recording');

%%

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

function timeVec = getTime()
Now = clock;
timeVec = (Now(1,4)*3600)+(Now(1,5)*60)+(Now(1,6));
end

function showVideo(input)
global videoOut %#ok<GVMIS>

if videoOut.QueueLength < 2

    tic
    frames = input{1,1};
    fps = input{1,2};

    numCam = size(frames,4);
    clf
    for ii = 1:numCam
        subplot(1,numCam,ii)
        image(frames(:,:,:,ii));
    end
    drawnow
    while toc < 1/fps
        % Pause for frame
    end
end
end

function plotTD(data)
global dataFig %#ok<GVMIS>

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

% Function to remap coordinates when image size is changed
function [Out] = RemapPoint(Points,SizeIn,SizeOut,RotAngle,FixedPoint)

Out = zeros(size(Points,1),2);

for i = 1:size(Points,1)

    point = Points(i,:);
    if all(point > 0) && all(point ~= FixedPoint)
        Loc = zeros(SizeIn(1,1),SizeIn(1,2));
        Loc(point(1,2),point(1,1)) = 1;

        if any(SizeIn ~= SizeOut)
            Loc = imresize(Loc,SizeOut);
            [M,mLoc] = max(Loc,[],1);
            [~,mLoc2] = max(M);
            Loc(mLoc(1,mLoc2),mLoc2) = 1;
            Loc(Loc < 1) = 0;
        end
        if RotAngle > 0
            Loc = imrotate(Loc,RotAngle);
        end
        [Out(i,2),Out(i,1)] = find(Loc == 1);
    else
        Out(i,:) = 0;
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