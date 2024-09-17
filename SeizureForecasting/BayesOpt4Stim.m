delete(gcp('nocreate'))
clear all %#ok<CLALL>

% Testing schedule:

StartTime = [7,00]; % Time to begin recording in 24hr time, format = [hr,min]
nextDay = 0; % Start recording tomorrow?

optSched(1,1) = 30; % Pre stim baseline.
optSched(2,1) = 30; % Post stim observation.
optSched(3,1) = 0; % Post stim wash-out.
repsDay = 25; % Number of tests to do per day.
numDay = 5; % Number of days to run.

% Box assignments:

% Subject in box, 'N/A' if empty.
boxSubject{1,1} = 'EDS 2.3'; % Box 1, port A
boxSubject{2,1} = 'EDS 3.0'; % Box 2, port B
boxSubject{3,1} = 'EDS 4.1'; % Box 3, port C
boxSubject{4,1} = 'EDS 4.2'; % Box 4, port D

% Active TCP channels to stream, [] = unused.
dataChan{1,1} = 0:7; % Port A
dataChan{2,1} = 0:7; % Port B
dataChan{3,1} = 0:7; % Port C
dataChan{4,1} = 0:7; % Port D

% Bayesian optimization settings:

% Optimization parameters:
optParam(1).name = 'Amplitude';
optParam(1).min = 5;
optParam(1).step = 1;
optParam(1).max = 5;
optParam(2).name = 'Frequency';
optParam(2).min = 0.1;
optParam(2).step = 0.1;
optParam(2).max = 200;
optParam(3).name = 'Duration';
optParam(3).min = 1;
optParam(3).step = 1;
optParam(3).max = 10;

gridStep = 3; % Divisions of grid search

% Target for optimization, 1 = T/D, 2 = Model.
targetOpt(1,1) = 1; % Box 1
targetOpt(2,1) = 1; % Box 2
targetOpt(3,1) = 1; % Box 3
targetOpt(4,1) = 1; % Box 4

% Theta/Delta:

showData = 1; % Display T/D plots
inputTD = 5000; % Amount of time to input into T/D filter in ms
TD2Plot = 60; % Amount of time to plot T/D in s

% Input channel for T/D.
tdChan(1,1) = 5; % Port A
tdChan(2,1) = 5; % Port B
tdChan(3,1) = 4; % Port C
tdChan(4,1) = 4; % Port D

% Prediction:

usePrediction = 1;
showLabels = 1; % Display T/D plots
inputMP = 2000; % Amount of time to input into classifier in ms
MP2Plot = 60; % Amount of time to plot predictions in s

% Model for each box
modelPath{1,1} = 'R:\McKenzieLab\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';
modelPath{2,1} = 'R:\McKenzieLab\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';
modelPath{3,1} = 'R:\McKenzieLab\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';
modelPath{4,1} = 'R:\McKenzieLab\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';

% Intan TCP:
intanIP_RHS = '127.0.0.1'; % IP for RHS.
intanPort1_RHS = 5000; % Port for RHS commands.
intanPort2_RHS = 5001; % Port for RHS data output.

% Video Settings:
fps = 25;
RecVid = 1; % Record from webcams, 1 = yes, 0 = no
showVideo = 1; % Show video feed from cameras, 1 = yes, 0 = no

% Arduino Settings:
ArdPin{1,1} = 'D13'; % Arduino pin for LED.
ArdPin{2,1} = 'D8'; % Arduino pin for Timestamp.
LEDtime = 10; % Time between LED blinks is seconds

% Mis settings:
sF = 20000;
dSF = 1000;
saveData = 1;

% Initial stim settings:

% Ports
stimConfig(1).port = 'a';
stimConfig(2).port = 'b';
stimConfig(3).port = 'c';
stimConfig(4).port = 'd';

% F Key for trigger
stimConfig(1).fKey = 1;
stimConfig(2).fKey = 2;
stimConfig(3).fKey = 3;
stimConfig(4).fKey = 4;

% Channel, [nan,nan] = unused, [c1,nan] = monopolar, [c1,c2] = bipolar.
stimConfig(1).chan = [4,5];
stimConfig(2).chan = [4,5];
stimConfig(3).chan = [4,5];
stimConfig(4).chan = [4,5];

% Amplitudes in uA.
stimConfig(1).amp = optParam(1).min;
stimConfig(2).amp = optParam(1).min;
stimConfig(3).amp = optParam(1).min;
stimConfig(4).amp = optParam(1).min;

% Stim polarity, -1 = Cathodic first, 1 = Anodic first. Ex: [-1,1].
stimConfig(1).polarity = [-1,-1];
stimConfig(2).polarity = [-1,-1];
stimConfig(3).polarity = [-1,-1];
stimConfig(4).polarity = [-1,-1];

% Phase duration in uS.
stimConfig(1).duration = 200;
stimConfig(2).duration = 200;
stimConfig(3).duration = 200;
stimConfig(4).duration = 200;

% Number stim train of pulses
stimConfig(1).pulseNum = 1;
stimConfig(2).pulseNum = 1; 
stimConfig(3).pulseNum = 1; 
stimConfig(4).pulseNum = 1; 

% Train period in mS
stimConfig(1).period = 1000;
stimConfig(2).period = 1000;
stimConfig(3).period = 1000;
stimConfig(4).period = 1000;

% Active channels for RHS
activeChan(1,1:2,1) = [{'a-000'},{'PrL(L)'}];
activeChan(2,1:2,1) = [{'a-001'},{'PrL(R)'}];
activeChan(3,1:2,1) = [{'a-002'},{'AVT(L)'}];
activeChan(4,1:2,1) = [{'a-003'},{'BLA(R)'}];
activeChan(5,1:2,1) = [{'a-004'},{'CA1(L)'}];
activeChan(6,1:2,1) = [{'a-005'},{'CA1(R)'}];
activeChan(7,1:2,1) = [{'a-006'},{'LDT(L1)'}];
activeChan(8,1:2,1) = [{'a-007'},{'LDT(L2)'}];

activeChan(1,1:2,2) = [{'b-000'},{'PrL(L)'}];
activeChan(2,1:2,2) = [{'b-001'},{'PrL(R)'}];
activeChan(3,1:2,2) = [{'b-002'},{'AVT(L)'}];
activeChan(4,1:2,2) = [{'b-003'},{'BLA(R)'}];
activeChan(5,1:2,2) = [{'b-004'},{'CA1(L)'}];
activeChan(6,1:2,2) = [{'b-005'},{'CA1(R)'}];
activeChan(7,1:2,2) = [{'b-006'},{'LDT(L1)'}];
activeChan(8,1:2,2) = [{'b-007'},{'LDT(L2)'}];

activeChan(1,1:2,3) = [{'c-000'},{'PrL(L)'}];
activeChan(2,1:2,3) = [{'c-001'},{'PrL(R)'}];
activeChan(3,1:2,3) = [{'c-002'},{'AVT(L)'}];
activeChan(4,1:2,3) = [{'c-003'},{'BLA(R)'}];
activeChan(5,1:2,3) = [{'c-004'},{'CA1(L)'}];
activeChan(6,1:2,3) = [{'c-005'},{'CA1(R)'}];
activeChan(7,1:2,3) = [{'c-006'},{'LDT(L1)'}];
activeChan(8,1:2,3) = [{'c-007'},{'LDT(L2)'}];

activeChan(1,1:2,4) = [{'d-000'},{'PrL(L)'}];
activeChan(2,1:2,4) = [{'d-001'},{'PrL(R)'}];
activeChan(3,1:2,4) = [{'d-002'},{'AVT(L)'}];
activeChan(4,1:2,4) = [{'d-003'},{'BLA(R)'}];
activeChan(5,1:2,4) = [{'d-004'},{'CA1(L)'}];
activeChan(6,1:2,4) = [{'d-005'},{'CA1(R)'}];
activeChan(7,1:2,4) = [{'d-006'},{'LDT(L1)'}];
activeChan(8,1:2,4) = [{'d-007'},{'LDT(L2)'}];
%%

totBox = length(boxSubject);
idxBox = 1:totBox;
for i = 1:length(boxSubject)
    if strcmp(boxSubject{i,1},'N/A')
        idxBox(1,i) = 0;
    end
end
idxBox(idxBox == 0) = [];
numBox = length(idxBox);
numChan = sum(cellfun(@length,dataChan));

numChan_Active = 0;
for i = 1:size(activeChan,3)
    for ii = 1:size(activeChan,1)
        if activeChan{ii,3,i} == 1
            numChan_Active = numChan_Active+1;
        end
    end
end

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
    while 1

        prompt = 'Label WebCam, hit Cancel to skip';
        dlgtitle = 'Input';
        camLabel = inputdlg(prompt,dlgtitle,[1,40],WL(i,1));
        if ~isempty(camLabel)
            camLabel = str2double(camLabel);
            if ~isnan(camLabel)
                WL{i,2} = camLabel;
                break
            else
                waitfor(warndlg('Enter the number of the sticker on the top-left of the cage','Warning'));
            end
        else
            WL{i,2} = 'N/A';
            break
        end
    end
    closePreview(cam)
end

if size(WL,1) > 0

    for i = 1:size(WL,1)

        if strcmp(WL{i,2},'N/A') == 0
            WL{i,3} = char(boxSubject{WL{i,2},1}+" ("+WL{i,2}+")");
        else
            WL{i,3} = 'N/A';
        end
    end
    
    [camChoice,~] = listdlg('PromptString','Select cameras.','ListString',WL(:,3),'SelectionMode','multiple');
    numCam = length(camChoice);
    WLc = WL(camChoice,:);
    [~,iWL] = sortrows(WLc(:,2));
    iWL = iWL';
    
    Res = cam.AvailableResolutions;
    [resChoice,~] = listdlg('PromptString','Select resolution.','ListString',Res,'SelectionMode','single');
    Dim = symsepchar(Res{1,resChoice},'x');
    Dim = [str2double(Dim{1,2}),str2double(Dim{1,1})];
    clear cam
end

port = {'a','b','c','d'};

for i = 1:size(activeChan,3)
    activeChan(:,3,i) = {0};
    [chanChoice,~] = listdlg('PromptString',{'Select active channels for';['port ',port{1,i},'.']},'ListString',activeChan(:,1,i),'SelectionMode','multiple');
    activeChan(chanChoice,3,i) = {1};
end

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

delete(gcp('nocreate'))
clear all
load('G:\RecordingData\test\baysOpt\config.mat')
clc

% CPU core roles
readTimer = 1;
intanWriter = 2;
thetaCalc = 3;
deltaCalc = 4;
predictModel = 5:8;
bayesOptimizer = 9:12;
dataSaver = 13;
camReader = 14;
camSaver = 15;

% Send tags
coreCom = 1;
startSig = 2;
timerSig = 3;
recPhase = 4;
bayesSig = 5;
ampData = 6;
videoData = 7;

if showVideo == 1
    global videoOut %#ok<GVMIS,TLEV>
    videoOut = parallel.pool.DataQueue;
    afterEach(videoOut,@plotFrame);
end

textOut = parallel.pool.DataQueue;
afterEach(textOut,@disp);

waitShow = parallel.pool.DataQueue;
afterEach(waitShow,@textWaitBar);

dataOut = parallel.pool.DataQueue;
afterEach(dataOut,@plotData);

labelOut = parallel.pool.DataQueue;
afterEach(labelOut,@plotLabel);

fileCount = 0;
vStartTime = (StartTime(1,1)*3600)+(StartTime(1,2)*60);
for Iter = 1:numDay    
    
    Now = clock; %#ok<CLOCK>
    dateToday = [num2str(Now(1,2)),'-',num2str(Now(1,3)),'-',num2str(Now(1,1))];

    dayPathName = [PathName,'\Day',num2str(Iter),'(',dateToday,')'];
    eval(['mkdir ',dayPathName]);    

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

    for I = 1:repsDay

        Now = clock; %#ok<CLOCK>
        Time = ['(',num2str(Now(1,4)),'.',num2str(Now(1,5)),')'];
        
        fileCount = fileCount+1;
        seshPathName = [dayPathName,'\Stim',num2str(fileCount),Time];

        eval(['mkdir ',seshPathName]);
        send(textOut,['Directory set to: ',seshPathName])

        if isempty(gcp('nocreate')) == 1
            parpool('local',15);
            disp('Connecting to hardware...');
        end
        
        spmd(15)
    
            mpiSettings('DeadlockDetection','off');    
     
            if spmdIndex == readTimer

                if I == 1
    
                    comIssue = cell(6,2);
                    comIssue{1,1} = 'RHS Data';
                    comIssue{2,1} = 'Arduino';
                    comIssue{3,1} = 'intanWrite';
                    comIssue{4,1} = 'dataSaver';
                    comIssue{5,1} = 'camReader';
                    comIssue{6,1} = 'camSaver';
                    comIssue(:,2) = {0};
                    
                    send(textOut,'Checking hardware');
    
                    % Connect to RHS data output and test connection
                    t = 0;
                    while 1
                        try
                            intanRead_RHS = tcpclient(intanIP_RHS,intanPort2_RHS);
                            flush(intanRead_RHS)
                            send(textOut,'Connected to Intan RHS data output')
                            break
                        catch
                            try
                                flush(intanRead_RHS)
                                send(textOut,'Connected to Intan RHS data output')
                                break
                            catch
                                send(textOut,'Issue connecting to Intan data output')
                                pause(1)
                                t = t+1;
                                if t == 60
                                    send(textOut,'Connection to Intan data output timed out')
                                    comIssue{1,2} = 1;
                                    break
                                end
                            end
                        end
                    end                    
    
                    % Connect to Arduino and test connection
                    t = 0;
                    while 1
                        try
                            Ard = arduino('com3','Uno');
                            send(textOut,'Connected to Arduino')
                            break
                        catch
                            try
                                writeDigitalPin(Ard,ArdPin{1,1},1);
                                writeDigitalPin(Ard,ArdPin{1,1},0);
                                send(textOut,'Connected to Arduino')
                                break
                            catch
    
                                send(textOut,'Issue connecting to Arduino')
                                pause(1)
                                t = t+1;
                                if t == 60
                                    send(textOut,'Connection to Arduino timed out')
                                    comIssue{2,2} = 1;
                                    break
                                end
                            end
                        end
                    end
                    
                    % Make sure all cores started correctly
                    comIssue{3,2} = spmdReceive(intanWriter,coreCom);
                    comIssue{4,2} = spmdReceive(dataSaver,coreCom);
                    comIssue{5,2} = spmdReceive(camReader,coreCom);
                    comIssue{6,2} = spmdReceive(camSaver,coreCom);
                    allClear = sum(cell2mat(comIssue(:,2)));

                    % Report results of hardware check
                    if allClear == 0
                        send(textOut,'Hardware check passed');
                    else
                        send(textOut,'Problem detected with:');
                        for i = 1:size(comIssue,1)
                            if comIssue{i,2} == 1
                                send(textOut,comIssue{i,1});
                            end
                        end
                    end
                    send(textOut,'----------')
                    spmdSend(allClear,[intanWriter,dataSaver,camReader,camSaver],coreCom);
                end
    
                if allClear == 0                    

                    % Create cell to hold all block titles
                    blockText = cell(4,1);
                    blockText{1,1} = 'Baseline';
                    blockText{2,1} = 'Stim';
                    blockText{3,1} = 'Observation';
                    blockText{4,1} = 'Wash Out';
            
                    % Calculations for accurate parsing
                    waveformBytesPerFrame = 4+2*numChan_Active;
                    waveformBytesPerBlock = 128*waveformBytesPerFrame+4;

                    % Pre-allocate counters and data stores
                    timeWindow = zeros(size(activeChan,1),inputTD,numBox);
                    tdSend = zeros(length(tdChan),inputTD);
                    labels = zeros(numBox,1);
                    plotTD = zeros(length(tdChan),TD2Plot);
                    plotLabels = zeros(numBox,MP2Plot);
                    stimOver = zeros(numBox,1);
                    obsTD = zeros(length(tdChan),optSched(2,1));
                    obsLabels = zeros(numBox,optSched(2,1));
                    obsCount = ones(numBox,1);
                    obsSent = zeros(numBox,1);
                    
                    % Warm up calc cores and open plots
                    spmdSend(tdSend,[thetaCalc,deltaCalc],ampData);                    

                    if showData == 1
                        send(dataOut,{plotTD,WLc,iWL})
                    end

                    if showLabels == 1
                        send(labelOut,{plotLabels,WLc,iWL})
                    end

                    theta = spmdReceive(thetaCalc,ampData);
                    delta = spmdReceive(deltaCalc,ampData);

                    flush(intanRead_RHS);

                    % Wait for stim upload on intanWriter
                    allDur = spmdReceive(intanWriter,startSig);
                    obsDur = optSched(2,1);
                    timerShed = optSched;
                    timerShed(2,1) = obsDur+ceil(max(allDur));
                    
                    % sync with system clock
                    send(textOut,'Syncing hardware');
                    if I == 1
                        while getTime < vStartTime
                            pause(0.001);
                        end
                    else
                        while getTime/round(getTime) ~= 1
                            pause(0.001);
                        end
                    end                    
                    spmdSend(1,[intanWriter,camReader],startSig)
                    
                    for i = 1:length(timerShed)
                        
                        % Tell intanWriter current recoring phase
                        if i == 1 || i == 3
                            spmdSend(1,intanWriter,recPhase)
                        elseif i == 2
                            spmdSend(2,intanWriter,recPhase)
                        end
                            
                        for ii = 1:timerShed(i,1)

                            % Wait for data
                            while intanRead_RHS.NumBytesAvailable < waveformBytesPerBlock
                                pause(0.0001)
                            end
                            
                            t = tic;

                            % Send sync signal to camReader
                            spmdSend(t,camReader,timerSig);

                            if rem(ii,LEDtime) == 1

                                t2 = tic;

                                % Turn on block marker
                                if ii == 1
                                    writeDigitalPin(Ard,ArdPin{2,1},1);
                                end

                                % Turn on LED
                                writeDigitalPin(Ard,ArdPin{1,1},1);

                                if ii == 1
                                    % Display current block
                                    send(textOut,'----------')
                                    send(textOut,blockText{i,1})
                                end
    
                                % Pause for 0.2s
                                while toc(t2) < 0.2
                                    pause(0.0001)
                                end
                            
                                % Turn off LED
                                writeDigitalPin(Ard,ArdPin{1,1},0);
                                
                                % Turn off block marker
                                if ii == 1
                                    writeDigitalPin(Ard,ArdPin{2,1},0);
                                end
                            end

                            % Wait till 0.5s has passed
                            while toc(t) < 0.4999
                                pause(0.0001)
                            end
                            
                            % Find total number of blocks ready to read
                            blocksToRead = floor(intanRead_RHS.NumBytesAvailable/waveformBytesPerBlock);
                            bytesToRead = waveformBytesPerBlock*blocksToRead;
                            
                            % Read data and decode
                            waveformArray = read(intanRead_RHS,bytesToRead);
                            [amplifierData,offset] = byte2double2(waveformArray,numChan_Active,blocksToRead,sF,dSF);

                            % If read cycle gets offset find the offset and correct
                            if offset > 0
                                send(textOut,['Read offset of ',num2str(offset),' detected.'])
                                waveformArray(1:offset) = [];
                                waveformArray = [waveformArray,read(intanRead_RHS,offset)];  %#ok<AGROW>
                                [amplifierData,offset] = byte2double2(waveformArray,numChan_Active,blocksToRead,sF,dSF);
                                send(textOut,'Corrected')
                            end
                            
                            % Parse channels into groups by box
                            readWindow = zeros(size(activeChan,1),size(amplifierData,2),numBox);
                            iChan = 1;
                            for i3 = 1:numBox
                                for i4 = 1:size(activeChan,1)
                                    if activeChan{i4,3,i3} == 1 && any(dataChan{i3,1} == i4-1)
                                        readWindow(i4,:,i3) = amplifierData(iChan,:);
                                        iChan = iChan+1;
                                    end
                                end
                            end
                            
                            % Update time window
                            timeWindow = [timeWindow,readWindow]; %#ok<AGROW>
                            timeWindow(:,1:size(amplifierData,2),:) = [];
                            
                            % Pull channels for T/D
                            for i3 = 1:numBox
                                tdSend(idxBox(1,i3),:) = timeWindow(tdChan(i3,1)+1,:,i3);
                            end                            

                            % Send timeWindow to calc cores for parallel computation
                            spmdSend(tdSend,[thetaCalc,deltaCalc],ampData);                            
                            for i3 = 1:numBox
                                spmdSend(timeWindow(:,end-inputMP+1:end,i3),predictModel(1,idxBox(1,i3)),ampData);
                            end                            
                            
                            % Get T/D ratio
                            theta = spmdReceive(thetaCalc,ampData);
                            delta = spmdReceive(deltaCalc,ampData);
                            ToD = theta./delta;
                            
                            % Get labels                           
                            for i3 = 1:numBox
                                modelOut = spmdReceive('any',ampData);
                                labels(modelOut(1,1),1) = modelOut(1,2);
                            end

                            % Update T/D plots
                            if showData == 1
                                plotTD(:,end+1) = ToD; %#ok<SAGROW>
                                plotTD(:,1) = [];
                                send(dataOut,{plotTD,WLc,iWL})
                            end
                            
                            % Update label plots
                            if showLabels == 1
                                plotLabels = [plotLabels,labels]; %#ok<AGROW>
                                plotLabels(:,1) = [];
                                send(labelOut,{plotLabels,WLc,iWL})
                            end
                            
                            % Add data to observations once stim finishes
                            if i == 2
                                for i3 = 1:numBox
                                    if ii > allDur(i3,1) && obsCount(i3,1) < obsDur

                                        if stimOver(i3,1) == 0
                                            send(textOut,['Box ',num2str(idxBox(1,i3)),' stim complete.'])
                                            stimOver(i3,1) = 1;
                                        end

                                        if targetOpt(i3,1) == 1
                                            obsTD(i3,obsCount(i3,1)) = ToD(i3,1);
                                        elseif targetOpt(i3,1) == 2
                                            obsLabels(i3,obsCount(i3,1)) = labels(i3,1);
                                        end
                                        obsCount(i3,1) = obsCount(i3,1)+1;
                                    elseif ii >= allDur(i3,1)+obsDur
                                        if obsSent(i3,1) == 0

                                            send(textOut,['Box ',num2str(idxBox(1,i3)),' observation complete.'])

                                            if targetOpt(i3,1) == 1
                                                spmdSend({obsTD(i3,:)},bayesOptimizer(1,idxBox(1,i3)),bayesSig)
                                            elseif targetOpt(i3,1) == 2
                                                spmdSend({obsLabels(i3,:)},bayesOptimizer(1,idxBox(1,i3)),bayesSig)
                                            end
                                            obsSent(i3,1) = 1;
                                        end
                                    end
                                end
                            end
                            
                            % Send to datasaver core if enabled
                            if saveData == 1
                                if usePrediction == 0
                                    spmdSend(ToD,dataSaver,ampData);
                                elseif usePrediction == 1
                                    spmdSend([ToD,labels],dataSaver,ampData);
                                end
                            end

                            % Wait till 1s has passed
                            while toc(t) < 0.9999
                                pause(0.0001)
                            end
                        end
                    end
                    
                    spmdSend([],camReader,timerSig)
                    spmdSend([],intanWriter,recPhase)
                    flush(intanRead_RHS)
                    send(textOut,'----------')
                    send(textOut,'Recording finished.')

                    % Send stop to datasaver core if enabled
                    if saveData == 1
                        spmdSend([],dataSaver,ampData);
                    end
                end
                spmdSend([],[thetaCalc,deltaCalc,predictModel],ampData);
                spmdSend([],bayesOptimizer,bayesSig);
            elseif spmdIndex == intanWriter
                
                if I == 1

                    comIssue = 0;
                    % Connect and test connection to RHS commands                
                    t = 0;
                    while 1
                        try
                            intanSend_RHS = tcpclient(intanIP_RHS,intanPort1_RHS);
                            write(intanSend_RHS,uint8('execute RescanPorts;'));
                            pause(0.1)
                            send(textOut,'Connected to Intan RHS commands')
                            break
                        catch
                            try
                                write(intanSend_RHS,uint8('execute RescanPorts;'));
                                send(textOut,'Connected to Intan RHS commands')
                                break
                            catch
                                send(textOut,'Issue connecting to Intan RHS commands')
                                pause(1)
                                t = t+1;
                                if t == 60
                                    send(textOut,'Connection to Intan timed out')
                                    comIssue = 1;
                                    break
                                end
                            end
                        end
                    end
                    
                    spmdSend(comIssue,readTimer,coreCom)
                    allClear = spmdReceive(readTimer,coreCom);
                end

                if allClear == 0

                    % Set up RHS channel parameters
                    if intanSend_RHS ~= 0
    
                        write(intanSend_RHS,uint8(['set Filename.Path ',seshPathName,';']));
                        write(intanSend_RHS,uint8('set FileFormat OneFilePerSignalType;'));
                        write(intanSend_RHS,uint8('set Filename.BaseFilename RHS;'));
    
                        % Enable TCP recording for data channels
                        for i = 1:size(dataChan,1)
                            if ~isempty(dataChan{i,1})
                                for ii = 1:length(dataChan{i,1})
                                    if dataChan{i,1}(1,ii) < 10
                                        write(intanSend_RHS,uint8(['set ',port{1,i},'-00',num2str(dataChan{i,1}(1,ii)),'.tcpdataoutputenabled true;']));
                                    elseif dataChan{i,1}(1,ii) > 9 
                                        write(intanSend_RHS,uint8(['set ',port{1,i},'-0',num2str(dataChan{i,1}(1,ii)),'.tcpdataoutputenabled true;']));
                                    end
                                end
                            end
                        end
                    end

                    configStruct = stimConfig;

                    % Ask bayesOptimizer cores to pick new stim settings
                    for i = idxBox
                        spmdSend(1,bayesOptimizer(1,i),bayesSig);
                    end
                    
                    phase = 1;
                    stimStart = 1;
                    ticSec = 10000000;
                    allDur = zeros(numBox,1);
                    ticStart = uint64(zeros(4,1));
                    ticStop = uint64(zeros(4,1));
                    clockStop = zeros(4,1);
                    numStim = zeros(4,1);
                    
                    % Recieve data from bayesOptimizer cores
                    nextParam = cell(totBox,1);
                    for i = 1:numBox
                        baysRec = spmdReceive('any',bayesSig);
                        nextParam{baysRec(1,1),1} = baysRec(1,2:end);
                    end
                    
                    for i = idxBox

                        configStruct(idxBox(1,i)).amp = nextParam{i,1}(1,1);
                             
                        if nextParam{i,1}(1,2) <= 5

                            configStruct(idxBox(1,i)).pulseNum = 1;
                            configStruct(idxBox(1,i)).period = 930;
                            clockStop(i,1) = 1/nextParam{i,1}(1,2);
                            numStim(i,1) = floor(nextParam{i,1}(1,3)/(clockStop(i,1)));
                            if numStim(i,1) == 0
                                clockStop(i,1) = 1/nextParam{i,1}(1,2);
                                numStim(i,1) = 1;
                            end
                        else
                            pulseNum = round(nextParam{i,1}(1,2));
                            period = (1/nextParam{i,1}(1,2))*1000;

                            while pulseNum*period > 930
                                pulseNum = pulseNum-1;
                            end

                            configStruct(idxBox(1,i)).pulseNum = round(nextParam{i,1}(1,2));
                            configStruct(idxBox(1,i)).period = (1/nextParam{i,1}(1,2))*1000;
                            clockStop(i,1) = 1;
                            numStim(i,1) = nextParam{i,1}(1,3);
                        end
                        allDur(i,1) = nextParam{i,1}(1,3);

                        send(textOut,['Box ',num2str(i),' -> Amplitude: ',num2str(nextParam{i,1}(1,1)),'uA, Frequency: ',num2str(nextParam{i,1}(1,2)),'Hz, Duration: ',num2str(nextParam{i,1}(1,3)),'s.'])
                    end
                    
                    % Set up stim parameters
                    setStim(intanSend_RHS,1,configStruct)
                    write(intanSend_RHS,uint8('execute uploadstimparameters;'));
                    tic

                    % Wait for stim upload to complete
                    send(textOut,'Uploading stim parameters')
                    send(waitShow,14)
                    pause(15)
                    
                    % Send durations to readTimer and let it know the
                    % upload is complete
                    spmdSend(allDur,readTimer,startSig)

                    % Wait for sync signal from readTimer
                    spmdReceive(readTimer,startSig);
    
                    % Start recording                    
                    write(intanSend_RHS,uint8('set runMode record;'));
                    
                    while 1

                        if spmdProbe(readTimer,recPhase) == 1
                            phase = spmdReceive(readTimer,recPhase);

                            if isempty(phase)
                                break
                            end
                        end

                        if phase == 2

                            commandArray = '';
                            if stimStart == 1

                                ticStart(:,1) = tic;
                                for i = idxBox
                                    ticStop(i,1) = ticStart(i,:)+(clockStop(i,1)*ticSec);
                                    commandArray = [commandArray,['execute ManualStimTriggerPulse F',num2str(configStruct(i).fKey),';']]; %#ok<AGROW>
                                    numStim(i,1) = numStim(i,1)-1;
                                end
                                stimStart = 0;
                            elseif stimStart == 0
                                for i = idxBox
                                    if tic >= ticStop(i,1) && numStim(i,1) > 0
                                        ticStart(i,1) = tic;
                                        ticStop(i,1) = ticStart(i,1)+(clockStop(i,1)*ticSec);
                                        commandArray = [commandArray,['execute ManualStimTriggerPulse F',num2str(configStruct(i).fKey),';']]; %#ok<AGROW>
                                        numStim(i,1) = numStim(i,1)-1;
                                    end
                                end
                            end   

                            % Execute f key presses to Intan
                            if ~isempty(commandArray)
                                write(intanSend_RHS,uint8(commandArray));
                            end
                        else
                            pause(0.001)
                        end
                    end    
                    
                    % Stop recording
                    write(intanSend_RHS,uint8('set runMode stop;'));
                    
                    % Disable stim and TCP channels after recording day finishes
                    if I == repsDay    
                        
                        setStim(intanSend_RHS,0,configStruct)                        
    
                        for i = 1:size(dataChan,1)
                            if ~isempty(dataChan{i,1})
                                for ii = 1:length(dataChan{i,1})
                                    if dataChan{i,1}(1,ii) < 10
                                        write(intanSend_RHS,uint8(['set ',port{1,i},'-00',num2str(dataChan{i,1}(1,ii)),'.tcpdataoutputenabled false;']));
                                    elseif dataChan{i,1}(1,ii) > 9 
                                        write(intanSend_RHS,uint8(['set ',port{1,i},'-0',num2str(dataChan{i,1}(1,ii)),'.tcpdataoutputenabled false;']));
                                    end
                                end
                            end
                        end
    
                        write(intanSend_RHS,uint8('execute uploadstimparameters;'));
                    end
                end
            elseif spmdIndex == thetaCalc                
                while 1
                    try
                        if spmdProbe(readTimer,ampData) == 1
    
                            timeWindow = spmdReceive(readTimer,ampData);
    
                            if isempty(timeWindow)
                                break
                            end
    
                            % Calculate mean theta amp
                            theta = abs(hilbert(bandpass(timeWindow',[5,12],dSF)));
                            theta = mean(theta(end-inputTD+1:end,:))';
        
                            spmdSend(theta,readTimer,ampData)
                        else
                            pause(0.001)
                        end
                    catch
                        spmdSend(0,readTimer,ampData)
                        send(textOut,'Unknown packet')
                    end
                end
            elseif spmdIndex == deltaCalc
                while 1
                    try
                        if spmdProbe(readTimer,ampData) == 1
    
                            timeWindow = spmdReceive(readTimer,ampData);
    
                            if isempty(timeWindow)
                                break
                            end
    
                            % Calculate mean delta amp
                            delta = abs(hilbert(bandpass(timeWindow',[1,4],dSF)));
                            delta = mean(delta(end-inputTD+1:end,:))';
        
                            spmdSend(delta,readTimer,ampData)
    
                        else
                            pause(0.001)
                        end
                    catch
                        spmdSend(0,readTimer,ampData)
                        send(textOut,'Unknown packet')
                    end
                end
            elseif any(spmdIndex == predictModel)

                boxID = find(predictModel == spmdIndex);

                model = load(modelPath{boxID,1});
                nChannel = model.ops.nCh_raw;
                freq = model.ops.freqs;
                nFreq = length(freq);
                id = spmdIndex-predictModel(1,1)+1;
                feat = nan(1,nFreq*nChannel);

                while 1                    
                    if spmdProbe(readTimer,ampData) == 1

                        data = spmdReceive(readTimer,ampData);

                        if isempty(data)
                            break
                        end

                        data = data/0.195;
                        if ~any(any(abs(data) > model.ops.art_thres))   
                            for i = 1:nChannel
                                
                                % wavelet decomposition
                                tmp = abs(sm_wavelet(data(i,:)',dSF,freq))';
                                
                                % loop over the frequencies
                                for ii = 1:nFreq
                                    
                                    % for the sampling window, take the mean power at each frequency.
                                    tmp = mean(tmp,2);
                                    feat((ii-1)*nChannel+i) = tmp(ii);
                                end
                            end
                            label = predict(model.rusTree,feat);
                            spmdSend([id,label],readTimer,ampData)
                        else
                            spmdSend([],readTimer,ampData)
                        end
                    else
                        pause(0.001)
                    end
                end                
            elseif any(spmdIndex == bayesOptimizer)

                boxID = find(bayesOptimizer == spmdIndex);

                if I == 1
                    
                    if targetOpt(boxID,1) == 1

                        fileName = [PathName,'\Box',num2str(boxID),'_gridData_TD.mat'];
                        bins = [0.5,1.5,2.5,3.5,4.5];
                        scoring = zeros(5,1);
                        scoring(1,1) = 1000;
                        scoring(2,1) = 100;
                        scoring(3,1) = 50;
                        scoring(4,1) = 10;
                        scoring(5,1) = 1;
                    elseif targetOpt(boxID,1) == 2

                        fileName = [PathName,'\Box',num2str(boxID),'_gridData_MP.mat'];
                        scoring = zeros(6,1);
                        scoring(1,1) = 1;
                        scoring(2,1) = 10;
                        scoring(3,1) = 50;
                        scoring(4,1) = 100;
                        scoring(5,1) = 1000;
                        scoring(6,1) = 1000;
                    end

                    gridData = matfile(fileName,'Writable',true);
                    if isprop(gridData,'parameters') == 0
                        
                        % Load initial parameters
                        temp = optParam;
                        numParam = size(temp,2);

                        % Get grid search coordinates
                        idx = cell(1,numParam);
                        for ii = 1:numParam
    
                            m = optParam(ii).min;
                            s = optParam(ii).step;
                            M = optParam(ii).max;
                            
                            temp(ii).range = m:s:M;

                            l = length(temp(ii).range);
                            idx{1,ii} = zeros(1,gridStep);
                            idx{1,ii}(1,1:gridStep-1) = round(1:l/(gridStep-1):l);
                            idx{1,ii}(1,end) = l;
                            idx{1,ii} = unique(temp(ii).range(1,idx{1,ii}));
                        end
                        gridData.parameters = temp;
                        
                        % Get initial grid search 
                        numCells = length(idx);
                        varNumel = cellfun(@length,idx);
                        numComb = prod(varNumel);
                        idc = zeros(numComb,numCells);
                        for ii = 1:numCells
                            var = idx{ii}(:);
                            r = repelem(var,prod(varNumel(ii+1:end)),1);
                            idc(:,ii) = repmat(r,prod(varNumel(1:ii-1)),1);
                        end
                        gridData.search = idc;
                        gridData.obs = zeros(0,1);
                    end
                end

                count = 1;
                useParm = 1:size(optParam,2);
                for i = 1:size(optParam,2)
                    if optParam(i).min < optParam(i).max
                        optVar(count).Var = optimizableVariable(optParam(i).name,[optParam(i).min,optParam(i).max]); %#ok<SAGROW>
                        count = count+1;
                    else
                        useParm(1,i) = 0;
                    end
                end
                useParm(useParm == 0) = [];

                while 1
                    if spmdProbe('any',bayesSig) == 1

                        stimOpt = spmdReceive('any',bayesSig);

                        if isempty(stimOpt)
                            break
                        end

                        if ~iscell(stimOpt)

                            numObs = length(gridData.obs);
                            if numObs >= size(gridData.search,1)

                                fun = defFun({gridData.search(1:numObs,useParm),gridData.obs(1:numObs,1)});
                                results = bayesopt(fun,[optVar.Var],...
                                    'Verbose',0,...
                                    'PlotFcn',[],...
                                    'AcquisitionFunctionName','expected-improvement-plus',...
                                    'NumSeedPoints',10,...
                                    'ExplorationRatio',0.5,...
                                    'IsObjectiveDeterministic',false,...
                                    'UseParallel',false,...
                                    'MaxObjectiveEvaluations',15);
                                nextParam = table2array(results.NextPoint);
                                
                                count = 1;
                                for i = 1:size(gridData.search,2)
                                    if any(i == useParm)
                                        gridData.search(numObs+1,i) = nextParam(1,count);
                                        count = count+1;
                                    else
                                        gridData.search(numObs+1,i) = gridData.search(numObs,i);
                                    end
                                end
                            end
                            nextParam = gridData.search(numObs+1,:);
                            spmdSend([boxID,nextParam],intanWriter,bayesSig)
                        else
                            data = stimOpt{1,1};
                            for i = 1:length(data)
                                if targetOpt(boxID,1) == 1
                                    [counts,~] = hist(data(1,i),bins); %#ok<HIST>
                                    data(1,i) = scoring(logical(counts),1);
                                elseif targetOpt(boxID,1) == 2
                                    data(1,i) = scoring(data(1,i),1);
                                end
                            end
                            gridData.obs = [gridData.obs;mean(data)];
                        end
                    else
                        pause(0.001)
                    end
                end               
            elseif spmdIndex == dataSaver
    
                % Send check to hardwareTimer
                comIssue = 0;
                try
                    if saveData == 1
                        fID = fopen([seshPathName,'\ToD_',num2str(length(tdChan)),'Ch_1Hz.dat'],'w');
                        if usePrediction == 1
                            fID2 = fopen([seshPathName,'\Labels_',num2str(length(tdChan)),'Ch_1Hz.dat'],'w');
                        end
                    end
                catch
                    comIssue = 1;
                end

                if I == 1
                    spmdSend(comIssue,readTimer,coreCom)
                    comIssue = spmdReceive(readTimer,coreCom);
                end
                
                if comIssue == 0
                    while 1
                        if spmdProbe(readTimer,ampData) == 1

                            data = spmdReceive(readTimer,ampData);
                            if ~isempty(data)                    
                                fwrite(fID,data(:,1)*100,'int16');
                                if usePrediction == 1
                                    fwrite(fID2,data(:,2)*100,'int16');
                                end
                            else
                                break
                            end
                        else
                            pause(0.001)
                        end
                    end
    
                    if saveData == 1
                        fclose(fID);
                    end
                end
            elseif spmdIndex == camReader
                
                if I == 1

                    comIssue = 0;
                    cam = cell(numCam,1); 
                    for ii = 1:numCam
                        try
                            % Create all webcam objects
                            cam{ii,1} = webcam(camChoice(1,ii));
                        catch
                            comIssue = 1;
                            send(textOut,['Could not connect to webcam ', num2str(ii)])
                        end
                    end
                    spmdSend(comIssue,readTimer,coreCom)
                    comIssue = spmdReceive(readTimer,coreCom);
                end
    
                if comIssue == 0
                    
                    ck = zeros(1,numCam);
                    frame = uint8(zeros(Dim(1,1),Dim(1,2),3,numCam)); 
                    errorCode = uint8([0,255,0;255,0,255;0,255,0]);

                    if showVideo == 1                        
                        send(videoOut,{frame+256,WLc,iWL}); %#ok<SPGV>                        
                    end
    
                    spmdReceive(readTimer,startSig);

                    while 1
                        if spmdProbe(readTimer,timerSig) == 1

                            % wait for timer signal
                            t = spmdReceive(readTimer,timerSig);

                            if ~isempty(t)
                                while toc(t) < 1
                                    for i3 = 1:numCam
                                        try
                                            frame(:,:,:,i3) = snapshot(cam{i3,1});
                                            ck(1,i3) = 0;
                                        catch
                                            try
                                                send(textOut,['Could not get frame from cam ',num2str(i3)])
                                                send(textOut,'Atempting to reconnect...')
                                                cam{i3,1} = webcam(camChoice(1,i3));
                                                frame(:,:,:,i3) = snapshot(cam{i3,1});
                                            catch
                                                if ck(1,i3) < 6
                                                    send(textOut,['Failed to connect to cam ',num2str(i3)])
                                                    ck(1,i3) = ck(1,i3)+1;
                                                else 
                                                    send(textOut,['Connection to cam ',num2str(i3),' lost.'])
                                                end
                                                frame(:,:,:,i3) = uint8(zeros(Dim(1,1),Dim(1,2),3)); 
                                                frame(1:3,1:3,1,i3) = errorCode;
                                            end
                                        end
                                    end
                                    spmdSend(frame,camSaver,videoData);
                                end
                                spmdSend(0,camSaver,videoData);
                            else
                                break
                            end
                        else
                            pause(0.0001)
                        end
                    end
                    spmdSend([],camSaver,videoData);
                end
            elseif spmdIndex == camSaver    
                
                comIssue = 0;
                if RecVid == 1
                    try
                        % Create all video writer objects
                        VidFileName = cell(numCam,1);
                        writerObj = cell(numCam,1);
                        errorCode = uint8([0,255,0;255,0,255;0,255,0]);
            
                        for ii = 1:numCam
                            VidFileName{ii,1} = [WL{camChoice(1,ii),3}(1,1:end-4),'.mp4'];
                            writerObj{ii,1} = VideoWriter([seshPathName,'\',VidFileName{ii,1}],'MPEG-4'); %#ok<TNMLP>
                            writerObj{ii,1}.FrameRate = fps;
                            open(writerObj{ii,1});
                        end
                    catch
                        send(textOut,'Could not create video files')
                        comIssue = 1;
                    end
                end

                if I == 1
                    spmdSend(comIssue,readTimer,coreCom)
                    comIssue = spmdReceive(readTimer,coreCom);
                end
                
                if comIssue == 0

                    frames = uint8(zeros(Dim(1,1),Dim(1,2),3,numCam,40));

                    while 1

                        numFrames = 0;                            
                        while 1
                            if spmdProbe(camReader,videoData) == 1

                                frame = spmdReceive(camReader,videoData);
                                if size(frame,1) > 1
    
                                    if showVideo == 1
                                        if videoOut.QueueLength < 10 %#ok<SPGV>
                                            send(videoOut,{frame,WLc,iWL}); %#ok<SPGV>
                                        end
                                    end
    
                                    numFrames = numFrames+1;
                                    frames(:,:,:,:,numFrames) = frame;
                                else
                                    break
                                end
                            else
                                pause(0.0001)
                            end
                        end

                        if isempty(frame)
                            break
                        end

                        syncFrames = frames(:,:,:,:,round(linspace(1,numFrames,fps)));

                        if RecVid == 1
                            for i3 = 1:fps
                                for i4 = 1:numCam
                                    if all(syncFrames(1:3,1:3,1,i4,i3) ~= errorCode,'all')
                                        writeVideo(writerObj{i4,1},syncFrames(:,:,:,i4,i3));
                                    end                   
                                end
                            end
                        end
                    end

                    send(textOut,'Finalizing video files.')
                    for ii = 1:numCam
                        close(writerObj{ii,1});
                    end
                    send(textOut,'----------')
                end
            end
        end
    end
end
send(textOut,'Complete.')

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
    Now = clock; %#ok<CLOCK>
    timeVec = (Now(1,4)*3600)+(Now(1,5)*60)+(Now(1,6));
end

function textWaitBar(dur)

    fprintf([repmat('.',1,dur),'\n\n']);
    for i = 1:dur
        fprintf('\b|\n');
        pause(1)
    end
end

function setStim(TCP,enable,config)

    % Set up a channel's stim parameters
    % TCP = Intan TCP object
    % enable = turn stim on or off, 0 = off, 1 = on
    % config = structure with fields:   
        % port = port on Intan: 'a','b','c','d'
        % fKey = F key to trigger stim
        % chan = channel number (base zero), [c1,c2] for bipolar
        % amp = Stim amplitude in uA
        % polarity = leading current, -1 = cathode first, 1 = anode first
        % duration = pulse duration in uS
        % pulseNum = number of pulses
        % period = pulse train period in mS
    
    for i = 1:size(config,2)
        for ii = 1:length(config(i).chan)
            if ~isnan(config(i).chan(ii))

                if config(i).chan(ii) < 10
                    chanStr = [config(i).port,'-00',num2str(config(i).chan(ii))];
                else 
                    chanStr = [config(i).port,'-0',num2str(config(i).chan(ii))];
                end

                if enable == 1              

                    write(TCP,uint8(['set ',chanStr,'.stimenabled true;']));
                    write(TCP,uint8(['set ',chanStr,'.source KeyPressF',num2str(config(i).fKey),';']));
                    write(TCP,uint8(['set ',chanStr,'.FirstPhaseAmplitudeMicroAmps ',num2str(config(i).amp),';']));
                    write(TCP,uint8(['set ',chanStr,'.SecondPhaseAmplitudeMicroAmps ',num2str(config(i).amp),';']));
    
                    if config(i).polarity(ii) == 1
                        write(TCP,uint8(['set ',chanStr,'.Polarity PositiveFirst;']));
                    elseif config(i).polarity(ii) == -1
                        write(TCP,uint8(['set ',chanStr,'.Polarity NegativeFirst;']));
                    end
    
                    write(TCP,uint8(['set ',chanStr,'.FirstPhaseDurationMicroseconds ',num2str(round(config(i).duration/2)),';']));
                    write(TCP,uint8(['set ',chanStr,'.SecondPhaseDurationMicroseconds ',num2str(round(config(i).duration/2)),';']));

                    if config(i).pulseNum == 1
                        write(TCP,uint8(['set ',chanStr,'.PulseOrTrain SinglePulse;']));
                    else
                        write(TCP,uint8(['set ',chanStr,'.PulseOrTrain PulseTrain;']));
                        write(TCP,uint8(['set ',chanStr,'.NumberOfStimPulses ',num2str(config(i).pulseNum),';']));
                        write(TCP,uint8(['set ',chanStr,'.PulseTrainPeriodMicroseconds ',num2str(config(i).period*1000),';']));
                    end
                elseif enable == 0
                    write(TCP,uint8(['set ',chanStr,'.stimenabled false;']));
                end
            end
        end
    end
end

function [amplifierData,offset] = byte2double(waveformArray,numChan,blocksPerRead,sampleFreq,downsampleFreq)
    
    % Pre-allocate memory
    amplifierData = 32768*ones(numChan,blocksPerRead*128);
    
    rawIndex = 1;
    offset = 0;
    while 1

        %Expect 4 bytes to be TCP Magic Number as uint32. 
        % If not what's expected, print that there was an error.
        Bytes = waveformArray(1,rawIndex:rawIndex+3);
        magicNumber = typecast(uint8(Bytes),'uint32');

        if magicNumber == 0x2ef07a08
            break
        else
            offset = offset+1;
            
            rawIndex = rawIndex+1;
        end
    end

    if offset == 0
    
        amplifierTimestampsIndex = 1;
        for iBlock = 1:blocksPerRead

            % Skip over magicNumber
            rawIndex = rawIndex+4;
    
            % Each block should contain 128 frames of data, process each one-by-one
            for iFrame = 1:128
    
                % Expect 4 bytes to be timestamp and skip these
                rawIndex = rawIndex+4;
    
                % Parse all bands of amplifier channels
                for iChannel = 1:numChan
    
                    % 2 bytes of wide            
                    Bytes = waveformArray(rawIndex:rawIndex+1);
                    amplifierData(iChannel,amplifierTimestampsIndex) = typecast(uint8(Bytes),'uint16');
                    rawIndex = rawIndex+2;
                end
                amplifierTimestampsIndex = amplifierTimestampsIndex+1;
            end
        end
    
        % Scale these data blocks and downsample
        amplifierData = amplifierData(:,1:sampleFreq/downsampleFreq:end);
        amplifierData = 0.195*(amplifierData'-32768);
    else
        amplifierData = [];
    end
end

function plotFrame(input)   
    global videoOut vidFig vidAx HvidAx %#ok<GVMIS>

    frame = input{1,1};
    numCam = size(frame,4);
    WL = input{1,2};
    iWL = input{1,3};
    
    if videoOut.QueueLength < 10
        if isempty(vidFig) == 1 || ishandle(vidFig) == 0
            vidFig = figure;
            vidFig.Name = 'Video';
            vidFig.Position = [0,805,960,190];
            vidAx = cell(1,numCam);
            HvidAx = cell(1,numCam);
    
            for i = 1:numCam
                vidAx{1,i} = subplot(1,numCam,i,'Parent',vidFig);
                HvidAx{1,i} = image(frame(:,:,:,iWL(1,i)),'Parent',vidAx{1,i});
                title(vidAx{1,i},WL{iWL(1,i),3})
            end
        else        
            for i = 1:numCam
                set(HvidAx{1,i},'CData',frame(:,:,:,iWL(1,i)))
            end  
        end
        drawnow
        pause(0.05)
    end
end

function plotData(input)
    global datFig datAx HdatAx %#ok<GVMIS>

    data = input{1,1};
    numChan = size(data,1);
    numTics = size(data,2);
    WL = input{1,2};
    iWL = input{1,3};
    
    if isempty(datFig) == 1 || ishandle(datFig) == 0

        datFig = figure;
        datFig.Name = 'Theta/Delta';
        datFig.Position = [0,530,960,190];
        datAx = cell(1,numChan);
        HdatAx = cell(numChan,1);
        for i = 1:numChan
             
            datAx{1,i} = subplot(1,numChan,i,'Parent',datFig);
            HdatAx{i,1} = plot(datAx{1,i},data(i,:));
            hold(datAx{1,i},'off')
            ylim(datAx{1,i},[0,4])
            ylabel(datAx{1,i},'T/D')
            xticks(0:round(numTics/6):numTics)
            xticklabels(numTics:-(round(numTics/6)):0)
            xlabel(datAx{1,i},'seconds')
            title(datAx{1,i},WL{iWL(1,i),3})
        end
    else
        for i = 1:numChan                   
            set(HdatAx{i,1},'YData',data(i,:))            
        end
    end
    drawnow
end

function plotLabel(input)
    global labelFig labelAx HlabelAx %#ok<GVMIS>

    data = input{1,1};
    numChan = size(data,1);
    numTics = size(data,2);
    WL = input{1,2};
    iWL = input{1,3};
    
    if isempty(labelFig) == 1 || ishandle(labelFig) == 0

        labelFig = figure;
        labelFig.Name = 'Classifier Output';
        labelFig.Position = [0,255,960,190];
        labelAx = cell(1,numChan);
        HlabelAx = cell(numChan,1);
        for i = 1:numChan
             
            labelAx{1,i} = subplot(1,numChan,i,'Parent',labelFig);         
            HlabelAx{i,1} = plot(labelAx{1,i},data(i,:));
            ylim(labelAx{1,i},[0,7])
            yticks(0:7)
            yticklabels(0:7)
            ylabel(labelAx{1,i},'Label')
            xticks(0:round(numTics/6):numTics)
            xticklabels(numTics:-(round(numTics/6)):0)
            xlabel(labelAx{1,i},'seconds')
            title(labelAx{1,i},WL{iWL(1,i),3})
        end
    else
        for i = 1:numChan                    
            set(HlabelAx{i,1},'YData',data(i,:))
        end
    end
    drawnow
end

function fun = defFun(input)
    
    xData = input{1,1};
    yData = input{1,2};
    fun = @(x) myPred1(xData,yData,table2array(x));
end