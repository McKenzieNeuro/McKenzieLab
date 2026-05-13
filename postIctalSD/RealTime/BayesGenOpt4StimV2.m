delete(gcp('nocreate'))
clear all %#ok<CLALL>

%%

% Testing schedule:

StartTime = [5,20]; % Time to begin recording in 24hr time, format = [hr,min]
StopTime = [17,5]; % Time to begin recording in 24hr time, format = [hr,min]
numDay = 5; % Number of days to run.

optSched(1,1) = 240; % Pre stim wash-out.
optSched(2,1) = 60; % Pre stim baseline.
optSched(3,1) = 60; % Post stim observation.
optSched(4,1) = 240; % Post stim wash-out.

% Box assignments:

% Subject in box, '' if empty.
boxSubject{1,1} = 'LPP 1.0'; % Box 1, port A
boxSubject{2,1} = 'LPP 2.0'; % Box 2, port B
boxSubject{3,1} = ''; % Box 3, port C
boxSubject{4,1} = ''; % Box 4, port D

% Optimization settings:

% Optimization parameters:
optParam(1).name = 'Amplitude';
optParam(1).range = [10,200];
optParam(1).step = 1;
optParam(1).exclude = [100,200];
optParam(2).name = 'Frequency';
optParam(2).range = [0.1,200];
optParam(2).step = 0.1;
optParam(2).exclude = [100,200];
optParam(3).name = 'Duration';
optParam(3).range = [1,300];
optParam(3).step = 1;
optParam(3).exclude = nan;

% Grid search settings

gridStep = 4; % Divisions of grid search
gridRep = 1; % Reps of initial grid search

% Optimizer order settings

whichOpt = [1,2]; % [1st,2nd], 1 = Bayesian, 2 = Genetic
switchOpt = 200; % Trials, not including grid, search to change optimizer

% Bayesian settings

boLP = 150; % When to start late phase mode
ER_early = 0.75; % early phase exploration ratio
ER_late = 0.5; % early phase exploration ratio

% Genetic settings

gaLP = 400; % When to start late phase mode
numOffspring_early = 10; % number of offspring from genOpt
numOffspring_late = 5; % number of offspring from genOpt

crossRate = 0.98; % Rate of trait crossover

mutationRate_early = 0.98; % Rate of mutation
mutationRate_late = 0.96; % Rate of mutation

dataStore = 5000; % Amount of data to hold in time window in ms

% Target for optimization, 1 = T/D, 2 = Model.
targetOpt(1,1) = 1; % Box 1
targetOpt(2,1) = 1; % Box 2
targetOpt(3,1) = 1; % Box 3
targetOpt(4,1) = 1; % Box 4

% Theta/Delta:
inputTD = 1000; % Amount of time to input into T/D filter in ms
thresholdTD = 1.4; % T/D threshold for stim, nan to disable.
thresholdDur = 60; % Duration that mean T/D must be under thresholdTD.
showData = 1; % Display T/D plots
tdRow = 2; % Which row to display plot, 1,2,3.
TD2Plot = 60; % Amount of time to plot T/D in sec

% Input channel for T/D, nan for unused
tdChan(1,1) = 4; % Port A
tdChan(2,1) = 4; % Port B
tdChan(3,1) = nan; % Port C
tdChan(4,1) = nan; % Port D

% Seizure dectection
useDetector = 1;
detectPath = 'D:\Data\ClassifierData\SeizureDetect\net1dCNN.mat';

% Prediction:
usePrediction = 1;
inputMP = 2000; % Amount of time to input into classifier in ms
showLabels = 0; % Display label plots
predRow = 3; % Which row to display plot, 1,2,3.
MP2Plot = 60; % Amount of time to plot predictions in s

% Model for each box
modelPath{1,1} = 'R:\McKenzieLab\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';
modelPath{2,1} = 'R:\McKenzieLab\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';
modelPath{3,1} = 'R:\McKenzieLab\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';
modelPath{4,1} = 'R:\McKenzieLab\Analysis\SeizureForecasting\IHKA_rat_RF\classification.mat';

% Input channel for Model, nan for unused
mpChan{1,:} = 0:7; % Port A
mpChan{2,:} = 0:7; % Port B
mpChan{3,:} = nan; % Port C
mpChan{4,:} = nan; % Port D

% Intan Settings:
intanIP_RHS = '127.0.0.1'; % IP for RHS.
intanPort1_RHS = 5000; % Port for RHS commands.
intanPort2_RHS = 5001; % Port for RHS data output.
sF = 20000; % Sampling frequency
dSF = 1000; % Downsample frequency
saveData = 1; % Save T/D and model data. 1 = yes, 0 = no.
saveLFP = 0; % Save downsampled LFP data through Matlab instead of intan.

% Active TCP channels to stream, [] = unused.
dataChan{1,1} = 0:7; % Port A
dataChan{2,1} = 0:7; % Port B
dataChan{3,1} = []; % Port C
dataChan{4,1} = []; % Port D

% Sim channel, [nan,nan] = unused, [c1,nan] = monopolar, [c1,c2] = bipolar.
stimChan(1,:) = [6,nan];
stimChan(2,:) = [6,nan];
stimChan(3,:) = [nan,nan];
stimChan(4,:) = [nan,nan];

% Stim polarity, -1 = Cathodic first, 1 = Anodic first. Ex: [-1,1].
stimPolarity(1,:) = [-1,1];
stimPolarity(2,:) = [-1,1];
stimPolarity(3,:) = [-1,1];
stimPolarity(4,:) = [-1,1];

% Phase duration in uS.
phaseDur(1,1) = 200;
phaseDur(2,1) = 200;
phaseDur(3,1) = 200;
phaseDur(4,1) = 200;

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

% Arduino Settings:
ArdPin{1,1} = 'D13'; % Arduino pin for LED.
ArdPin{2,1} = 'D8'; % Arduino pin for Timestamp.
LEDtime = 10; % Time between LED blinks is seconds

% Video Settings:
fps = 10;
useMotion = 1; % Use motion decector to control for movement
thresholdM = 20; % pixel change threshold to start stim.
thresholdDurM = 60; % Duration that max pixel change must be under thresholdM to start stim.
stimKillM = 20; % pixel change threshold to end stim if too strong.
RecVid = 1; % Record from webcams, 1 = yes, 0 = no
showVideo = 1; % Show video feed from cameras, 1 = yes, 0 = no
vidRow = 1; % Which row to display plot, 1,2,3.
showMot = 1; % Show motion plots from cameras, 1 = yes, 0 = no
MV2Plot = 60; % Amount of time to plot motion in s
motRow = 3; % Which row to display plot, 1,2,3.

%%

answer = questdlg('Would you like to start a new test or resume?', ...
	'Start Options','Start New','Resume','Start New');

switch answer
    case 'Start New'

        totBox = length(boxSubject);
        idxBox = 1:totBox;
        for i = 1:length(boxSubject)
            if isempty(boxSubject{i,1})
                idxBox(1,i) = 0;
            end
        end
        idxBox(idxBox == 0) = [];
        numBox = length(idxBox);
        
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
        activeBox = ~cellfun(@isempty,boxSubject);
        for i = 1:size(activeChan,3)

            activeChan(:,3,i) = {0};
            if activeBox(i,1) == 1
                [chanChoice,~] = listdlg('PromptString',{'Select active channels for';['port ',port{1,i},'.']},'ListString',activeChan(:,1,i),'SelectionMode','multiple');
                activeChan(chanChoice,3,i) = {1};
            end
        end
        
        numChan = 0;
        for i = 1:size(activeChan,3)
            for ii = 1:size(activeChan,1)
                if activeChan{ii,3,i} == 1
                    numChan = numChan+1;
                end
            end
        end
        
        PathName = uigetdir(cd,'Select folder for Datastore');
        cd(PathName);

        recPathName = [PathName,'\rec'];
        eval(['mkdir ',recPathName]);

        promptMessage = sprintf('Would you like to copy to another directory?');
        titleBarCaption = 'settings';
        copyAction = questdlg(promptMessage, titleBarCaption, 'Yes','No','Yes');
        
        if strcmp(copyAction,'Yes') == 1
            copyPathName = uigetdir(cd,'Select directory');
            copyRecPathName = [copyPathName,'\rec'];
            eval(['mkdir ',copyRecPathName]);
        end

        codeUsed = mfilename("fullpath");

        save([PathName,'\config.mat'])
        if strcmp(copyAction,'Yes') == 1
            save([copyPathName,'\config.mat'])
        end
    case 'Resume'

        PathName = uigetdir(cd,'Select folder for Datastore');
        cd(PathName);
        load([PathName,'\config.mat'])
end

%%

clc

debugMode = 0;

answer = questdlg('Would you like to start today or tomorrow?', ...
	'Start Options','Today','Tomorrow','Today');
switch answer
    case 'Today'
        nextDay = 0;
    case 'Tomorrow'
        nextDay = 1;
end

% CPU core roles
readTimer = 1;
intanWriter = 2;
thetaCalc = 3;
deltaCalc = 4;
seizureDetect = 5;
predictModel = 6:5+sum(activeBox);
bayesOptimizer = predictModel(1,end)+1:predictModel(1,end)+sum(activeBox);
dataSaver = bayesOptimizer(1,end)+1;
camReader = dataSaver+1;
camSaver = camReader+1;
moveCalc = camSaver+1:camSaver+sum(activeBox);

% Send tags
coreCom = 1;
testNoise = 2;
startSig = 3;
timerSig = 4;
recPhase = 5;
bayesSig = 6;
ampData = 7;
videoData = 8;
moveData = 9;

textOut = parallel.pool.DataQueue;
afterEach(textOut,@disp);

global boxPos %#ok<GVMIS>
boxPos = [];
stopBox = parallel.pool.DataQueue;
afterEach(stopBox,@makeBox);

if showVideo == 1
    global videoOut %#ok<GVMIS,TLEV>
    videoOut = parallel.pool.DataQueue;
    afterEach(videoOut,@plotFrame);
end

moveOut = parallel.pool.DataQueue;
afterEach(moveOut,@plotMotion);

waitShow = parallel.pool.DataQueue;
afterEach(waitShow,@textWaitBar);

dataOut = parallel.pool.DataQueue;
afterEach(dataOut,@plotData);

labelOut = parallel.pool.DataQueue;
afterEach(labelOut,@plotLabel);

global figPos %#ok<GVMIS>
figPos = [];
updateInfo = parallel.pool.DataQueue;
afterEach(updateInfo,@infoFig);

vStartTime = (StartTime(1,1)*3600)+(StartTime(1,2)*60);
vStopTime = (StopTime(1,1)*3600)+(StopTime(1,2)*60);
for Iter = 1:numDay

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

    fileCount = size(dir(recPathName),1)-2;
    
    I = 1;
    while getTime < vStopTime

        Now = clock; %#ok<CLOCK>
        Time = [num2str(Now(1,2)),'-',num2str(Now(1,3)),'-',num2str(Now(1,1)),'(',num2str(Now(1,4)),'.',num2str(Now(1,5)),')'];
        
        fileCount = fileCount+1;

        zeroPadStr = repmat('0',[1,4-length(num2str(fileCount))]);
        seshName = ['Stim',zeroPadStr,num2str(fileCount),'_',Time];
        seshPathName = [recPathName,'\',seshName];

        eval(['mkdir ',seshPathName]);
        send(textOut,['Starting trial ',num2str(fileCount)])
        send(textOut,['Directory set to: ',seshPathName])

        if isempty(gcp('nocreate')) == 1
            parpool('local',moveCalc(1,end));
            disp('Connecting to hardware...');
        end
        
        catchProb = 0;
        try
            spmd(moveCalc(1,end))
        
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
                    else
                        allClear = spmdReceive(intanWriter,coreCom);
                    end
        
                    if allClear == 0
                        
                        % Update channel values to relative index
                        tdChanRel = tdChan;
                        dataChanRel = dataChan;
                        for i = idxBox
    
                            if tdChan(i,1) < 10
                                chanStrTD = [port{1,i},'-00',num2str(tdChan(i,1))];
                            else 
                                chanStrTD = [port{1,i},'-0',num2str(tdChan(i,1))];
                            end
                            tdChanRel(i,1) = find(contains(activeChan(:,1,i),chanStrTD));
    
                            for ii = 1:size(activeChan(:,1,i),1)
    
                                chanStr = activeChan(ii,1,i);
                                dataChanRel{i,1}(1,ii) = find(contains(activeChan(:,1,i),chanStr));
                            end
                        end
    
                        % Make pollable queue for stopBox
                        qW = parallel.pool.PollableDataQueue;
    
                        % Create cell to hold all block titles
                        blockText = cell(5,1);
                        blockText{1,1} = 'Pre-Wash';
                        blockText{2,1} = 'Baseline';
                        blockText{3,1} = 'Stim';
                        blockText{4,1} = 'Observation';
                        blockText{5,1} = 'Post-Wash';
                        blockText{6,1} = 'Complete';
                
                        % Calculations for accurate parsing
                        waveformBytesPerFrame = 4+2*numChan;
                        waveformBytesPerBlock = 128*waveformBytesPerFrame+4;
    
                        % Pre-allocate counters and data stores
                        countLED = 1;
                        currentBlock = ones(1,totBox);
                        blockTime = zeros(1,totBox);
                        timeWindow = zeros(size(activeChan,1),dataStore,totBox);
                        lostSig = zeros(totBox,1)+10;
                        tdSend = zeros(totBox,dataStore);
                        labels = zeros(totBox,1);
                        plotTD = zeros(totBox,TD2Plot);
                        plotDect = zeros(totBox,TD2Plot);
                        plotLabels = zeros(totBox,MP2Plot);
                        obsTD_1 = nan(totBox,optSched(2,1));
                        mean_obsTD_1 = zeros(totBox,1);
                        obsTD_2 = nan(totBox,optSched(3,1));
                        obsLabels_1 = nan(totBox,optSched(2,1));
                        obsLabels_2 = nan(totBox,optSched(3,1));
                        
                        % Warm up calc cores and 
                        spmdSend(tdSend,[thetaCalc,deltaCalc],ampData);                    
                        
                        % Open plot figures
                        if showData == 1
                            if useDetector == 0
                                send(dataOut,{plotTD,tdRow,WLc,iWL})
                            elseif useDetector == 1
                                 send(dataOut,{cat(3,plotTD,plotDect),tdRow,WLc,iWL})
                            end
                        end
                        if showMot == 1
                            send(moveOut,{numBox,1:numBox,zeros(numBox,MV2Plot),3,WLc,iWL});
                        end
                        if showLabels == 1
                            send(labelOut,{plotLabels,predRow,WLc,iWL})
                        end
    
                        theta = spmdReceive(thetaCalc,ampData);
                        delta = spmdReceive(deltaCalc,ampData);
    
                        flush(intanRead_RHS);
                        
                        timerShed = optSched;
                        timerShed = [timerShed(1:2,1);0;timerShed(3:4,1)];
                        timerShed = repmat(timerShed,1,numBox);
    
                        for i = 1:4
                            if any(idxBox == i)
                                send(updateInfo,{i,0,'Pre-Wash'})
                            else
                                send(updateInfo,{i,0,'Empty'})
                            end
                        end
                        
                        % Wait for stim upload on intanWriter
                        allDur = spmdReceive(intanWriter,startSig);
                        timerShed(3,:) = ceil(allDur)';              
    
                        % sync with system clock
                        send(textOut,'----------')
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
    
                        % Open stop box
                        stopExp = 0;
                        send(stopBox,qW);
    
                        send(textOut,'----------')
                        send(textOut,'Running')              
    
                        % Main loop
                        while ~all(currentBlock(1,idxBox) == 6)
    
                            % Wait for data
                            while intanRead_RHS.NumBytesAvailable < waveformBytesPerBlock
                                pause(0.0001)
                            end
                            
                            % Sync signal
                            t = tic;
    
                            % Send sync signal to camReader
                            spmdSend(t,camReader,timerSig);
                            
                            % Trigger LED
                            if countLED == 1
    
                                t2 = tic;
    
                                % Turn on LED
                                writeDigitalPin(Ard,ArdPin{1,1},1);
    
                                % Pause for 0.2s
                                while toc(t2) < 0.2
                                    pause(0.0001)
                                end
                            
                                % Turn off LED
                                writeDigitalPin(Ard,ArdPin{1,1},0);
    
                                countLED = 2;
                            elseif countLED > 1 && countLED < LEDtime
                                countLED = countLED+1;
                            else
                                countLED = 1;
                            end
                            
                            % Send currentBlock to moveCalc
                            if useMotion == 1
                                for i = idxBox
                                    if blockTime(1,i) == 0

                                        spmdSend(currentBlock(1,i),moveCalc(1,i),recPhase)

                                        if currentBlock(1,i) == 6
                                            blockTime(1,i) = -1;
                                        end
                                    end
                                end
                            end
    
                            % Check if stop box has been closed
                            send(stopBox,1);
                            boxClose = poll(qW);
                            if ~isempty(boxClose)
                                stopExp = 1;
                                break
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
                            [amplifierData,offset] = byte2double2(waveformArray,numChan,blocksToRead,sF,dSF);
    
                            % If read cycle gets offset find the offset and correct
                            if offset > 0
                                send(textOut,['Read offset of ',num2str(offset),' detected.'])
                                waveformArray(1:offset) = [];
                                waveformArray = [waveformArray,read(intanRead_RHS,offset)];  %#ok<AGROW>
                                [amplifierData,offset] = byte2double2(waveformArray,numChan,blocksToRead,sF,dSF);
                                send(textOut,'Corrected')
                            end

                            % Parse channels into groups by box and check if signal is dropped
                            numSample = size(amplifierData,2);
                            readWindow = zeros(size(activeChan,1),numSample,totBox);
                            iChan = 0;
                            for i = idxBox 
    
                                for ii = 1:size(activeChan,1)
                                    if activeChan{ii,3,i} == 1 && any(dataChanRel{i,1} == ii)
                                        iChan = iChan+1;
                                        readWindow(ii,:,i) = amplifierData(iChan,:);
                                    end
                                end
                                
                                if lostSig(i,1) > 0 && currentBlock(1,i) < 6 && currentBlock(1,i) ~= 3
                                    
                                    meanChan = mean(readWindow(:,:,i),2);
                                    zcVolt = readWindow(:,:,i)-meanChan;
                                    stdVolt = std(zcVolt,0,2);
                                    maxVolt = max(abs(zcVolt),[],2);
                                    if max(abs(meanChan)) > 250 || sum(stdVolt < 20) >= 4 || sum(maxVolt > 1000) >= 4
                                        if lostSig(i,1) > 1
                                            lostSig(i,1) = lostSig(i,1)-1;
                                        elseif lostSig(i,1) == 1
                                            lostSig(i,1) = 0;
                                            blockTime(1,i) = -1;
                                            currentBlock(1,i) = 6;
                                            spmdSend({nan,[]},bayesOptimizer(1,i),bayesSig)
                                            send(updateInfo,{i,0,'Ended'})
                                        end
                                        send(updateInfo,{i,1,lostSig(i,1)})
                                    else
                                        if lostSig(i,1) > 0 && lostSig(i,1) < 10
                                            lostSig(i,1) = lostSig(i,1)+1;
                                            send(updateInfo,{i,1,lostSig(i,1)})
                                        end
                                    end
                                end
                            end
    
                            % Send to LFP                        
                            
                            % Reset timewindow after the stim
                            for i = idxBox
                                if currentBlock(1,i) == 4 && blockTime(1,i) == 0
                                    timeWindow(i,:) = 0;
                                end
                            end
                            
                            % Update time window
                            timeWindow = circshift(timeWindow,numSample,2);
                            timeWindow(:,1:numSample,:) = readWindow;
                            
                            % Pull channels for T/D
                            for i = idxBox
                                if currentBlock(1,i) ~= 3 && currentBlock(1,i) ~= 6
                                    tdSend(i,:) = timeWindow(tdChanRel(i,1),:,i);
                                else
                                    tdSend(i,:) = 0;
                                end
                            end
    
                            % Send timeWindow to calc cores for parallel computation
                            spmdSend(tdSend,[thetaCalc,deltaCalc],ampData);
    
                            % Send timeWindow to dectecor core
                            if useDetector == 1
                                spmdSend(timeWindow(:,1:1000,:),seizureDetect,ampData);
                            end
    
                            % Send timeWindow to prediction cores for parallel computation
                            coreCount = 1;
                            if usePrediction == 1
                                for i = idxBox
                                    spmdSend(timeWindow(:,1:inputMP,i),predictModel(1,coreCount),ampData);
                                    coreCount = coreCount+1;
                                end
                            end
                            
                            % Get T/D ratio
                            theta = spmdReceive(thetaCalc,ampData);
                            delta = spmdReceive(deltaCalc,ampData);
                            ToD = theta./delta;
                            
                            % Get detection
                            if useDetector == 1
                                detectOut = spmdReceive(seizureDetect,ampData);
                                plotDect(:,end+1) = detectOut*4; %#ok<SAGROW>
                                plotDect(:,1) = [];
                            end
                            
                            % Get labels
                            if usePrediction == 1
                                for i3 = 1:numBox
                                    modelOut = spmdReceive('any',ampData);
                                    labels(modelOut(1,1),1) = modelOut(1,2);
                                end
                            end
    
                            % Update T/D plots
                            if showData == 1
                                plotTD(:,end+1) = ToD; %#ok<SAGROW>
                                plotTD(:,1) = [];
    
                                if useDetector == 0
                                    data = plotTD;
                                elseif useDetector == 1
                                    data = cat(3,plotTD,plotDect);
                                end
                                send(dataOut,{data,tdRow,WLc,iWL})
                            end
                            
                            % Update label plots
                            if usePrediction == 1 && showLabels == 1
                                plotLabels = [plotLabels,labels]; %#ok<AGROW>
                                plotLabels(:,1) = [];
                                send(labelOut,{plotLabels,predRow,WLc,iWL})
                            end
    
                            % Add data to observations
                            for i = 1:numBox
                                if currentBlock(1,i) == 1 || currentBlock(1,i) == 5 % Wash
                                    blockTime(1,i) = blockTime(1,i)+1;
                                elseif currentBlock(1,i) == 2 % Baseline
    
                                    % Get motion data from moveCalc
                                    if useMotion == 1                                   
                                        if spmdProbe(moveCalc(1,i),moveData) == 1
                                            maxMove = spmdReceive(moveCalc(1,i),moveData);
                                        else
                                            maxMove = nan;
                                        end
                                    else
                                        maxMove = 0;
                                    end
    
                                    if targetOpt(i,1) == 1 % T/D
                                        
                                        obsTD_1(i,:) = circshift(obsTD_1(i,:),1);
                                        obsTD_1(i,1) = ToD(i,1);
                                        
                                        if ~isnan(thresholdTD)
                                            if all(~isnan(obsTD_1(i,1:thresholdDur)))
                                                mean_obsTD_1(i,1) = mean(obsTD_1(i,1:thresholdDur));
                                                if mean_obsTD_1(i,1) < thresholdTD && maxMove < thresholdM && lostSig(i,1) > 6
                                                    blockTime(1,i) = timerShed(2,i);
                                                end
                                            end
                                        else
                                            blockTime(1,i) = blockTime(1,i)+1;
                                        end
                                    elseif targetOpt(i,1) == 2 % Model
                                        
                                        obsLabels_1(i,:) = circshift(obsTD_1(i,:),1);
                                        obsLabels_1(i,1) = labels(i,1);
                                        
                                        if ~isnan(thresholdTD)
                                            if all(~isnan(obsLabels_1(i,:)))
                                                mean_obsTD_1(i,1) = mean(obsLabels_1(i,1:thresholdDur));
                                                if mean_obsTD_1(i,1) < thresholdTD && maxMove < thresholdM && lostSig(i,1) > 6
                                                    blockTime(1,i) = timerShed(2,i);
                                                end
                                            end
                                        else
                                            blockTime(1,i) = blockTime(1,i)+1;
                                        end
                                    end
    
                                    if blockTime(1,i) >= timerShed(2,i)
                                        if targetOpt(i,1) == 1
                                            spmdSend({2,obsTD_1(i,:)},bayesOptimizer(1,i),bayesSig)
                                        elseif targetOpt(i,1) == 2
                                            spmdSend({2,obsLabels_1(i,:)},bayesOptimizer(1,i),bayesSig)
                                        end
                                    elseif blockTime(1,i) == 0
                                        blockTime(1,i) = 1;
                                    end
                                elseif currentBlock(1,i) == 3 % Stim
                                    
                                    if blockTime(1,i) == 0                              
                                        spmdSend([i,1],intanWriter,recPhase)
                                    elseif blockTime(1,i) >= timerShed(3,i)-1
                                        spmdSend([i,0],intanWriter,recPhase)
                                    end
                                    blockTime(1,i) = blockTime(1,i)+1;
    
                                    if spmdProbe(moveCalc(1,i),moveData) == 1
                                        if spmdReceive(moveCalc(1,i),moveData) == -1
                                            blockTime(1,i) = 0;
                                            currentBlock(1,i) = 5;
                                            send(updateInfo,{i,0,'Post-Wash'})
                                            spmdSend({2,0},bayesOptimizer(1,i),bayesSig)
                                            spmdSend({4,0},bayesOptimizer(1,i),bayesSig)
                                        end
                                    end
                                elseif currentBlock(1,i) == 4 % Observation
    
                                    if spmdProbe(moveCalc(1,i),moveData) == 1
                                        if spmdReceive(moveCalc(1,i),moveData) == -1
                                            blockTime(1,i) = 0;
                                            currentBlock(1,i) = 5;
                                            send(updateInfo,{i,0,'Post-Wash'})
                                            spmdSend({2,0},bayesOptimizer(1,i),bayesSig)
                                            spmdSend({4,0},bayesOptimizer(1,i),bayesSig)
                                        end
                                    end
    
                                    if targetOpt(i,1) == 1
                                        obsTD_2(i,:) = circshift(obsTD_2(i,:),1);
                                        obsTD_2(i,1) = ToD(i,1);
                                    elseif targetOpt(i,1) == 2
                                         obsLabels_2(i,:) = circshift(obsLabels_2(i,:),1);
                                        obsLabels_2(i,1) = labels(i,1);
                                       
                                    end
                                    blockTime(1,i) = blockTime(1,i)+1;
    
                                    if blockTime(1,i) >= timerShed(4,i)
                                        if targetOpt(i,1) == 1
                                            spmdSend({4,obsTD_2(i,:)},bayesOptimizer(1,i),bayesSig)
                                        elseif targetOpt(i,1) == 2
                                            spmdSend({4,obsLabels_2(i,:)},bayesOptimizer(1,i),bayesSig)
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
                            
                            % Update current block
                            for i = idxBox
                                if currentBlock(1,i) < 6
                                    if blockTime(1,i) >= timerShed(currentBlock(1,i),i)
                                        blockTime(1,i) = 0;
                                        currentBlock(1,i) = currentBlock(1,i)+1;
                                        send(updateInfo,{i,0,blockText{currentBlock(1,i),1}})
                                    end
                                end
                            end
    
                            % Wait till 1s has passed
                            while toc(t) < 0.9999
                                pause(0.0001)
                            end

                            if getTime >= vStopTime
                                break
                            end
                        end                    
                        
                        spmdSend([],camReader,timerSig)
                        spmdSend([],intanWriter,recPhase)
                        % Send stop to datasaver core if enabled
                        if saveData == 1
                            spmdSend([],dataSaver,ampData);
                        end
    
                        flush(intanRead_RHS)
                        send(updateInfo,[])
    
                        if stopExp == 0
                            send(stopBox,[]);
                        end
    
                        send(textOut,'----------')
                        send(textOut,'Recording finished.')                 
                    end
                    spmdSend([],[thetaCalc,deltaCalc,seizureDetect,predictModel],ampData);
                    spmdSend([],bayesOptimizer,bayesSig);

                    if debugMode == 1
                        send(textOut,spmdIndex)
                    end
                elseif spmdIndex == intanWriter
                    
                    % Connect and test connection to RHS commands
                    if I == 1
    
                        t = 0;
                        comIssue = 0;
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
                    else
                        try
                            write(intanSend_RHS,uint8('execute RescanPorts;'));
                            spmdSend(0,readTimer,coreCom)
                        catch
                            spmdSend(1,readTimer,coreCom)
                            send(textOut,'Connection to Intan lost')
                        end
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
    
                        % Set up config structure
                        stimConfig(1).port = 'a';
                        stimConfig(2).port = 'b';
                        stimConfig(3).port = 'c';
                        stimConfig(4).port = 'd';
                        stimConfig(1).fKey = 1;
                        stimConfig(2).fKey = 2;
                        stimConfig(3).fKey = 3;
                        stimConfig(4).fKey = 4;
                        stimConfig(1).chan = stimChan(1,:);
                        stimConfig(2).chan = stimChan(2,:);
                        stimConfig(3).chan = stimChan(3,:);
                        stimConfig(4).chan = stimChan(4,:);
                        stimConfig(1).amp = [];
                        stimConfig(2).amp = [];
                        stimConfig(3).amp = [];
                        stimConfig(4).amp = [];
                        stimConfig(1).polarity = stimPolarity(1,:);
                        stimConfig(2).polarity = stimPolarity(2,:);
                        stimConfig(3).polarity = stimPolarity(3,:);
                        stimConfig(4).polarity = stimPolarity(4,:);
                        stimConfig(1).duration = phaseDur(1,1);
                        stimConfig(2).duration = phaseDur(2,1);
                        stimConfig(3).duration = phaseDur(3,1);
                        stimConfig(4).duration = phaseDur(4,1);
                        stimConfig(1).pulseNum = [];
                        stimConfig(2).pulseNum = []; 
                        stimConfig(3).pulseNum = []; 
                        stimConfig(4).pulseNum = [];
                        stimConfig(1).period = [];
                        stimConfig(2).period = [];
                        stimConfig(3).period = [];
                        stimConfig(4).period = [];
    
                        % Ask bayesOptimizer cores to pick new stim settings
                        spmdSend(1,bayesOptimizer(1,idxBox),bayesSig);
                        
                        phase = zeros(1,numBox);
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
    
                            stimConfig(idxBox(1,i)).amp = nextParam{i,1}(1,1);
                                 
                            if nextParam{i,1}(1,2) <= 5
    
                                stimConfig(idxBox(1,i)).pulseNum = 1;
                                stimConfig(idxBox(1,i)).period = 975;
                                clockStop(i,1) = 1/nextParam{i,1}(1,2);
                                numStim(i,1) = floor(nextParam{i,1}(1,3)/(clockStop(i,1)));
                                if numStim(i,1) == 0
                                    clockStop(i,1) = 1/nextParam{i,1}(1,2);
                                    numStim(i,1) = 1;
                                end
                            else
                                pulseNum = round(nextParam{i,1}(1,2));
                                period = (1/nextParam{i,1}(1,2))*1000;
    
                                while pulseNum*period > 975
                                    pulseNum = pulseNum-1;
                                end
    
                                stimConfig(idxBox(1,i)).pulseNum = pulseNum;
                                stimConfig(idxBox(1,i)).period = period;
                                clockStop(i,1) = 1;
                                numStim(i,1) = nextParam{i,1}(1,3);
                            end
                            allDur(i,1) = nextParam{i,1}(1,3);
                        end
                        clockStop = uint64(clockStop*10000000);
                        
                        % Set up stim parameters
                        write(intanSend_RHS,uint8('execute RescanPorts;'));
                        pause(0.1)
                        setStim(intanSend_RHS,1,stimConfig)
                        write(intanSend_RHS,uint8('execute uploadstimparameters;'));
                        
                        % Send stim duration to moveCalc cores 
                        if useMotion == 1
                            for i = 1:numBox
                                spmdSend(allDur(i,1),moveCalc(1,i),startSig)
                            end
                        end
    
                        % Wait for stim upload to complete
                        send(textOut,'Uploading stim parameters')
                        send(waitShow,13)
                        pause(16)
                        send(textOut,'----------')
                        
                        % Print stim params in command window
                        for i = idxBox
                            send(textOut,['Box ',num2str(i),' -> Amplitude: ',num2str(nextParam{i,1}(1,1)),'uA, Frequency: ',num2str(nextParam{i,1}(1,2)),'Hz, Duration: ',num2str(nextParam{i,1}(1,3)),'s.'])
                        end
                        
                        % Send durations to readTimer and let it know the upload is complete
                        spmdSend(allDur,readTimer,startSig)
    
                        % Wait for sync signal from readTimer
                        spmdReceive(readTimer,startSig);
        
                        % Start recording
                        if saveLFP == 0
                            write(intanSend_RHS,uint8('set runMode record;'));
                        elseif saveLFP == 1
                            write(intanSend_RHS,uint8('set runMode run;'));
                        end
                        
                        while 1
    
                            if spmdProbe('any',recPhase) == 1 
                                
                                in = spmdReceive('any',recPhase);
    
                                if isempty(in)
                                    break
                                end
    
                                phase(1,in(1,1)) = in(1,2);
                            end
    
                            if any(phase == 1)
                                
                                % Create string array of commands to send to intan
                                commandArray = '';
                                for i = idxBox
                                    if phase(1,i) == 1
                                        if ticStart(i,1) == 0
                                            ticStart(i,1) = tic;
                                            ticStop(i,1) = ticStart(i,1)+clockStop(i,1);
                                            commandArray = [commandArray,['execute ManualStimTriggerPulse F',num2str(stimConfig(i).fKey),';']]; %#ok<AGROW>
                                            numStim(i,1) = numStim(i,1)-1;
                                        else
                                            if tic >= ticStop(i,1) && numStim(i,1) > 0
                                                ticStart(i,1) = tic;
                                                ticStop(i,1) = ticStart(i,1)+clockStop(i,1);
                                                commandArray = [commandArray,['execute ManualStimTriggerPulse F',num2str(stimConfig(i).fKey),';']]; %#ok<AGROW>
                                                numStim(i,1) = numStim(i,1)-1;
                                            end
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
                        if getTime >= vStopTime
                            
                            setStim(intanSend_RHS,0,stimConfig)
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

                    if debugMode == 1
                        send(textOut,spmdIndex)
                    end
                elseif spmdIndex == thetaCalc                
                    while 1
                        try
                            if spmdProbe(readTimer,ampData) == 1
        
                                tdSend = spmdReceive(readTimer,ampData);
        
                                if isempty(tdSend)
                                    break
                                end
        
                                % Calculate mean theta amp
                                theta = abs(hilbert(bandpass(tdSend',[5,12],dSF)));
                                thetaMean = mean(theta(1:inputTD,:))';
            
                                spmdSend(thetaMean,readTimer,ampData)
                            else
                                pause(0.001)
                            end
                        catch
                            spmdSend(0,readTimer,ampData)
                            send(textOut,'Unknown packet')
                        end
                    end

                    if debugMode == 1
                        send(textOut,spmdIndex)
                    end
                elseif spmdIndex == deltaCalc
                    while 1
                        try
                            if spmdProbe(readTimer,ampData) == 1
        
                                tdSend = spmdReceive(readTimer,ampData);
        
                                if isempty(tdSend)
                                    break
                                end
        
                                % Calculate mean delta amp
                                delta = abs(hilbert(bandpass(tdSend',[1,4],dSF)));
                                deltaMean = mean(delta(1:inputTD,:))';
            
                                spmdSend(deltaMean,readTimer,ampData)    
                            else
                                pause(0.001)
                            end
                        catch
                            spmdSend(0,readTimer,ampData)
                            send(textOut,'Unknown packet')
                        end
                    end

                    if debugMode == 1
                        send(textOut,spmdIndex)
                    end
                elseif any(spmdIndex == seizureDetect)
    
                    if useDetector == 1
                        detectOut = zeros(totBox,1);
                        net1dCNN = load(detectPath);
                        net1dCNN = net1dCNN.net1dCNN;
                        predict(net1dCNN,zeros(size(activeChan,1),dSF));
                    end
    
                    while 1
                        if spmdProbe(readTimer,ampData) == 1
                            
                            try
                                data = spmdReceive(readTimer,ampData);
                            catch
                                continue
                            end
    
                            if isempty(data)
                                break
                            end
    
                            for i = idxBox
                                detectOut(i,1) = mean(predict(net1dCNN,data(:,:,i)/10000));
                            end
    
                            spmdSend(detectOut,readTimer,ampData)
                        end
                    end

                    if debugMode == 1
                        send(textOut,spmdIndex)
                    end
                elseif any(spmdIndex == predictModel)
    
                    boxID = find(predictModel == spmdIndex);
                    
                    if ~isempty(boxSubject{boxID,1})
                        model = load(modelPath{boxID,1});
                        nChannel = model.ops.nCh_raw;
                        freq = model.ops.freqs;
                        nFreq = length(freq);
                        id = spmdIndex-predictModel(1,1)+1;
                        feat = nan(1,nFreq*nChannel);
                    end
    
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

                    if debugMode == 1
                        send(textOut,spmdIndex)
                    end
                elseif any(spmdIndex == bayesOptimizer)
    
                    boxID = find(bayesOptimizer == spmdIndex);
                    
                    if any(idxBox == boxID)
                        if I == 1
                            
                            if targetOpt(boxID,1) == 1
    
                                fileName = [PathName,'\Box',num2str(boxID),'_gridData_TD.mat'];
                                if strcmp(copyAction,'Yes') == 1
                                    fileNameC = [copyPathName,'\Box',num2str(boxID),'_gridData_TD.mat'];
                                end
                            elseif targetOpt(boxID,1) == 2
        
                                fileName = [PathName,'\Box',num2str(boxID),'_gridData_MP.mat'];
                                if strcmp(copyAction,'Yes') == 1
                                    fileNameC = [copyPathName,'\Box',num2str(boxID),'_gridData_MP.mat'];
                                end
    
                                scoring = zeros(6,1);
                                scoring(1,1) = 1;
                                scoring(2,1) = 10;
                                scoring(3,1) = 50;
                                scoring(4,1) = 100;
                                scoring(5,1) = 1000;
                                scoring(6,1) = 1000;
                            end
        
                            gridData = matfile(fileName,'Writable',true);
                            if strcmp(copyAction,'Yes') == 1
                                gridDataC = matfile(fileNameC,'Writable',true);
                            end
    
                            if isprop(gridData,'parameters') == 0
                                
                                % Load initial parameters
                                param = optParam;
                                numParam = size(param,2);
                                gridData.parameters = param;
                                
                                gridData.folders = cell(0,1);
                                gridData.obs = nan(0,7);
                                gridData.obsHeader = {'Pre-Stim BL mean','Post-Stim BL mean','Post-Stim max','Post-Stim max location','Post-Stim num over threshold','Movement','Score'};
        
                                % Get grid search coordinates
                                idx = cell(1,numParam);
                                for ii = 1:numParam
            
                                    m = optParam(ii).range(1);
                                    s = optParam(ii).step;
                                    M = optParam(ii).range(2);
                                    
                                    range = m:s:M;
        
                                    l = length(range);
                                    idx{1,ii} = zeros(1,gridStep);
                                    idx{1,ii}(1,1:gridStep-1) = round(1:l/(gridStep-1):l);
                                    idx{1,ii}(1,end) = l;
                                    idx{1,ii} = unique(range(1,idx{1,ii}));
                                end
                                
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
                                
                                % Find points that fall into exclusion zone
                                exclude = zeros(numCells,2);
                                for i = 1:size(optParam,2)                           
                                    exclude(i,:) = optParam(i).exclude;
                                end
                                useCon = find(~isnan(exclude(:,1)))';
                                exPoint = xConstraint(idc,exclude);
                                
                                % Move excluded points to boundary of exclusion zone
                                for i = 1:length(exPoint)
                                    if exPoint(i,1) == 0
                                        rmvPoint = 0;
                                        for ii = 1:numCells
                                            if idc(i,ii) >= exclude(ii,1) && idc(i,ii) <= exclude(ii,2)*0.6
                                                idc(i,ii) = exclude(ii,1)-1;
                                            elseif idc(i,ii) >= exclude(ii,2)*0.6 && idc(i,ii) <= exclude(ii,2)
                                                rmvPoint = rmvPoint+1;
                                            end
                                        end
                                        if rmvPoint == length(useCon)
                                            idc(i,:) = nan;
                                        end
                                    end
                                end
                                idc(isnan(idc(:,1)),:) = [];
        
                                idc = repmat(idc,[gridRep,1]);
                                idc = datasample(idc,size(idc,1),'Replace',false);
                                gridData.search = idc;
                                gridData.intGrid = idc;
                                gridData.tracker = [mutationRate_early,0];
    
                                if strcmp(copyAction,'Yes') == 1
                                    gridDataC.folders = gridData.folders;
                                    gridDataC.intGrid = gridData.intGrid; 
                                    gridDataC.obs = gridData.obs; 
                                    gridDataC.obsHeader = gridData.obsHeader; 
                                    gridDataC.parameters = gridData.parameters; 
                                    gridDataC.search = gridData.search; 
                                    gridDataC.tracker = gridData.tracker; 
                                end
                            else
                                param = gridData.parameters;
                            end
                        else
                            gridData = matfile(fileName,'Writable',true);
                            if strcmp(copyAction,'Yes') == 1
                                gridDataC = matfile(fileNameC,'Writable',true);
                            end
                        end
    
                        % See if any entry is incomplete
                        numObs = size(gridData.obs,1);
                        numObs_search = numObs-size(gridData.intGrid,1);
                        
                        if numObs == 0
                            numObs = 1;
                            gridData.obs = nan(1,7);
                        else
                            numObs = numObs+1;
                            missData = find(isnan(gridData.obs(:,2)),1,'first');
                            if isempty(missData)
                                gridData.obs(numObs,:) = nan;
                            else
                                numObs = missData;
                            end
                        end
                        gridData.folders(numObs,1) = {seshName};
                        
                        % Decide whitch opt to use
                        if numObs_search < switchOpt
                            useOpt = whichOpt(1,1);
                        else
                            useOpt = whichOpt(1,2);
                        end

                        % Set up optimizable variables
                        if useOpt == 1
                        
                            count = 1;
                            exclude = zeros(size(optParam,2),2);
                            useParm = 1:size(optParam,2);
                            for i = 1:size(optParam,2)
                                if length(optParam(i).range) == 2
                                    if rem(optParam(i).step,1) == 0
                                        optVar(count).Var = optimizableVariable(optParam(i).name,[optParam(i).range(1),optParam(i).range(2)],'Type','integer'); %#ok<SAGROW>
                                    else
                                        optVar(count).Var = optimizableVariable(optParam(i).name,[optParam(i).range(1),optParam(i).range(2)],'Type','real'); %#ok<SAGROW>
                                    end
                                    exclude(i,:) = optParam(i).exclude;
                                    count = count+1;
                                else
                                    useParm(1,i) = 0;
                                end
                            end
                            useParm(useParm == 0) = [];
                        else
    
                            useParm = 1:size(optParam,2);
                            for i = 1:size(optParam,2)
                                if length(optParam(i).range) < 2
                                    useParm(1,i) = 0;
                                end
                            end
    
                            range = {optParam.range};
                            steps = {optParam.step};
                            exclude = {optParam.exclude};
                            
                            newMin = 0;
                            countStale = 0;
                            minFitBest = min(gridData.obs(:,7));
                        end
                    end
        
                    pM = 0;
                    gotObs = 0;
                    while 1
                        if spmdProbe('any',bayesSig) == 1
    
                            stimOpt = spmdReceive('any',bayesSig);
    
                            if isempty(stimOpt)
                                break
                            end                      
                            
                            if ~iscell(stimOpt)                            
                                if numObs > size(gridData.search,1) % use bayesOpt
                                    if useOpt == 1
    
                                        if numObs_search < boLP
                                            explorRatio = ER_early;
                                        else
                                            explorRatio = ER_late;
                                        end
    
                                        nextParam = bayesOpt(gridData.search(1:numObs-1,useParm),gridData.obs(1:numObs-1,7),optVar,explorRatio,exclude);
        
                                        count = 1;
                                        for i = 1:size(gridData.search,2)
        
                                            x = param(i).step;
                                            n = 0;
                                            while floor(x*10^n) ~= x*10^n
                                                n = n+1;
                                            end
                                            nextParam(1,i) = round(nextParam(1,i),n);
                                            
                                            if any(i == useParm)
                                                gridData.search(numObs,i) = nextParam(1,count);
                                                count = count+1;
                                            else
                                                gridData.search(numObs,i) = gridData.search(numObs,i);
                                            end
                                        end
                                    elseif useOpt == 2 % use genOpt
    
                                        if numObs < gaLP
                                            numOffspring = numOffspring_early;
                                            mutationRateS = mutationRate_early;
                                        else
                                            numOffspring = numOffspring_late;
                                            mutationRateS = mutationRate_late;
                                        end
    
                                        population = gridData.search(1:numObs-1,:);
                                        fitness = gridData.obs(1:numObs,7);
                                        mutationRate = gridData.tracker(1,1);
    
                                        if min(fitness) < minFitBest
    
                                            minFitBest = min(fitness);
                                            gridData.tracker(1,2) = 0;
                                            if newMin == 0
                                                newMin = 1;
                                                mutationRate = 0.9;
                                            else
                                                mutationRate = mutationRate-0.02;
                                            end
                                        else
                                            if gridData.tracker(1,2) < numOffspring
                                                gridData.tracker(1,2) = gridData.tracker(1,2)+1;
                                            elseif gridData.tracker(1,2) >= numOffspring
                                                newMin = 0;
                                                gridData.tracker(1,2) = 0;
                                                if round(mutationRate,2) <= mutationRateS-0.02
                                                    mutationRate = mutationRate+0.02;
                                                end
                                            end
                                        end
    
                                        nextParam = genOpt(population,fitness(1:size(population,1),1),numOffspring,range,steps,exclude,0,crossRate,mutationRate,[]);
                                        nextParam = datasample(nextParam,size(nextParam,1),'Replace',false);
                                        gridData.search(numObs:numObs+numOffspring-1,:) = nextParam;
                                        gridData.tracker(1,1) = mutationRate;
                                    end
                                end
                                nextParam = gridData.search(numObs,:);
                                spmdSend([boxID,nextParam],intanWriter,bayesSig)
                            else
                                
                                data = stimOpt{1,2};
                                if ~isempty(data)
                                    for i = 1:length(data)
                                        if targetOpt(boxID,1) == 2
                                            data(1,i) = scoring(data(1,i),1);
                                        end
                                    end
                                end
                                
                                if ~isnan(stimOpt{1,1})
                                    if stimOpt{1,1} == 2
                                        gotObs = 1;
                                        gridData.obs(numObs,1) = round(mean(data),3);
                                    elseif stimOpt{1,1} == 4
    
                                        gotObs = 2;
                                        gridData.obs(numObs,2) = round(mean(data),3);
                                        [gridData.obs(numObs,3),gridData.obs(numObs,4)] = max(data);
                                        gridData.obs(numObs,5) = sum(data > thresholdTD);
                                        
                                        if useMotion == 1
                                            pM = spmdReceive(moveCalc(1,boxID),moveData);
                                            gridData.obs(numObs,6) = pM;
                                        end
                                        
                                        % Scoring
                                        preTD = gridData.obs(numObs,1);
                                        postTD = gridData.obs(numObs,2);
                                        peakTD = gridData.obs(numObs,3);
                                        peakLoc = gridData.obs(numObs,4);
                                        numOT = gridData.obs(numObs,5);
                                        gridData.obs(numObs,7) = -((postTD-preTD)+((peakTD/peakLoc)/2))*numOT;
                                        send(textOut,['Box ',num2str(boxID),' Score = ',num2str(gridData.obs(numObs,7))])
                                    end
                                else
                                    if gotObs == 1
                                        gridData.obs(numObs,:) = NaN;
                                    end
                                end
                            end
                        else
                            pause(0.001)
                        end
                    end
                    if strcmp(copyAction,'Yes') == 1
                        gridDataC.folders = gridData.folders;
                        gridDataC.intGrid = gridData.intGrid; 
                        gridDataC.obs = gridData.obs; 
                        gridDataC.obsHeader = gridData.obsHeader; 
                        gridDataC.parameters = gridData.parameters; 
                        gridDataC.search = gridData.search; 
                        gridDataC.tracker = gridData.tracker; 
                        gridDataC = [];
                    end
                    gridData = [];

                    if debugMode == 1
                        send(textOut,spmdIndex)
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
                        if usePrediction == 1
                            fclose(fID2);
                        end
                    end

                    if debugMode == 1
                        send(textOut,spmdIndex)
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
                        
                        if useMotion == 1

                            frameTest = uint8(zeros(Dim(1,1),Dim(1,2),3,numCam,fps*5));
                            for i = 1:2
                                numFrames = 0;
                                t = tic;
                                while toc(t) < 0.99
                                    numFrames = numFrames+1;
                                    for ii = 1:numCam
                                        frameTest(:,:,:,ii,numFrames) = snapshot(cam{ii,1});
                                    end
                                end  
                                t = 0;
                            end
                            frameTest = frameTest(:,:,:,:,round(linspace(1,numFrames,fps)));

                            for i = 1:numCam
                                spmdSend(frameTest(:,:,:,iWL(1,i),:),moveCalc(1,i),testNoise);
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
                            send(videoOut,{frame+256,vidRow,WLc,iWL}); %#ok<SPGV>                        
                        end
                        
                        try
                            spmdReceive(readTimer,startSig);
                        catch
                        end
    
                        while 1
                            if spmdProbe(readTimer,timerSig) == 1
    
                                % wait for timer signal
                                try
                                    t = spmdReceive(readTimer,timerSig);
                                catch
                                    t = tic;
                                end
    
                                if ~isempty(t)
                                    while toc(t) < 0.99
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

                    if debugMode == 1
                        send(textOut,spmdIndex)
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
    
                        frames = uint8(zeros(Dim(1,1),Dim(1,2),3,numCam,fps*5));
    
                        while 1
    
                            numFrames = 0;                            
                            while 1
                                if spmdProbe(camReader,videoData) == 1
                                    
                                    try
                                        frame = spmdReceive(camReader,videoData);
                                    catch
                                        if numFrames < 20
                                            continue
                                        else
                                            frame = 0;
                                        end
                                    end

                                    if size(frame,1) > 1
        
                                        if showVideo == 1
                                            if videoOut.QueueLength < 10 %#ok<SPGV>
                                                send(videoOut,{frame,vidRow,WLc,iWL}); %#ok<SPGV>
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
                            
                            if numFrames > 0
                                syncFrames = frames(:,:,:,:,round(linspace(1,numFrames,fps)));
                                
                                if useMotion == 1
                                    for i = 1:numCam
                                        spmdSend(syncFrames(:,:,:,iWL(1,i),:),moveCalc(1,i),videoData);
                                    end
                                end
        
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
                        end                    
                        spmdSend([],moveCalc,videoData);                    
                        
                        if RecVid == 1
                            send(textOut,'Finalizing video files.')
                            for ii = 1:numCam
                                close(writerObj{ii,1});
                            end
                            writerObj = [];
                            send(textOut,'----------')
                        end
                    end

                    if debugMode == 1
                        send(textOut,spmdIndex)
                    end
                elseif any(spmdIndex == moveCalc)
    
                    boxID = find(moveCalc == spmdIndex);
                    
                    pM = 0;
                    phase = 1;
                    obsComplete = 0;
                    
                    if useMotion == 1

                        baseMotion = nan(1,thresholdDurM);
                        strike = 0;
                        strikeLim = [0,0,3,5];
                        motionPlot = zeros(1,MV2Plot);

                        if I == 1

                            frameTest = spmdReceive(camReader,testNoise);
                            frameTest = squeeze(frameTest);
                            frameTest = frameTest(1:100,1:100,:,:);

                            meanDiff = zeros(1,fps-1);
                            firstFrame = rgb2gray(frameTest(:,:,:,1));
                            for i = 2:fps
                                frameCurrent = rgb2gray(frameTest(:,:,:,i));
                                d = double(abs(frameCurrent-firstFrame))/255;
                                d = d(:);
                                d = d(d > 0);
                                meanDiff(1,i-1) = max(d);
                                firstFrame = frameCurrent;
                            end
                            camNoise = mean(meanDiff);                            
                        end
                        stimDur = spmdReceive(intanWriter,startSig);
                        stimMotion = nan(1,stimDur);
                    end                
    
                    while 1
    
                        if spmdProbe(readTimer,recPhase) == 1
                            phase = spmdReceive(readTimer,recPhase);
                        end
                        
                        if spmdProbe(camSaver,videoData) == 1
    
                            frames = spmdReceive(camSaver,videoData);
        
                            if ~isempty(frames)
        
                                numFrames = size(frames,5);
                                sumDiff = zeros(1,numFrames-1);
                                camFrames = squeeze(frames);
                                firstFrame = rgb2gray(camFrames(:,:,:,1));
                                for i = 2:numFrames
                                    frameCurrent = rgb2gray(camFrames(:,:,:,i));
                                    d = double(abs(frameCurrent-firstFrame))/255;
                                    d(d > 0) = d(d > 0)-camNoise;
                                    d(d < 0) = 0; 
                                    sumDiff(1,i-1) = sum(d,'all');
                                    firstFrame = frameCurrent;
                                end                           
                                meanD = mean(sumDiff/prod(Dim))*20000;

                                if showMot == 1
                                    motionPlot = circshift(motionPlot,-1);
                                    motionPlot(1,end) = meanD;
                                    send(moveOut,{1,boxID,motionPlot,3,WLc,iWL});
                                end
                                
                                if obsComplete == 0
                                    if phase == 2
                                        baseMotion(1,end) = meanD;
                                        baseMotion = circshift(baseMotion,1);
                                        spmdSend(max(baseMotion),readTimer,moveData)
                                    elseif phase == 3 || phase == 4
    
                                        if any(isnan(stimMotion))
                                            stimMotion(1,end) = meanD;
                                            stimMotion = circshift(stimMotion,1);
                                        end
    
                                        % if motion exceeds stimKillM end session
                                        if meanD > stimKillM
    
                                            strike = strike+1;
                                            if strike == strikeLim(1,phase)
                                                spmdSend([boxID,0],intanWriter,recPhase)
                                                spmdSend(-1,readTimer,moveData)
                                                spmdSend(1000,bayesOptimizer(1,boxID),moveData)
                                                phase = 5;
                                                obsComplete = 1;
                                            end
                                        else
                                            if strike > 0
                                                strike = strike-1;
                                            end
                                        end
                                    elseif phase == 5
                                        pM = (length(find(stimMotion > thresholdM))/length(stimMotion));
                                        spmdSend(pM,bayesOptimizer(1,boxID),moveData)
                                        obsComplete = 1;
                                    end
                                end
                            else
                                break
                            end
                        else
                            pause(0.001)
                        end
                    end

                    if debugMode == 1
                        send(textOut,spmdIndex)
                    end
                end
            end
        catch ME
            catchProb = 1;
            spmd(moveCalc(1,end))
                if spmdIndex == intanWriter
                    write(intanSend_RHS,uint8('set runMode stop;'));
                elseif spmdIndex == dataSaver
                    if saveData == 1
                        fclose(fID);
                    end
                    if usePrediction == 1
                        fclose(fID2);
                    end
                elseif spmdIdex == camSaver
                    if RecVid == 1
                        for ii = 1:numCam
                            close(writerObj{ii,1});
                        end
                        writerObj = [];
                    end
                end
            end
        end

        if catchProb == 0 && any(currentBlock{1}(idxBox) > 4)
            if strcmp(copyAction,'Yes') == 1
                disp(['Copying data to ',copyRecPathName]);
                copyfile(seshPathName,[copyRecPathName,'\',seshName]);
                send(textOut,'----------')
            end
        elseif all(currentBlock{1}(idxBox) < 5)
            rmdir(seshPathName,'s')
        end
        
        if stopExp{1} == 1
            answer = questdlg('Would you like to stop the experiment?', ...
	        'Start Options','No','Yes','No');
            if strcmp(answer,'Yes')
                break
            end
        end
        I = I+1;
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

function makeBox(input)
    global closeBox boxPos qSend %#ok<GVMIS>

    if ~isempty(input)
        if input ~= 1
            
            qSend = input;

            if isempty(boxPos)
                boxPos = [1278,728,161,60];
            end

            closeBox = errordlg('Press ok to stop.','Stop','non-modal');
            closeBox.Position = boxPos;
            drawnow
        else
            if ~ishandle(closeBox)
                send(qSend,1)
            else
                boxPos = closeBox.Position;
            end
        end
    else
        close(closeBox)
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
    rowNum = input{1,2};
    numCam = size(frame,4);
    WL = input{1,3};
    iWL = input{1,4};
    
    if videoOut.QueueLength < 10
        if isempty(vidFig) == 1 || ishandle(vidFig) == 0

            vidFig = figure;
            vidFig.Name = 'Video';
            if rowNum == 1
                vidFig.Position = [0,805,960,190];
            elseif rowNum == 2
                vidFig.Position = [0,530,960,190];
            elseif rowNum == 3
                vidFig.Position = [0,255,960,190];
            end

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

function plotMotion(input)
    global movFig movAx HmovAx %#ok<GVMIS>
    
    totChan = input{1,1};
    curChan = input{1,2};
    data = input{1,3};
    rowNum = input{1,4};
    WL = input{1,5};
    iWL = input{1,6};
    
    if isempty(movFig) == 1 || ishandle(movFig) == 0

        movFig = figure;
        movFig.Name = 'Motion';
        if rowNum == 1
            movFig.Position = [0,805,960,190];
        elseif rowNum == 2
            movFig.Position = [0,530,960,190];
        elseif rowNum == 3
            movFig.Position = [0,255,960,190];
        end 

        movAx = cell(1,totChan);
        HmovAx = cell(totChan,1);
        for i = 1:length(curChan)
             
            movAx{1,curChan(1,i)} = subplot(1,totChan,curChan(1,i),'Parent',movFig);
            HmovAx{curChan(1,i),1} = plot(movAx{1,curChan(1,i)},data(i,:));
            hold(movAx{1,curChan(1,i)},'off')
            ylim(movAx{1,curChan(1,i)},[0,100])
            ylabel(movAx{1,curChan(1,i)},'Movement')
            xlabel(movAx{1,curChan(1,i)},'seconds')
            title(movAx{1,curChan(1,i)},WL{iWL(1,curChan(1,i)),3})
        end
    else
        for i = 1:length(curChan)
            set(HmovAx{curChan(1,i),1},'YData',data(i,:))
        end
    end
    drawnow
end

function plotData(input)
    global datFig datAx HdatAx %#ok<GVMIS>

    data = input{1,1};
    rowNum = input{1,2};
    numTics = size(data,2);
    WL = input{1,3};
    iWL = input{1,4};
    numChan = length(iWL);
    
    if isempty(datFig) == 1 || ishandle(datFig) == 0

        datFig = figure;
        datFig.Name = 'Theta/Delta';
        if rowNum == 1
            datFig.Position = [0,805,960,190];
        elseif rowNum == 2
            datFig.Position = [0,530,960,190];
        elseif rowNum == 3
            datFig.Position = [0,255,960,190];
        end 

        datAx = cell(1,numChan);
        HdatAx = cell(numChan,2);
        for i = 1:numChan
             
            datAx{1,i} = subplot(1,numChan,i,'Parent',datFig);
            HdatAx{i,1} = plot(datAx{1,i},data(i,:,1));
            HdatAx{i,2} = yline(mean(data(i,:,1)),'LineWidth',2);
            if size(data,3) == 2
                hold(datAx{1,i},'on')
                HdatAx{i,3} = plot(datAx{1,i},data(i,:,2));
            end           
            hold(datAx{1,i},'off')
            ylim(datAx{1,i},[0,4])
            ylabel(datAx{1,i},'T/D')
            xlim(datAx{1,i},[0,60])
            xticks(0:round(numTics/6):numTics)
            xticklabels(numTics:-(round(numTics/6)):0)
            xlabel(datAx{1,i},'seconds')
            title(datAx{1,i},WL{iWL(1,i),3})
            if size(data,3) == 1
                legend('T/D','Mean T/D')
            else
                legend('T/D','Mean T/D','Seizure (4x)')
            end
        end
    else
        for i = 1:numChan                   
            set(HdatAx{i,1},'YData',data(i,:,1))
            set(HdatAx{i,2},'Value',mean(data(i,:,1)))
            if size(data,3) == 2
                set(HdatAx{i,3},'YData',data(i,:,2))
            end
        end
    end
    drawnow
end

function plotLabel(input)
    global labelFig labelAx HlabelAx %#ok<GVMIS>

    data = input{1,1};
    rowNum = input{1,2};
    numTics = size(data,2);
    WL = input{1,3};
    iWL = input{1,4};
    numChan = length(iWL);
    
    if isempty(labelFig) == 1 || ishandle(labelFig) == 0

        labelFig = figure;
        labelFig.Name = 'Classifier Output';
        if rowNum == 1
            labelFig.Position = [0,805,960,190];
        elseif rowNum == 2
            labelFig.Position = [0,530,960,190];
        elseif rowNum == 3
            labelFig.Position = [0,255,960,190];
        end 

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

function infoFig(input)
    global fig figPos gl pnl lbl_phase lbl_sig %#ok<GVMIS>

    if isempty(fig) == 1 || ishandle(fig) == 0 

        fig = uifigure;

        if isempty(figPos)
            figPos = [962,42,959,90];
        end

        fig.Position = figPos;
        gl = uigridlayout(fig,[1,4]);
        gl.RowHeight = {'1x'};

        pnl = cell(1,4);
        lbl_phase = cell(1,4);
        lbl_sig = cell(1,4);
        for i = 1:4

            pnl{1,i} = uipanel(gl);
            pnl{1,i}.Layout.Column = i;
            pnl{1,i}.Title = ['Box ',num2str(i)];
            pnl{1,i}.TitlePosition = 'centertop';
            pnl{1,i}.FontSize = 18;
            pnl{1,i}.FontWeight = 'bold';

            drawnow

            lbl_sig{1,i} = uilabel(pnl{1,i});
            lbl_sig{1,i}.InnerPosition = [1,20,200,30];
            lbl_sig{1,i}.Text = 'Signal strength: ||||||||||';
            lbl_sig{1,i}.FontSize = 16;

            lbl_phase{1,i} = uilabel(pnl{1,i});
            lbl_phase{1,i}.Text = '';
            lbl_phase{1,i}.FontSize = 20;
            lbl_phase{1,i}.InnerPosition = [25,0,200,30];
            lbl_phase{1,i}.HorizontalAlignment = 'center';
        end
    end

    if ~isempty(input)

        box = input{1,1};
        sig = input{1,2};
        data = input{1,3};

        if sig == 0
            lbl_phase{1,box}.Text = data;
        elseif sig == 1
            if data > 0

                sigText = 'Signal strength: ';
                for i = 1:data
                    sigText = [sigText,'|']; %#ok<AGROW>
                end
                lbl_sig{1,box}.Text = sigText;
            else
                lbl_sig{1,box}.Text = 'Signal strength: lost';
            end
        end
        drawnow

        figPos = fig.Position;
    else
        close(fig)
    end
end

function nextParam = bayesOpt(population,score,optVar,explorRatio,exclude)

    funObj = defObjFun({population,score});

    if isempty(exclude)

        results = bayesopt(funObj,[optVar.Var],...
            'Verbose',0,...
            'PlotFcn',[],...
            'AcquisitionFunctionName','expected-improvement-plus',...
            'NumSeedPoints',10,...
            'ExplorationRatio',explorRatio,...
            'GPActiveSetSize',200,...
            'MaxObjectiveEvaluations',50,...
            'IsObjectiveDeterministic',false,...
            'UseParallel',false,...
            'MaxObjectiveEvaluations',50);
            
    else

        funCon = defConFun(exclude);
        results = bayesopt(funObj,[optVar.Var],...
            'Verbose',0,...
            'PlotFcn',[],...
            'AcquisitionFunctionName','expected-improvement-plus',...
            'NumSeedPoints',10,...
            'ExplorationRatio',explorRatio,...
            'GPActiveSetSize',200,...
            'MaxObjectiveEvaluations',50,...
            'IsObjectiveDeterministic',false,...
            'UseParallel',false,...
            'MaxObjectiveEvaluations',50,...
            'XConstraintFcn',funCon);
    end
    nextParam = table2array(results.NextPoint);

    function fun = defObjFun(input)
    
        xData = input{1,1};
        yData = input{1,2};
        fun = @(x) predictGauss(xData,yData,table2array(x));
    
        function pred = predictGauss(x0,y0,x)
            pred = predict(fitrgp(x0,y0),x);
        end
    end

    function fun = defConFun(exclude)
        fun = @(x) xConstraint(table2array(x),exclude);

        function tf = xConstraint(x,exclude)

            useCon = find(~isnan(exclude(:,1)))';
            tf = true(size(x,1),length(useCon));
            for i = useCon
                tf(:,i) = x(:,i) >= exclude(i,1) & x(:,i) <= exclude(i,2);
            end
            tf = ~all(tf == 1,2);
        end
    end
end

function offSpring = genOpt(population,fitness,numOffspring,range,steps,exclude,dieOff,crossRate,mutationRate,immigrant)
    
    crossRateN = 10*(1-crossRate);
    mutationRateN = 10*(mutationRate-1);

    % Rank the fitness of the population
    [rankFit,idx] = sort(fitness);
    rankPop = population(idx,:);
    
    % Remove individuals above the dieOff threshold
    rmv = find(rankFit >= dieOff);
    rankFit(rmv,:) = [];
    rankPop(rmv,:) = [];
    numLeft = length(rankFit);

    % Get probability of reproduction using half normal pdf
    prob = normpdf(-numLeft:numLeft,0,5);
    prob = prob(1,1:numLeft+1);
    prob = prob/max(prob);
    prob = prob(1,1:numLeft);
    prob = flip(prob);
    
    % Pick successful breeders, individuals can breed more than once
    breaders = zeros(numOffspring,1);
    while all(breaders == mean(breaders))
        for i = 1:numOffspring    
            breaders(i,1) = randsample(1:numLeft,1,true,prob);
        end
    end
    
    % Find all possible pairings, breeders cannot pair with themselves
    allPairs = cell(numOffspring,numOffspring);
    for i = 1:numOffspring
        for ii = 1:numOffspring
            if breaders(i,1) ~= breaders(ii,1)
                allPairs{i,ii} = sort([breaders(i,1),breaders(ii,1)]);
            end
        end
    end
    allPairs = allPairs(:);
    allPairs = allPairs(~cellfun('isempty',allPairs),1);
    
    % Pick the final pairings
    pairsFin = cell(1,max(breaders));
    for i = unique(sort(breaders))'
    
        count = 1;
        pairPool = zeros(1,1);
        for ii = 1:length(allPairs)
            if allPairs{ii,1}(1,1) == i || allPairs{ii,1}(1,2) == i
                pairPool(count,1) = ii;
                count = count+1;
            end
        end
    
        numOff = sum(breaders == i);
        idx = randsample(pairPool,numOff);
        pairsFin(1:numOff,i) = allPairs(idx,1);
    end
    pairsFin = pairsFin(:);
    pairsFin = pairsFin(~cellfun('isempty',pairsFin),1);
    
    % Mate immigrants to pool of selected breeders
    if ~isempty(immigrant)
        numImmigrant = size(immigrant,1);
        for i = 1:numImmigrant
            pairsFin{end+1,1} = [randsample(breaders,1),-i]; %#ok<AGROW>
        end
    else
        numImmigrant = 0;
    end
    
    % Recombination of traits
    numTraits = size(population,2);
    offSpring = zeros(numOffspring+numImmigrant,numTraits);
    for i = 1:numOffspring+numImmigrant
        
        % Genes from each parent
        G1 = rankPop(pairsFin{i,1}(1,1),:);
        if pairsFin{i,1}(1,2) > 0
            G2 = rankPop(pairsFin{i,1}(1,2),:);
        else
            G2 = immigrant(-pairsFin{i,1}(1,2),:);
        end

        while 1

            isExclude = zeros(1,numTraits);
            for ii = 1:numTraits
                while 1
                    
                    % If the traits are different, set up crossover range
                    if G1(1,ii) < G2(1,ii)
                        crossRange = G1(1,ii):steps{1,ii}:G2(1,ii);
                    elseif G1(1,ii) > G2(1,ii)
                        crossRange = G2(1,ii):steps{1,ii}:G1(1,ii);
                    else
                        crossRange = [];
                    end
                    
                    % Crosserover of traits
                    if ~isempty(crossRange)

                        crLh = floor(length(crossRange)/2);
                        isOdd = (length(crossRange)/2)-crLh > 0;
    
                        % Get probability for crossover using exponential decay
                        % Each gene has a 50% chance to select unchanged
                        % Blending chance drops off exponentialy
                        if isOdd == 0
                            crossProb = exp(-crossRateN*(1:crLh));
                            crossProb = [crossProb,flip(crossProb)]; %#ok<AGROW>
                        else
                            crossProb = [exp(-crossRateN*(1:crLh)),flip(exp(-crossRateN*(1:crLh+1)))];
                        end
                        crossProb = crossProb/max(crossProb)/2;
                        newTrait = randsample(crossRange,1,true,crossProb);
                    else
                        newTrait = G1(1,ii);
                    end

                    % Get probability for mutation using exponential decay
                    % Each gene has a 50% chance to remain unchanged
                    % mutation amount drops off exponentialy
                    mutationRange{1,1} = 1:newTrait;
                    mutationRange{1,2} = 1:range{1,ii}(1,2)-newTrait-1;
                    muteRangeFull = [flip(-mutationRange{1,1}),0,mutationRange{1,2}];
                   
                    mutationProb{1,1} = exp(mutationRateN*(1:length(mutationRange{1,1})+1));
                    mutationProb{1,1} = (mutationProb{1,1}/max(mutationProb{1,1}))/2;
                    mutationProb{1,2} = (exp(mutationRateN*(1:length(mutationRange{1,2})+1)));
                    mutationProb{1,2} = (mutationProb{1,2}/max(mutationProb{1,2}))/2;
                    muteProbFull = [flip(mutationProb{1,1}(1,2:end)),mutationProb{1,2}];

                    mutation = randsample(muteRangeFull,1,true,muteProbFull);
                    offSpring(i,ii) = newTrait+mutation;
    
                    if offSpring(i,ii) >= range{1,ii}(1,1) && offSpring(i,ii) <= range{1,ii}(1,2)
                        break
                    end
                end
                
                % Check if new gene falls within exclusion zone
                if ~isempty(exclude)
                    if offSpring(i,ii) >= exclude{1,ii}(1,1) && offSpring(i,ii) <= exclude{1,ii}(1,2)
                        isExclude(1,ii) = 1;
                    end
                end
            end

            % If the new offspring does not fall into exclusion zone exit loop
            if sum(isExclude) < numTraits-1
                break
            end
        end
    end
end