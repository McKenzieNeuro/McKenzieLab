delete(gcp('nocreate'))
clear all %#ok<CLALL>

StartTime = [9,0]; % Time to begin recording in 24hr time, format = [hr,min]
nextDay = 0; % Start recording tomorrow?

startIter = 1; % Start on which day of the schedule
startBlock = 1; % Start on which block of the schedule

stimParam(1).name = 'Info';

% only non-empty strings will be active
stimParam(1).config(1).subject = 'Test1'; % Subject in box 1
stimParam(1).config(2).subject = 'Test2'; % Subject in box 2
stimParam(1).config(3).subject = ''; % Subject in box 3
stimParam(1).config(4).subject = ''; % Subject in box 4

% Program for experiment
% 1 = baseline
% 2 = open loop
% 3 = kindling
% 4 = closed loop
% 5 = induction peak
% 6 = induction decay
program = 4;

% Block types: 

% BL = Baseline
% stimNum_stimTrigger = Stim experiment.
%   stimNum = S1 or S2.
%   S1X or S2X = sham of respective stim
%   stimTrigger = O, tdC, ptC, dtC, bOp, ipC
%       O = Open loop
%       tdC = theta-delta trigger closed loop 
%       ptC = peak theta-delta trigger closed loop
%       dtC = decay theta-delta trigger closed loop
%       ipC = intelligent trigger

% Block schedule:
% [{Day1 block1},{Day2 block1},...
%  {Day1 block2},{Day2 block2},...];

for i = 1:4
    if program == 1 % Baseline

        schedule = [{'BL'},{'BL'},{'BL'},{'BL'},{'BL'},{'BL'},{'BL'}];
        duration = 3600*8;

    elseif program == 2 % Open loop

        schedule = [{'BL'  },{'BL'   },{'BL'  },{'BL'   },{'BL'  },{'BL'   };...
                    {'S1_O'},{'S1X_O'},{'S1_O'},{'S1X_O'},{'S1_O'},{'S1X_O'};...
                    {'BL'  },{'BL'   },{'BL'  },{'BL'   },{'BL'  },{'BL'   }];

        duration = [120;120;120];

    elseif program == 3 % Kindling

        schedule = [{'BL'  },{'BL'   },{'BL'  },{'BL'   },{'BL'  },{'BL'   };...
                    {'S2_O'},{'S2X_O'},{'S2_O'},{'S2X_O'},{'S2_O'},{'S2X_O'};...
                    {'BL'  },{'BL'   },{'BL'  },{'BL'   },{'BL'  },{'BL'   }];

        duration = [3600;3600*6;3600];

    elseif program == 4 % Closed loop

        schedule = [{'S1X_ipC'},{'S1X_ipC'},{'S1X_ipC'},{'S1X_ipC'},{'S1X_ipC'},{'S1X_ipC'};...
                    {'S1_ipC' },{'S1X_ipC'},{'S1_ipC' },{'S1X_ipC'},{'S1_ipC' },{'S1X_ipC'};...
                    {'S1X_ipC'},{'S1X_ipC'},{'S1X_ipC'},{'S1X_ipC'},{'S1X_ipC'},{'S1X_ipC'}];

        duration = [3600*1;3600*2;3600*1];

    elseif program == 5 % Induction peak

        schedule = [{'BL'    },{'BL'},{'BL'    },{'BL'},{'BL'    },{'BL'};...
                    {'S1_O'  },{'SX'},{'S1_O'  },{'SX'},{'S1_O'  },{'SX'};...
                    {'S2_ptC'},{'SX'},{'S2_ptC'},{'SX'},{'S2_ptC'},{'SX'};...
                    {'BL'    },{'BL'},{'BL'    },{'BL'},{'BL'    },{'BL'}];

        duration = [3600;3600;10;3600];

    elseif program == 6 % Induction decay

        schedule = [{'BL'    },{'BL'},{'BL'    },{'BL'},{'BL'    },{'BL'};...
                    {'S1_O'  },{'SX'},{'S1_O'  },{'SX'},{'S1_O'  },{'SX'};...
                    {'S2_dtC'},{'SX'},{'S2_dtC'},{'SX'},{'S2_dtC'},{'SX'};...
                    {'BL'    },{'BL'},{'BL'    },{'BL'},{'BL'    },{'BL'}];

        duration = [3600;3600;10;3600];
    end

    stimParam(1).config(i).schedule = schedule;
    stimParam(1).config(i).duration = duration;
end

% Post Recording baseline
usePBL = 0; % 1 = yes, 0 = no
durPBL = 3600*22; % Duration in seconds

% Video Settings:
fps = 10; % Frame rate of video recording
RecVid = 1; % Record from webcams, 1 = yes, 0 = no
showVideo = 1; % Show video feed from cameras, 1 = yes, 0 = no
labelFrames = 0; % Use CNN to label frames, 1 = yes, 0 = no

% Arduino Settings:
ArdPin{1,1} = 'D13'; % Arduino pin for LED.
ArdPin{2,1} = 'D8'; % Arduino pin for Timestamp.
LEDtime = 8:12; % Time between LED blinks is seconds

% Intan TCP:
useRHS = 1;
intanIP_RHS = '127.0.0.1'; % IP for RHS.
intanPort1_RHS = 5000; % Port for RHS commands.
intanPort2_RHS = 5001; % Port for RHS data output.
recIntan = 1; % Record data from Intan, 1 = yes, 0 = no
saveDS = 1; % Save downsampled data, 1 = yes, 0 = no
saveTS = 1; % Save downsampled timestamps, 1 = yes, 0 = no
saveStim = 1; % Save downsampled stim timeSpamps, 1 = yes, 0 = no

useRHD = 0;
intanIP_RHD = '127.0.0.2'; % IP for RHD.
intanPort1_RHD = 5000; % Port for RHD commands.
intanPort2_RHD = 5001; % Port for RHD data output.

boxRHS(1,1) = 1; % Box for RHS port A
boxRHS(2,1) = 2; % Box for RHS port B
boxRHS(3,1) = 3; % Box for RHS port C
boxRHS(4,1) = 4; % Box for RHS port D

boxRHD(1,1) = 5; % Box for RHD port A
boxRHD(2,1) = 6; % Box for RHD port B
boxRHD(3,1) = 7; % Box for RHD port C
boxRHD(4,1) = 8; % Box for RHD port D

% Stim parameters:

% Stim 1 (prophylactic):
stimParam(2).name = 'Stim 1';

% Ports for S1
stimParam(2).config(1).port = 'a';
stimParam(2).config(2).port = 'b';
stimParam(2).config(3).port = 'c';
stimParam(2).config(4).port = 'd';

% F Key for trigger
stimParam(2).config(1).fKey = 1;
stimParam(2).config(2).fKey = 2;
stimParam(2).config(3).fKey = 3;
stimParam(2).config(4).fKey = 4;

% Channel, [nan,nan] = unused, [c1,nan] = monopolar, [c1,c2] = bipolar.
% -1 = Digital-Out 1, -2 = Digital-Out 2
stimParam(2).config(1).chan = [7,nan];
stimParam(2).config(2).chan = [7,nan];
stimParam(2).config(3).chan = [7,nan];
stimParam(2).config(4).chan = [-1,nan];

% Amplitudes in uA.
stimParam(2).config(1).amp = 100;
stimParam(2).config(2).amp = 100;
stimParam(2).config(3).amp = 100;
stimParam(2).config(4).amp = 100;

% Stim polarity, -1 = Cathodic first, 1 = Anodic first. Ex: [-1,1].
stimParam(2).config(1).polarity = [-1,-1];
stimParam(2).config(2).polarity = [-1,-1];
stimParam(2).config(3).polarity = [-1,-1];
stimParam(2).config(4).polarity = [-1,-1];

% Phase duration in uS.
stimParam(2).config(1).duration = 200;
stimParam(2).config(2).duration = 200;
stimParam(2).config(3).duration = 200;
stimParam(2).config(4).duration = 50000;

% use freq and duration
stimFreq = zeros(4,2);
stimFreq(1,1) = 0;
stimFreq(2,1) = 0;
stimFreq(3,1) = 0;
stimFreq(4,1) = 0;

stimDur = zeros(4,2);
stimDur(1,1) = 0;
stimDur(2,1) = 0;
stimDur(3,1) = 0;
stimDur(4,1) = 0;

pulsNum = zeros(4,1,2);
pulsNum(:,2,:) = 1;
period = zeros(4,1);
for i = 1:4

    pulsNum(i,1,1) = round(stimFreq(i,1)*stimDur(i,1));
    period(i,1) = 1000/stimFreq(i,1);
    if pulsNum(i,1,1) > 256

        pulsNum(i,1,1) = round(stimFreq(i,1));
        pulsNum(i,2,1) = stimDur(i,1);
        period(i,1) = 1000/stimFreq(i,1);
        while pulsNum(i,1,1)*period(i,1) > 950
            pulsNum(i,1,1) = pulsNum(i,1,1)-1;
        end
    end
end

% Number stim train of pulses
stimParam(2).config(1).pulseNum = 4;
stimParam(2).config(2).pulseNum = 4;
stimParam(2).config(3).pulseNum = 4; 
stimParam(2).config(4).pulseNum = 10; 

% Train period in mS
stimParam(2).config(1).period = 10;
stimParam(2).config(2).period = 10;
stimParam(2).config(3).period = 10;
stimParam(2).config(4).period = 100;

% Inter stim intervals in seconds, randomly selected.
ISI{1,1} = 10:15; % Open loop
ISI{2,1} = 10:15; % Theta/delt trig closed loop 
ISI{3,1} = 10:15; % Peak theta trig closed loop 
ISI{4,1} = 10:15; % Decay theta trig closed loop 
ISI{5,1} = 300; % Intelligent predictor closed loop

% Total number of stims applied in block. Inf = limited only by block time.
stimNum(1,1) = Inf; 

% Stim 2 (induction):
stimParam(3).name = 'Stim 2';

% Ports for S1
stimParam(3).config(1).port = 'a';
stimParam(3).config(2).port = 'b';
stimParam(3).config(3).port = 'c';
stimParam(3).config(4).port = 'd';

% F Key for trigger
stimParam(3).config(1).fKey = 5;
stimParam(3).config(2).fKey = 6;
stimParam(3).config(3).fKey = 7;
stimParam(3).config(4).fKey = 8;

% Channel, [nan,nan] = unused, [c1,nan] = monopolar, [c1,c2] = bipolar.
% -1 = Digital-Out 1, -2 = Digital-Out 2
stimParam(3).config(1).chan = [4,nan];
stimParam(3).config(2).chan = [4,nan];
stimParam(3).config(3).chan = [4,nan];
stimParam(3).config(4).chan = [4,nan];

% Amplitudes in uA.
stimParam(3).config(1).amp = 1000;
stimParam(3).config(2).amp = 0;
stimParam(3).config(3).amp = 0;
stimParam(3).config(4).amp = 0;

% Stim polarity, -1 = Cathodic first, 1 = Anodic first. Ex: [-1,1].
stimParam(3).config(1).polarity = [-1,-1];
stimParam(3).config(2).polarity = [-1,-1];
stimParam(3).config(3).polarity = [-1,-1];
stimParam(3).config(4).polarity = [-1,-1];

% Phase duration in uS.
stimParam(3).config(1).duration = 200;
stimParam(3).config(2).duration = 200;
stimParam(3).config(3).duration = 200;
stimParam(3).config(4).duration = 200;

% use freq and duration
stimFreq(1,2) = 10;
stimFreq(2,2) = 10;
stimFreq(3,2) = 0;
stimFreq(4,2) = 0;

stimDur(1,2) = 10;
stimDur(2,2) = 10;
stimDur(3,2) = 0;
stimDur(4,2) = 0;

period = zeros(4,1);
for i = 1:4

    pulsNum(i,1,2) = round(stimFreq(i,2)*stimDur(i,2));
    period(i,1) = 1000/stimFreq(i,2);
    if pulsNum(i,1,2) > 256

        pulsNum(i,1,2) = round(stimFreq(i,2));
        pulsNum(i,2,2) = stimDur(i,2);
        period(i,1) = 1000/stimFreq(i,21);
        while pulsNum(i,1,2)*period(i,1) > 950
            pulsNum(i,1,2) = pulsNum(i,1,2)-1;
        end
    end
end

% Number stim train of pulses
stimParam(3).config(1).pulseNum = 4;
stimParam(3).config(2).pulseNum = 4; 
stimParam(3).config(3).pulseNum = 4;
stimParam(3).config(4).pulseNum = 4; 

% Train period in mS
stimParam(3).config(1).period = 100;
stimParam(3).config(2).period = 100;
stimParam(3).config(3).period = 100;
stimParam(3).config(4).period = 100;

% Inter stim intervals in seconds, randomly selected.
ISI{1,2} = 1790; % Open loop
ISI{2,2} = 10; % Theta/delt trig closed loop 
ISI{3,2} = 10; % Peak theta trig closed loop 
ISI{4,2} = 10; % Decay theta trig closed loop 
ISI{5,2} = 10; % Intelligent predictor closed loop

% Total number of stims applied in block. Inf = limited only by block time.
stimNum(1,2) = 12; 

% Input settings:

% Input channel for model, [] = unused.
dataChan{1,1} = 0:7; % Port A
dataChan{2,1} = 0:7; % Port B
dataChan{3,1} = 0:7; % Port C
dataChan{4,1} = 0:7; % Port D

% Input channel for T/D closed-loop, nan = unused.
tdChan(1,1) = 5; % Port A
tdChan(2,1) = 5; % Port B
tdChan(3,1) = 5; % Port C
tdChan(4,1) = 5; % Port D

lengthTW_Filt = 5000; % Amount of data held for filter calculations in ms.
lengthTW_Mean = 1000; % Amount of data held for averaging in ms. Must be <= lengthTW_Filt
showData = 1; % Plot theta/dela data.
lengthData = 60; % seconds of data to store for plot.
saveTD = 1; % Save theda/delta to .dat.

threshold_tdC = 1; % Theda/delta threshold for closed loop.
threshold_ptC = 2; % Peak Theda/delta threshold for closed loop.
sustain_ptC = 4; % Peak Theda/delta sustain for closed loop in seconds.
threshold_dtC = 1; % Decay theda/delta threshold for closed loop.
sustain_dtC = 4; % Decay Theda/delta sustain for closed loop in seconds.

% Prediction:
usePrediction = 1;
modelPath{1,1} = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\model_pilo_bilateral.mat';
modelPath{2,1} = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\model_pilo_bilateral.mat';
modelPath{3,1} = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\model_pilo_bilateral.mat';
modelPath{4,1} = 'R:\Analysis\SeizureForecasting\IHKA_rat_RF\model_pilo_bilateral.mat';

 



% info for full feature space
ops.Fs  = 1250; % sampling rate of lfp files

% frequency bands for coherence. must match sm_getPowerPerChannel

% frequency selection for phase/amplitude (must match index in getPowerPerChannel)dSF




inputDur = 4000; % Amount of time to input into classifier in ms
stimClass = 2; % Label to trigger stim
numConsec = 5; % Threshold of stimClass needed to trigger stim

% Analysis Settings
runAnalysis = 1;
sF = 20000; % Recording sample frequency
dSF = 1000; % Downsample frequency
use60HzFilt = 0;
ACfiltPath = 'G:\RecordingData\EDS\CNN60Hz.mat';
replaceFile = 1;

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

answer = questdlg('Would you like to start a new test or resume?', ...
	'Start Options','Start New','Resume','Start New');

switch answer
    case 'Start New'

        WL = webcamlist;

        useCam = true(1,length(WL));
        for i = 1:length(WL)
             if strcmp(WL{i,1},'HD Pro Webcam C920') == 0
                 useCam(1,i) = false;
             end
        end
        
        for i = 1:length(useCam)
            if useCam(1,i) == 1
        
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
            else
                WL{i,2} = 'N/A';
            end
        end
        
        if size(WL,1) > 0
        
            for i = 1:length(useCam)
                if strcmp(WL{i,2},'N/A') == 0
                    WL{i,3} = char(stimParam(1).config(WL{i,2}).subject+" ("+WL{i,2}+")");
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
            if ~isempty(stimParam(1).config(i).subject)
                [chanChoice,~] = listdlg('PromptString',{'Select active channels for';['port ',port{1,i},'.']},'ListString',activeChan(:,1,i),'SelectionMode','multiple');
                activeChan(chanChoice,3,i) = {1};
            end
        end
        numChan = sum(cell2mat(activeChan(:,3,:)),'all');
        numChanStream = 2; %hardcoded for now
        if labelFrames == 1
            CNNPathName = uigetdir(cd,'Select folder for CNN');
        end
        
        PathName = uigetdir(cd,'Select folder for Datastore');
        cd(PathName);
        
        promptMessage = sprintf('Would you like to copy to another directory?');
        titleBarCaption = 'settings';
        copyAction = questdlg(promptMessage, titleBarCaption, 'Yes','No','Yes');
        
        codeUsed = mfilename("fullpath");
        
        if strcmp(copyAction,'Yes') == 1
            copyPathName = uigetdir(cd,'Select directory');
            save([copyPathName,'\config.mat'])
        end
        
        save([PathName,'\config.mat'])
        save([PathName,'\iter.mat'],'startIter')
            
    case 'Resume'

        delete(gcp('nocreate'))
        clear all %#ok<CLALL>

        PathName = uigetdir(cd,'Select folder for Datastore');
        cd(PathName);
        load([PathName,'\config.mat'])
        load([PathName,'\iter.mat'])
end
%%

clc

% some quick checks
if lengthTW_Mean > lengthTW_Filt
    error('lengthTW_Mean must be <= lengthTW_Filt.')
end

if any(extractfield(stimParam(2).config,'chan') == extractfield(stimParam(3).config,'chan'),'all')
    sameS1S2 = 1;
else
    sameS1S2 = 0;
end

activeBox = ~cellfun(@isempty,{stimParam(1).config.subject});
for i = 1:4
    a = -(~any(cellfun(@(x) x == 1,activeChan(:,3,i)) == 1));
    if activeBox(1,i) == 1
        activeBox(1,i) = activeBox(1,i)+a;
    end
end

for i = 2:3
    for ii = 1:4
        if stimParam(i).config(ii).chan(1,1) > 0
            if isempty(stimParam(1).config(ii).subject) || activeChan{stimParam(i).config(ii).chan(1,1)+1,3,ii} == 0
                stimParam(i).config(ii).chan = [nan,nan];
            end
        end
    end
end

% CPU core roles
intanIO = 1;
thetaCalc = 2;
deltaCalc = 3;
predictModel = 4:4+sum(activeBox)-1;
dataSaver = predictModel(1,end)+1;
camReader = dataSaver+1;
camSaver = camReader+1;
frameLabeler = camSaver+1;

% Send tags
coreCom = 1;
startSig = 2;
timerSig = 3;
ampData = 4;
modelOut = 5;
videoData = 6;
diskFull = 7;

% Stim F keys
fKeyStim = [1:4;5:8]';

% Color bank for target labels
CB = [255,0,0]; % Red
CB(2,:) = [0,255,0]; %Green
CB(3,:) = [0,0,255]; % Blue
CB(4,:) = [255,255,0]; % Yellow
CB(5,:) = [255,0,255]; % Magenta
CB(6,:) = [0,255,255]; % Cyan
CB(7,:) = [255,102,0]; % Orange
CB(8,:) = [138,43,226]; % Violet
CB(9,:) = [64,224,208]; % Turquoise
CB(10,:) = [93,60,33]; % Chocolate

textOut = parallel.pool.DataQueue;
afterEach(textOut,@disp);

stopBox = parallel.pool.DataQueue;
afterEach(stopBox,@makeBox);

if showVideo == 1
    global videoOut %#ok<GVMIS,TLEV>
    videoOut = parallel.pool.DataQueue;
    afterEach(videoOut,@plotFrame);
end

if showData == 1
    dataOut = parallel.pool.DataQueue;
    afterEach(dataOut,@plotData);
end

if usePrediction == 1
    labelOut = parallel.pool.DataQueue;
    afterEach(labelOut,@plotLabel);
end

if labelFrames == 1
    load([CNNPathName,'\netCNN.mat'])
    inputSize = netCNN.Layers(1,1).InputSize(1,1:2);
    numTarget = netCNN.Layers(end-1,1).NumFilters;
end

vStartTime = (StartTime(1,1)*3600)+(StartTime(1,2)*60);
for Iter = startIter:size(stimParam(1).config(1).schedule,2)

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

    for I = 1:1+usePBL

        Now = clock; %#ok<CLOCK>
        Time = [num2str(Now(1,2)),'-',num2str(Now(1,3)),'-',num2str(Now(1,1)),'(',num2str(Now(1,4)),'.',num2str(Now(1,5)),')'];
        
        currentBlock = cell(1,1,4);
        currentTimes = nan(1,4);
        for i = 1:4
            if activeBox(1,i) == 1
                if I == 1
                    numBlock = length(stimParam(1).config(i).schedule(:,Iter));
                    currentBlock(1:numBlock,1,i) = stimParam(1).config(i).schedule(:,Iter);
                    for ii = 1:length(stimParam(1).config(i).duration)
                        currentTimes(ii,i) = stimParam(1).config(i).duration(ii,1);
                    end
                else
                    currentBlock{1,1,i} = 'BL';
                    currentTimes(1,:) = durPBL;
                end
            else
                
                for ii = 1:length(stimParam(1).config(i).duration)
                    currentBlock(ii,1,i) = {'N/A'};
                    currentTimes(ii,i) = nan;
                end
            end
        end
        currentTimes = mean(currentTimes,2,'omitnan');
        
        if I == 1
            SB = startBlock;
            seshPathName = [PathName,'\',Time];
        else
            SB = 1;
            seshPathName = [PathName,'\',Time,'_BL'];
        end            
    
        seshStim = zeros(size(currentBlock,3),3,size(currentBlock,1));
        for i = 1:size(currentBlock,1)
    
            % Check stimNum
            seshStim(squeeze(contains(currentBlock(i,1,:),'S1')),1,i) = 1;
            seshStim(squeeze(contains(currentBlock(i,1,:),'S2')),1,i) = 2;
            
            % Check stimTrig
            seshStim(squeeze(contains(currentBlock(i,1,:),'O')),2,i) = 1;
            seshStim(squeeze(contains(currentBlock(i,1,:),'tdC')),2,i) = 2;
            seshStim(squeeze(contains(currentBlock(i,1,:),'ptC')),2,i) = 3;
            seshStim(squeeze(contains(currentBlock(i,1,:),'dtC')),2,i) = 4;
            seshStim(squeeze(contains(currentBlock(i,1,:),'ipC')),2,i) = 5;

            % Check if sham
            seshStim(squeeze(contains(currentBlock(i,1,:),'X')),3,i) = 1;
        end
    
        if useRHS == 1 && sameS1S2 == 1
            [~,~,reStart] = ind2sub(size(seshStim(:,1,:)),find(seshStim(:,1,:) == 2,1,'first'));
            if isempty(reStart)
                reStart = 0;
            end
        else
            reStart = 0;
        end
        
        eval(['mkdir ',seshPathName]);
        send(textOut,['Directory set to: ',seshPathName])
    
        if isempty(gcp('nocreate')) == 1
            parpool('local',frameLabeler);
            disp('Connecting to hardware...');
        end
        
        spmd(frameLabeler)
    
            mpiSettings('DeadlockDetection','off');    
     
            if spmdIndex == intanIO
    
                comIssue = cell(8,2);
                comIssue{1,1} = 'RHS Commands';
                comIssue{2,1} = 'RHS Data';
                comIssue{3,1} = 'RHD Commands';
                comIssue{4,1} = 'Arduino';
                comIssue{5,1} = 'dataSaver';
                comIssue{6,1} = 'camReader';
                comIssue{7,1} = 'camSaver';
                comIssue{8,1} = 'frameLabeler';
                comIssue(:,4) = {0};
                
                % Connect and test connection to RHS commands
                
                if useRHS == 1
                    t = 0;
                    while 1
                        try
                        
                            intanSend_RHS = tcpclient(intanIP_RHS,intanPort1_RHS);
                            write(intanSend_RHS,uint8('execute RescanPorts;'));
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
                                    comIssue{1,2} = 1;
                                    break
                                end
                            end
                        end
                    end
                else
                    intanSend_RHS = 0;
                end
    
                % Connect to RHS data output and test connection
                if useRHS == 1
    
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
                                    comIssue{2,2} = 1;
                                    break
                                end
                            end
                        end
                    end
                else
                    intanRead_RHS = 0;
                end
                
                % Connect and test connection to RHD commands
                if useRHD == 1
                    t = 0;
                    while 1
                        try                        
                            intanSend_RHD = tcpclient(intanIP_RHD,intanPort1_RHD);
                            write(intanSend_RHD,uint8('execute RescanPorts;'));
                            send(textOut,'Connected to Intan RHD commands')
                            break
                        catch
                            try                               
                                write(intanSend_RHD,uint8('execute RescanPorts;'));
                                send(textOut,'Connected to Intan RHD commands')
                                break
                            catch
                                send(textOut,'Issue connecting to Intan RHS commands')
                                pause(1)
                                t = t+1;
                                if t == 60
                                    send(textOut,'Connection to Intan timed out')
                                    comIssue{3,2} = 1;
                                    break
                                end
                            end
                        end
                    end
                else
                    intanSend_RHD = 0;
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
                                comIssue{4,2} = 1;
                                break
                            end
                        end
                    end
                end
                
                % Make sure all cores started correctly
                send(textOut,'Checking hardware');
                comIssue{5,2} = spmdReceive(dataSaver,coreCom);
                comIssue{6,2} = spmdReceive(camReader,coreCom);
                comIssue{7,2} = spmdReceive(camSaver,coreCom);
                comIssue{8,2} = spmdReceive(frameLabeler,coreCom);
                allClear = sum(cell2mat(comIssue(:,2)));
                if allClear == 0
                    send(textOut,'Hardware check passed');
                end
               
                spmdSend(allClear,[dataSaver,camReader,camSaver,frameLabeler],coreCom);
         
                if allClear == 0
                     
                    % Set up RHS channel parameters
                    if intanSend_RHS ~= 0
                
                        write(intanSend_RHS,uint8(['set Filename.Path ',seshPathName,';']));
                        write(intanSend_RHS,uint8('set FileFormat OneFilePerSignalType;'));
                        write(intanSend_RHS,uint8('set Filename.BaseFilename RHS;'));
      
                        % Enable TCP recording for active channels
                        for i = 1:size(activeChan,3)
                            for ii = 1:size(activeChan,1)
                                 if activeChan{ii,3,i} == 1
                                    %write(intanSend_RHS,uint8(['set ',activeChan{ii,1,i},'.tcpdataoutputenabled true;']));
                                end
                            end
                        end
                         disp('9')
                        % Set up stim parameters
                        if I == 1                                            
                            for i = 2:size(stimParam,2)
                                if i == 2
                                    setStim(intanSend_RHS,1,stimParam(2).config)
                                else
                                    if sameS1S2 == 0
                                        setStim(intanSend_RHS,1,stimParam(3).config)
                                    end
                                end
                            end
                        end
                        send(textOut,'Uploading stim parameters');
                        write(intanSend_RHS,uint8('execute uploadstimparameters;'));
                          
                        tic
                    end
                    
                    % Set file info for RHD
                    if intanSend_RHD ~= 0
                        write(intanSend_RHD,uint8(['set Filename.Path ',seshPathName,';']));
                        write(intanSend_RHD,uint8('set FileFormat OneFilePerSignalType;'));
                        write(intanSend_RHD,uint8('set Filename.BaseFilename RHD;'));
                    end
                    
                    % Make pollable queue for stopBox
                    qW = parallel.pool.PollableDataQueue;
    
                    % Create cell to hold all block titles
                    blockText = cell(size(currentBlock,1),size(currentBlock,3));
                    for i = 1:size(currentBlock,1)
                        for ii = 1:size(currentBlock,3)
                            if strcmp(currentBlock{i,ii},'N/A') == 0
                                if strcmp(currentBlock{i,1,ii},'BL') == 1
                                    blockText{i,ii} = ['Box ',num2str(boxRHS(ii,1)),' Baseline'];
                                elseif contains(currentBlock{i,1,ii},'O') == 1
                                    blockText{i,ii} = ['Box ',num2str(boxRHS(ii,1)),' Open Loop'];
                                elseif contains(currentBlock{i,1,ii},'tdC') == 1
                                    blockText{i,ii} = ['Box ',num2str(boxRHS(ii,1)),' Theta/delta trigger closed loop'];
                                elseif contains(currentBlock{i,1,ii},'ptC') == 1
                                    blockText{i,ii} = ['Box ',num2str(boxRHS(ii,1)),' Peak theta/delta trigger closed loop'];
                                elseif contains(currentBlock{i,1,ii},'dtC') == 1
                                    blockText{i,ii} = ['Box ',num2str(boxRHS(ii,1)),' Decay theta/delta trigger closed loop'];
                                elseif contains(currentBlock{i,1,ii},'ipC') == 1
                                    blockText{i,ii} = ['Box ',num2str(boxRHS(ii,1)),' Intelligent prediction closed loop'];
                                elseif strcmp(currentBlock{i,1,ii},'SX') == 1
                                    blockText{i,ii} = ['Box ',num2str(boxRHS(ii,1)),' Sham Stim'];
                                end
                            else
                                 blockText{i,ii} = ['Box ',num2str(boxRHS(ii,1)),' Empty'];
                            end
                        end
                    end
            
                    % Calculations for accurate parsing
                    waveformBytesPerFrame = 4+2*numChanStream;
                    waveformBytesPerBlock = 128*waveformBytesPerFrame+4;
                    blocksToRead_hs1 = 78;
                    blocksToRead_hs2 = 78;
                    blocksToRead_hs2(2,1) = 78;
                    blocksToRead_hs2(3,1) = 78;
                    blocksToRead_hs2(4,1) = 79;
                    bytesToRead_hs1 = waveformBytesPerBlock*blocksToRead_hs1;
                    bytesToRead_hs2 = waveformBytesPerBlock*blocksToRead_hs2;
                    blockCycle = 1;
    
                    % Pre-allocate counters and data stores
                    countP = 1;
                    sendTry = 0;
                    timeWindow = zeros(size(activeChan,1),lengthTW_Filt,size(activeChan,3));
                    tdSend = zeros(size(tdChan,1),lengthTW_Filt);
                    dataSec = zeros(5,1,length(tdChan));
                    allData = zeros(5,lengthData,length(activeBox));
                    labels = zeros(size(activeChan,3),1);
                    allLabels = zeros(size(activeChan,3),lengthData);
                    stimDat = uint8(zeros(8,dSF));
                    saveCell = cell(1,5);
                    
                    % Warm up calc cores
                    spmdSend(timeWindow(:,:,1),[thetaCalc,deltaCalc],ampData);
                    theta = spmdReceive(thetaCalc,ampData);
                    delta = spmdReceive(deltaCalc,ampData);
    
                    if showData == 1
                        send(dataOut,{allData,activeBox,WLc,iWL})
                    end
    
                    if usePrediction == 1
                        send(labelOut,{allLabels,activeBox,WLc,iWL})
                    end
                    
                    % Set index of TCP data
                    count = 1;
                    idxRead = zeros(128,2);
                    for i3 = 1:size(activeChan,3)
                        for i4 = 1:size(activeChan,1)
                            if activeChan{i4,3,i3} == 1
                                idxRead(count,1) = i3;
                                idxRead(count,2) = str2double(activeChan{i4,1,i3}(3:end))+1;
                                count = count+1;
                            end
                        end
                    end
    
                    % Wait for stim upload to complete
                    pause(15-toc)
                    
                    % sync with system clock
                    send(textOut,'Syncing hardware');
                    while getTime < vStartTime
                        pause(0.001);
                    end
    
                    % Open stop box
                    send(stopBox,qW);

                    stopExp = 0;
                    diskWarn = 0;
                    
                    % Start camReader core
                    spmdSend(1,camReader,startSig)
                    
                    % Start recording
                    if intanSend_RHS ~= 0
                        if recIntan == 0
                            write(intanSend_RHS,uint8('set runMode run;'));
                        else
                            write(intanSend_RHS,uint8('set runMode record;'));
                        end
                    end
                    if intanSend_RHD ~= 0
                        write(intanSend_RHD,uint8('set runMode record;'));
                    end
                    
                    LED = 0;
                    countLED = 0;
                    for i = SB:size(currentBlock,1)

                        stimCounter = zeros(size(currentBlock,3),2);
                        for ii = 1:currentTimes(i,1)
    
                            % Pause if the recording needs to be stopped to change the parameters
                            if i == reStart && ii == 1
                                
                                countP = countP+1;
                                send(textOut,'----------')
                                send(textOut,'Pausing recording to update parameters...')
                                write(intanSend_RHS,uint8('set runMode stop;'));
                                pause(0.1)
                                write(intanSend_RHS,uint8('execute RescanPorts;'));
                                write(intanSend_RHS,uint8(['set Filename.BaseFilename RHS_P',num2str(countP),';']));
                                setStim(intanSend_RHS,1,stimParam(3).config);
                                write(intanSend_RHS,uint8('execute uploadstimparameters;'));
                                flush(intanRead_RHS)
                                pause(15)
                                write(intanSend_RHS,uint8('set runMode record;'));
                            end
    
                            % Wait for data
                            while intanRead_RHS.NumBytesAvailable < waveformBytesPerBlock
                                pause(0.0001)
                            end
    
                            % Send sync signal to camReader
                            sendTry = 0;
                            while 1
                                try
                                    spmdSend(1,camReader,timerSig);
                                    break
                                catch ME
                                    pause(0.01)
                                    sendTry = sendTry+1;
                                    warning('Worker %dIntanIO failed to send start signal to CamReader: %s\n', spmdIndex, ME.message);                                    
                                end
                                if sendTry == 3
                                    break
                                end
                            end
                            if sendTry == 3
                                stopExp = 1;
                                break
                            end
                            
                            % Turn on block marker on first iteration
                            if ii == 1
                                writeDigitalPin(Ard,ArdPin{2,1},1);
                            end
                            
                            % Turn on LED if interval is reached
                            if countLED == LED
                                writeDigitalPin(Ard,ArdPin{1,1},1);
                            end

                            ticLED = tic;
                            
                            % First iteration set up, done during the LED
                            % pause to save time
                            if ii == 1
                                
                                % Prealocate Stim trigger counters
                                nextISI = zeros(size(currentBlock,3),1);
                                for i3 = 1:length(nextISI)
                                    if seshStim(i3,1,i) > 0
                                        nextISI(i3,1) = datasample(ISI{seshStim(i3,2,i),seshStim(i3,1,i)},1);
                                    end
                                end
                                countISI = zeros(size(currentBlock,3),2);
                                stimCount = zeros(size(currentBlock,3),2);
                                stimTrig = zeros(size(currentBlock,3),2);
                                ptC_Count = zeros(size(currentBlock,3),2);
                                dtC_Count = zeros(size(currentBlock,3),2);
                                pastPeak = zeros(size(currentBlock,3),2);
                                labelCount = zeros(size(currentBlock,3),2);
                                
                                % Display current block for each box
                                send(textOut,'----------')
                                for i3 = 1:size(currentBlock,3)
                                    send(textOut,blockText{i,i3})
                                end
                            end

                            % Pause for 0.2s
                            while toc(ticLED) < 0.2
                                pause(0.0001)
                            end

                            % Turn off block marker
                            if ii == 1
                                writeDigitalPin(Ard,ArdPin{2,1},0);
                            end

                            % Turn off LED and reset counter or increment
                            % counter
                            if countLED == LED
                                writeDigitalPin(Ard,ArdPin{1,1},0);
                                LED = randsample(LEDtime,1);
                                countLED = 0;
                            else
                                countLED = countLED+1;
                            end                              
    
                            % Check if stop box has been closed
                            send(stopBox,1);
                            boxClose = poll(qW);
                            if ~isempty(boxClose)
                                stopExp = 1;
                                break
                            end
                            
                            % If the disk is full stop
                            if spmdProbe(camSaver,diskFull)
                                spmdReceive(camSaver,diskFull);
                                diskWarn = 1;
                                break
                            end

                            % Wait till 0.5s on Intan has passed and read                          
                            while intanRead_RHS.NumBytesAvailable < bytesToRead_hs1
                                pause(0.0001)
                            end

                            % Read data and decode
                            waveformArray = read(intanRead_RHS,bytesToRead_hs1);
                            [amplifierData,timeStamps,offset] = byte2doubleConv(waveformArray,numChanStream,blocksToRead_hs1,sF,dSF,false);

                            % If read cycle gets offset find the offset and correct
                            if offset > 0
                                send(textOut,['Read offset of ',num2str(offset),' detected.'])
                                waveformArray(1:offset) = [];
                                waveformArray = [waveformArray,read(intanRead_RHS,offset)];  %#ok<AGROW>
                                [amplifierData,timeStamps,offset] = byte2doubleConv(waveformArray,numChanStream,blocksToRead_hs1,sF,dSF,false);
                                send(textOut,'Corrected')
                            end
                            numSample = size(amplifierData,2);
                            
                            % Sort amp data by subject into an array [channel,data,subject]
                            readWindow = zeros(size(activeChan,1),numSample,size(activeChan,3));
                            for i3 = 1:numChanStream
                                readWindow(idxRead(i3,2),:,idxRead(i3,1)) = amplifierData(i3,:);
                            end
                            
                            % Update time window
                            timeWindow = circshift(timeWindow,numSample,2);
                            timeWindow(:,1:numSample,:) = readWindow;
                            
                            % Pull channels for T/D
                            for i3 = 1:size(tdChan,1)
                                if ~isnan(tdChan(i3,1))
                                    tdSend(i3,:) = timeWindow(tdChan(i3,1)+1,:,i3);
                                end
                            end
    
                            % Send timeWindow to calc cores for parallel computation
                            spmdSend(tdSend,[thetaCalc,deltaCalc],ampData);
                            if usePrediction == 1

                                iCore = 1;
                                for i3 = 1:length(activeBox)
                                    if activeBox(1,i3) == 1
                                        spmdSend({i3,timeWindow(:,1:inputDur,i3)},predictModel(1,iCore),ampData);
                                        iCore = iCore+1;
                                    end
                                end
                            end
                            
                            % Receive T/D data from calc cores
                            try
                                theta = spmdReceive(thetaCalc,ampData);
                            catch
                                theta = 0;
                            end
                            try
                                delta = spmdReceive(deltaCalc,ampData);
                            catch
                                delta = nan;
                            end                            

                            % Get T/D ratio
                            ToD = theta./delta;
                            
                            % Receive labels from predict cores
                            if usePrediction == 1
                                for i3 = predictModel
                                    try
                                        modelData = spmdReceive('any',modelOut);
                                        labels(modelData(1,1),1) = modelData(1,2);
                                    catch
                                        labels(modelData(1,1),1) = 0;
                                    end
                                end
                                allLabels = circshift(allLabels,1,2);
                                allLabels(:,1) = labels;
                                send(labelOut,{allLabels,activeBox,WLc,iWL})
                            end
    
                            % Send to datasaver core if enabled
                            saveCell = cell(1,5);
                            if saveDS == 1
                                saveCell{1,1} = amplifierData;
                            end
                            if saveTS == 1
                                saveCell{1,2} = timeStamps;
                            end
                            if saveTD == 1
                                saveCell{1,3} = ToD;
                            end
                            if usePrediction == 1
                                saveCell{1,4} = labels;
                            end
                            spmdSend(saveCell,dataSaver,ampData);
                                
                            % Decide which channels to stim and when
                            stim = zeros(size(currentBlock,3),2);
                            for i3 = 1:size(seshStim,1)
                                if activeBox(1,i3) == 1
    
                                    dataSec(1,1,i3) = ToD(i3,1);
                                    if seshStim(i3,2,i) == 0 % BL
                                        dataSec(2:5,1,i3) = 0;
                                        continue
                                    else
                                        if countISI(i3,seshStim(i3,1,i)) == 0
                                            if stimCount(i3,seshStim(i3,1,i)) < stimNum(1,seshStim(i3,1,i))
                                                if seshStim(i3,2,i) == 1 % O
                                                    stim(i3,seshStim(i3,1,i)) = 1;
                                                    stimTrig(i3,seshStim(i3,1,i)) = 1;
                                                    nextISI(i3,1) = datasample(ISI{1,seshStim(i3,1,i)},1)+stimDur(i3,seshStim(i3,1,i));
                                                elseif seshStim(i3,2,i) == 2 % tdC
        
                                                    stim(i3,seshStim(i3,1,i)) = (ToD(i3,1) > 0 && ToD(i3,1) < threshold_tdC);
                                                    if stim(i3,seshStim(i3,1,i)) == 1
                                                        stimTrig(i3,seshStim(i3,1,i)) = 1;
                                                        nextISI(i3,1) = datasample(ISI{2,seshStim(i3,1,i)},1)+stimDur(i3,seshStim(i3,1,i));
                                                    end
                                                elseif seshStim(i3,2,i) == 3 % ptC
                                                    if ToD(i3,1) >= threshold_ptC
                                    
                                                        ptC_Count(i3,seshStim(i3,1,i)) = ptC_Count(i3,seshStim(i3,1,i))+1;
                                                        if ptC_Count(i3,seshStim(i3,1,i)) >= sustain_ptC
                                                            stim(i3,seshStim(i3,1,i)) = 1;
                                                            stimTrig(i3,seshStim(i3,1,i)) = 1;
                                                            nextISI(i3,1) = datasample(ISI{3,seshStim(i3,1,i)},1)+stimDur(i3,seshStim(i3,1,i));
                                                            ptC_Count(i3,seshStim(i3,1,i)) = 0;
                                                        end
                                                    else
                                                        ptC_Count(i3,seshStim(i3,1,i)) = 0;
                                                    end
                                                elseif seshStim(i3,2,i) == 4 % dtC
                                                    if ToD(i3,1) >= threshold_ptC && pastPeak(i3,seshStim(i3,1,i)) == 0
                                                        pastPeak(i3,seshStim(i3,1,i)) = 1;                    
                                                    elseif ToD(i3,1) > 0 && ToD(i3,1) < threshold_dtC && pastPeak(i3,seshStim(i3,1,i)) == 1
                                                        
                                                        dtC_Count(i3,seshStim(i3,1,i)) = dtC_Count(i3,seshStim(i3,1,i))+1;
                                                        if dtC_Count(i3,seshStim(i3,1,i)) >= sustain_dtC                            
                                                            stim(i3,seshStim(i3,1,i)) = 1;
                                                            stimTrig(i3,seshStim(i3,1,i)) = 1;
                                                            nextISI(i3,1) = datasample(ISI{4,seshStim(i3,1,i)},1)+stimDur(i3,1);
                                                            pastPeak(i3,seshStim(i3,1,i)) = 0;
                                                            dtC_Count(i3,seshStim(i3,1,i)) = 0;
                                                        end
                                                    else
                                                        dtC_Count(i3,seshStim(i3,1,i)) = 0;
                                                    end
                                                elseif seshStim(i3,2,i) == 5 % ipC
                                                    if labels(i3,1) == stimClass
    
                                                        labelCount(i3,1) = labelCount(i3,1)+1;
                                                        if labelCount(i3,1) >= numConsec
                                                            stim(i3,seshStim(i3,1,i)) = 1;
                                                            stimTrig(i3,seshStim(i3,1,i)) = 1;
                                                            nextISI(i3,1) = datasample(ISI{5,seshStim(i3,1,i)},1)+stimDur(i3,1);
                                                            labelCount(i3,1) = 0;
                                                        end
                                                    else
                                                        if labelCount(i3,1) > 0
                                                            labelCount(i3,1) = labelCount(i3,1)-1;
                                                        end
                                                    end
                                                end
    
                                                if stim(i3,seshStim(i3,1,i)) == 1
                                                    stimCounter(i3,seshStim(i3,1,i)) = pulsNum(i3,2,seshStim(i3,1,i));
                                                end
                                            end
                                        end
                                        if stimTrig(i3,seshStim(i3,1,i)) == 1
                                            countISI(i3,seshStim(i3,1,i)) = countISI(i3,seshStim(i3,1,i))+1;
                                        end
                                        if countISI(i3,seshStim(i3,1,i)) >= nextISI(i3,1)
                                            countISI(i3,seshStim(i3,1,i)) = 0;
                                            stimTrig(i3,seshStim(i3,1,i)) = 0;
                                        end
                                    end                            
                                    dataSec(2,1,i3) = (nextISI(i3,1)-countISI(i3,seshStim(i3,1,i)))/nextISI(i3,1);
                                    if isnan(dataSec(2,1,i3)) || dataSec(2,1,i3) == 1
                                        dataSec(2,1,i3) = 0;
                                    end
                                    dataSec(3,1,i3) = ptC_Count(i3,seshStim(i3,1,i))/sustain_ptC;
                                    dataSec(4,1,i3) = pastPeak(i3,seshStim(i3,1,i));
                                    dataSec(5,1,i3) = dtC_Count(i3,seshStim(i3,1,i))/sustain_dtC;
                                end
                            end
                            stimCount = stimCount+stim;
                            
                            % Update data plots
                            if showData == 1
                                allData = circshift(allData,1,2);
                                allData(:,1,:) = dataSec;
                                send(dataOut,{allData,activeBox,WLc,iWL})
                            end
    
                            % Make command array
                            stimOut = '';
                            idxStim = 0;
                            for i3 = 1:size(stim,2)
                                for i4 = 1:size(stim,1)
                                    idxStim = idxStim+1;
                                    if seshStim(i4,3,i) == 0
                                        if stimCounter(i4,i3) > 0
                                            stimOut = [stimOut,['execute ManualStimTriggerPulse F',num2str(stimParam(i3+1).config(i4).fKey),';']]; %#ok<AGROW>
                                            stimCounter(i4,i3) = stimCounter(i4,i3)-1;
                                        end
                                        stimDat(idxStim,round(dSF*0.85)) = 1;
                                    end
                                end
                            end
                            
                           % Wait till 0.75s on Intan has passed
                            while intanRead_RHS.NumBytesAvailable < waveformBytesPerBlock*blocksToRead_hs1/2
                                pause(0.0001)
                            end
    
                            % Execute f key presses to Intan
                            if ~isempty(stimOut)
                                write(intanSend_RHS,uint8(stimOut));
                            end
    
                            % Wait till 1.0s on Intan has passed
                            while intanRead_RHS.NumBytesAvailable < bytesToRead_hs2(blockCycle,1)
                                pause(0.0001)
                            end

                            % Send sync signal to camReader
                            sendTry = 0;
                            while 1
                                try
                                    spmdSend(1,camReader,timerSig);
                                    break
                                catch ME
                                    pause(0.01)
                                    sendTry = sendTry+1;
                                    warning('Worker %dIntanIO failed to send stop signal to CamReader: %s\n', spmdIndex, ME.message);                                    
                                end
                                if sendTry == 3
                                    break
                                end
                            end
                            if sendTry == 3
                                stopExp = 1;
                                break
                            end

                            % Read data and decode
                            waveformArray = read(intanRead_RHS,bytesToRead_hs2(blockCycle,1));
                            [amplifierData,timeStamps,~] = byte2doubleConv(waveformArray,numChanStream,blocksToRead_hs2(blockCycle,1),sF,dSF,false);
                            numSample = size(amplifierData,2);
                            
                            % Send to datasaver core if enabled
                            saveCell = cell(1,5);
                            if saveDS == 1
                                saveCell{1,1} = amplifierData;
                            end
                            if saveTS == 1
                                saveCell{1,2} = timeStamps;
                            end
                            if saveStim == 1
                                saveCell{1,5} = stimDat;
                                stimDat(:) = 0;
                            end
                            spmdSend(saveCell,dataSaver,ampData);
                            
                            %fix this for multi subject
                            %
                            % % Sort amp data by subject into an array [channel,data,subject]
                            % readWindow = zeros(size(activeChan,1),numSample,size(activeChan,3));
                            % for i3 = 1:numChanStream
                            %     readWindow(idxRead(i3,2),:,idxRead(i3,1)) = amplifierData(i3,:);
                            % end
                            

                            readWindow = amplifierData;
                            % Update time window
                            timeWindow = circshift(timeWindow,numSample,2);
                            timeWindow(:,1:numSample,:) = readWindow;

                            blockCycle = blockCycle+1;
                            if blockCycle >= 5
                                blockCycle = 1;
                            end
                        end

                        if stopExp == 1 || diskWarn == 1
                            break
                        end
                    end                    
                    
                    if intanSend_RHS ~= 0
                        write(intanSend_RHS,uint8('set runMode stop;'));
                    end
                    if intanSend_RHD ~= 0
                        write(intanSend_RHD,uint8('set runMode stop;'));
                    end

                    if stopExp == 0 || diskWarn == 1
                        send(stopBox,[]);
                    end

                    spmdSend([],camReader,timerSig)
                    flush(intanRead_RHS)

                    send(textOut,'----------')
                    send(textOut,'Recording finished.')
    
                    % Send stop to datasaver core if enabled
                    if saveTD == 1
                        spmdSend([],dataSaver,ampData);
                    end
                
                    % Disable stim and TCP channels after recording finishes
                    if intanSend_RHS ~= 0
    
                        for i = 2:3
                            setStim(intanSend_RHS,0,stimParam(i).config)
                        end
    
                        % Disable TCP recording for active channels
                        for i = 1:size(activeChan,3)
                            for ii = 1:size(activeChan,1)
                                 if activeChan{ii,3,i} == 1
                                    write(intanSend_RHS,uint8(['set ',activeChan{ii,1,i},'.tcpdataoutputenabled false;']));
                                end
                            end
                        end
    
                        write(intanSend_RHS,uint8('execute uploadstimparameters;'));
                        pause(15)
                    end
                end
                spmdSend([],[thetaCalc,deltaCalc,predictModel],ampData);
            elseif spmdIndex == thetaCalc                
                while 1
                    try
                        if spmdProbe(intanIO,ampData) == 1
                            
                            try
                                timeWindow = spmdReceive(intanIO,ampData);
                            catch
                                % Do nothing
                            end
    
                            if ~isempty(timeWindow)
                                % Calculate mean theta amp
                                theta = abs(hilbert(bandpass(timeWindow',[5,12],dSF)));
                                thetaMean = mean(theta(1:lengthTW_Mean,:))';
                                spmdSend(thetaMean,intanIO,ampData)
                            else
                                break
                            end
                        else
                            pause(0.001)
                        end
                    catch
                        spmdSend(0,intanIO,ampData)
                        send(textOut,'Unknown packet')
                    end
                end
            elseif spmdIndex == deltaCalc
                while 1
                    try
                        if spmdProbe(intanIO,ampData) == 1
                            
                            try
                                timeWindow = spmdReceive(intanIO,ampData);
                            catch
                                % Do nothing
                            end
    
                            if ~isempty(timeWindow)
                                % Calculate mean delta amp
                                delta = abs(hilbert(bandpass(timeWindow',[1,4],dSF)));
                                deltaMean = mean(delta(1:lengthTW_Mean,:))';
                                spmdSend(deltaMean,intanIO,ampData)
                            else
                                break
                            end
                        else
                            pause(0.001)
                        end
                    catch
                        spmdSend(0,intanIO,ampData)
                        send(textOut,'Unknown packet')
                    end
                end
            elseif any(spmdIndex == predictModel)
                
                boxID = find(predictModel == spmdIndex);
                
                if usePrediction == 1
                    model = load(modelPath{boxID,1});
                    nChannel = length(model.ops.ch_subj);
                    freq = model.ops.freqs;
                    nFreq = length(freq);
                    feat = nan(1,nFreq*nChannel);

                    modelOPs.freqs = freq;
                    modelOPs.durFeat = model.ops.durFeat; % 4s feature bins
                    modelOPs.art_thres = 5e4; % NOTE probably needs to be decreased
                    modelOPs.features = model.ops.features;
                    modelOPs.ch_subj = model.ops.ch_subj;
                    modelOPs.Fs = dSF;
                end
    
                while 1                    
                    if spmdProbe(intanIO,ampData) == 1
    
                        data = spmdReceive(intanIO,ampData);                        
    
                        if ~isempty(data)

                            boxNum = data{1,1};
                            data = data{1,2};
        
                            data = data/0.195;

                          

                            %this will need to be changed if feature
                            %changes
                            if ~any(any(abs(data) > model.ops.art_thres))
                                 feat = sm_GetDataFeature_rat(data,[],modelOPs);
                                 size(feat)
                                %run through model
                                label = predict(model.rusTree,feat);
                                spmdSend([boxNum,label],intanIO,modelOut)
                            else
                                spmdSend([],intanIO,ampData)
                            end
                        else
                            break
                        end
                    else
                        pause(0.001)
                    end
                end
            elseif spmdIndex == dataSaver
    
                % Send check to hardwareTimer
                comIssue = 0;
                try
                    if saveDS == 1
                        fID_amp = fopen([seshPathName,'\amplifier_',num2str(numChan),'Ch_',num2str(dSF/1000),'KHz_int16.dat'],'w');
                    end
                    if saveTS == 1
                        fID_ts = fopen([seshPathName,'\timestamps_4Ch_',num2str(dSF/1000),'KHz_uint8.dat'],'w');
                    end
                    if saveTD == 1
                        fID_ToD = fopen([seshPathName,'\ToD_4Ch_1Hz_uint8.dat'],'w');
                    end
                    if usePrediction == 1
                        fID_labels = fopen([seshPathName,'\Labels_4Ch_1Hz_uint8.dat'],'w');
                    end
                    if saveStim == 1
                        fID_stim = fopen([seshPathName,'\Stim_8Ch_',num2str(dSF/1000),'KHz_uint8.dat'],'w');
                    end
                catch
                    comIssue = 1;
                end
                spmdSend(comIssue,intanIO,coreCom)
                comIssue = spmdReceive(intanIO,coreCom);
                
                if comIssue == 0
                    while 1
                        if spmdProbe(intanIO,ampData) == 1
                            
                            try
                                saveCell = spmdReceive(intanIO,ampData);
                            catch
                                % Do nothing
                            end

                            if ~isempty(saveCell)
                                if saveDS == 1
                                   fwrite(fID_amp,saveCell{1,1},'int16');
                                end
                                if saveTS == 1
                                   fwrite(fID_ts,saveCell{1,2},'uint8');
                                end
                                if saveTD == 1
                                    fwrite(fID_ToD,saveCell{1,3}*25,'uint8');
                                end
                                if usePrediction == 1
                                    fwrite(fID_labels,saveCell{1,4},'uint8');
                                end
                                if saveStim == 1
                                    fwrite(fID_stim,saveCell{1,5},'uint8');
                                end
                            else
                                break
                            end
                        else
                            pause(0.001)
                        end
                    end

                    if saveDS == 1
                       fclose(fID_amp);
                    end
                    if saveTS == 1
                        fclose(fID_ts);
                    end
                    if saveTD == 1
                        fclose(fID_ToD);
                    end
                    if usePrediction == 1
                        fclose(fID_labels);
                    end
                    if saveStim == 1
                        fclose(fID_stim);
                    end
                end
            elseif spmdIndex == camReader
                
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
                spmdSend(comIssue,intanIO,coreCom)
                comIssue = spmdReceive(intanIO,coreCom);                
    
                if comIssue == 0
                    
                    ck = zeros(1,numCam);
                    frame = uint8(zeros(Dim(1,1),Dim(1,2),3,numCam)); 
                    errorCode = uint8([0,255,0;255,0,255;0,255,0]);
    
                    if showVideo == 1                        
                        send(videoOut,{frame+255,WLc,iWL}); %#ok<SPGV>                        
                    end
                    
                    while spmdProbe(intanIO,timerSig) == 0
                        try
                            spmdReceive(intanIO,startSig);
                            break
                        catch
                        end
                    end
    
                    while 1
                        if spmdProbe(intanIO,timerSig) == 1
    
                            % wait for timer signal
                            try
                                go = spmdReceive(intanIO,timerSig);
                            catch
                                go = 1;
                            end
    
                            if ~isempty(go)
                                while spmdProbe(intanIO,timerSig) == 0
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

                                go = spmdReceive(intanIO,timerSig);
                                if isempty(go)
                                    break
                                end
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
                            writerObj{ii,1} = VideoWriter([seshPathName,'\',VidFileName{ii,1}],'MPEG-4');
                            writerObj{ii,1}.FrameRate = fps;
                            open(writerObj{ii,1});
                        end
            
                        if labelFrames == 1
                            
                            writerObj_Pose = cell(numCam,1);
                            for ii = 1:numCam    
                                VidFileName{ii,1} = [WL{camChoice(1,ii),3}(1,1:end-4),'_Pose.mp4'];
                                writerObj_Pose{ii,1} = VideoWriter([seshPathName,'\',VidFileName{ii,1}],'MPEG-4');
                                writerObj_Pose{ii,1}.FrameRate = fps;
                                open(writerObj_Pose{ii,1});
                            end
                        end
                    catch
                        send(textOut,'Could not create video files')
                        comIssue = 1;
                    end
                end
                spmdSend(comIssue,intanIO,coreCom)
                comIssue = spmdReceive(intanIO,coreCom);
                
                if comIssue == 0
    
                    frames = uint8(zeros(Dim(1,1),Dim(1,2),3,numCam,40));
    
                    while 1

                        numFrames = 0;                            
                        while 1
                            if spmdProbe(camReader,videoData) == 1
                                
                                try
                                    frame = spmdReceive(camReader,videoData);
                                catch
                                    continue
                                end

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
                                pause(0.001)
                            end
                        end
    
                        if isempty(frame)
                            break
                        end
                        
                        if numFrames > 0
                            syncFrames = frames(:,:,:,:,round(linspace(1,numFrames,fps)));
                        end
    
                        if labelFrames == 1
                            spmdSend(syncFrames,frameLabeler,videoData);
                        end                    
    
                        if RecVid == 1
                            for i3 = 1:fps
                                for i4 = 1:numCam    
                                    if all(syncFrames(1:3,1:3,1,i4,i3) ~= errorCode,'all')
                                        try
                                            writeVideo(writerObj{i4,1},syncFrames(:,:,:,i4,i3));
                                        catch
                                            spmdSend(1,intanIO,diskFull);
                                        end
                                    end
                                    if labelFrames == 1
                                        frames = spmdReceive(frameLabeler,videoData);
                                        writeVideo(writerObj_Pose{i4,1},frames(:,:,:,i4,i3));
                                    end                       
                                end
                            end
                        end
                    end

                    spmdSend([],frameLabeler,videoData);
    
                    send(textOut,'Finalizing video files.')
                    for ii = 1:numCam
                        close(writerObj{ii,1});
                    end  
                end
            elseif spmdIndex == frameLabeler                
                
                comIssue = 0;
                blobAnalyserCNN = vision.BlobAnalysis('BoundingBoxOutputPort',false,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',10);
                spmdSend(comIssue,intanIO,coreCom)
                comIssue = spmdReceive(intanIO,coreCom);
                
                if comIssue == 0
                    
                    if labelFrames == 1

                        QP = zeros(prod(Dim),2);
                        countISI = 1;
                        for i4 = 1:Dim(1,2)
                            for i5 = 1:Dim(1,1)
                                QP(countISI,:) = [i4,i5]; 
                                countISI = countISI+1;
                            end
                        end
    
                        allCamXY = cell(1,numCam);
                        allCamXY(1,:) = {zeros(1,4,numTarget)};
                    end
                   
                    count = 1;
                    while 1

                        % Wait for frames from camSaver
                        while 1
                            try
                                if spmdProbe(camSaver,videoData)
                                    break
                                else
                                    pause(0.1)
                                end
                            catch
                                continue
                            end
                        end

                        try
                            frames = spmdReceive(camSaver,videoData);
                        catch
                            continue
                        end

                        if ~isempty(frames)
        
                            framesR = double(imresize(frames,inputSize))/255;
                            YPred = predict(netCNN,framesR);
                            for ii = 1:numCam
                                for i3 = 1:numTarget
                                        
                                    P = YPred(:,:,i3,ii);
                                    P(P < 0.2) = 0;
                                    P(P > 1) = 1;
                                    mask = logical(P);
                                    [area,centroids] = blobAnalyserCNN(mask);
                                    if isempty(centroids) == 0
                
                                        [~,aMax] = max(area);
                                        xy = round(centroids(aMax,:)); 
                                        xy = RemapPoint(xy,[inputSize(1,1),inputSize(1,2)],Dim,0,[1,1]); 
                
                                        [counts,centers] = hist(P(P > 0),linspace(0.01,1,20)); %#ok<HIST> 
                                        [~,idx] = max(counts);
                                        Con = centers(idx);
                                    else
                                        xy = [nan,nan];
                                        Con = 0;
                                    end
        
                                    allCamXY{1,ii}(count,1:3,i3) = [xy(1,1),xy(1,2),Con];
                                    if count == 1
                                        allCamXY{1,ii}(count,4,i3) = 0;
                                    else
                                        allCamXY{1,ii}(count,4,i3) = pdist([allCamXY{1,ii}(count-1,1:2,i3);allCamXY{1,ii}(count,1:2,i3)],'euclidean');
                                    end
                
                                    if ~any(isnan(xy))
                                        
                                        p = nsidedpoly(1000,'Center',xy,'Radius',6);
                                        shp = alphaShape(p.Vertices);
                                        a = criticalAlpha(shp,'all-points')*2;
                                        shp.Alpha = a; 
                    
                                        IS = inShape(shp,QP);
                                        IS = find(IS == 1);
                                        
                                        for i4 = 1:length(IS) 
                                            
                                            col = QP(IS(i4,1),1);
                                            row = QP(IS(i4,1),2);                        
                                            frames(row,col,:,ii) = CB(i3,:);
                                        end
                                    end     
                                end
                            end
                            count = count+1;
                            spmdSend(camSaver,frames);
                        else
                            break
                        end
                    end
                end
            end
        end

        send(textOut,'Complete.')

        if diskWarn{1} == 1
            f = errordlg('Disk is full, free more space to continue');
            waitfor(f);
            diskWarn = 0;
        end
    
        if labelFrames == 1
            allCamXY = allCamXY{frameLabeler};
            for i = 1:numCam
    
                allXY = cell(2,numTarget+1);
                allXY(1,1:numTarget) = Targets';
                allXY{1,end} = 'Behavior';
                for ii = 1:numTarget
                    allXY{2,ii} = allCamXY{1,i}(:,:,ii);
                end
                save([seshPathName,'\',WLc{i,3}(1,1:end-4),'_PoseData.mat'],'allXY')
            end
        end
    
        tic
    
        if runAnalysis == 1
            
            send(textOut,'----------')
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
    
                            recInfo(i).fileData(ii).subject = stimParam(1).config(ii).subject;
                            
                            timestamps = cell(sum(activeBox),2);
                            x = 1;
                            y = 0;
                            for i3 = 1:size(currentTimes,1)
    
                                timestamps{i3,1} = currentBlock{i3,1,ii};
                                y = y+currentTimes(i3,1);
                                timestamps{i3,2} = [x,y];
                                x = x+currentTimes(i3,1);
                            end
                            recInfo(i).fileData(ii).blockDuration = timestamps;
    
                            recInfo(i).fileData(ii).stim1_uA = stimParam(2).config(ii).amp; 
                            recInfo(i).fileData(ii).stim1_chan = stimParam(2).config(ii).chan;
                            if any(isnan(stimParam(2).config(ii).chan))
                                recInfo(i).fileData(ii).stim1_type = 'monopole';
                            else
                                recInfo(i).fileData(ii).stim1_type = 'bipole';
                            end
    
                            recInfo(i).fileData(ii).stim2_uA = stimParam(3).config(ii).amp;
                            recInfo(i).fileData(ii).stim2_chan = stimParam(3).config(ii).chan;
                            if any(isnan(stimParam(3).config(ii).chan))
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
                L = ceil(numSamples/(sF*600)); 
                
                fileID = fopen([recInfo(i).file,'\digitalin.dat']);
                if replaceFile == 1
                    fID = fopen(recInfo(i).file+"\digitalin_1Ch_"+dSF/1000+"KHz_.dat",'w');
                end
            
                fseek(fileID,0,'bof');
                for ii = 1:L
                    temp = fread(fileID,[1,sF*600],'int16');
                    temp = temp(:,1:sF/dSF:end);
                    if replaceFile == 1
                        fwrite(fID,temp,'uint16');
                    end
                end
                fclose(fileID);
    
                if replaceFile == 1                
                    fclose(fID);
                    eval(['delete ',recInfo(i).file,'\digitalin.dat']);
                end
            
                % Stim TS
                disp('Working on stim.dat')
                s = dir([recInfo(i).file,'\stim.dat']);
                numSamples = s.bytes/(2*numChan);
                L = ceil(numSamples/(sF*600));
                stimData = nan(numChan,L*600*dSF);

                fileID = fopen([recInfo(i).file,'\stim.dat']);
                if replaceFile == 1
                    fID = fopen([recInfo(i).file,'\stim_',num2str(numChan),'Ch_',num2str(dSF/1000),'KHz_.dat'],'w');
                end
                
                countISI = 1;
                fseek(fileID,0,'bof');
                wb = waitbar(0,'Working on stim.dat');
                for ii = 1:L
                    
                    temp = fread(fileID,[numChan,sF*600],'int16');
                    temp = temp(:,1:sF/dSF:end);
                    temp(temp > 0) = 1;
                    temp(temp < 0) = 0;

                    if replaceFile == 1
                        fwrite(fID,temp,'int16');
                    end
                    
                    stimData(:,countISI:countISI+size(temp,2)-1) = temp;
                    countISI = countISI+size(temp,2);
                    waitbar(ii/L,wb,'Working on stim.dat');
                end
                idx = find(~isnan(stimData(1,:)));
                stimData = stimData(:,idx);
                close(wb);
                
                fclose(fileID);
                if replaceFile == 1
                    fclose(fID);
                    eval(['delete ',recInfo(i).file,'\stim.dat']);
                end
                
                stimTS = cell(size(activeChan,3),2);
                count = 1;
                for ii = 1:size(activeChan,3)
                    for i3 = 1:size(activeChan,1)
                        if ~isempty(recInfo(i).fileData(ii).subject) && activeChan{i3,3,ii} == 1
                            for i4 = 1:2
                                if any(stimParam(i4+1).config(ii).chan == i3-1) 
                                    temp = stimData(count,:);
                                    idx = find(diff(temp) > 0)+1;
                                    stimTS{ii,i4} = idx';
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
            end
            save([seshPathName,'\recInfo.mat'],'recInfo')
            disp('Analysis complete');
            send(textOut,'----------')
        end                                     
    
        if strcmp(copyAction,'Yes') == 1
            
            if I == 1
                copyDir = [copyPathName,'\',Time];
            else
                copyDir = [copyPathName,'\',Time,'_BL'];
            end
            
            disp(['Copying data to ',copyPathName]);
            eval(['mkdir ',copyDir]);
            copyfile(seshPathName,copyDir);
            send(textOut,'----------')
        end

        if stopExp{1} == 1
            break
        end
    end
    startIter = Iter+1;
    save([PathName,'\iter.mat'],'startIter')
end
close all
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
    Now = clock; %#ok<CLOCK>
    timeVec = (Now(1,4)*3600)+(Now(1,5)*60)+(Now(1,6));
end

function makeBox(input)
    global closeBox qSend %#ok<GVMIS>

    if ~isempty(input)
        if input ~= 1
            
            qSend = input;
            closeBox = errordlg('Press ok to stop.','Stop','non-modal');
            closeBox.Position = [1278,728,161,60];
            drawnow
        else
            if ~ishandle(closeBox)
                send(qSend,0)
            end
        end
    else
        close(closeBox)
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

                if config(i).chan(ii) > 0
                    if config(i).chan(ii) < 10
                        chanStr = [config(i).port,'-00',num2str(config(i).chan(ii))];
                    else 
                        chanStr = [config(i).port,'-0',num2str(config(i).chan(ii))];
                    end
                elseif config(i).chan(ii) < 0
                    chanStr = ['DIGITAL-OUT','-0',num2str(abs(config(i).chan(ii)))];
                end

                if enable == 1

                    write(TCP,uint8(['set ',chanStr,'.stimenabled true;']));
                    write(TCP,uint8(['set ',chanStr,'.source KeyPressF',num2str(config(i).fKey),';']));

                    if config(i).chan(ii) > 0

                        write(TCP,uint8(['set ',chanStr,'.FirstPhaseAmplitudeMicroAmps ',num2str(config(i).amp),';']));
                        write(TCP,uint8(['set ',chanStr,'.SecondPhaseAmplitudeMicroAmps ',num2str(config(i).amp),';']));

                        if config(i).polarity(ii) == 1
                            write(TCP,uint8(['set ',chanStr,'.Polarity PositiveFirst;']));
                        elseif config(i).polarity(ii) == -1
                            write(TCP,uint8(['set ',chanStr,'.Polarity NegativeFirst;']));
                        end

                        write(TCP,uint8(['set ',chanStr,'.FirstPhaseDurationMicroseconds ',num2str(round(config(i).duration/2)),';']));
                        write(TCP,uint8(['set ',chanStr,'.SecondPhaseDurationMicroseconds ',num2str(round(config(i).duration/2)),';']));
                    elseif config(i).chan(ii) < 0
                         write(TCP,uint8(['set ',chanStr,'.FirstPhaseDurationMicroseconds ',num2str(config(i).duration),';']));
                    end       
                    
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

function [amplifierData,timeStamps,offset] = byte2doubleConv(waveformArray,numChan,blocksToRead,sampleFreq,downsampleFreq,useConv)
    
    % Expect first 4 bytes to be TCP Magic Number as uint32, if not, determine the offset.
    for i = 1:length(waveformArray)        
        Bytes = waveformArray(1,i:i+3);
        magicNumber = typecast(uint8(Bytes),'uint32');
        if magicNumber == 0x2ef07a08
            break
        end
    end
    offset = i-1;

    if offset == 0    
        
        % Reshape into 1 block per column
        waveformArray = reshape(waveformArray,[length(waveformArray)/blocksToRead,blocksToRead]);

        % Remove magic numbers
        waveformArray(1:4,:) = [];

        % Reshape into 1 frame per column
        waveformArray = reshape(waveformArray,[size(waveformArray,1)/128,blocksToRead*128]);

        % Downsample
        if useConv == 0
            waveformArray = waveformArray(:,1:round(sampleFreq/downsampleFreq):end);
        elseif useConv == 1
            waveformArray = convLPFilt(waveformArray',downsampleFreq/2,sampleFreq,[],sampleFreq/downsampleFreq)';
        end

        % Separate timestamps
        timeStamps = waveformArray(1:4,:);
        waveformArray(1:4,:) = [];
        
        % Pre-allocate memory
        amplifierData = 32768*ones(numChan,size(waveformArray,2));
        
        % Decode each byte pair into a single measurement
        iAmp = 1;
        for i = 1:2:numChan*2
            for ii = 1:size(waveformArray,2)
                amplifierData(iAmp,ii) = typecast(uint8(waveformArray(i:i+1,ii)'),'uint16');
            end
            iAmp = iAmp+1;
        end
        amplifierData = 0.195*(amplifierData-32768);
    else
        amplifierData = [];
        timeStamps = [];
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
    
        for i3 = 1:size(data,2)
    
            P = numel(data(:,i3));
            Q = numel(B);
            filtLength = P+Q-1;
            K = 2^nextpow2(filtLength);
        
            % do the convolution
            afft = fft(data(:,i3),K);
            bfft = fft(B,K);
            c = ifft(afft(:).*bfft(:));
        
            range = [floor(length(B)/2)+1,filtLength-ceil(length(B)/2)+1];
            data(:,i3) = c(range(1):range(2));
        end
    
        if ~isempty(downSampleRatio)
            data = downsample(data,downSampleRatio);
        end
    
        out = data;
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
    activeBox = input{1,2};
    numBox = sum(activeBox);
    numChan = size(data,3);
    WL = input{1,3};
    iWL = input{1,4};
    
    if isempty(datFig) == 1 || ishandle(datFig) == 0

        datFig = figure;
        datAx = cell(1,numBox);
        HdatAx = cell(size(data,1),numBox);
        count = 1;
        for i = 1:numChan
            if activeBox(1,i) == 1
             
                datAx{1,count} = subplot(1,numBox,count,'Parent',datFig);
                hold(datAx{1,count},'on')
                for ii = 1:size(data,1)            
                    HdatAx{ii,count} = plot(datAx{1,count},data(ii,:,i));
                end
                hold(datAx{1,count},'off')
                ylim(datAx{1,count},[0,4])
                ylabel(datAx{1,count},'T/D')
                xticks(0:round(size(data,2)/6):size(data,2))
                xticklabels(0:round(size(data,2)/6):size(data,2))
                xlabel(datAx{1,count},'seconds')
                title(datAx{1,count},WL{iWL(1,count),3})

                count = count+1;
            end
        end
    else
        count = 1;
        for i = 1:numChan
            if activeBox(1,i) == 1
                for ii = 1:size(data,1)
                    set(HdatAx{ii,count},'YData',data(ii,:,i))
                end
                count = count+1;
            end
        end
    end
    drawnow
end

function plotLabel(input)
    global labelFig labelAx HlabelAx %#ok<GVMIS>

    data = input{1,1};
    activeBox = input{1,2};
    numBox = sum(activeBox);
    numChan = size(data,1);
    WL = input{1,3};
    iWL = input{1,4};
    
    if isempty(labelFig) == 1 || ishandle(labelFig) == 0

        labelFig = figure;
        labelAx = cell(1,numBox);
        HlabelAx = cell(numBox,1);
        count = 1;
        for i = 1:numChan
            if activeBox(1,i) == 1
             
                labelAx{1,count} = subplot(1,numBox,count,'Parent',labelFig);         
                HlabelAx{count,1} = plot(labelAx{1,count},data(i,:));
                hold(labelAx{1,count},'off')
                ylim(labelAx{1,count},[0,7])
                yticks(0:7)
                yticklabels(0:7)
                ylabel(labelAx{1,count},'Label')
                xticks(0:round(size(data,2)/6):size(data,2))
                xticklabels(0:round(size(data,2)/6):size(data,2))
                xlabel(labelAx{1,count},'seconds')
                title(labelAx{1,count},WL{iWL(1,count),3})

                count = count+1;
            end
        end
    else
        count = 1;
        for i = 1:numChan
            if activeBox(1,i) == 1
                set(HlabelAx{count,1},'YData',data(i,:))
                count = count+1;
            end
        end
    end
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

function [wt,freqlist,psi_array] = sm_wavelet(x,Fs,freqlist)

    % Inputs:
    % x = signal to be analyzed
    % Fs = sampling frequency
    % freqlist = list of frequencies at which to compute wt (optional) (or set to [] for automatic definition)
    % Outputs:
    % wt = time-frequency image
    % freqlist = useful in case of automatic definition of frequencies
    % psi_array = array of analysis functions (complex values)
    
    x = x(:);
    n = length(x);
    sigma2 = 1;
    xi = 5;
    tolerance = 0.5;
    omega = [(0:n/2) (-ceil(n/2)+1:-1)].*Fs/n;
    omega = omega(:);
    fftx = fft(x);   

    mincenterfreq = 2*tolerance*sqrt(sigma2)*Fs*xi./n;
    maxcenterfreq = Fs*xi/(xi+tolerance/sqrt(sigma2));

    s_array = xi./freqlist;
    minscale = xi./maxcenterfreq;
    maxscale = xi./mincenterfreq;

    freq =  (omega*s_array-xi);
    Psi = bsxfun(@times,realpow(4.*pi.*sigma2,1/4)*sqrt(s_array),exp(-sigma2/2*realpow(freq,2)));

    wt = ifft(bsxfun(@times,Psi,fftx),[],1);
    wt(:,s_array(:)'<minscale | s_array(:)' > maxscale) = 0;
    psi_array = ifft(Psi);
end