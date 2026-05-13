clear all %#ok<CLALL>

boxSubject{1,1} = 'EDS 3.0'; % Subject in box 1
boxSubject{2,1} = 'EDS 4.0'; % Subject in box 2
boxSubject{3,1} = 'EDS 4.1'; % Subject in box 3
boxSubject{4,1} = 'EDS 4.2'; % Subject in box 4

boxRHS(1,1) = 1; % Box for RHS port A
boxRHS(2,1) = 2; % Box for RHS port B
boxRHS(3,1) = 3; % Box for RHS port C
boxRHS(4,1) = 4; % Box for RHS port D

stimI = 50:50:1000; % Stim amp in uA.
stimP = 100; % Stim train period in ms.
stimNum = 100; % Stim train pulse number.
BL = 600; % Baseline interval in seconds.

% Channel, [nan,nan] = unused, [c1,nan] = monopolar, [c1,c2] = bipolar.
stimChan(1,:) = [5,6]; % Port A
stimChan(2,:) = [5,6]; % Port B
stimChan(3,:) = [5,6]; % Port C
stimChan(4,:) = [5,6]; % Port D

% Current direction for bipolar stim. P = positive first, N = negative
% first
currentFlow = {'N','N'};

useStim = logical([1,1,1,1]);

intanIP_RHS = '127.0.0.1'; % IP for RHS.
intanPort1_RHS = 5000; % Port for RHS commands.
ArdPin{1,1} = 'D13'; % Arduino pin for LED.

% Intan channel ID
chanID(1,1:2) = [{'000'},{'PrL(L)'}];
chanID(2,1:2) = [{'001'},{'PrL(R)'}];
chanID(3,1:2) = [{'002'},{'AVT(L)'}];
chanID(4,1:2) = [{'003'},{'BLA(R)'}];
chanID(5,1:2) = [{'004'},{'CA1(L)'}];
chanID(6,1:2) = [{'005'},{'CA1(R)'}];
chanID(7,1:2) = [{'006'},{'LDT(L1)'}];
chanID(8,1:2) = [{'007'},{'LDT(L2)'}];

% Video Settings:
fps = 25;
RecVid = 1; % Record from webcams, 1 = yes, 0 = no
showVideo = 1; % Show video feed from cameras, 1 = yes, 0 = no

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

PathName = uigetdir(cd,'Select folder for Datastore');
cd(PathName);
save([PathName,'\config.mat'])

%%
clc
port = {'a','b','c','d'};

hardwareTimer = 1;
camReader = 2;
camSaver = 3;

textOut = parallel.pool.DataQueue;
afterEach(textOut,@disp);

notification = parallel.pool.DataQueue;
afterEach(notification,@parBeep);

global videoOut %#ok<GVMIS>
videoOut = parallel.pool.DataQueue;
afterEach(videoOut,@plotFrame);

if isempty(gcp('nocreate')) == 1
    parpool('local',3);
end 

for Iter = 1:length(stimI)

    seshPathName = [PathName,'\',num2str(stimI(1,Iter)),'uA'];
    eval(['mkdir ',seshPathName]);

    spmd (3)
        if spmdIndex == hardwareTimer
            
            if Iter == 1
                send(textOut,'Connecting to hardware...')
                intanSend = tcpclient(intanIP_RHS,intanPort1_RHS);
                Ard = arduino('com3','uno');
            end

            write(intanSend,uint8(['set Filename.Path ',seshPathName]));
            write(intanSend,uint8('set FileFormat OneFilePerSignalType'));
            write(intanSend,uint8('set Filename.BaseFilename RHS'));
            write(intanSend,uint8('execute RescanPorts;'));

            for i = 1:size(stimChan,1)
                for ii = 1:size(stimChan(i,:),2)
                    if ~isnan(stimChan(i,ii)) && useStim(1,i) == 1
            
                        write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.stimenabled true;']));
                        write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.source KeyPressF1;']));
                        write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.FirstPhaseAmplitudeMicroAmps ',num2str(stimI(1,Iter)),';']));
                        write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.SecondPhaseAmplitudeMicroAmps ',num2str(stimI(1,Iter)),';']));
                        write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.PulseOrTrain PulseTrain;']));
                        write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.NumberOfStimPulses ',num2str(stimNum),';']));
                        write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.PulseTrainPeriodMicroseconds ',num2str(stimP*1000),';']));
                        
                        if strcmp(currentFlow{1,ii},'N') == 1
                            write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.Polarity NegativeFirst;']));
                        elseif strcmp(currentFlow{1,ii},'P') == 1
                            write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.Polarity PositiveFirst;']));
                        end                        
                    
                        send(textOut,['Port',num2str(i),' amp set to ',num2str(stimI(1,Iter)),'uA.'])
                        
                    elseif ~isnan(stimChan(i,ii)) && useStim(1,i) == 0
                        write(intanSend,uint8(['set ',port{1,i},'-',chanID{stimChan(i,ii),1},'.stimenabled false;']));
                    end
                end
            end
            write(intanSend,uint8('execute uploadstimparameters'));
            pause(15)
            
            spmdBarrier   

            send(textOut,['Starting ',num2str(stimI(1,Iter)),'uA test.'])
            write(intanSend,uint8('set runMode record;'));
            pause(BL)
            send(notification,1);
            writeDigitalPin(Ard,ArdPin{1,1},1);
            write(intanSend,uint8('execute ManualStimTriggerPulse F1'));
            pause((stimP/1000)*stimNum)
            writeDigitalPin(Ard,ArdPin{1,1},0);
            pause(BL)
            write(intanSend,uint8('set runMode stop;'));

            spmdSend([],camReader);

        elseif spmdIndex == camReader

            % Create all webcam objects           
            cam = cell(numCam,1); 
            for ii = 1:numCam                
                cam{ii,1} = webcam(camChoice(1,ii));
            end           
            frames = uint8(zeros(Dim(1,1),Dim(1,2),3,numCam,40)); 
            errorCode = uint8([0,255,0;255,0,255;0,255,0]);
            
            spmdBarrier

            ck = zeros(1,numCam);
            while 1

                Count = 1;
                tic
                while toc < 0.999
                    for ii = 1:numCam
                        try
                            frames(:,:,:,ii,Count) = snapshot(cam{ii,1});
                            ck(1,ii) = 0;
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
                            frames(1:3,1:3,1,ii,Count) = errorCode;                                
                        end
                    end
                    if showVideo == 1
                        if videoOut.QueueLength < 2 %#ok<SPGV>
                            send(videoOut,{frames(:,:,:,:,Count),fps,WLc,iWL}); %#ok<SPGV>
                        end
                    end
                    Count = Count+1;
                end
                spmdSend({frames(:,:,:,:,1:Count),Count},camSaver);

                if spmdProbe(hardwareTimer) == 1
                    spmdReceive(hardwareTimer);
                    break
                end
            end 
            spmdSend([],camSaver);            

        elseif spmdIndex == camSaver

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
            
            spmdBarrier            

            ck = zeros(1,numCam);
            while 1
    
                data = spmdReceive(camReader);
                
                if isempty(data)
                    break
                end 

                Count = data{1,2}-1;
                frames = data{1,1};

                idx = round(linspace(1,Count,fps));
                frames = frames(:,:,:,:,idx);

                for ii = 1:numCam
                    for iii = 1:fps
                        if ~all(frames(1:3,1:3,1,ii,iii) == errorCode,'all')
                            if mean(sum(frames(:,:,:,ii,iii),3),'all') > 250
                                writeVideo(writerObj{ii,1},frames(:,:,:,ii,iii));
                                ck = zeros(1,numCam);
                            else
                                if all(ck < 5)
                                    ck(1,ii) = ck(1,ii)+1;
                                else
                                    spmdSend(0,intanWriter);
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
    
    if Iter < length(stimI)
        promptMessage = sprintf(['Would you like to continue to ',num2str(stimI(1,Iter+1)),'uA?']);
        titleBarCaption = 'settings';
        nextAmp = questdlg(promptMessage,titleBarCaption,'Yes','No','Yes');
    
        if strcmp(nextAmp,'Yes')
            
            idx = find(useStim == 1);
            [boxChoice,~] = listdlg('PromptString',[{'Did any subjects'};{'have a seizure?'}],'ListString',boxSubject(boxRHS(useStim > 0,1),1),'SelectionMode','multiple');
            useStim(1,idx(1,boxChoice)) = false;

            if any(useStim == 1)
                continue
            else
                break
            end
        elseif strcmp(nextAmp,'No')
            break
        end
    end
end
disp('Complete')

%%

function parBeep(in)
    if in == 1
        beep
    end
end

function plotFrame(input)   
    global videoOut vidFig vidAx HvidAx %#ok<GVMIS>

    frame = input{1,1};
    numCam = size(frame,4);
    WL = input{1,3};
    iWL = input{1,4};
    
    if videoOut.QueueLength < 3
        
        tic
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
    
        while toc < 0.05
            pause(0.001)
        end
    end
end