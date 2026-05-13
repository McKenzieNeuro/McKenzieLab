%% *** NOTE this version only works on R2024a and up ***

checkMatprefs()
delete(gcp('nocreate'))
clear all %#ok<CLALL>
reset(gpuDevice)

StartTime = [6,30];
recTime = 3600*23.5;
numDays = 7;

% Subject in cage, '' if empty.
boxSubject{1,1} = 'PTP_7.0'; % Cage 1
boxSubject{2,1} = 'PTP_7.2'; % Cage 2
boxSubject{3,1} = 'PTP_8.0'; % Cage 3
boxSubject{4,1} = 'PTP_8.1'; % Cage 4
boxSubject{5,1} = ''; % Cage 5
boxSubject{6,1} = ''; % Cage 6
boxSubject{7,1} = ''; % Cage 7
boxSubject{8,1} = ''; % Cage 8

% IP adress for raspberry pi's, '' if empty.
ipCell{1,1} = '169.254.97.194'; % Cage 1
ipCell{2,1} = '169.254.69.249'; % Cage 2 
ipCell{3,1} = '169.254.101.138'; % Cage 3 
ipCell{4,1} = '169.254.137.129'; % Cage 4 
ipCell{5,1} = '169.254.159.223'; % Cage 5
ipCell{6,1} = '169.254.171.252'; % Cage 6
ipCell{7,1} = ''; % Cage 7
ipCell{8,1} = ''; % Cage 8

usePose = 0;
execution = 'gpu'; % Execution environment: cpu,gpu, multi-gpu.
%netPath = 'D:\Data\ClassifierData\PoseEstimation_EDS+HC\netCNN.mat';
netPath = 'D:\Data\ClassifierData\PoseEstimation_EDS+HC_uint8\netCNN.mat';

useClassify = 0;
classRes = 0.5;
execution = 'gpu'; % Execution environment: cpu,gpu, multi-gpu.
net3dPath = 'D:\Data\ClassifierData\PoseEstimation_EDS+HC_uint8\net3dCNN.mat';

labelDisplay = 1;

numCam = 8;%sum(~cellfun(@isempty,ipCell));
useCam = ~cellfun(@isempty,boxSubject) & ~cellfun(@isempty,ipCell);
useCore = sum(useCam);

PathName = uigetdir(cd,'Select folder for Datastore');
cd(PathName);

vidPaths = cell(8,1);
for i = 1:8
    if useCam(i,1) == 1
        vidPaths{i,1} = [PathName,'\',boxSubject{i,1}];
        if ~isfolder(vidPaths{i,1})
            eval(['mkdir ',vidPaths{i,1}]);
        end
    end
end

global videoOut %#ok<GVMIS>
videoOut = parallel.pool.DataQueue;
afterEach(videoOut,@plotFrame);

textOut = parallel.pool.DataQueue;
afterEach(textOut,@disp);

stopBox = parallel.pool.DataQueue;
afterEach(stopBox,@makeBox);

% CPU core roles
camTimer = 1;
camReader = 2:1+useCore;
frameHandeler = camReader(1,end)+1:camReader(1,end)+useCore;
frameSync = frameHandeler(1,end)+1;
camSaver = frameSync+1:frameSync+2;
numCore = camSaver(1,end);

% Send tags
camCheck = 1;
useSave = 2;
startSig = 3;
qSend = 4;
timerSig = 5;
frameSend = 6;

Targets{1,1} = 'RightEar';
Targets{2,1} = 'LeftEar';
Targets{3,1} = 'TailBase';

Classes{1,1} = 'Active-Light';
Classes{2,1} = 'Active-Dark';
Classes{3,1} = 'Rear-Light';
Classes{4,1} = 'Rear-Dark';
Classes{5,1} = 'Sleep-Light';
Classes{6,1} = 'Sleep-Dark';
Classes{7,1} = 'Seizure';

%%
clc

if usePose == 1

    load(netPath)
    numTarget = length(Targets);
    inputSize = netCNN.Layers(1,1).InputSize(1,1:2);

    % Create object detector to detect CNN parts
    blobAnalyserCNN = vision.BlobAnalysis('BoundingBoxOutputPort',false,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',10);    
end

if useClassify == 1
    load(net3dPath)
    numClass = length(Classes);
    inputSize = net3dCNN.Layers(1,1).InputSize(1,1:2);
    classRes = round(10*classRes);
end

nextDay = 0;
standBy = 2;
whichSaver = 1;
for I = 1:numDays % Day loop
    checkMatprefs()
    % Wait for start time minus two minutes to begin startup checks
    vStartTime = (StartTime(1,1)*3600)+(StartTime(1,2)*60);
    send(textOut,['Will start at ',num2str(StartTime(1,1)),':',num2str(StartTime(1,2))]);
    while 1
        if nextDay == 0
            if getTime >= vStartTime-120
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

    % Get formatted date and time
    DT = getDT;
    
    % Start parallel pool
    if isempty(gcp('nocreate')) == 1
        parpool('local',camSaver(1,end));
    end
    
    skipStart = 0;
    errorCount = 0;
    while 1
        try
            spmd(numCore)            
                mpiSettings('DeadlockDetection','off');

                if spmdIndex == camTimer

                    if skipStart == 0
        
                        % Data queue for stop box
                        qW = parallel.pool.PollableDataQueue;
                        
                        % Check to see if each core connected to it's camera
                        send(textOut,'Connecting to cameras...')
                        camConnect = ~useCam;
                        for i = 1:useCore
                            r = spmdReceive('any',camCheck);
                            camConnect(r(1,1),1) = r(1,2);
                        end
                        
                        finTime = 1;
                        timeCount = [0,0,0];

                        startTime = 1;
                    else
                        startTime = i;
                    end

                    if all(camConnect(useCam,1))                        
                        
                        % Start to frameHandeler cores
                        spmdSend(1,frameHandeler,startSig)
            
                        % Start frameSync core
                        spmdSend(1,frameSync,startSig)
                        
                        % Send signal to the core that will save the frames
                        spmdSend(1,camSaver(1,whichSaver),startSig)

                        if skipStart == 0
                            
                            % Wait for start time
                            while getTime < vStartTime
                                pause(0.001);
                            end
                            
                            % Open video plot
                            send(videoOut,cell(3,numCam)) %#ok<SPGV>
            
                            % Open stop box
                            send(stopBox,qW);
                        end

                        % Start camReader cores
                        spmdSend(1,camReader,startSig)
                        
                        send(textOut,'Starting recording')
                        for i = startTime:recTime % Session loop
            
                            % Send sync signal to camReaders
                            t = tic;
                            spmdSend({t,timeCount},camReader,timerSig);
                            
                            % Keep time in hr:min:sec
                            timeCount(1,3) = timeCount(1,3)+1;
                            if timeCount(1,3) > 59
                                timeCount(1,3) = 0;
                                timeCount(1,2) = timeCount(1,2)+1;
                                if timeCount(1,2) > 59
                                    timeCount(1,2) = 0;
                                    timeCount(1,1) = timeCount(1,1)+1;
                                end
                            end
                            
                            % Check if stop box has been closed
                            send(stopBox,1);
                            boxClose = poll(qW);
                            if isempty(boxClose)    
                                while toc(t) < 0.9999
                                    pause(0.0001)
                                end
                            else
                                finTime = 0;
                                break
                            end
                        end
            
                        % Stop camera core's session loop
                        spmdSend([],camReader,timerSig);
                        send(textOut,'Recording complete')
                        send(textOut,'Finializing videos')
                    
                        % Close stop box
                        if finTime == 1
                            send(stopBox,[]);
                        end
                    else            
                        send(textOut,'Connection to cameras failed')
                        spmdSend([],camReader,startSig) 
                        spmdSend([],frameHandeler,startSig)
                        spmdSend([],frameSync,startSig)
                        spmdSend([],camSaver(1,whichSaver),startSig)
                    end
                    spmdSend([],camSaver(1,standBy),startSig)

                    send(textOut,['camTimer ',num2str(spmdIndex),' Done'])
                elseif any(spmdIndex == camReader)
                    
                    % Find the camera this core uses
                    camID = find(useCam == 1);
                    coreID = find(camReader == spmdIndex);
                    camID = camID(coreID,1);

                    if skipStart == 0
                    
                        % Connect to camera
                        frame = [];
                        for i = 1:3
                            try
                                raspObj = raspi(ipCell{camID,1},'pi','raspberry');
                                camObj = cameraboard(raspObj,'Resolution','640x480');
            
                                for ii = 1:10
                                    frame = snapshot(camObj);
                                end
            
                                break
                            catch
                                try
                                    frame = snapshot(camObj);
                                    break
                                catch
                                    camObj = [];
                                    send(textOut,['Cam ',num2str(camID),' not connected, will try again (',num2str(i),'/3)'])
                                end
                            end
                        end
            
                        frameQ = spmdReceive(frameHandeler(1,coreID),qSend);
                        
                        % Preallocate and inform camTimer of camera connection status
                        if ~isempty(frame)   
                            
                            secCount = 0;
                            lastFrame = 1;
                            frameSec = nan(10,1);
                            frameSize = size(frame);
                            fpsTrack = zeros(recTime,2);
                            frames1Sec = uint8(zeros([frameSize,10]));
                            spmdSend([camID,1],camTimer,camCheck)
                            send(textOut,['Cam ',num2str(camID),' connected'])
                        else
                            spmdSend([camID,0],camTimer,camCheck)
                        end
                    end                    
                    
                    go = spmdReceive(camTimer,startSig);
                    if ~isempty(go)
                        while 1 % Session loop
                            
                            try
                                timeData = spmdReceive(camTimer,timerSig);
                            catch
                                continue
                            end

                            if ~isempty(timeData)
                                
                                fps = 0;
                                t = tic;%timeData{1,1};
                                timeCount = timeData{1,2};
                                while toc(t) < 1 && fps < 10
        
                                    fs = tic;
                                    try
                                        frame = snapshot(camObj);
                                        fps = fps+1;
                                        frames1Sec(:,:,:,fps) = frame;
                                        frameSec(fps,1) = toc(fs);
                                    catch
                                        % Do nothing
                                    end
                                end
                                
                                secCount = secCount+1;                            
                                
                                fs = tic;
                                if fps > 0
                                    % Send one second batch of frames to frameHandeler
                                    send(frameQ,{frames1Sec(:,:,:,1:fps),timeCount})
                                    fpsTrack(secCount,1) = mean(frameSec(1:fps,1),'omitnan');
                                    lastFrame = fps;
                                else
                                    send(frameQ,{frames1Sec(:,:,:,lastFrame),timeCount})
                                end
                                fpsTrack(secCount,2) = toc(fs);
                            else                        
                                break
                            end
                        end
                    end
        
                    % Stop frameHandeler's session loop
                    send(frameQ,[])

                    send(textOut,['camReader ',num2str(spmdIndex),' Done'])
                elseif any(spmdIndex == frameHandeler)
            
                    % Find camera this core uses
                    camID = find(useCam == 1);
                    coreID = find(frameHandeler == spmdIndex);
                    camID = camID(coreID,1);
                    
                    if skipStart == 0
                    
                        % Initialize netCNN
                        if usePose == 1                
                            for i = 1:5
                                predict(netCNN,zeros([inputSize,3]),'ExecutionEnvironment',execution,'Acceleration','auto');
                            end
                        end
            
                        % Initialize net3dCNN
                        if useClassify == 1                
                            for i = 1:5
                                predict(net3dCNN,zeros([inputSize,10,2]),'ExecutionEnvironment',execution,'Acceleration','auto');
                            end
                        end
            
                        frameQ = parallel.pool.PollableDataQueue;
                        spmdSend(frameQ,camReader(1,coreID),qSend)
                    end
                    
                    go = spmdReceive(camTimer,startSig);
                    if ~isempty(go)
        
                        % Preallocate
                        frames1Sec = uint8(zeros(480,640,3,10));
                        frames1SecR = uint8(zeros(256,256,3,10));
        
                        if usePose == 1 || useClassify == 1
                            % Create array to store data to disk
                            poseData = matfile([vidPaths{camID,1},'\',DT,'_poseData.mat'],'Writable',true);
                        end
                        
                        if usePose == 1
                
                            % Create array to store pose coordinates
                            poseData.targets = Targets;
                            poseData.headers = [{'X'},{'Y'},{'Confidence'},{'Distance'}];
                            poseData.XY = zeros([],4,numTarget);
                            poseData.diffPix = zeros([],1);
                            
                            % Preallocate
                            diff1Sec = zeros(10,1);
                            XY = zeros(10,4,numTarget);
                            lastPoint = nan(1,2,numTarget);
                        end
                    
                        if useClassify == 1
                            
                            % Create array to store pose coordinates
                            poseData.classes = Classes;
                            poseData.confidence = zeros([],length(Classes));
            
                            % Preallocate
                            Yclass = 0;
                            behaveCon = single(zeros(10,numClass));
                            data2Sec = uint8(zeros([inputSize,20,2])*127);
                        end
                        
                        secCount = 0;
                        fpsTrack = zeros(recTime,3);
                        while 1 % Session loop
                            if frameQ.QueueLength > 0
                                
                                % Receive 1s load of frames from camReader
                                frameData = poll(frameQ);
        
                                if ~isempty(frameData)
        
                                    secCount = secCount+1;
        
                                    fs = tic;
                
                                    frames1Sec = frameData{1,1};
                                    timeCount = frameData{1,2};
                                    
                                    % Fill in missing frames if fps < 10
                                    fps = size(frames1Sec,4);
                                    frames1Sec = frames1Sec(:,:,:,round(linspace(1,fps,10)));
                                    
                                    % Send to camSaver
                                    spmdSend({frames1Sec,camID},camSaver(1,whichSaver),frameSend)
        
                                    fpsTrack(secCount,1) = toc(fs);
                                    fs = tic;
                                    
                                    % Resize frames for neural networks.
                                    if usePose == 1 || useClassify == 1
                                        frames1SecR = imresize(frames1Sec,inputSize);
                                        frameSecOut = frames1SecR;
                                    else
                                        frameSecOut = frames1Sec;
                                    end
                                    
                                    % Do pose estimation
                                    if usePose == 1
                
                                        YPred = predict(netCNN,frames1SecR,'ExecutionEnvironment',execution);
                                        YPred(YPred < 0.25) = 0;
                                        YPred(YPred > 1) = 1;
                                        YPred = uint8(YPred*255);
                
                                        for i = 1:10
                                            data2Sec(:,:,i+10,3:5) = YPred(:,:,:,i);
                                        end
                
                                        if labelDisplay == 1
                                            frameSecOut = frameSecOut+uint8(sum(YPred/2,3));
                                        end
                                    end
                                    
                                    % Do behavior classification
                                    if useClassify == 1
        
                                        count = 1;
                                        for i = 11:20
                                            data2Sec(:,:,i,1) = rgb2gray(frames1SecR(:,:,:,count));
                                            d = abs(data2Sec(:,:,i,1)-data2Sec(:,:,i-1,1));
                                            d(d < 8) = 0;
                                            data2Sec(:,:,i,2) = d;
                                            count = count+1;
                                        end
        
                                        for i = 1:classRes:10    
                
                                            y = predict(net3dCNN,data2Sec(:,:,i:i+9,1:2),'ExecutionEnvironment',execution,'Acceleration','auto');
                                            [~,Yclass] = max(y);
                                            
                                            idx = i:i+classRes-1;
                                            idx(idx > fps) = [];
                                            for ii = idx
                                                behaveCon(ii,:) = y;
                                                if labelDisplay == 1
                                                    frameSecOut(:,:,:,ii) = insertText(frameSecOut(:,:,:,ii),[inputSize(1,2),0],Classes(Yclass,1),AnchorPoint="RightTop",FontSize=18,BoxColor="blue",BoxOpacity=0.5,TextColor="white");
                                                end
                                            end
                                        end
                                        
                                        data2Sec(:,:,1:10,:) = data2Sec(:,:,11:20,:);
                                        poseData.confidence = [poseData.confidence;behaveCon];
                                    end
        
                                    fpsTrack(secCount,2) = toc(fs);
                                    fs = tic;
                                    
                                    % Send frames to frameSync for output to display
                                    spmdSend({camID,{frameSecOut;fps;timeCount}},frameSync,frameSend);
                                    
                                    % Find XY coordinates from
                                    if usePose == 1
                                        for i = 1:10
                    
                                            diff1Sec(i,1) = mean(data2Sec(:,:,i+10,2),'all');
                                            for ii = 1:numTarget
                    
                                                P = data2Sec(:,:,i+10,2+ii);
                                                mask = logical(imresize(P,inputSize));
                                                [area,centroids] = blobAnalyserCNN(mask);
                                                if isempty(centroids) == 0
                    
                                                    [~,aMax] = max(area);
                                                    xy = round(centroids(aMax,:));
                    
                                                    [counts,centers] = hist(P(P > 0),linspace(0.01,1,20)); %#ok<HIST> 
                                                    [~,idx] = max(counts);
                                                    Con = centers(idx);
                                                else
                                                    xy = [nan,nan];
                                                    Con = 0;
                                                end
                                                XY(i,1:3,ii) = [xy(1,1),xy(1,2),Con];
                    
                                                if isnan(lastPoint(1,1,ii))
                                                    XY(i,4,ii) = 0;
                                                else
                                                    if isnan(XY(i,1,ii))
                                                        XY(i,4,ii) = 0;
                                                    else
                                                        XY(i,4,ii) = pdist([lastPoint(1,1:2,ii);XY(i,1:2,ii)],'euclidean');
                                                    end
                                                end
                                                lastPoint(1,1:2,ii) = XY(i,1:2,ii);
                                            end
                                        end                    
                                        poseData.diffPix = [poseData.diffPix;diff1Sec];
                                        poseData.XY = cat(1,poseData.XY,XY);
                                    end
                                    fpsTrack(secCount,3) = toc(fs);
                                else
                                    break
                                end
                            else
                                pause(0.001)
                            end            
                        end
                    end
        
                    % Stop frameSync's session loop
                    spmdSend([],frameSync,frameSend);
                    % Stop camSaver's session loop
                    spmdSend({[],camID},camSaver(1,whichSaver),frameSend)
                    
                    send(textOut,['frameHandeler ',num2str(spmdIndex),' Done'])
                elseif spmdIndex == frameSync
                    
                    isRunning = useCam;
                    camID = find(useCam == 1)';                    
                   
                    go = spmdReceive(camTimer,startSig);
                    if ~isempty(go)
                        while 1 % Session loop
                            
                            allCam = cell(3,numCam);
                            while sum(cellfun(@isempty,allCam(1,camID))) > 0
                                for i = 1:useCore
                                    if isempty(allCam{1,camID(1,i)})
                                        if spmdProbe(frameHandeler(1,i),frameSend)
                                            
                                            try
                                                frameData = spmdReceive(frameHandeler(1,i),frameSend);
                                            catch
                                                frameData = 0;
                                            end
            
                                            if ~isempty(frameData)
                                                if iscell(frameData)
                                                    allCam(:,camID(1,i)) = frameData{1,2};
                                                end
                                            else
                                                allCam(:,camID(1,i)) = {0};
                                                isRunning(camID(1,i),1) = 0;
                                            end
                                        end
                                    end
                                end
                            end
        
                            % Output data to screen
                            if sum(isRunning) > 0
                                if labelDisplay == 1
                                    send(videoOut,allCam) %#ok<SPGV>
                                end
                            else
                                break
                            end
                        end
                    end
                    send(textOut,['frameSync ',num2str(spmdIndex),' Done'])
                elseif any(spmdIndex == camSaver)    
                   
                    while 1
                        if spmdProbe(camTimer,startSig)
                            
                            try
                                go = spmdReceive(camTimer,startSig);
                            catch
                                go = 1;
                            end

                            if ~isempty(go)
                
                                % Create all video writer objects
                                if skipStart == 0
                                    writerObj = cell(8,1);
                                    recFinished = ~useCam;
                                    VidFileName = [DT,'.mp4'];
                                    for i = 1:8
                                        if useCam(i,1) == 1
                                            writerObj{i,1} = VideoWriter([vidPaths{i,1},'\',VidFileName],'MPEG-4');
                                            writerObj{i,1}.FrameRate = 10;
                                            open(writerObj{i,1});
                                        else
                                             writerObj{i,1} = [];
                                        end
                                    end
                                end
                                
                                while ~all(recFinished) % Session loop                        
                                    try
                                        camData = spmdReceive('any',frameSend);        
                                        frames1Sec = camData{1,1};
                                        camID = camData{1,2};
                                        
                                        % Write frames to file
                                        if ~isempty(frames1Sec)
                                            for i = 1:10                    
                                                writeVideo(writerObj{camID,1},frames1Sec(:,:,:,i));
                                            end
                                        else
                                            if recFinished(camID,1) == 0
                                                close(writerObj{camID,1});
                                                recFinished(camID,1) = 1;
                                                send(textOut,'Cam ',num2str(camID), 'complete.');
                                            end
                                        end
                                    catch
                                        send(textOut,'camSaver received packet of unknown type')
                                    end                        
                                end
                                send(videoOut,[]) %#ok<SPGV>
                            end
                            break
                        else
                            pause(0.01)
                        end
                    end

                    send(textOut,['camSaver ',num2str(spmdIndex),' Done'])
                end
            end
            break
        catch ME

            fprintf('MPI error on worker %d: %s\n', spmdIndex, ME.message);
            disp(ME.getReport());
            skipStart = 1;

            errorCount = errorCount+1;
            if errorCount == 3
                delete(gcp('nocreate'))
                break
            end
        end
    end

    % Switch camSaver cores
    if whichSaver == 1
        whichSaver = 2;
        standBy = 1;
    elseif whichSaver == 2
        whichSaver = 1;
        standBy = 2;
    end

    l = i{1};
    i{1}
    allTrack_Mean1 = zeros(l,2,sum(useCam));
    allTrack_Mean2 = zeros(l,3,sum(useCam));
    count = [0;0];
    for i = 1:size(fpsTrack,2)
        if any(i == camReader)
            temp = fpsTrack{i};
            count(1,1) = count(1,1)+1;
            allTrack_Mean1(:,:,count(1,1)) = temp(1:l,:);
        elseif any(i == frameHandeler)
            temp = fpsTrack{i};
            count(2,1) = count(2,1)+1;
            allTrack_Mean2(:,:,count(2,1)) =  temp(1:l,:);
        end
    end

    allTrack_Mean1 = mean(allTrack_Mean1,3);
    allTrack_Mean2 = mean(allTrack_Mean2,3);
    
    ledg = cell(5,1);
    ledg{1,1} = 'Get Frame';
    ledg{2,1} = 'Send Frame';
    ledg{3,1} = 'Send Frame 2';
    ledg{4,1} = 'ANN';
    ledg{5,1} = 'Pose';
    
    core = 3;
    clf
    hold on
    for i = 1:2
        plot(allTrack_Mean1(:,i))
    end
    for i = 1:3
        plot(allTrack_Mean2(:,i))
    end
    legend(ledg)
end
send(textOut,'Done')

%%

function DT = getDT
    now = clock; %#ok<CLOCK>
    DT = [num2str(now(1,2)),'-',num2str(now(1,3)),'-',num2str(now(1,1)),'(',num2str(now(1,4)),'.',num2str(now(1,5)),')'];
end

function plotFrame(input)
    global videoOut vidFig vidAx HvidAx %#ok<GVMIS>

    if ~isempty(input)
        numCam = size(input,2);
    
        if isempty(vidFig) == 1 || ishandle(vidFig) == 0
            vidFig = figure;
            vidFig.Name = 'Video';
            vidFig.Position = [1921,267,1920,730];
            vidAx = cell(1,numCam);
            HvidAx = cell(1,numCam);
        end
        
        for i = 1:10

            tic
            for ii = 1:8
                if ~isempty(input{1,ii})
                    
                    frame = insertText(input{1,ii}(:,:,:,i),[0,0],['FPS: ',num2str(input{2,ii})],AnchorPoint="LeftTop",FontSize=18,BoxColor="blue",BoxOpacity=0,TextColor="black");
                    if isempty(HvidAx{1,ii})
                        vidAx{1,ii} = subplot(2,4,ii,'Parent',vidFig);
                        HvidAx{1,ii} = image(frame,'Parent',vidAx{1,ii});
                        title(vidAx{1,ii},['Cam ',num2str(ii),' ',num2str(input{3,ii}(1,1)),':',num2str(input{3,ii}(1,2)),':',num2str(input{3,ii}(1,3))])
                    else
                        set(HvidAx{1,ii},'CData',frame)
                        title(vidAx{1,ii},['Cam ',num2str(ii),' ',num2str(input{3,ii}(1,1)),':',num2str(input{3,ii}(1,2)),':',num2str(input{3,ii}(1,3))])
                    end
                end
            end
            drawnow

            if videoOut.QueueLength <= 1
                pause(0.095-toc)
            else
                pause(0.025)
            end            
        end
    else
        close(vidFig)
    end
end

function makeBox(input)
    global closeBox qSend %#ok<GVMIS>

    if ~isempty(input)
        if input ~= 1
            
            qSend = input;
            closeBox = errordlg('Press ok to stop.','Stop','non-modal');
            closeBox.Position = [2719,728,161,60];
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