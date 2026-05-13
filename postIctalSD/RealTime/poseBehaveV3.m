
% Tracking targets for pose estimation
Targets{1,1} = 'RightEar';
Targets{2,1} = 'LeftEar';
Targets{3,1} = 'TailBase';
SymPair = [1,2]; % Targets that are symetrical pairs. Format as [p1,p2;p3,p4;...]
 
% Pose
execution = 'gpu'; % Execution environment: cpu,gpu, multi-gpu.
labelVidPose = 1; % Create labeled video, will still save coordinate data if off
xSmooth = 1; % Use smoothing to reduce target jitter. 0 = off, 1 = on.
numWorkers = 3; % Number of CPU cores to use for labeling.
pixCm = 13; % Number of pixels per cm
overWrite = 0; % Over write previously processed data

% Image input size for CNN, larger gives better results but uses more memory. 
% Must be divisible by 16.
inputSize = [256,256];

% Training Data
numTrainingFrames = 60; % Number of frames to pull from selected video.
ClusterPool = 5000; % Number of frames from video to use.
numClust = 10; % Number of clusters to sort frames.
BaseSize = [480,640]; % Display size for frames in UI.

% Net training
numIter = 50000; % Number of training iterations.
miniBatchSizeCNN = 10; % Number of images in training batches, more gives better results but uses more memory. Best as a factor of 100.
miniBatchSize3dCNN = 5;
actLayer = 'reLuRegBN_6.1';
actLayer_numFilt = 512;

% Behavioral classes
Classes{1,1} = 'Active-Light';
Classes{2,1} = 'Active-Dark';
Classes{3,1} = 'Rear-Light';
Classes{4,1} = 'Rear-Dark';
Classes{5,1} = 'Sleep-Light';
Classes{6,1} = 'Sleep-Dark';
Classes{7,1} = 'Seizure';

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

isLoad = 0;
%%

promptMessage = sprintf('Would you like to make a new model, or load an existing?');
titleBarCaption = 'settings';
loadAction = questdlg(promptMessage, titleBarCaption, 'New','Load','New');

if exist('ProjectPathName','var') == 1 && ischar(ProjectPathName) == 1
    ProjectPathName = uigetdir(ProjectPathName,'Select folder for Datastore');
else
    ProjectPathName = uigetdir(cd,'Select folder for Datastore');
end
cd(ProjectPathName)

if strcmp(loadAction,'New') == 1

    % make data dir folders
    eval(['mkdir ',ProjectPathName,'\Training']);
    eval(['mkdir ',ProjectPathName,'\Training\Data']);
    eval(['mkdir ',ProjectPathName,'\Training\LED']);

    numTarget = length(Targets);
    numClass = length(Classes);
    totalFrames = 0;
    save([ProjectPathName,'\ProjectData.mat'],'inputSize','Targets','totalFrames');

    % make data dir folders
    eval(['mkdir ',ProjectPathName,'\Training\Sequences']);
    for i = 1:numClass
        eval(['mkdir ',ProjectPathName,'\Training\Sequences\',num2str(i),'_',Classes{i,1}]);
    end

    sequenceData = Classes;
    sequenceData(:,2) = {0};
    save([ProjectPathName,'\Training\SequenceData.mat'],'sequenceData');

    layersLED = makeLED;
    save([ProjectPathName,'\layersLED.mat'],'layersLED');

    layersCNN = makeCNN(numTarget,inputSize);
    save([ProjectPathName,'\NetworkLayers.mat'],'layersCNN');
    
    promptMessage = sprintf('Would you like to label some training data?');
    titleBarCaption = 'settings';
    labelAction_CNN = questdlg(promptMessage,titleBarCaption,'Yes','Later','Yes');
    labelAction_LSTM = 'NA';
    updateAction = 'NA';

elseif strcmp(loadAction,'Load') == 1

    if isLoad == 0

        findDir = dir([ProjectPathName,'\Training']);
        findDir = {findDir.name};
        hasVidDir = 0;
        for i = 1:length(findDir)
            if strcmp(findDir{1,i},'Sequences') == 1
                hasVidDir = 1;
                break
            end
        end
    
        load('netLED.mat')
        layersLED = layerGraph(netLED);
    
        load('netCNN.mat')
        layersCNN = layerGraph(netCNN);
        isLoad = 1;
    end
    
    promptMessage = sprintf('What would you like to do');
    titleBarCaption = 'settings';
    loadAction = questdlg(promptMessage,titleBarCaption,'Analyse video','Update networks','Analyse video');

    if strcmp(loadAction,'Update networks') == 1
        
        promptMessage = sprintf('What would you like to do');
        titleBarCaption = 'settings';
        if hasVidDir == 0
            loadAction = questdlg(promptMessage,titleBarCaption,'Update LED network','Update pose network','Create behavior network','Update LED network');
        else
            loadAction = questdlg(promptMessage,titleBarCaption,'Update LED network','Update pose network','Update behavior network','Update LED network');
        end

        promptMessage = sprintf('What would you like to do');
        titleBarCaption = 'settings';
        updateAction = questdlg(promptMessage, titleBarCaption,'Add training data','Refine network','Retrain network','Add training data');
    else
        updateAction = 'N/A';
    end
    
    load('ProjectData.mat')
end

if strcmp(loadAction,'Analyse video') == 1

    promptMessage = sprintf('Select manually or batch?');
    titleBarCaption = 'settings';
    selFile = questdlg(promptMessage,titleBarCaption, 'Manual','Batch','Manual');

    if strcmp(selFile,'Manual') == 1
        
        count = 0;
        allFiles = cell(0,5);
        while 1
    
            [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select video(s)','MultiSelect','on');
            if length(VidFilePath) > 1

                cd(VidFilePath)
                [DigDFileName,DigFilePath] = uigetfile('*.dat','Sync to Intan?','MultiSelect','off');
                if length(DigFilePath) > 1
                    prompt = {'Enter sample rate:'};
                    definput = {'20000'};
                    Sfreq = str2double(cell2mat(inputdlg(prompt,'Input',[1,25],definput)));
                else
                    Sfreq = NaN;
                end

                if ~iscell(VidFileName)
                    count = count+1;
                    if isnan(Sfreq)
                        allFiles(count,1:2) = {VidFilePath,VidFileName};
                    else
                        allFiles(count,:) = {VidFilePath,VidFileName,DigFilePath(length(VidFilePath)+1:end),DigDFileName,Sfreq};
                    end
                else
                    for i = 1:length(VidFileName)
                        count = count+1;
                        if isnan(Sfreq)
                            allFiles(count,1:2) = {VidFilePath,VidFileName};
                        else
                            allFiles(count,:) = {VidFilePath,VidFileName{1,i},DigFilePath(length(VidFilePath)+1:end),DigDFileName,Sfreq};
                        end
                    end
                end               

                promptMessage = sprintf('Add more files?');
                addFile = questdlg(promptMessage,'Pick files','Yes','No','Yes');
    
                if strcmp(addFile,'No') == 1
                    break
                end
            else
                break
            end
        end
    else

        [file,path] = uigetfile('*.xlsx');
        t = readtable([path,file],'VariableNamingRule','preserve');
        t = table2cell(t);

        count = 0;
        allFiles = cell(0,5);
        for i = 1:size(t,1)
            if ~isempty(t{i,3}) && ~isnan(t{i,4})
                count = count+1;
                allFiles{count,1} = [t{i,2},'\'];
                allFiles{count,2} = t{i,3};
                allFiles{count,3} = [t{i,5},'\'];
                allFiles{count,4} = t{i,6};
                allFiles{count,5} = t{i,7};
            end
        end
    end

    if ~isempty(allFiles)
        %%
        clc
        numFile = size(allFiles,1);
        for i = 1:numFile

            disp(['Working on file ',num2str(i),'/',num2str(numFile)])

            % Create and pre-allocate datastore
            vidFile = [allFiles{i,1},allFiles{i,2}];
            poseData = matfile([vidFile(1,1:end-4),'_poseData.mat'],'Writable',true);
            
            Sfreq = allFiles{i,5};
            if ~isprop(poseData,'digData1KHz') || overWrite == 1
                if ~isempty(allFiles{i,5})
                    digInFile = [allFiles{i,1},allFiles{i,3},allFiles{i,4}];
                    [digData,tsN] = scanDigitalin(digInFile,Sfreq);
                else
                    digData = 'N/A';
                    tsN = 'N/A';
                end
                poseData.digData1KHz = digData;
                poseData.tsPulse_neural = tsN;
            else
                digData = poseData.digData1KHz;
                tsN = poseData.tsPulse_neural;
            end
            
            if  ~isprop(poseData,'frameRate') || ~isprop(poseData,'diffPix') || overWrite == 1
                [fps,diffPix,conLED,XY] = scanVideo(vidFile,netLED,{netCNN,length(Targets)},[]);
                if ~isempty(fps)
                    poseData.frameRate = fps;
                    poseData.diffPix = diffPix;
                    poseData.conLED = conLED;
                    poseData.XY = XY;
                else
                    continue
                end
            else
                fps = poseData.frameRate;
                diffPix = poseData.diffPix;
                conLED = poseData.conLED;
                XY = poseData.XY;
            end
            numFrames = length(poseData.diffPix);
            
            if ~ischar(conLED) && sum(conLED) > 0
                conLED(conLED < mode(conLED(conLED > 0.8))*0.95) = 0;
                conLED(conLED > 0) = 1;
                tsV = findVideoPulse(conLED,fps);
                poseData.tsPulse_video = tsV;
                [tsNsync,tsVsync] = syncPulse(tsN,tsV,digData,conLED,fps,numFrames);
            else
                 poseData.tsPulse_video = 'N/A';
                 tsNsync = tsN;
                 tsVsync = ((tsN*Sfreq)/(Sfreq/fps))/fps;
            end

            poseData.tsPulse_neuralSync = tsNsync;
            poseData.tsPulse_videoSync = tsVsync;

            % Interpolate video timestamps to neural   
            try
                neuralTS = interp1(tsVsync(:,1),tsNsync(:,1),(0:numFrames-1)/fps,'linear','extrap')';
            catch
                neuralTS = 'N/A';
                disp('Could not interp')
            end            
            poseData.neuralTS = neuralTS;
            
            if ~isnan(pixCm) && ~isempty(XY)

                if ~ischar(neuralTS)
                    ts = neuralTS;
                else
                    ts = (0:numFrames-1)*(1/fps);
                end
                velocity = XY2Vel(XY,ts,pixCm);
            else
                velocity = 'N/A';
            end
            poseData.velocity = velocity;
        end
        disp('Complete!')
        %%
    end
elseif strcmp(loadAction,'Update LED network') == 1
    if strcmp(updateAction,'Add training data') == 1
        [frameDataN,frameDataY] = pickLEDFrames(ProjectPathName,[],0);
    elseif strcmp(updateAction,'Refine network') == 1
        [frameDataN,frameDataY] = pickLEDFrames(ProjectPathName,netLED,1);
        promptMessage = sprintf('Train now?');
        trainNow = questdlg(promptMessage,'Training','Yes','No','Yes');    
        if strcmp(trainNow,'Yes') == 1
            netLED = trainLED(ProjectPathName,frameDataN,frameDataY,netLED,0.001);
        end
    elseif strcmp(updateAction,'Retrain network') == 1
        layersLED = makeLED;
        [frameDataN,frameDataY] = loadFrameData(ProjectPathName);
        netLED = trainLED(ProjectPathName,frameDataN,frameDataY,layersLED,0.001);
    end
elseif strcmp(loadAction,'Update pose network') == 1

    if strcmp(updateAction,'Add training data') == 1
        frameData = pullFrames(numTrainingFrames,ClusterPool,numClust);
    elseif strcmp(updateAction,'Refine network') == 1
        frameData = pickFrames(netCNN,CB);
    end

    xyLabels = labelFrames(frameData,Targets,CB);
    [trainingData,validationData] = makeTrain(frameData,xyLabels,inputSize,SymPair);
    
    if strcmp(updateAction,'Refine network') == 1

        learnRate = 0.00064;
        for i = 1:5
            [layersCNN,netCNN] = trainCNN(ProjectPathName,Targets,{trainingData,validationData},layersCNN,miniBatchSizeCNN,learnRate);
            learnRate = learnRate/2;
        end
    end

    saveTrain(ProjectPathName,{trainingData,validationData});

elseif strcmp(loadAction,'Create behavior network') == 1

    numFrame = 10;
    layers3dCNN = make3dCNN(numClass,[inputSize,numFrame,2]);

    analyzeNetwork(layers3dCNN)
elseif strcmp(loadAction,'Update behavior network') == 1
    sequenceData = pickVideo(ProjectPathName,[256,256,10,2],net3dCNN,CB);
    [layers3dCNN,net3dCNN] = train3dCNN(ProjectPathName,Classes,layers3dCNN,20,0.0001);
end

function layersLED = makeLED

    layersLED = [
    imageInputLayer([100,100,3],'Normalization','none')
    
    convolution2dLayer(7,8,'Padding',[1,1,1,1])
    batchNormalizationLayer
    reluLayer   

    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer 

    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer

    maxPooling2dLayer(2,'Stride',2)

    convolution2dLayer(3,64,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];
    %analyzeNetwork(layersLED)
end

function [frameDataN,frameDataY] = pickLEDFrames(ProjectPathName,netLED,refine)

    try
        load([ProjectPathName,'\','Training\LED\trainingData.mat']) %#ok<LOAD>
    catch
        frameDataN = uint8(zeros(100,100,3,0));
        frameDataY = uint8(zeros(100,100,3,0));
        Count = [0,0];
    end    
    
    while 1
    
        [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select training video','MultiSelect','off');
        if ischar(VidFilePath)

            cd(VidFilePath)
            FullName = fullfile(VidFilePath,VidFileName);
    
            if refine == 1
                try
                    poseData = load([FullName(1,1:end-4),'_poseData.mat']);
                catch
                    poseData = [];
                end
            else
                poseData = [];
            end
            
            % Create video reader
            disp('Reading frames...')
            vr = VideoReader(FullName); %#ok<TNMLP>
            fps = round(vr.FrameRate);
            numFrames = floor(fps*vr.Duration);
            
            figure
            definput = {'0','0','1','0','10'};
            while 1
    
                if ~isempty(poseData)
                    promptMessage = sprintf('Type in time or select frame from confindence plot?');
                    titleBarCaption = 'settings';
                    selAction = questdlg(promptMessage,titleBarCaption,'Type','Select','Type');
                else
                    selAction = [];
                end
    
                if isempty(selAction) || strcmp(selAction,'Type')
            
                    dlgtitle = '';
                    prompt = {'Hours:','Minutes:','Seconds:','Or frame number:','Sample length(s):'};
                    dims = [1,10];
                    t = inputdlg(prompt,dlgtitle,dims,definput);
                    if isempty(t) == 1
                        break
                    end
                
                    tNum(1,1) = str2double(t{1,1}); 
                    tNum(2,1) = str2double(t{2,1}); 
                    tNum(3,1) = str2double(t{3,1}); 
                    tNum(4,1) = str2double(t{4,1});
                    tNum(5,1) = str2double(t{5,1});
                else
                    
                    f = figure;
                    ax = axes(f); %#ok<LAXES>
                    plot(ax,poseData.conLED);
                    axis tight
                    typeIn = input('enter frame number here: ');
                    close(f)
    
                    tNum(1,1) = str2double(definput{1,1}); 
                    tNum(2,1) = str2double(definput{1,2}); 
                    tNum(3,1) = str2double(definput{1,3}); 
                    tNum(4,1) = typeIn;
                    tNum(5,1) = str2double(definput{1,5});
                end
                
                if tNum(4,1) == 0
                    s = floor(((tNum(1,1)*3600)+(tNum(2,1)*60)+tNum(3,1))*fps);
                else
        
                    s = tNum(4,1);
                    tNum(3,1) = round(s/fps);
                    for i = 0:1
                        d = floor(tNum(3-i,1)/60);
                        r = (tNum(3-i,1)/60)-d;
                        tNum(3-i,1) = round(r*60);
                        tNum(3-i-1,1) = d;
                    end
                    tNum(4,1) = 0;
                end
        
                L = tNum(5,1);
                if s+(L*fps) > numFrames
                    wBox = warndlg('Not enough frames in file.','Warning');
                    waitfor(wBox)
                    continue
                end
               
                tNum(3,1) = tNum(3,1)+L;
                if tNum(3,1) >= 60
        
                    tNum(2,1) = tNum(2,1)+1;
                    tNum(3,1) = tNum(3,1)-60;
            
                    if tNum(2,1) >= 60
            
                        tNum(1,1) = tNum(1,1)+1;
                        tNum(2,1) = tNum(2,1)-60;
                    end
                end
            
                for i = 1:5
                    definput{1,i} = num2str(tNum(i,1));
                end
                
                % Read frames, resize
                frame = read(vr,[s,s+(fps*L)-1]);
                numFrame = size(frame,4);
                frameR = imresize(frame,[100,100]);
        
                try                
                    con = predict(netLED,frameR);
                catch
                    con = [];
                end
        
                RGB = zeros(4,(fps*L));
                CB = [1,0,0;0,1,0;0,0,1];
                for i = 1:numFrame
                    
                    clf
                    subplot(2,1,1)
                    title('RGB Values')
                    hold on
                    for ii = 1:3
                        temp = double(frameR(:,:,ii,i))-200;
                        temp(temp < 1) = nan;
                        RGB(ii,i) = mean(temp,'all','omitnan');
                        plot(RGB(ii,:),'Color',CB(ii,:))
                    end
        
                    if ~isempty(con)
                        RGB(4,i) = con(i,2);
                        plot(RGB(4,:)*50,'Color','k')
                    end        
        
                    hold off
                    xlim([1,numFrame])
                    ylim([0,50])
            
                    subplot(2,1,2)
                    f = frame(:,:,:,i);
                    imshow(f)
                    drawnow
                end
            
                while 1
                
                    [sel(1,1),sel(1,2),~] = ginput(1);
                    x = round(sel(1,1));
                
                    if x > numFrame
                        break
                    elseif x < 1
                        x = 1;
                    end           
                    
                    subplot(2,1,2)
                    f = frame(:,:,:,x);
                    imshow(f)
                
                    promptMessage = sprintf('LED?');
                    titleBarCaption = 'settings';
                    answer = questdlg(promptMessage,titleBarCaption,'No','Yes','Skip','Skip');
                    
                    if ~isempty(answer)
                        if strcmp(answer,'No') == 1
                            Count(1,1) = Count(1,1)+1;
                            frameDataN(:,:,:,Count(1,1)) = frameR(:,:,:,x);
                        elseif strcmp(answer,'Yes') == 1
                            Count(1,2) = Count(1,2)+1;
                            frameDataY(:,:,:,Count(1,2)) = frameR(:,:,:,x);
                        end
                        clc
                        disp(['Off: ',num2str(Count(1,1))])
                        disp(['On:  ',num2str(Count(1,2))])
                    end
                end
            
                promptMessage = sprintf('Select another time?');
                titleBarCaption = 'settings';
                answer = questdlg(promptMessage,titleBarCaption,'Yes','No','Yes');
                
                if strcmp(answer,'No') == 1
                    break
                end
            end
            close(gcf)
            clear vr
    
            save([ProjectPathName,'\','Training\LED\trainingData.mat'],'frameDataN','frameDataY','Count')
        
            promptMessage = sprintf('Select another video?');
            titleBarCaption = 'settings';
            answer = questdlg(promptMessage,titleBarCaption,'Yes','No','Yes');
            
            if strcmp(answer,'No') == 1
                break
            end
        else
            break
        end
    end
end

function [frameDataN,frameDataY] = loadFrameData(ProjectPathName)
    d = load([ProjectPathName,'\','Training\LED\trainingData.mat']);
    frameDataN = d.frameDataN;
    frameDataY = d.frameDataY;
end

function netLED = trainLED(ProjectPathName,frameDataN,frameDataY,netLED,learnRate)

    X = cat(4,frameDataN,frameDataY);
    Y = [false(size(frameDataN,4),1);true(size(frameDataY,4),1)];
    numObs = length(Y);
    
    X_Aug = uint8(zeros(100,100,3,numObs*8));
    Y_Aug = false(numObs,1);
    
    count = 0;
    for i = 1:numObs
    
        frame = X(:,:,:,i);
        for ii = 1:4
    
            count = count+1;
            X_Aug(:,:,:,count) = frame;
            Y_Aug(count,1) = Y(i,1);
    
            count = count+1;
            X_Aug(:,:,:,count) = flip(frame,1);
            Y_Aug(count,1) = Y(i,1);
    
            frame = imrotate(frame,90);
        end
    end
    
    idx = randperm(numObs*8);
    X_Aug = X_Aug(:,:,:,idx);
    Y_Aug = Y_Aug(idx,1);
    Y_Aug = categorical(Y_Aug);
    
    mB = 20;
    options = trainingOptions('sgdm', ...
                'MiniBatchSize',mB, ...
                'MaxEpochs',50, ...
                'InitialLearnRate',learnRate, ...
                'LearnRateSchedule','piecewise', ...
                'LearnRateDropPeriod',10, ...
                'LearnRateDropFactor',0.5, ...
                'GradientThreshold',2, ...
                'Shuffle','every-epoch', ...
                'Plots','training-progress', ...
                'Verbose',true, ...
                'VerboseFrequency',floor(((numObs*8)/mB)/2), ...
                'ExecutionEnvironment','gpu', ...
                'OutputNetwork','last-iteration');
    
    layers = layerGraph(netLED);
    netLED = trainNetwork(X_Aug,Y_Aug,layers,options);
    delete(findall(0));
    reset(gpuDevice(1));
    save([ProjectPathName,'\netLED.mat'],'netLED')
end

function layersCNN = makeCNN(numTarget,inputSize)
    
    % CNN network architecture, the brains of the operation
    numUnits = 20;
    depth = 5;
    branchDepth = 4;
    scale = 0.1;

    layersCNN = [
        imageInputLayer([inputSize,3],'Normalization','none','Name','Input')
        functionLayer(@(X) (single(X)/255),'Name','zero2one')
        convolution2dLayer(7,10,'Padding','same',Name="conv7_Stem")
        batchNormalizationLayer(Name="BN_Stem_1")
        reluLayer(Name="reLu_Stem_1")
        convolution2dLayer(5,15,'Padding','same',Name="conv5_Stem")
        batchNormalizationLayer(Name="BN_Stem_2")
        reluLayer(Name="reLu_Stem_2")
        convolution2dLayer(3,20,'Padding','same',Name="conv3_Stem")
        batchNormalizationLayer(Name="BN_Stem_3")
        reluLayer(Name="reLu_Stem_3")];

    outputName = 'reLu_Stem_3';
    layersCNN = layerGraph(layersCNN);   
    
    % Encoder
    for i = 1:depth-1

        numUnits = numUnits*2;
    
        hiddenLayers = [
            %convolution2dLayer(1,numUnits,'Padding','same',Name="conv_FC_"+i)
            maxPooling2dLayer(2,'Stride',2,Name="maxpool_"+i)
            depthConcatenationLayer(2,Name="downCat_"+i)];

        layersCNN = addLayers(layersCNN,hiddenLayers);
        layersCNN = connectLayers(layersCNN,outputName,"maxpool_"+i);

        layersCNN = addLayers(layersCNN,convolution2dLayer(2,numUnits,'Stride',2,'Padding',[0,0],Name="downConv_"+i));
        layersCNN = connectLayers(layersCNN,outputName,"downConv_"+i);
        layersCNN = connectLayers(layersCNN,"downConv_"+i,"downCat_"+i+"/in2");
        outputName = "downCat_"+i;

        for ii = 1:branchDepth

            branchLayers = [
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_1_"+i+"_"+ii+".0")
                batchNormalizationLayer(Name="BN_1_"+i+"_"+ii+".0")
                reluLayer(Name="reLu_BN_1_"+i+"_"+ii+".0")
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_1_"+i+"_"+ii+".1")
                batchNormalizationLayer(Name="BN_1_"+i+"_"+ii+".1")
                reluLayer(Name="reLu_BN_1_"+i+"_"+ii+".1")
                depthConcatenationLayer(2,Name="depthCat_"+i+"_"+ii)
                convolution2dLayer(1,numUnits,'Padding','same',Name="convFC_BN_"+i+"_"+ii+".0")
                scalingLayer(Scale=scale,Name="scale_"+i+"_"+ii)
                additionLayer(2,Name="add_BN_"+i+"_"+ii+".0")
                reluLayer(Name="reLu_BN_Out_"+i+"_"+ii)];

            layersCNN = addLayers(layersCNN,branchLayers);
            layersCNN = connectLayers(layersCNN,outputName,"conv_BN_1_"+i+"_"+ii+".0");

            branchLayers2 = [
                convolution2dLayer(1,numUnits,'Padding','same',Name="conv_BN_2_"+i+"_"+ii+".0")
                batchNormalizationLayer(Name="BN_2_"+i+"_"+ii+".0")
                reluLayer(Name="reLu_BN_2_"+i+"_"+ii+".0")];            

            layersCNN = addLayers(layersCNN,branchLayers2);
            layersCNN = connectLayers(layersCNN,outputName,"conv_BN_2_"+i+"_"+ii+".0");
            layersCNN = connectLayers(layersCNN,"reLu_BN_2_"+i+"_"+ii+".0","depthCat_"+i+"_"+ii+"/in2");

            if ii == 1
                layersCNN = addLayers(layersCNN,convolution2dLayer(1,numUnits,'Padding','same',Name="conv_BN_3_"+i));
                layersCNN = connectLayers(layersCNN,outputName,"conv_BN_3_"+i);
                layersCNN = connectLayers(layersCNN,"conv_BN_3_"+i,"add_BN_"+i+"_"+ii+".0/in2");
            else
                layersCNN = connectLayers(layersCNN,outputName,"add_BN_"+i+"_"+ii+".0/in2");
            end
            outputName = "reLu_BN_Out_"+i+"_"+ii;
        end
    end
   
    % Bridge
    numUnits = numUnits*2;
    bridgeLayers = [
        convolution2dLayer(1,numUnits,'Padding','same','Name','conv_FC_Bridge')
        maxPooling2dLayer(2,'Stride',2,Name="maxpool_"+(i+1))
        depthConcatenationLayer(2,Name="downCat_Bridge")];

    layersCNN = addLayers(layersCNN,bridgeLayers);
    layersCNN = connectLayers(layersCNN,outputName,"conv_FC_Bridge");

    layersCNN = addLayers(layersCNN,convolution2dLayer(2,numUnits,'Stride',2,'Padding',[0,0],Name="downConv_Bridge"));
    layersCNN = connectLayers(layersCNN,'conv_FC_Bridge',"downConv_Bridge");
    layersCNN = connectLayers(layersCNN,"downConv_Bridge","downCat_Bridge/in2");
    outputName = "downCat_Bridge";

    for ii = 1:branchDepth

        branchLayers = [
            convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_1_Bridge_"+ii+".0")
            batchNormalizationLayer(Name="BN_Bridge_1_"+ii+".0")
            reluLayer(Name="reLu_Bridge_1_"+ii+".0")
            convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_1_Bridge_"+ii+".1")
            batchNormalizationLayer(Name="BN_Bridge_1_"+ii+".1")
            reluLayer(Name="reLu_Bridge_1_"+ii+".1")
            depthConcatenationLayer(2,Name="depthCat_Bridge_"+ii)
            convolution2dLayer(1,numUnits,'Padding','same',Name="convFC_BN_1_Bridge_"+ii+".0")
            scalingLayer(Scale=scale,Name="scale_Bridge_"+ii)
            additionLayer(2,Name="add_BN_Bridge_"+ii+".0")
            reluLayer(Name="reLu_BN_Bridge_Out_"+i+"_"+ii)];
        
        layersCNN = addLayers(layersCNN,branchLayers);
        layersCNN = connectLayers(layersCNN,outputName,"conv_BN_1_Bridge_"+ii+".0");

        branchLayers2 = [
            convolution2dLayer(1,numUnits,'Padding','same',Name="conv_BN_2_Bridge_"+ii+".0")
            batchNormalizationLayer(Name="BN_Bridge_2_"+ii+".0")
            reluLayer(Name="reLu_Bridge_2_"+ii+".0")];        

        layersCNN = addLayers(layersCNN,branchLayers2);
        layersCNN = connectLayers(layersCNN,outputName,"conv_BN_2_Bridge_"+ii+".0");
        layersCNN = connectLayers(layersCNN,"reLu_Bridge_2_"+ii+".0","depthCat_Bridge_"+ii+"/in2");

        if ii == 1
            layersCNN = addLayers(layersCNN,convolution2dLayer(1,numUnits,'Padding','same',Name="conv_BN_3_Bridge"));
            layersCNN = connectLayers(layersCNN,outputName,"conv_BN_3_Bridge");
            layersCNN = connectLayers(layersCNN,"conv_BN_3_Bridge","add_BN_Bridge_"+ii+".0/in2");
        else
            layersCNN = connectLayers(layersCNN,outputName,"add_BN_Bridge_"+ii+".0/in2");
        end
        outputName = "reLu_BN_Bridge_Out_"+i+"_"+ii;
    end
    
    % Decoder
    Count = 1;
    for i = depth-1:-1:1
    
        numUnits = numUnits/2;
        hiddenLayers = [
            resize2dLayer('Scale',2,Name="resize_"+i)
            depthConcatenationLayer(2,Name="upCat_"+i)
            depthConcatenationLayer(2,Name="Fuse_"+Count)
            convolution2dLayer(1,numUnits,'Padding','same',Name="conv_up_FC_"+i)]; 
            
        layersCNN = addLayers(layersCNN,hiddenLayers);
        layersCNN = connectLayers(layersCNN,outputName,"resize_"+i);

        layersCNN = addLayers(layersCNN,transposedConv2dLayer(2,numUnits*2,'Stride',2,Name="upConv_"+i));
        layersCNN = connectLayers(layersCNN,outputName,"upConv_"+i);
        layersCNN = connectLayers(layersCNN,"upConv_"+i,"upCat_"+i+"/in2");
        outputName = "conv_up_FC_"+i;
        
        for ii = 1:branchDepth

            branchLayers = [
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_up_"+i+"_"+ii+".0")
                batchNormalizationLayer(Name="BN_BN_up_"+i+"_"+ii+".0")
                reluLayer(Name="reLu_BN_up_"+i+"_"+ii+".0")
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_up_"+i+"_"+ii+".1")
                batchNormalizationLayer(Name="BN_BN_up_"+i+"_"+ii+".1")
                reluLayer(Name="reLu_BN_up_"+i+"_"+ii+".1")
                depthConcatenationLayer(2,Name="depthCat_up_"+i+"_"+ii)
                convolution2dLayer(1,numUnits,'Padding','same',Name="conv_BN_FC_up_"+i+"_"+ii)
                scalingLayer(Scale=scale,Name="scale_up_"+i+"_"+ii)
                additionLayer(2,Name="add_BN_up_"+i+"_"+ii+".0")
                reluLayer(Name="reLu_BN_Up_Out_"+i+"_"+ii)];

            branchLayers2 = [
                convolution2dLayer(1,numUnits,'Padding','same',Name="conv_BN_up_1_"+i+"_"+ii+".0")
                batchNormalizationLayer(Name="BN_BN_up_1"+i+"_"+ii+".0")
                reluLayer(Name="reLu_BN_up_1_"+i+"_"+ii+".0")];

            layersCNN = addLayers(layersCNN,branchLayers);
            layersCNN = connectLayers(layersCNN,outputName,"conv_BN_up_"+i+"_"+ii+".0");

            layersCNN = addLayers(layersCNN,branchLayers2);
            layersCNN = connectLayers(layersCNN,outputName,"conv_BN_up_1_"+i+"_"+ii+".0");

            layersCNN = connectLayers(layersCNN,"reLu_BN_up_1_"+i+"_"+ii+".0","depthCat_up_"+i+"_"+ii+"/in2");
            layersCNN = connectLayers(layersCNN,outputName,"add_BN_up_"+i+"_"+ii+".0/in2");
            outputName = "reLu_BN_Up_Out_"+i+"_"+ii;
        end
        Count = Count+1;
    end

    outputLayers = [
        resize2dLayer('Scale',2,Name="resize_"+0)
        depthConcatenationLayer(2,Name="upCat_"+0)
        depthConcatenationLayer(2,Name="Fuse_"+Count)
        convolution2dLayer(3,numUnits/2,'Padding','same',Name="conv_Output")
        batchNormalizationLayer(Name="BN_Output")
        reluLayer(Name="reLu_Output")
        convolution2dLayer(1,numTarget,'Padding','same',Name="conv_FC_Target_1.0")
        depthConcatenationLayer(2,Name="Fuse_"+(Count+1))
        convolution2dLayer(1,numTarget,'Padding','same',Name="conv_FC_Target_1.1")];
    layersCNN = addLayers(layersCNN,outputLayers);
    layersCNN = connectLayers(layersCNN,outputName,"resize_"+0);

    layersCNN = addLayers(layersCNN,transposedConv2dLayer(2,numUnits,'Stride',2,Name="upConv_"+0));
    layersCNN = connectLayers(layersCNN,outputName,"upConv_"+0);
    layersCNN = connectLayers(layersCNN,"upConv_"+0,"upCat_"+0+"/in2");

    Count = depth-1;
    for i = 1:depth+1
        if i < depth
            layersCNN = connectLayers(layersCNN,"reLu_BN_Out_"+Count+"_"+branchDepth,"Fuse_"+i+"/in2");
        elseif i == depth
            layersCNN = connectLayers(layersCNN,"reLu_Stem_3","Fuse_"+i+"/in2");
        elseif i == depth+1
            layersCNN = connectLayers(layersCNN,"zero2one","Fuse_"+i+"/in2");
        end
        Count = Count-1;
    end

    layersCNN = addLayers(layersCNN,regressionLayer('Name','RegOut'));
    layersCNN = connectLayers(layersCNN,'conv_FC_Target_1.1','RegOut');
end

function frameData = pullFrames(numTrainingFrames,ClusterPool,numClust)

    [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select training video','MultiSelect','off');
    cd(VidFilePath)
    FullName = {fullfile(VidFilePath,VidFileName)};

    % Create video reader
    disp('Reading frames.')
    vr = VideoReader(FullName{1,1});
    H = vr.Height;
    W = vr.Width;
    numFrames = vr.NumFrames; 
    
    % Pull frames from video and cluster
    idx = round(linspace(1,numFrames-1,ClusterPool));   
    frameData = uint8(zeros(H,W,3,ClusterPool));
    frameVect = zeros(2500,ClusterPool);
    iFrame = 1;
    i = 0;
    while iFrame <= ClusterPool 
        
        try
            frame = read(vr,idx(1,iFrame));
        catch
            break
        end
        iFrame = iFrame+1;
        if mean(frame,'all') > 10
            i = i+1;
            frameData(:,:,:,i) = frame;
            GS = rgb2gray(frame);
            GS = imresize(GS,[50,50]);
            frameVect(:,i) = double(GS(:));
        end
    end
    frameData = frameData(:,:,:,1:i);
    frameVect = frameVect(:,1:i);    
    
    disp('Clustering frames.')
    clust = kmeans(frameVect',numClust);
    clear frameVect
        
    % Group clustered frames for labeling
    idx = zeros(0,1);
    Add = 0;
    for i = 1:numClust
        
        loc = find(clust == i);
        if length(loc) >= numTrainingFrames/numClust
            idx = [idx;datasample(loc,numTrainingFrames/numClust,'Replace',false)]; %#ok<AGROW> 
        else
            idx = [idx;datasample(loc,length(loc),'Replace',false)]; %#ok<AGROW> 
            Add = Add+(numTrainingFrames/numClust)-length(loc);
        end
    end   

    if Add > 0
        
        pull = (1:ClusterPool)';
        Count = 1;
        while Count <= Add

            p = datasample(pull,1);
            if isempty(find(idx == p,1))
                idx = [idx;p]; %#ok<AGROW> 
                Count = Count+1;
            end
        end            
    end
    frameData = frameData(:,:,:,idx);
end

function frameData = pickFrames(netCNN,CB)

    [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select training video','MultiSelect','off');
    cd(VidFilePath)
    FullName = {fullfile(VidFilePath,VidFileName)};

    inputSize = netCNN.Layers(1,1).InputSize(1,1:2);
    numTarget = netCNN.Layers(end-2,1).NumFilters;

    % Create video reader
    disp('Reading frames.')
    vr = VideoReader(FullName{1,1});
    H = vr.Height;
    W = vr.Width;
    fps = vr.FrameRate;
    numFrames = floor(fps*vr.Duration);

    % Create object detector to detect CNN parts
    blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort',false,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',10);

    definput = {'0','0','1','5'};
    frameData = uint8(zeros(H,W,3,0));
    Count = 1;
    figure
    set(gcf,'Position',get(0,'Screensize'));
    while 1
    
        dlgtitle = '';
        prompt = {'Hours:','Minutes:','Seconds','Sample length(s):'};
        dims = [1,10];
        t = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(t) == 1
            break
        end

        tNum(1,1) = str2double(t{1,1}); 
        tNum(2,1) = str2double(t{2,1}); 
        tNum(3,1) = str2double(t{3,1}); 
        tNum(4,1) = str2double(t{4,1}); 
                        
        s = floor(((tNum(1,1)*3600)+(tNum(2,1)*60)+tNum(3,1))*fps);
        L = tNum(4,1);
        if s+(L*fps) > numFrames
            wBox = warndlg('Not enough frames in file.','Warning');
            waitfor(wBox)
            continue
        end
       
        tNum(3,1) = tNum(3,1)+L;
        if tNum(3,1) >= 60

            tNum(2,1) = tNum(2,1)+1;
            tNum(3,1) = tNum(3,1)-60;

            if tNum(2,1) >= 60

                tNum(1,1) = tNum(1,1)+1;
                tNum(2,1) = tNum(2,1)-60;
            end
        end

        for i = 1:4
            definput{1,i} = num2str(tNum(i,1));
        end
        
        % Read frames, resize and predict
        frame = read(vr,[s,s+(fps*L)-1]);
        frameSize = size(frame,[1,2]);
        numFrame = size(frame,4);
        frameR = imresize(frame,inputSize);
        reg = predict(netCNN,frameR);
        
        xy = zeros(numFrame,2,numTarget);
        vel = zeros(1,numFrame);
        lastPoint = [0,0];
        Con = zeros(numTarget,(fps*L));
        for i = 1:size(reg,4)
            
            clf
            subplot(2,1,1)
            title('CNN Confidence')
            hold on
            for ii = 1:numTarget
                plot(Con(ii,:),'Color',CB(ii,:)/255)
            end
            plot(vel(1,:),'Color','k')
            hold off
            xlim([1,numFrame])
            ylim([0,1])
    
            subplot(2,1,2)
            f = frame(:,:,:,i);
            imshow(f)
            hold on
            for ii = 1:numTarget
            
                P = reg(:,:,ii,i);
                P(P < 0.25) = 0;
                P(P > 1) = 1;

                mask = logical(imresize(P,frameSize));
                [area,centroids] = blobAnalyser(mask);
                if isempty(centroids) == 0

                    [~,aMax] = max(area);
                    xy(i,:,ii) = round(centroids(aMax,:));
                    plot(xy(i,1,ii),xy(i,2,ii),'o-','MarkerFaceColor',CB(ii,:)/255,'MarkerEdgeColor',CB(ii,:)/255);

                    [counts,centers] = hist(P(P > 0),linspace(0,1,20)); %#ok<HIST> 
                    [~,idx] = max(counts);
                    Con(ii,i) = centers(idx);
                else
                    xy(i,:,ii) = [nan,nan];
                    Con(ii,i) = 0;
                end
            end

            mXY = [mean(xy(i,1,:),'omitnan'),mean(xy(i,2,:),'omitnan')];
            pd = pdist([lastPoint;mXY],'euclidean');
            vel(1,i) = (pd/(1/fps))/10000;
            lastPoint = mXY;

            drawnow
        end
    
        while 1
    
            [sel(1,1),sel(1,2),~] = ginput(1);
            x = round(sel(1,1));

            if x > numFrame
                break
            end           
            
            subplot(2,1,2)
            f = frame(:,:,:,x);
            imshow(f)
            hold on
            for ii = 1:numTarget     
                plot(xy(x,1,ii),xy(x,2,ii),'o-','MarkerFaceColor',CB(ii,:)/256,'MarkerEdgeColor',CB(ii,:)/256);
            end
            hold off
    
            promptMessage = sprintf('Add frame to training data?');
            titleBarCaption = 'settings';
            answer = questdlg(promptMessage,titleBarCaption,'Yes','No','Exit','Yes');
            
            if strcmp(answer,'Yes') == 1
                frameData(:,:,:,Count) = f;
                Count = Count+1;
            elseif strcmp(answer,'Exit') == 1
                break
            end
        end
    
        promptMessage = sprintf('Select another time?');
        titleBarCaption = 'settings';
        answer = questdlg(promptMessage,titleBarCaption,'Yes','No','Yes');
    
        if strcmp(answer,'No') == 1
            break
        end
    end
    close(gcf)            
end

function xyLabels = labelFrames(frameData,Targets,CB)

    % Begin UI labeling
    
    numTrainingFrames = size(frameData,4);
    numTarget = length(Targets);
    xyLabels = zeros(numTrainingFrames,2,numTarget);
    goBack = 0;
    Stop = 0;
    BaseSize = [480,640];
    f = figure('Resize','on');
    f.Position = [0,41,1920,963];
    i = 0;            
    while i < numTrainingFrames

        if goBack == 0
            i = i+1;
        else
            i = i-1;  
        end
        
        frame = frameData(:,:,:,i);
        if any(size(frame) ~= [BaseSize,3])
            frame = imresize(frame,BaseSize);
        end
        imshow(frame);
        hold on
        f.Position = [0,41,1920,963];

        if goBack == 1
            for ii = 1:numTarget-1
                 PlotPoint(1,ii) = plot(xyLabels(i,1,ii),xyLabels(i,2,ii),'o-','MarkerFaceColor',CB(ii,:)/255,'MarkerEdgeColor',CB(ii,:)/255); %#ok<AGROW>
            end
        end
        
        sel = ones(1,2);
        if goBack == 0
            ii = 1;
        else
            ii = numTarget;
            goBack = 0;
        end

        while ii <= numTarget
            
            title(['Selecting: ',Targets{ii,1},'    Frame ',num2str(i),'/',num2str(numTrainingFrames),'    Left-click to select, right-click to skip, middle to go back'])
            try
                [sel(1,1),sel(1,2),button] = ginput(1);
            catch
                Stop = 1;
                break
            end
            if button == 1
                if all(sel > 0) && sel(1,1) <= BaseSize(1,2) && sel(1,2) <= BaseSize(1,1)

                    sel = round(sel);
                    xyLabels(i,:,ii) = sel;
                    PlotPoint(1,ii) = plot(sel(1,1),sel(1,2),'o-','MarkerFaceColor',CB(ii,:)/255,'MarkerEdgeColor',CB(ii,:)/255);
                    ii = ii+1;
                else
                    Stop = 2;
                    break
                end
            elseif button == 2
                if ii > 1
                    ii = ii-1;
                    delete(PlotPoint(1,ii));
                else
                    if i > 1
                        goBack = 1;
                        break
                    end
                end                        
            elseif button == 3
                PlotPoint(1,ii) = plot(0);
                ii = ii+1;
            end
        end
        hold off                
        if Stop == 1
            break
        elseif Stop == 2
            Stop = 0;
            continue
        end                
    end
    if ishandle(f) == 1
        close(f)
    end  
end

function [trainingData,validationData] = makeTrain(frameData,xyLabels,inputSize,SymPair)
    
    numObs = size(xyLabels,1);
    numTarget = size(xyLabels,3);
    frameSize = size(frameData,[1,2]);
    tX = uint8(zeros([inputSize,3,numObs*8]));
    tY = single(zeros([inputSize,numTarget,numObs*8]));

    % Create vector of all points in image for alphashape query
    QP = zeros(prod(frameSize),2);
    Count = 1;
    for i = 1:frameSize(1,2)
        for ii = 1:frameSize(1,1)
            QP(Count,:) = [i,ii]; 
            Count = Count+1;
        end
    end
    
    % Generate training data
    disp('Generating training data')
    Count = 1;
    for i = 1:numObs
        
        frame = frameData(:,:,:,i);
        xy = xyLabels(i,:,:);
        Y = single(zeros([frameSize,numTarget]));
        for ii = 1:numTarget
        
            x = xy(1,1,ii);
            y = xy(1,2,ii);
    
            if x == 0 && y == 0
                continue
            end
            
            % Label ground truth frames
            p = nsidedpoly(1000,'Center',[x,y],'Radius',8);        
            shp = alphaShape(p.Vertices);
            a = criticalAlpha(shp,'all-points')*2;
            shp.Alpha = a;
            
            IS = inShape(shp,QP);
            IS = find(IS == 1);
           
            for i3 = 1:length(IS)
                row = QP(IS(i3,1),2); 
                col = QP(IS(i3,1),1);
                Y(row,col,ii) = 1;
            end               
        end

        frame = imresize(frame,inputSize);
        Y = imresize(Y,inputSize);
        
        for ii = 1:4           

            tX(:,:,:,Count) = frame;
            tY(:,:,:,Count) = Y;
            Count = Count+1;        
                    
            % Reflect
            frame_R = flip(frame,1);
            Y_R = flip(Y,1);
            if size(SymPair,2) > 1
                for i3 = 1:size(SymPair,1)                                
                    y = Y_R(:,:,SymPair(1,1));
                    Y_R(:,:,SymPair(1,1)) = Y_R(:,:,SymPair(1,2));
                    Y_R(:,:,SymPair(1,2)) = y;
                end
            end
            tX(:,:,:,Count) = frame_R;
            tY(:,:,:,Count) = Y_R;
            Count = Count+1;
            
            % Rotate 90
            frame = imrotate(frame,90);
            Y = imrotate(Y,90);
        end 
    end

    numObs = size(tX,4);
    idx = randperm(numObs);
    tX = tX(:,:,:,idx);
    tY = tY(:,:,:,idx);

    numVal = rem(numObs,10);
    if numObs >= 100
        numVal = numVal+50;
    elseif numObs >= 20 && numObs < 100
        numVal = numVal+10;
    end

    vX = tX(:,:,:,1:numVal);
    vY = tY(:,:,:,1:numVal);
    tX(:,:,:,1:numVal) = [];
    tY(:,:,:,1:numVal) = [];

    trainingData = {tX,tY};
    validationData = {vX,vY};
end

function [layersCNN,netCNN] = trainCNN(ProjectPathName,Targets,Data,layersCNN,miniBatchSizeCNN,learnRate)

    imDirTrain = [ProjectPathName,'\Training\Data'];
    listing = dir(imDirTrain);
    numFile = size(listing,1)-2;
    if numFile > 0
        imds = fileDatastore(imDirTrain,'ReadFcn',@load2VarLabeled);
        imds = shuffle(imds);
    end

    if ~isempty(Data)

        trainingData = Data{1,1};
        validationData = Data{1,2};
        numData = size(trainingData{1,1},4);
    else
        numData = 0;
    end
       
    allObs = numFile+numData;

    inputSize = layersCNN.Layers(1,1).InputSize;
    outputSize = [inputSize(1,1:2),layersCNN.Layers(end-2,1).NumFilters];

    if allObs > 200
        if numData == 0

            numVal = rem(numFile,miniBatchSizeCNN);
            if numVal < 50
                numVal = numVal+(miniBatchSizeCNN*round(50/miniBatchSizeCNN));
            end

            vX = uint8(zeros([inputSize,numVal]));
            vY = single(zeros([outputSize,numVal]));
            for i = 1:numVal
                readDS = read(imds);
                vX(:,:,:,i) = readDS{1,1};
                vY(:,:,:,i) = readDS{1,2};
            end
            imds.Files(1:numVal,:) = [];
            
            if size(imds.Files,1) <= 1000
                usePiece = 0;
            else

                usePiece = 1;
                tX = uint8(zeros([inputSize,1000]));
                tY = single(zeros([outputSize,1000]));
                for i = 1:1000
                    readDS = read(imds);
                    tX(:,:,:,i) = readDS{1,1};
                    tY(:,:,:,i) = readDS{1,2};
                end
            end
            S = [1,2];
        else

            idx = randperm(numData);
            tX = trainingData{1,1}(:,:,:,idx);
            tY = trainingData{1,2}(:,:,:,idx);

            vX = validationData{1,1};
            vY = validationData{1,2};
            
            for i = size(tX,4)+1:size(tX,4)*2
                readDS = read(imds);
                tX(:,:,:,i) = readDS{1,1};
                tY(:,:,:,i) = readDS{1,2};
            end
            S = 3;
        end

        for i = S
        
            % Set options for network training
            if i == 1

                if usePiece == 1
                    numObs = 1000;
                    numEpoch = 5;
                else
                    numObs = size(imds.Files,1);
                    numEpoch = round(500/(numObs/miniBatchSizeCNN));
                end
                learnRateDrop = 5;
                valPat = 5;
                outPut = 'last-iteration';
            elseif i == 2

                numObs = size(imds.Files,1);
                numEpoch = round(50000/(numObs/miniBatchSizeCNN));
                learnRateDrop = 5;
                valPat = 10;
                outPut = 'best-validation-loss';
            elseif i == 3

                numObs = size(tX,4);
                numEpoch = 30;
                learnRateDrop = 30;
                valPat = 5;
                outPut = 'best-validation-loss';
            end

            options = trainingOptions('sgdm', ...
            'MiniBatchSize',miniBatchSizeCNN, ...
            'MaxEpochs',numEpoch, ...
            'InitialLearnRate',learnRate, ...
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropPeriod',learnRateDrop, ...
            'LearnRateDropFactor',0.5, ...
            'GradientThreshold',2, ...
            'Shuffle','every-epoch', ...
            'ValidationData',{vX,vY}, ... 
            'ValidationFrequency',floor(numObs/miniBatchSizeCNN), ...
            'ValidationPatience',valPat, ...
            'Plots','training-progress', ...
            'Verbose',true, ...
            'VerboseFrequency',floor((numObs/miniBatchSizeCNN)/4), ...
            'ExecutionEnvironment','gpu', ...
            'OutputNetwork',outPut);

            if i == 1 || i == 2
                if usePiece == 0
                    netCNN = trainNetwork(imds,layersCNN,options);
                else
                    netCNN = trainNetwork(tX,tY,layersCNN,options);
                    usePiece = 0;
                end
            else
                netCNN = trainNetwork(tX,tY,layersCNN,options);
            end
            delete(findall(0));
            layersCNN = layerGraph(netCNN);
        end
        save([ProjectPathName,'\netCNN.mat'],'Targets','netCNN')
        disp('Your network is now trained =D')
    else
        disp('Please add more training frames.')
    end

    % Datastore read function
    function Out = load2VarLabeled(In)
        data = load(In);
        name = fieldnames(data);
        Out = data.(name{1,1});   
    end
end

function saveTrain(ProjectPathName,Data)

    imDirTrain = [ProjectPathName,'\Training\Data'];
    listing = dir(imDirTrain);
    numFile = size(listing,1)-2;
    
    disp('Saving data')
    for i = 1:length(Data)

        numObs = size(Data{1,i}{1,1},4);
        for ii = 1:numObs
            numFile = numFile+1;
            XY = {Data{1,i}{1,1}(:,:,:,ii),Data{1,i}{1,2}(:,:,:,ii)};
            save([imDirTrain,'\image_',num2str(numFile),'.mat'],'XY')
        end
    end
    disp('Done')
end

function layers3dCNN = make3dCNN(numClass,inputSize)
    
    % 3dCNN network architecture, the brains of the operation

    numFilt = 5;

    inputLayers = [
    image3dInputLayer(inputSize,'Normalization','none','Name','Input')
    functionLayer(@(X) (single(X)/255),'Name','zero2one')
    convolution3dLayer([7,7,1],numFilt,'Padding','same','Name',"conv3d_space_int")
    batchNormalizationLayer(Name="BN_space_int")
    reluLayer('Name',"reLu_space_int")
    convolution3dLayer([1,1,5],numFilt,'Padding',[0,0,1],'Name',"conv3d_time_int")
    batchNormalizationLayer(Name="BN_time_int")
    reluLayer('Name',"reLu_time_int")];

    outputName = 'reLu_time_int';
    layers3dCNN = layerGraph(inputLayers);
    numFilt = numFilt*2;

    for i = 1:6
    
        if rem(i,2) == 0
            tMP = 2;
        else
            tMP = 1;
        end
        
        hiddenLayers = [
        convolution3dLayer([3,3,1],numFilt,'Padding','same','Name',"conv3d_space"+i)
        batchNormalizationLayer(Name="BN_"+i+".0")
        reluLayer('Name',"reLu_"+i+".0")
        convolution3dLayer([1,1,3],numFilt,'Padding','same','Name',"conv3d_time"+i)
        batchNormalizationLayer(Name="BN_"+i+".1")
        reluLayer('Name',"reLu_"+i+".1")
        additionLayer(2,'Name',"add_"+i)
        maxPooling3dLayer([2,2,tMP],'Stride',[2,2,tMP],'Name',"maxPool_"+i)];
        
        layers3dCNN = addLayers(layers3dCNN,hiddenLayers);
        layers3dCNN = connectLayers(layers3dCNN,outputName,"conv3d_space"+i);
        
        skipLayers = [
            convolution3dLayer(1,numFilt,'Padding','same','Name',"conv3d_skip_"+i)
            batchNormalizationLayer(Name="BN_skip_"+i)
            reluLayer('Name',"reLu_skip_"+i)];
    
        layers3dCNN = addLayers(layers3dCNN,skipLayers);
        layers3dCNN = connectLayers(layers3dCNN,outputName,"conv3d_skip_"+i);
        layers3dCNN = connectLayers(layers3dCNN,"reLu_skip_"+i,"add_"+i+"/in2");
        outputName = "maxPool_"+i;
        numFilt = numFilt*2;
    end        
    
    outputLayers = [
        flattenLayer('Name','flatten')
        fullyConnectedLayer(numClass,'Name','fc_Out')
        softmaxLayer('Name','softMax')
        classificationLayer('Name','Output')];
    
    layers3dCNN = addLayers(layers3dCNN,outputLayers);
    layers3dCNN = connectLayers(layers3dCNN,outputName,"flatten");
end

function sequenceData = pickVideo(ProjectPathName,inputSize,net3dCNN,CB)

    [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select training video','MultiSelect','off');
    cd(VidFilePath)
    FullName = fullfile(VidFilePath,VidFileName);

    try
        load([FullName(1,1:end-4),'_poseData.mat'],'confidence');

        for i = 1:size(confidence,2)
            confidence(:,i) = smooth(confidence(:,i));
        end

        idx = find(confidence(:,7) > 0.8);
        idx = unique(round(idx*0.1));
        times = zeros(length(idx),3);
        
        for i = 1:length(idx)
        
            d = idx(i,1)/60;
            if d >= 1
                r = rem(idx(i,1),60);
                times(i,3) = r;
                times(i,2) = floor(d);
        
                d2 = times(i,2)/60;
                if d2 >= 1
                    r2 = rem(times(i,2),60);
                    times(i,1) = r2;
                    times(i,2) = floor(d2);
                end
            else
                 times(i,3) = idx(i,1);
            end
        end
    catch

    end

    seqDirTrain = [ProjectPathName,'\Training\Sequences'];
    listing = dir(seqDirTrain);
    listing(1:2) = [];
    numClass = size(listing,1);
    classCat = categorical(1:numClass);

    sequenceData = {listing.name}';
    for i = 1:numClass
        dataDir = [listing(i).folder,'\',listing(i).name];
        listingClass = dir(dataDir);
        listingClass(1:2) = [];
        sequenceData{i,2} = size(listingClass,1);
    end

    % Create video reader
    disp('Reading frames.')
    vr = VideoReader(FullName);
    fps = vr.FrameRate;
    frames = read(vr,[1,10]); %#ok<NASGU>

    definput = {'0','0','1','10'};
    vidFig = figure;
    set(vidFig,'Position',get(0,'Screensize'));
    labelAx = cell(numClass,1);
    while 1
    
        dlgtitle = '';
        prompt = {'Hours:','Minutes:','Seconds','Sample length(s):'};
        dims = [1,10];
        t = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(t) == 1
            break
        end

        tNum(1,1) = str2double(t{1,1}); 
        tNum(2,1) = str2double(t{2,1}); 
        tNum(3,1) = str2double(t{3,1}); 
        tNum(4,1) = str2double(t{4,1}); 
                        
        s = floor(((tNum(1,1)*3600)+(tNum(2,1)*60)+tNum(3,1))*fps);
        L = tNum(4,1);
        if s+(L*fps) > vr.NumFrames
            wBox = warndlg('Not enough frames in file.','Warning');
            waitfor(wBox)
            continue
        end
       
        tNum(3,1) = tNum(3,1)+L;
        if tNum(3,1) >= 60

            tNum(2,1) = tNum(2,1)+1;
            tNum(3,1) = tNum(3,1)-60;
            if tNum(2,1) >= 60
                tNum(1,1) = tNum(1,1)+1;
                tNum(2,1) = tNum(2,1)-60;
            end
        end

        for i = 1:4
            definput{1,i} = num2str(tNum(i,1));
        end
        
        % Read frames, resize and predict
        frames = read(vr,[s,s+(fps*L)]);
        frameR = imresize(frames,inputSize(1,1:2));
        Con = zeros(numClass,(fps*L));
        
        plotAx = subplot(2,1,1,'Parent',vidFig);
        title('CNN Confidence')
        hold on
        for i = 1:numClass
            labelAx{i,1} = plot(Con(i,:),'Color',CB(i,:)/255,'Parent',plotAx);
        end
        hold off
        xlim([1,fps*L])
        ylim([0,1])
        legend(sequenceData(:,1))
        
        inputData = uint8(zeros([inputSize(1,1:2),size(frameR,4)-1,2]));
        for i = 2:size(frames,4)

            frame = rgb2gray(frameR(:,:,:,i));
            inputData(:,:,i-1,1) = frame;

            diffFrame = frame-rgb2gray(frameR(:,:,:,i-1));
            diffFrame(diffFrame < 15) = 0;
            inputData(:,:,i-1,2) = diffFrame;
            
            if ~isempty(net3dCNN)
                if i > inputSize(1,3)
                    Con(:,i-1) = predict(net3dCNN,inputData(:,:,i-inputSize(1,3):i-1,:));
                end
            end
            
            for ii = 1:numClass
                set(labelAx{ii,1},'YData',Con(ii,:))
            end
    
            subplot(2,1,2)
            f = frames(:,:,:,i);
            imshow(f)
            drawnow
        end
    
        while 1

            ClassCount = cell(numClass,1);
            for ii = 1:numClass
                ClassCount{ii,1} = [sequenceData{ii,1},' (',num2str((sequenceData{ii,2})),')'];
            end  
    
            [sel(1,1),sel(1,2),~] = ginput(1);
            x = round(sel(1,1));
            
            if x+inputSize(1,3)-1 < size(frames,4)

                for i = x:x+inputSize(1,3)-1
                    subplot(2,1,2)
                    f = frames(:,:,:,i);
                    imshow(f)
                    pause(0.1)
                    drawnow
                end            
               
                promptMessage = sprintf('Add sequence to training data?');
                titleBarCaption = 'settings';
                answer = questdlg(promptMessage,titleBarCaption,'Yes','No','Exit','Yes');
                
                if strcmp(answer,'Yes') == 1
                    
                    while 1
                        pick = listdlg('PromptString',{'Select a class,';'Press cancel to exit.'}, ... 
                        'ListString',[{'Play Again'};ClassCount;],'SelectionMode','single','InitialValue',1);
        
                        if isempty(pick)
                            break
                        elseif pick == 1
        
                            for i = x:x+inputSize(1,3)-1
                                subplot(2,1,2)
                                f = frames(:,:,:,i);
                                imshow(f)
                                pause(0.1)
                                drawnow
                            end
                        elseif pick > 1

                            X = inputData(:,:,x:x+inputSize(1,3)-1,:);
                            for i = 1:4

                                XY = {X,classCat(1,pick-1)};
                                sequenceData{pick-1,2} = sequenceData{pick-1,2}+1;
                                save([ProjectPathName,'\Training\Sequences\',sequenceData{pick-1,1},'\seq_',num2str(sequenceData{pick-1,2}),'.mat'],'XY')
                                X = imrotate(X,90);
                            end
                            break
                        end
                    end
                elseif strcmp(answer,'Exit') == 1
                    break
                end
            else
                break
            end
        end
    
        promptMessage = sprintf('Select another time?');
        titleBarCaption = 'settings';
        answer = questdlg(promptMessage,titleBarCaption,'Yes','No','Yes');
        
        if strcmp(answer,'Yes') == 1
            cla(plotAx)
        elseif strcmp(answer,'No') == 1
            break
        end
    end
    close(gcf)
end

function [layers3dCNN,net3dCNN] = train3dCNN(ProjectPathName,Classes,layers3dCNN,miniBatchSize3dCNN,learnRate)

    seqDirTrain = [ProjectPathName,'\Training\Sequences'];
    listing = dir(seqDirTrain);
    listing(1:2) = [];
    numClass = size(listing,1);

    inputSize = layers3dCNN.Layers(1,1).InputSize;

    if numClass > 0
        
        numVal = 5;
        vX = uint8(zeros([inputSize,numVal*numClass]));
        vY = zeros([numVal*numClass,1]);

        targetDS = cell(numClass,1);
        numObs = zeros(numClass,1);
        countVal = 0;
        for i = 1:numClass

            dataDir = [listing(i).folder,'\',listing(i).name];
            targetDS{i,1} = fileDatastore(dataDir,'ReadFcn',@load2VarLabeled,'IncludeSubfolders',false);
            targetDS{i,1} = shuffle(targetDS{i,1});
            
            for ii = 1:numVal
                countVal = countVal+1;
                readDS = read(targetDS{i,1});
                vX(:,:,:,:,countVal) = readDS{1,1};
                vY(countVal,1) = readDS{1,2};
            end
            targetDS{i,1}.Files(1:numVal,:) = [];

            numObs(i,1) = length(targetDS{i,1}.Files);
        end
        numObs = min(numObs);
        vY = categorical(vY);

        for i = 1:numClass
            targetDS{i,1}.Files(numObs+1:end) = [];
        end
        
        imds = targetDS{1,1};
        for i = 2:numClass
            imds.Files = [imds.Files;targetDS{i,1}.Files];
        end
        imds = shuffle(imds);
    end

    numObs = numObs*numClass;

    outPut = 'best-validation-loss';

    options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize3dCNN, ...
    'MaxEpochs',250, ...
    'InitialLearnRate',learnRate, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',50, ...
    'LearnRateDropFactor',0.5, ...
    'GradientThreshold',2, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{vX,vY}, ... 
    'ValidationFrequency',floor(numObs/miniBatchSize3dCNN), ...
    'ValidationPatience',inf, ...
    'Plots','training-progress', ...
    'Verbose',true, ...
    'VerboseFrequency',10, ...
    'ExecutionEnvironment','gpu', ...
    'OutputNetwork',outPut);

    net3dCNN = trainNetwork(imds,layers3dCNN,options);
    layers3dCNN = layerGraph(net3dCNN);

    save([ProjectPathName,'\net3dCNN.mat'],'Classes','net3dCNN')
    disp('Your network is now trained =D')

    % Datastore read function
    function Out = load2VarLabeled(In)
        data = load(In);
        name = fieldnames(data);
        Out = data.(name{1,1});
    end
end

function [digData,tsN] = scanDigitalin(digInFile,Sfreq)
                
    disp('Reading digitalin file:')
    disp(digInFile)
    
    s = dir(digInFile);
    numSamples = s.bytes/2;
    numSamplesDS = floor(numSamples/(Sfreq/1000));
    digData = zeros(numSamplesDS,1);
        
    % Read data from digitalin.dat and downsample to 1kHz
    sampleCount = 0;
    fileID = fopen(digInFile);
    fseek(fileID,0,'bof');
    wb = waitbar(sampleCount/numSamplesDS,{'Reading file...';['Sample ',num2str(sampleCount),'/',num2str(numSamples)]});
    while 1
        temp = fread(fileID,[1,Sfreq*100],'int16');
        if isempty(temp) == 0
            tempDS = temp(1,1:Sfreq/1000:end);
    
            for ii = 1:length(tempDS)
                sampleCount = sampleCount+1;
                digData(sampleCount,1) = tempDS(1,ii);
            end
            waitbar(sampleCount/numSamplesDS,wb,{'Reading file...';['Sample ',num2str(sampleCount),'/',num2str(numSamples)]});
        else
            close(wb)
            fclose(fileID);
            break
        end
    end
    
    % Find incorrectly coded values
    idx = find(abs(diff(digData)) > 3,1,'last');
    if ~isempty(idx)
        digData(digData > 3) = digData(digData > 3)-digData(idx,1);
    end

    % Find the channel with LED pulses
    LEDchannel = mode(digData(digData > 0));
    if LEDchannel == 1
        phaseChannel = 2;
    else
        phaseChannel = 1;
    end

    if LEDchannel == 1
        digData(digData > phaseChannel) = LEDchannel;
    else
        digData(digData > LEDchannel) = phaseChannel;
    end

    % Analyse data to find all pulses and count samples between
    Ck = 2;
    iPulse = 1;
    tsN = [0,0,0];
    sampleCount = 0;
    for ii = 1:length(digData)

        sampleCount = sampleCount+1;
        if digData(ii,1) ~= LEDchannel
            if Ck == 3
                Ck = 2;
            end
            tsN(iPulse,3) = tsN(iPulse,3)+1;
        elseif digData(ii,1) == LEDchannel
            if Ck == 2
                iPulse = iPulse+1;
                tsN(iPulse,:) = 0;
                tsN(iPulse,1) = sampleCount;
                Ck = 3;
            end
            tsN(iPulse,2) = tsN(iPulse,2)+1;
        end
    end

    digData(digData ~= LEDchannel) = 0;
    digData(digData ~= 0) = 1;

    tsN = tsN/1000;
end

function [fps,diffPix,conLED,XY] = scanVideo(vidFile,netLED,netCNN,net3dCNN)

    disp('Reading video file:')
    disp(vidFile)

    LEDPoseClass = [0,0,0];
    if ~isempty(netLED)
        LEDPoseClass(1,1) = 1;
    end
    if ~isempty(netCNN)
        LEDPoseClass(1,2) = 1;
        numTarget = netCNN{1,2};
        netCNN = netCNN{1,1};       
        inputSize = netCNN.Layers(1,1).InputSize(1,1:2);
    end
    if ~isempty(net3dCNN)
        LEDPoseClass(1,3) = 1;
    end

    try
        vr = VideoReader(vidFile);
    catch
        disp('Could not read video.')
        fps = [];
        diffPix = [];
        conLED = [];
        XY = [];
        return
    end
    
    fps = vr.FrameRate;
    vr.CurrentTime = vr.Duration-1;
    vr.CurrentTime = 0;
    numFrames = floor(fps*vr.Duration);
    frameDim = [vr.Height,vr.Width];

    % Detect if LED is visible
    if LEDPoseClass(1,1) == 1

        disp('Finding LED')
        
        hasLED = 0;
        checkTime = [0,round((numFrames/2)/fps),round((numFrames/fps)-15)];
        for ii = 1:3

            vr.CurrentTime = checkTime(1,ii);
            frames = uint8(zeros(frameDim(1,1),frameDim(1,2),3,round(fps)*15));
            for i3 = 1:round(fps)*15
                if hasFrame(vr)
                    frames(:,:,:,i3) = readFrame(vr);
                else
                    break
                end
            end                

            frames = imresize(frames,[100,100]);
            con = predict(netLED,frames);

            if ~all(con(:,2) < 0.8)
                hasLED = 1;
                break
            end
        end
        vr.CurrentTime = 0;

        if hasLED == 1
            disp('LED dectected')
        else
            LEDPoseClass(1,1) = 0;
            disp('LED not dectected, switching to linear sync')
        end

        if mean(con(:,2)) > 0.2
            disp('Net needs to be trained on this video')
        end
    end    

    if isempty(gcp('nocreate')) == 1
        parpool('local',3);
    else
        pool = gcp;
        if pool.NumWorkers ~= 3
            delete(gcp('nocreate'))
            parpool('local',3);
        end
    end

    parWB = parallel.pool.DataQueue;
    afterEach(parWB,@makeWB);
          
    % Scan video and get net confidence for each frame
    while 1
        try
            numRead = 50;
            spmd(3)
                mpiSettings('DeadlockDetection','off');
        
                if spmdIndex == 1 % Video reader core
                    
                    vr = VideoReader(vidFile); %#ok<TNMLP>
                    vr.CurrentTime = vr.Duration;  % trigger internal indexing?
                    vr.CurrentTime = 0;
                    
                    % Pre-allocate
                    preFrame = [];
                    countFrame = 0;
                    diffPix = zeros(numFrames,1);
                    frames = uint8(zeros(frameDim(1,1),frameDim(1,2),3,numRead));
                    
                    % Open spmd progress bar
                    send(parWB,[0,countFrame,numFrames])

                    % Main loop
                    tic
                    while hasFrame(vr)
                        
                        % Read batch of frames from video reader
                        ii = 0;
                        while ii < numRead && hasFrame(vr)
                            
                            % Read a frame from video reader and store
                            ii = ii+1;
                            frames(:,:,:,ii) = readFrame(vr);
                            countFrame = countFrame+1;
                            
                            % Calculate mean pixel difference between current
                            % and last frame, subtracting noise
                            if isempty(preFrame)
                                preFrame = double(rgb2gray(frames(:,:,:,ii)));
                            else
                                curFrame = double(rgb2gray(frames(:,:,:,ii)));
                                diffFrame = abs(curFrame-preFrame);
                                diffFrame(diffFrame < 10) = 0;
                                diffPix(countFrame,1) = mean(diffFrame,'all');
                                preFrame = curFrame;
                            end
                        end
                        
                        % Send to predictor core
                        spmdSend(frames(:,:,:,1:ii),2,1)

                        % Update spmd progress bar
                        send(parWB,[toc,countFrame,numFrames])
                    end

                    % Send stop signal to predictor core
                    spmdSend([],2,1)
                    
                    tic
                    fin = 0;
                    while toc < 60
                        if spmdProbe(3,3)
                            try
                                spmdReceive(3,3);
                            catch
                            end
                            fin = 1;
                        else
                            pause(0.01)
                        end
                    end
                    if fin == 0
                        spmdSend([],2,1)
                    end

                    % Close spmd progress bar
                    send(parWB,[])
                elseif spmdIndex == 2 % CNN predictor core
                    
                    % Pre-allocate
                    s = 1;
                    countFrame = 0;
                    con = zeros(numFrames,2);

                    % Main loop
                    while 1
                        try                            
                            % Receive batch of frames from reader core
                            while 1
                                if spmdProbe(1,1)
                                    frames = spmdReceive(1,1);
                                    break
                                else
                                    pause(0.01)
                                end
                            end

                            if ~isempty(frames)
                            
                                % Get confidence of netLED if enabled
                                if LEDPoseClass(1,1) == 1
                                    countFrame = countFrame+size(frames,4);
                                    con(s:countFrame,:) = predict(netLED,imresize(frames,[100,100]));
                                    s = s+size(frames,4);
                                end
                                
                                % Get confidence of netCNN if enabled
                                if LEDPoseClass(1,2) == 1
                                    yPose = predict(netCNN,imresize(frames,inputSize));
                                    spmdSend(yPose,3,2)
                                end
                            else
                                % Send stop signal to pose estimation core and
                                % exit main loop
                                spmdSend([],3,2)
                                break
                            end
                        catch
                            continue
                        end
                    end
                    
                elseif spmdIndex == 3 % Pose estimation core

                    % Create object detector to detect CNN parts
                    blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort',false,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',10);
                    
                    % Pre-allocate
                    s = 1;
                    countFrame = 0;
                    XY = zeros(numFrames,3,numTarget);

                    % Main loop
                    while 1
                        try
                            % Receive batch of frames from predictor core
                            while 1
                                if spmdProbe(2,2)
                                    yPose = spmdReceive(2,2);
                                    break
                                else
                                    pause(0.01)
                                end
                            end

                            if ~isempty(yPose)
                                
                                % Remove noise and set to uint8
                                yPose(yPose < 0.25) = 0;
                                yPose(yPose > 1) = 1;
                                yPose = uint8(yPose*255);
                                
                                % Pre-allocate
                                numY = size(yPose,4);
                                xy = zeros(numY,3,numTarget);

                                for ii = 1:numY % Loop through frames
                                    for i3 = 1:numTarget % Loop through targets
                                        
                                        % Make logical map from frame and
                                        % find area and centroids
                                        y = yPose(:,:,i3,ii);
                                        mask = logical(imresize(y,frameDim));
                                        [area,centroids] = blobAnalyser(mask);

                                        if isempty(centroids) == 0
                                            
                                            % Get target x and y coordinates
                                            [~,aMax] = max(area);
                                            xy(ii,1:2,i3) = round(centroids(aMax,:),1);
                                            
                                            % Get confidence by finding the
                                            % mode of values in the target
                                            % frame
                                            [counts,centers] = hist(double(y(y > 0))/255,linspace(0.01,1,20)); %#ok<HIST> 
                                            [~,idx] = max(counts);
                                            xy(ii,3,i3) = centers(idx);
                                        else
                                            % Target is not visible
                                            xy(ii,:,i3) = [nan,nan,0];
                                        end
                                    end
                                end
                                
                                % Add coordinate data to main data store
                                countFrame = countFrame+numY;
                                XY(s:countFrame,:,:) = xy;
                                s = s+numY;
                            else
                                break
                            end
                        catch
                            continue
                        end
                    end
                    spmdSend(1,1,3)
                end
            end
            reset(gpuDevice(1));
            break
        catch
            send(parWB,[])
            continue
        end
    end
    countFrame = countFrame{1};

    diffPix = diffPix{1};
    diffPix = smooth(diffPix(1:countFrame,1));
    
    if LEDPoseClass(1,1) == 1

        conLED = con{2};
        conLED = conLED(1:countFrame,2);

        idx = find(conLED > 0.1 & conLED < 0.8);
        idxD = diff(idx);
        
        lowCon = 0;
        if length(idxD) > 3 && mode(idxD) == 1
    
            count = 1;
            for ii = 3:length(idxD)
                if all(idxD(ii-2:ii) == 1)
                    lowCon(count,1) = idx(ii-1); %#ok<AGROW>
                    count = count+1;
                end
            end
        end
    
        if length(lowCon) >= 2
        
            h = [0,0];
            m = [0,0];
            s = [0,0];
            for ii = 1:2
    
                if ii == 1
                    f = lowCon(1,1);
                else
                    f = lowCon(end,1);
                end
                h(1,ii) = (f/fps)/3600;
                if h(1,ii) >= 1
                    m(1,ii) = (h(1,ii)-floor(h(1,ii)))*60;
                    h(1,ii) = floor(h(1,ii));
                else
                    m(1,ii) = h(1,ii)*60;
                    h(1,ii) = 0;
                end
                s(1,ii) = round((m(1,ii)-floor(m(1,ii)))*60);
                m(1,ii) = floor(m(1,ii));
            end               
            
            t1 = [num2str(h(1,1)),':',num2str(m(1,1)),':',num2str(s(1,1))];
            t2 = [num2str(h(1,2)),':',num2str(m(1,2)),':',num2str(s(1,2))];
            disp(['Area of low confidence found at: ',t1,' - ',t2])
            disp('Please add this to the net training data')
        end
    else
        conLED = 'N/A';
    end

    if LEDPoseClass(1,2) == 1
        XY = XY{3};
        XY = XY(1:countFrame,:);
    else
        XY = 'N/A';
    end

    function makeWB(input)
        global wb_global %#ok<GVMIS>
        
        if ~isempty(input)
    
            t = input(1,1);
            numTic = input(1,2);
            total = input(1,3);
            
            t2c = round((t/numTic)*(total-numTic));
            minRem = floor(t2c/60);
            secRem = round(((t2c/60)-minRem)*60);
    
            if isempty(wb_global) || ~ishandle(wb_global)
                 wb_global = waitbar(0/total,{'Scanning frames...';['Frame ',num2str(0),'/',num2str(total)];['~ ',num2str(Inf),' min ',num2str(Inf),' sec remaining']});
            else
                 waitbar(numTic/total,wb_global,{'Scanning frames...';['Frame ',num2str(numTic),'/',num2str(total)];['~ ',num2str(minRem),' min ',num2str(secRem),' sec remaining']});
            end
            drawnow
        else
            close(wb_global)
        end
    end
end

function tsV = findVideoPulse(conLED,fps)

    % Find all LED pulses
    Ck = 0;
    tsV = [0,0,0];
    iPulse = 1;
    for ii = 1:length(conLED)
        if conLED(ii,1) == 0
            if Ck == 1
                Ck = 0;
            end
            tsV(iPulse,3) = tsV(iPulse,3)+1;
        else
            if Ck == 0
                iPulse = iPulse+1;
                tsV(iPulse,:) = 0;
                tsV(iPulse,1) = ii;
                Ck = 1;
            end
            tsV(iPulse,2) = tsV(iPulse,2)+1;
        end
    end
    tsV = tsV/fps;
end

function [tsNsync,tsVsync] = syncPulse(tsN,tsV,digData,conLED,fps,numFrames)

    tsNsync = tsN;
    tsVsync = tsV;

    idx = find(tsVsync(2:end,2) < 0.1)+1;
    tsVsync(idx,:) = [];
    
    numDiff = size(tsNsync,1)-size(tsVsync,1);
    if numDiff ~= 0

        tsNsync(end+1:size(tsVsync,1),:) = nan;
        tsVsync(end+1:size(tsNsync,1),:) = nan;
        for ii = 1:abs(numDiff)

            tsNdiff_norm = [nan;diff(tsNsync(:,1))];
            tsVdiff_norm = [nan;diff(tsVsync(:,1))];
            tsNdiff_norm = tsNdiff_norm/mode(tsNdiff_norm(~isnan(tsNdiff_norm)));
            tsVdiff_norm = tsVdiff_norm/mode(tsVdiff_norm(~isnan(tsVdiff_norm)));
            tsDiff = tsNdiff_norm-tsVdiff_norm;
            
            if numDiff > 0
                idx = find(tsDiff < -0.5,1,'first');
                if isempty(idx)
                    [M,idx] = max(tsDiff);
                    if M < 1
                        idx = find(~isnan(tsNsync(:,1)),1,'last');
                    end
                end
                tsNsync(idx,:) = [];
                tsNsync(end+1,:) = nan; %#ok<AGROW>
            elseif numDiff < 0
                [M,idx] = max(tsDiff);
                if M > 0.5
                    tsVsync(idx,:) = [];
                    tsVsync(end+1,:) = nan; %#ok<AGROW>
                else
                    idx = find(~isnan(tsVsync(:,1)),1,'last');
                    tsVsync(idx,:) = nan;
                end
            end
        end
        rmv = min([sum(~isnan(tsNsync)),sum(~isnan(tsVsync))]);
        tsNsync = tsNsync(1:rmv,:);
        tsVsync = tsVsync(1:rmv,:);
    end

    neuralTS = interp1(tsVsync(:,1),tsNsync(:,1),(0:numFrames-1)/fps,'linear','extrap')';
    videoTS = ((0:numFrames-1)/fps)';

    neuralTS_ms = round(neuralTS*1000)+1;
    neuralTS_ms(isnan(neuralTS_ms)) = [];

    [~,idx] = min(abs(neuralTS_ms-length(digData)));

    neuralFrame = [digData(neuralTS_ms(1:idx-1,1),1),conLED(1:length(digData(neuralTS_ms(1:idx-1,1))))];
    neuralFrame(:,3) = neuralFrame(:,1)+neuralFrame(:,2);
    neuralFrame(neuralFrame(:,3) < 2,3) = 0;
    idx = find(diff(neuralFrame(:,3)) == 2)+1;
    idx = [1;idx];
    
    keepTS = zeros(length(idx),2);
    keepTS(:,1) = neuralTS(idx);
    keepTS(:,2) = videoTS(idx);
    
    keepN = false(size(tsN,1),1);
    keepV = false(size(tsV,1),1);
    for ii = 1:size(keepTS,1)
    
        [~,idxN] = min(abs(keepTS(ii,1)-tsN(:,1)));
        [~,idxV] = min(abs(keepTS(ii,2)-tsV(:,1)));
    
        if ~keepN(idxN,1) && ~keepV(idxV,1)
            keepN(idxN,1) = true;
            keepV(idxV,1) = true;
        end
    end
    tsNsync = tsN(keepN,:);
    tsVsync = tsV(keepV,:);
end

function velocity = XY2Vel(XY,ts,pixCm)

    numPoint = size(XY,1);
    numTarget = size(XY,3);

    if isempty(gcp('nocreate')) == 1
        parpool('local',numTarget);
    end

    mag = zeros(numPoint,numTarget);
    parfor ii = 1:numTarget
        
        % Remove NaN's by filling with mean of first and
        % last observation
        xy = XY(:,1:2,ii); %#ok<PFBNS>
        for i3 = 1:2

            nanStart = 0;
            firstValue = 0;
            for i4 = 1:length(xy(:,i3))
                if isnan(xy(i4,i3))
                    if nanStart == 0

                        nanStart = i4;
                        if i4 > 1
                            firstValue = xy(i4-1,i3);
                        else
                            firstValue = -1;
                        end
                    end
                else
                    if nanStart > 0

                        nanFinish = i4-1;
                        lastValue = xy(i4,i3);
                        if firstValue ~= -1    
                            xy(nanStart:nanFinish,i3) = (firstValue+lastValue)/2;
                        else
                            xy(nanStart:nanFinish,i3) = lastValue;
                        end
                        nanStart = 0;
                    end
                end
            end
            if nanStart > 0
                xy(nanStart:end,i3) = firstValue;
            end
            xy(:,i3) = smooth(xy(:,i3));
            
            % Calculate velocity for each target
            [~,~,~,vx,vy,~,~] = KalmanVel(xy(:,1),xy(:,2),ts,1);
            v = sqrt((vx.^2)+(vy.^2));
            if length(v) < numPoint
                v = [zeros(numPoint-length(v),1);v];
            end
            mag(:,ii) = v;
        end
    end
    mag = mean(mag,2);

    if ~isempty(pixCm)
        velocity = mag/(pixCm/10);
    else
        velocity = mag;
    end
end

function [t,x,y,vx,vy,ax,ay] = KalmanVel(posx,posy,post,order,Q,R)

    % Copyright (C) 2004 Sturla Molden
    % Centre for the Biology of Memory
    % Norwegian University of Science and Technology
    
    % Default values. These are manually tuned for good performance
    if (nargin < 5)
        Q = 0.1 * eye(2 + 2*order);
    end
    if (nargin < 6)
        R = eye(2);
    end
    % Chop off any missing data in the start of session
    n = length(posx);
    lastmissing = find(isfinite(posx), 1 )-1;
    posx = posx(lastmissing+1:end);
    posy = posy(lastmissing+1:end);
    post = post(lastmissing+1:end);
    
    % Run Kalman filter on the remaining samples
    missing = zeros(n,1);
    missing(isnan(posx)) = 1;
    
    if (sum(missing)) 
        % Missing data, use EM algorithm to get MAP estimators (Dempster et al. 1977)
        % Find the missing data we want to augment, initially augment 
        % with Kalman predicted positions, then augment with Kalman filtered
        % positions. Iterate ten times to allow the augmented data to converge.
        missing_index = find(missing);
        [x,y,vx,vy,ax,ay] = kfilter(posx,posy,post,order,Q,R,1,missing);
        for i = 1:10 
            posx(missing_index) = x(missing_index);
            posy(missing_index) = y(missing_index);
            [x,y,vx,vy,ax,ay] = kfilter(posx,posy,post,order,Q,R,0);
        end 
    else
        % No missing data, get the MAP estimates from a single pass
        [x,y,vx,vy,ax,ay,~] = kfilter(posx,posy,post,order,Q,R,0);
    end
    t = post;

    % This is the actual Kalman filter
    function [x,y,vx,vy,ax,ay,mm] = kfilter(posx,posy,post,order,Q,R,firstrun,missing)
    
        % allocate memory for return values
        n2 = length(posx);
        x = zeros(n2,1); 
        y = zeros(n2,1); 
        vx = zeros(n2,1); 
        vy = zeros(n2,1); 
        ax = zeros(n2,1); 
        ay = zeros(n2,1);
        mm = zeros(n2,1);
    
        % initialise return values and filtered state estimate
        x(1) = posx(1);
        y(1) = posy(1);
        vx(1) = 0;
        vy(1) = 0;
        ax(1) = 0;
        ay(1) = 0;
    
        switch (order)
            case {0}
                cX = [posx(1) posy(1)]';
            case {1}
                cX = [posx(1) posy(1) 0 0]';
            case {2}
                cX = [posx(1) posy(1) 0 0 0 0]';
        end
    
        cP = 0.1*eye(2 + 2*order);
        I = eye(2+2*order);
        outlier = 0;
    
        for k = 2:n2
    
            % compute A and H from the time lag 
            T = post(k) - post(k-1);
            switch(order)
                case {0}
                  A = [1 0; ... 
                       0 1];
                  H = [1 0; ...
                       0 1];
                case {1}
                  A = [1 0 T 0; ... 
                       0 1 0 T; ...
                       0 0 1 0; ... 
                       0 0 0 1];
                  H = [1 0 0 0; ... 
                       0 1 0 0]; ...     
                case {2}
                  A = [1 0 T 0 0.5*T*T 0; ...
                       0 1 0 T 0 0.5*T*T; ... 
                       0 0 1 0 T 0; ...
                       0 0 0 1 0 T; ...
                       0 0 0 0 1 0; ...
                       0 0 0 0 0 1;];        
                  H = [1 0 0 0 0 0; ... 
                       0 1 0 0 0 0];
            end
            if (firstrun && missing(k))
        
                % missing data, only predict
                % in the next EM steps they will be augmented with their MAP estimates
                pX = A * cX;
                pP = A * cP * A' + Q;  
                % the equations are obtained by setting the Kalman gain to zero
                cX = pX;  
                cP = pP;  
            else % data not missing, predict and correct
        
                % prediction step  
                pX = A * cX;
                pP = A * cP * A' + Q;
                % observation
                z = [posx(k); posy(k)];
                residual = z - H*pX;
                % validation gate for robustifying against extreme outliers
                chisq = [1000, 1000, 1000]; % just use a ridicilous threshold
                invS = (H*pP*H' + R)^-1; % measurement prediction covariance
                mahalanobis = residual' * invS * residual;
                mm(k) = mahalanobis;
                if ((outlier>5) || (mahalanobis < chisq(order+1)))
                    % within validation gate -- perform correction step with
                    % the new measurement (i.e. non-zero Kalman gain)
                    K = pP * H' * invS; 
                    cX = pX + K*residual; 
                    cP = (I - K*H)*pP;
                    outlier = 0;
                else
                    % outlier -- ignore this measurement   
                    % the equations are obtained by setting the Kalman gain to zero
                    cX = pX;  
                    cP = pP;  
                    outlier = outlier + 1;
                end       
            end
    
            % save Kalman filtered states (cX)
            x(k) = cX(1);
            y(k) = cX(2);
            if (order > 0)
                vx(k) = cX(3);
                vy(k) = cX(4);
            end
            if (order > 1)
                ax(k) = cX(5);
                ay(k) = cX(6);
            end
        end
    end
end