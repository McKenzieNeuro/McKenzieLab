
% Tracking targets for pose estimation
Targets{1,1} = 'RightEar';
Targets{2,1} = 'LeftEar';
Targets{3,1} = 'TailBase';
SymPair = [1,2]; % Targets that are symetrical pairs. Format as [p1,p2;p3,p4;...]

% Pose
execution = 'gpu'; % Execution environment: cpu,gpu, multi-gpu.
usePar = 1; % Use multi-core to label videos
numWorkers = 2; % Number of CPU cores to use for labeling. Total CPU cores needed is 3+numWorkers
pixCm = 13; % Number of pixels per cm

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

if strcmp(loadAction,'New') == 1n

    % make data dir folders
    eval(['mkdir ',ProjectPathName,'\Training']);
    eval(['mkdir ',ProjectPathName,'\Training\Data']);

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

    layersCNN = makeCNN(numTarget,inputSize);
    save([ProjectPathName,'\NetworkLayers.mat'],'layersCNN');

    promptMessage = sprintf('Would you like to label some training data?');
    titleBarCaption = 'settings';
    labelAction_CNN = questdlg(promptMessage,titleBarCaption,'Yes','Later','Yes');
    labelAction_LSTM = 'NA';
    updateAction = 'NA';

elseif strcmp(loadAction,'Load') == 1

    findDir = dir([ProjectPathName,'\Training']);
    findDir = {findDir.name};
    hasVidDir = 0;
    for i = 1:length(findDir)
        if strcmp(findDir{1,i},'Sequences') == 1
            hasVidDir = 1;
            break
        end
    end

    load('netCNN.mat')
    layersCNN = layerGraph(netCNN);

    promptMessage = sprintf('What would you like to do');
    titleBarCaption = 'settings';
    if hasVidDir == 0
        loadAction = questdlg(promptMessage,titleBarCaption,'Analyse video','Update pose network','Create behavior network','Analyse video');
    else
        loadAction = questdlg(promptMessage,titleBarCaption,'Analyse video','Update pose network','Update behavior network','Analyse video');
    end

    if strcmp(loadAction,'Update pose network') == 1 || strcmp(loadAction,'Update behavior network') == 1
        promptMessage = sprintf('What would you like to do');
        titleBarCaption = 'settings';
        updateAction = questdlg(promptMessage, titleBarCaption,'Add training data','Refine network','Retrain network','Add training data');
    else
        updateAction = 'NA';
    end

    load('ProjectData.mat')
end

%%

if strcmp(loadAction,'Analyse video') == 1

    allFiles = cell(0,2);
    promptMessage = sprintf('Load from file or GUI?');
    loadAction2 = questdlg(promptMessage,titleBarCaption,'File','GUI','File');

    if strcmp(loadAction2,'GUI')
        while 1

            [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select training video','MultiSelect','on');
            cd(VidFilePath)

            if length(VidFilePath) > 1
                if ~iscell(VidFileName)
                    allFiles = [allFiles;{VidFilePath,VidFileName}]; %#ok<AGROW>
                else
                    for i = 1:length(VidFileName)
                        allFiles = [allFiles;{VidFilePath,VidFileName{1,i}}]; %#ok<AGROW>
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
        [MatFileName,MatFilePath] = uigetfile('*.mat','Select video file');
        v=  load([MatFilePath filesep MatFileName]);
        allFiles = v.allFiles;
    end

    if ~isempty(allFiles)
        if usePar == 0
            analyzeVideo(ProjectPathName,allFiles,pixCm)
        elseif usePar == 1
            analyzeVideo_Par(ProjectPathName,allFiles,pixCm,numWorkers)
        end
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

end

%%

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
    %convolution2dLayer(1,numUnits/2,'Padding','same',Name="conv_FC_Output_1.0")
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

% Count = depth;
% for i = 1:depth+1
%     if i == 1
%         layersCNN = connectLayers(layersCNN,'conv_FC_Bridge',"Fuse_"+i+"/in2");
%     elseif i > 1 && i < depth+1
%         layersCNN = connectLayers(layersCNN,"reLu_BN_Out_"+i+"_"+branchDepth,"Fuse_"+i+"/in2");
%     elseif i == depth+1
%         layersCNN = connectLayers(layersCNN,"zero2one","Fuse_"+i+"/in2");
%     end
%     Count = Count-1;
% end

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

numVal = rem(numObs,10)+50;

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
        %c = categorical(1:7);
        %Out{1,2} = c(Out{1,2});
    end
end

function analyzeVideo(ProjectPathName,allFiles,pixCm)

PoseClass = [0,0];
try
    load([ProjectPathName,'\netCNN.mat']); %#ok<LOAD>
    inputSize = netCNN.Layers(1,1).InputSize(1,1:2); %#ok<NODEF>
    numTarget = length(Targets);

    % Create object detector to detect CNN parts
    blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort',false,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',10);

    PoseClass(1,1) = 1;
    try
        load([ProjectPathName,'\net3dCNN.mat']); %#ok<LOAD>
        numFrames = net3dCNN.Layers(1,1).InputSize(1,3); %#ok<NODEF>
        numClass = length(Classes);
        PoseClass(1,2) = 1;
    catch
        net3dCNN = [];
        numClass = 0;
    end
catch
    netCNN = [];
end

numFile = size(allFiles,1);
if sum(PoseClass) > 0
    for i = 1:numFile

        disp(['Loading File ',num2str(i),'/',num2str(numFile)])
        outfil = [allFiles{i,1},allFiles{i,2}(1,1:end-4),'_poseData.mat'];
        if ~exist(outfil)

            % Create video reader
            vr = VideoReader([allFiles{i,1},allFiles{i,2}]); %#ok<TNMLP>
            numFrame = vr.NumFrames;
            fps = vr.FrameRate;
            vr.CurrentTime = 0;
            frameDim = [vr.Height,vr.Width];

            poseData = matfile([allFiles{i,1},allFiles{i,2}(1,1:end-4),'_poseData.mat'],'Writable',true);

            if PoseClass(1,1) == 1
                poseData.targets = Targets;
                poseData.headers = [{'X'},{'Y'},{'Confidence'},{'Distance'}];
                poseData.XY = zeros([],4,numTarget);
                poseData.diffPix = zeros([],1);
            end

            if PoseClass(1,2) == 1
                poseData.classes = Classes;
                poseData.confidence = zeros([],length(Classes));
            else
                numFrames = fps;
            end

            data1Sec = uint8(zeros([inputSize,numFrames,2]));
            lastPoint = nan(1,2,numTarget);
            frames = uint8(zeros([inputSize,3,100]));
            frameLast = uint8(ones([inputSize,3])*127);
            data1Frame = uint8(zeros([inputSize,1,2]));

            iFrame = 0;
            wb = waitbar(iFrame/numFrame,{'Reading frames...';['Frame ',num2str(0),'/',num2str(numFrame)]});
            while hasFrame(vr)

                waitbar(iFrame/numFrame,wb,{'Reading frames...';['Frame ',num2str(iFrame),'/',num2str(numFrame)]});

                for ii = 1:100
                    if hasFrame(vr)
                        frames(:,:,:,ii) = imresize(readFrame(vr),inputSize);
                    else
                        break
                    end
                end
                pulledFrames = ii;

                waitbar(iFrame/numFrame,wb,{'Detecting targets...';['Frame ',num2str(iFrame),'/',num2str(numFrame)]});

                if PoseClass(1,1) == 1

                    yPose = single(zeros([inputSize,numTarget,pulledFrames]));
                    for ii = 1:25:pulledFrames
                        try
                            yPose(:,:,:,ii:ii+24) = predict(netCNN,frames(:,:,:,ii:ii+24),'ExecutionEnvironment','auto','Acceleration','auto');
                        catch
                            yPose(:,:,:,ii:pulledFrames) = predict(netCNN,frames(:,:,:,ii:pulledFrames),'ExecutionEnvironment','auto','Acceleration','auto');
                        end
                    end

                    yPose(yPose < 0.25) = 0;
                    yPose(yPose > 1) = 1;
                    yPose = uint8(yPose*255);
                end

                diffPix = zeros(pulledFrames,1);
                XY = zeros(pulledFrames,4,numTarget);
                yBehave = single(zeros(pulledFrames,numClass));
                for ii = 1:pulledFrames

                    data1Sec = circshift(data1Sec,-1,3);
                    data1Frame(:,:,1,1) = rgb2gray(frames(:,:,:,ii));
                    data1Frame(:,:,1,2) = abs(data1Frame(:,:,1,1)-rgb2gray(frameLast));
                    diffPix(ii,1) = mean(data1Frame(:,:,1,2),'all');
                    frameLast = frames(:,:,:,ii);
                    data1Sec(:,:,end,:) = data1Frame;

                    if PoseClass(1,2) == 1
                        yBehave(ii,:) = predict(net3dCNN,data1Sec,'ExecutionEnvironment','auto','Acceleration','auto');
                    end

                    if PoseClass(1,1) == 1
                        for i3 = 1:numTarget

                            y = yPose(:,:,i3,ii);
                            mask = logical(imresize(y,frameDim));
                            [area,centroids] = blobAnalyser(mask);
                            if isempty(centroids) == 0

                                [~,aMax] = max(area);
                                xy = round(centroids(aMax,:),1);

                                [counts,centers] = hist(double(y(y > 0))/255,linspace(0.01,1,20)); %#ok<HIST>
                                [~,idx] = max(counts);
                                Con = centers(idx);
                            else
                                xy = [nan,nan];
                                Con = 0;
                            end
                            XY(ii,1:3,i3) = [xy(1,1),xy(1,2),Con];

                            if isnan(lastPoint(1,1,i3))
                                XY(ii,4,i3) = 0;
                            else
                                if isnan(XY(ii,1,i3))
                                    XY(ii,4,i3) = 0;
                                else
                                    XY(ii,4,i3) = pdist([lastPoint(1,1:2,i3);XY(ii,1:2,i3)],'euclidean');
                                end
                            end
                            lastPoint(1,1:2,i3) = XY(ii,1:2,i3);
                        end
                    end
                    waitbar(iFrame/numFrame,wb,{'Classifying data...';['Frame ',num2str(iFrame),'/',num2str(numFrame)]});
                    iFrame = iFrame+1;
                end

                if PoseClass(1,1) == 1
                    poseData.XY = cat(1,poseData.XY,XY);
                    poseData.diffPix = [poseData.diffPix;diffPix];
                end

                if PoseClass(1,2) == 1
                    poseData.confidence = [poseData.confidence;yBehave];
                end
            end

            if ~isnan(pixCm)
                vel = zeros(size(XY,1),1);
                for ii = 1:size(XY,1)
                    vel(ii,1) = (mean(XY(ii,4,:),'omitnan')/(1/fps))/(pixCm/10);
                end
                poseData.velocity = vel;
            end

            close(wb);
        end
    end
else
    disp('No networks found')
end
end

function analyzeVideo_Par(ProjectPathName,allFiles,pixCm,numWorkers)

PoseClass = [0,0];
try
    netCNN = load([ProjectPathName,'\netCNN.mat']);
    Targets = netCNN.Targets;
    netCNN = netCNN.netCNN;
    inputSize = netCNN.Layers(1,1).InputSize(1,1:2);
    numTarget = length(Targets);

    % Create object detector to detect CNN parts
    blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort',false,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',10);

    PoseClass(1,1) = 1;
    try
        load([ProjectPathName,'\net3dCNN.mat']); %#ok<LOAD>
        numFrames = net3dCNN.Layers(1,1).InputSize(1,3); %#ok<NODEF>
        numClass = length(Classes);
        PoseClass(1,2) = 1;
    catch
        net3dCNN = [];
        numClass = 0;
    end
catch
    netCNN = [];
end

numFile = size(allFiles,1);
if sum(PoseClass) > 0

    % CPU core roles
    reader = 1;
    pose = 2;
    classify = 3;
    labeler = 4:4+numWorkers-1;

    if isempty(gcp('nocreate')) == 1
        parpool('local',labeler(1,end));
    end

    for i = 1:numFile

        disp(['Loading File ',num2str(i),'/',num2str(numFile)])
        outfil = [allFiles{i,1},allFiles{i,2}(1,1:end-4),'_poseData.mat'];
        if ~exist(outfil)
            % Create video reader
            vr = VideoReader([allFiles{i,1} filesep allFiles{i,2}]); %#ok<TNMLP>
            numFrame = vr.NumFrames;
            fps = vr.FrameRate;
            vr.CurrentTime = 0;
            frameDim = [vr.Height,vr.Width];

            poseData = matfile([allFiles{i,1},allFiles{i,2}(1,1:end-4),'_poseData.mat'],'Writable',true);

            if PoseClass(1,1) == 1
                poseData.targets = Targets;
                poseData.headers = [{'X'},{'Y'},{'Confidence'},{'Distance'}];
                poseData.XY = zeros([],4,numTarget);
                poseData.diffPix = zeros([],1);
            end

            if PoseClass(1,2) == 1
                poseData.classes = Classes;
                poseData.confidence = zeros([],length(Classes));
            else
                numFrames = round(fps);
            end

            dataOut = parallel.pool.DataQueue;
            afterEach(dataOut,@makeWB);

            spmd(labeler(1,end))
                if spmdIndex == reader

                    vr.CurrentTime = 0;

                    qCell = cell(2,10);
                    qCell(1,:) = {uint8(zeros([inputSize,3,100]))};
                    qCell(2,:) = {0};

                    iFrame = 0;
                    doneRead = 0;
                    send(dataOut,[iFrame,numFrame])
                    while 1

                        if qCell{2,end} < 100 && doneRead == 0
                            if hasFrame(vr)

                                iFrame = iFrame+1;
                                qCell{2,end} = qCell{2,end}+1;
                                qCell{1,end}(:,:,:,qCell{2,end}) = imresize(readFrame(vr),inputSize);

                                if qCell{2,end} == 100 && qCell{2,1} == 0
                                    qCell = circshift(qCell,[2,-1]);
                                end
                            else
                                doneRead = 1;
                            end
                        end

                        if spmdProbe(pose)

                            for qNum = 1:10
                                if qCell{2,qNum} > 0
                                    break
                                end
                            end

                            if qCell{2,qNum} > 0
                                spmdReceive(pose);
                                spmdSend(qCell{1,qNum}(:,:,:,1:qCell{2,qNum}),[pose,classify]);
                                qCell{2,qNum} = 0;
                            else
                                if doneRead == 1
                                    break
                                end
                            end
                        end
                    end
                    spmdReceive(pose);
                    spmdSend([],[pose,classify]);

                elseif spmdIndex == pose

                    iRead = 0;
                    coreSend = 1;

                    predict(netCNN,uint8(zeros([inputSize,3,25])));

                    while 1

                        spmdSend(1,reader)

                        try
                            frames = spmdReceive(reader);
                        catch
                            continue
                        end

                        if isempty(frames)
                            break
                        else
                            iRead = iRead+1;
                            pulledFrames = size(frames,4);
                        end

                        if PoseClass(1,1) == 1

                            yPose = single(zeros([inputSize,numTarget,pulledFrames]));
                            for ii = 1:25:pulledFrames
                                try
                                    yPose(:,:,:,ii:ii+24) = predict(netCNN,frames(:,:,:,ii:ii+24),'ExecutionEnvironment','auto','Acceleration','auto');
                                catch
                                    yPose(:,:,:,ii:pulledFrames) = predict(netCNN,frames(:,:,:,ii:pulledFrames),'ExecutionEnvironment','auto','Acceleration','auto');
                                end
                            end

                            yPose(yPose < 0.25) = 0;
                            yPose(yPose > 1) = 1;
                            yPose = uint8(yPose*255);

                            spmdSend({yPose,iRead},labeler(1,coreSend))
                            coreSend = coreSend+1;
                            if coreSend > numWorkers
                                coreSend = 1;
                            end
                        end
                    end
                    spmdSend([],labeler)

                elseif spmdIndex == classify

                    iFrame = 0;
                    data1Sec = uint8(zeros([inputSize,numFrames,2]));
                    frameLast = uint8(ones([inputSize,3])*127);
                    data1Frame = uint8(zeros([inputSize,1,2]));
                    diffPix = zeros(numFrame,1);
                    yBehave = single(zeros(numFrame,numClass));

                    while 1

                        try
                            frames = spmdReceive(reader);
                        catch
                            continue
                        end

                        if isempty(frames)
                            break
                        else
                            pulledFrames = size(frames,4);
                        end

                        for ii = 1:pulledFrames

                            iFrame = iFrame+1;
                            data1Sec = circshift(data1Sec,-1,3);
                            data1Frame(:,:,1,1) = rgb2gray(frames(:,:,:,ii));
                            data1Frame(:,:,1,2) = abs(data1Frame(:,:,1,1)-rgb2gray(frameLast));
                            diffPix(iFrame,1) = mean(data1Frame(:,:,1,2),'all');
                            frameLast = frames(:,:,:,ii);
                            data1Sec(:,:,end,:) = data1Frame;

                            if PoseClass(1,2) == 1
                                yBehave(iFrame,:) = predict(net3dCNN,data1Sec,'ExecutionEnvironment','auto','Acceleration','auto');
                            end

                            send(dataOut,[iFrame,numFrame])
                        end
                    end
                elseif any(spmdIndex == labeler)

                    iRead = 0;
                    xyCell = cell(ceil((numFrame/100)/numWorkers),2);

                    blobAnalyser(false(frameDim));

                    while 1

                        try
                            data = spmdReceive(pose);
                        catch
                            continue
                        end

                        if isempty(data)
                            break
                        else
                            iRead = iRead+1;
                            yPose = data{1,1};
                            xyCell{iRead,2} = data{1,2};
                            pulledFrames = size(yPose,4);
                            xy = zeros(pulledFrames,3,numTarget);
                        end

                        for ii = 1:pulledFrames
                            if PoseClass(1,1) == 1
                                for i3 = 1:numTarget

                                    y = yPose(:,:,i3,ii);
                                    mask = logical(imresize(y,frameDim));
                                    [area,centroids] = blobAnalyser(mask);
                                    if isempty(centroids) == 0

                                        [~,aMax] = max(area);
                                        xy(ii,1:2,i3) = round(centroids(aMax,:),1);

                                        [counts,centers] = hist(double(y(y > 0))/255,linspace(0.01,1,20)); %#ok<HIST>
                                        [~,idx] = max(counts);
                                        xy(ii,3,i3) = centers(idx);
                                    else
                                        xy(ii,1:3,i3) = [nan,nan,0];
                                    end
                                end
                            end
                        end
                        xyCell{iRead,1} = xy;
                    end
                end
            end

            diffPix = diffPix{3};

            xyCore = cell(ceil((numFrame/100)/numWorkers),numWorkers);
            for ii = 1:numWorkers
                l = length(xyCell{labeler(1,ii)}(:,1));
                xyCore(1:l,ii) = xyCell{labeler(1,ii)}(:,1);
            end

            XY = zeros(numFrame,4,numTarget);
            iRow = 1;
            for ii = 1:size(xyCore,1)
                for i3 = 1:size(xyCore,2)
                    temp = xyCore{ii,i3};
                    l = size(temp,1);
                    XY(iRow:iRow+l-1,1:3,:) = temp;
                    iRow = iRow+l;
                end
            end

            for ii = 2:size(XY,1)
                for i3 = 1:numTarget
                    if isnan(XY(ii-1,1,i3)) || isnan(XY(ii,1,i3))
                        XY(ii,4,i3) = 0;
                    else
                        XY(ii,4,i3) = pdist([XY(ii-1,1:2,i3);XY(ii,1:2,i3)],'euclidean');
                    end
                end
            end

            if PoseClass(1,1) == 1
                poseData.diffPix = diffPix;
                poseData.XY = XY;
            end

            if PoseClass(1,2) == 1
                poseData.confidence = yBehave;
            end

            if ~isnan(pixCm)
                vel = zeros(size(XY,1),1);
                for ii = 1:size(XY,1)
                    vel(ii,1) = (mean(XY(ii,4,:),'omitnan')/(1/fps))/(pixCm/10);
                end
                poseData.velocity = vel;
            end

            send(dataOut,[])
        end
    end
else
    disp('No networks found')
end

reset(gpuDevice(1));

    function makeWB(input)
        global wb %#ok<GVMIS>

        if ~isempty(input)
            if isempty(wb) || ~ishandle(wb)
                wb = waitbar(0/input(1,2),{'Reading and processing frames...';['Frame ',num2str(0),'/',num2str(input(1,2))]});
            else
                waitbar(input(1,1)/input(1,2),wb,{'Reading and processing frames...';['Frame ',num2str(input(1,1)),'/',num2str(input(1,2))]});
            end
            drawnow
        else
            close(wb)
        end
    end
end
