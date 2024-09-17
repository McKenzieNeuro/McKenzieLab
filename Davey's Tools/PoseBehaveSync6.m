% Analysis
Sync = 1; % Enable video syncing to LED pulse. 0 = off, 1 = on.
Pose = 1; % Enable pose estimation. 0 = off, 1 = on.
Behavior = 0; %Enable behavior classification, Pose must be on. 0 = off, 1 = on.
 
% Pose
execution = 'gpu'; % Execution environment: cpu,gpu, multi-gpu.
labelVidPose = 1; % Create labeled video, will still save coordinate data if off
xSmooth = 1; % Use smoothing to reduce target jitter. 0 = off, 1 = on.
numWorkers = 32; % Number of CPU cores to use for labeling.

% Tracking targets for pose estimation
Targets{1,1} = 'RightEar';
Targets{2,1} = 'LeftEar';
Targets{3,1} = 'TailBase';
SymPair = [1,2]; % Targets that are symetrical pairs. Format as [p1,p2;p3,p4;...]

% Image input size for CNN, larger gives better results but uses more memory. 
% Must be divisible by 16.
inputSize = [256,256]; %320

% Training Data
numTrainingFrames = 60; % Number of frames to pull from selected video.
ClusterPool = 5000; % Number of frames from video to use.
numClust = 10; % Number of clustersxzfsdef to sort frames.
BaseSize = [480,640]; % Display size for frames in UI.

% Net training
numIter = 50000; % Number of training iterations.
miniBatchSizeCNN = 20; % Number of images in training batches, more gives better results but uses more memory. Best as a factor of 100.
miniBatchSizeLSTM = 100;
actLayer = 'reLuRegBN_6.1';
actLayer_numFilt = 512;

% Behavior
labelVidBehav = 1; % Create labeled video, will still save behavioral data if off
framesPerObs = 25;

% Behavioral classes
Classes{1,1} = 'Sleep';
Classes{2,1} = 'Awake';
Classes{3,1} = 'Seizure';

% Video Sync
LEDtime = 10;
LEDchannel = 2;
LEDloc = 'bottom-right';
Sfreq = 20000;

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

if strcmp(loadAction,'New') == 1

    % make data dir folders
    eval(['mkdir ',ProjectPathName,'\Training']);
    eval(['mkdir ',ProjectPathName,'\Training\Data']);

    numTarget = length(Targets);
    numClass = length(Classes);
    totalFrames = 0;
    save([ProjectPathName,'\ProjectData.mat'],'inputSize','Targets');
    
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

    imFiles = dir(fullfile([ProjectPathName,'\Training\Data'],'*.mat'));
    totalFrames = size(imFiles,1);
    clear imFiles
    numTarget = length(Targets);

    if strcmp(loadAction,'Update pose network') == 1

        if strcmp(updateAction,'Retrain network') == 0
            labelAction_CNN = 'Yes';
            trainAction_CNN = 'NA'; 
        else
            labelAction_CNN = 'NA';
            trainAction_CNN = 'Yes';
        end        
        labelAction_LSTM = 'NA';

    elseif strcmp(loadAction,'Create behavior network') == 1

        numClass = length(Classes);

        % make data dir folders
        eval(['mkdir ',ProjectPathName,'\Training\Sequences']);
        for i = 1:numClass
            eval(['mkdir ',ProjectPathName,'\Training\Sequences\',Classes{i,1}]);
        end

        sequenceData = cell(2,4);
        sequenceData{1,1} = 'Totals';
        sequenceData{1,2} = 'Files';
        sequenceData{1,3} = 'Sequences';
        sequenceData{1,4} = 'Labels';
        sequenceData{2,1} = Classes;
        sequenceData{2,1}(:,2:3) = {0};
        save([ProjectPathName,'\Training\SequenceData.mat'],'sequenceData');

        promptMessage = sprintf('Would you like to generate some training data?');
        titleBarCaption = 'settings';
        labelAction_LSTM = questdlg(promptMessage,titleBarCaption,'Yes','Later','Yes');

        labelAction_CNN = 'NA';
        trainAction_CNN = 'NA'; 
        
    elseif strcmp(loadAction,'Update behavior network') == 1
        
        if strcmp(updateAction,'Retrain network') == 0
           labelAction_LSTM = 'Yes';
           reactAction = 'No';
        else
            promptMessage = sprintf('Would you like to update CNN activations?');
            titleBarCaption = 'settings';
            reactAction = questdlg(promptMessage,titleBarCaption,'Yes','No','Yes');

            if strcmp(reactAction,'Yes') == 1
                labelAction_LSTM = 'Yes';
                iFile = 1;
            else
                labelAction_LSTM = 'NA';
            end            
            trainAction_LSTM = 'Yes'; 
        end
        labelAction_CNN = 'NA';
        trainAction_CNN = 'NA'; 
    elseif strcmp(loadAction,'Analyse video') == 1

        labelAction_CNN = 'NA';
        trainAction_CNN = 'NA'; 
        labelAction_LSTM = 'NA';
        trainAction_LSTM = 'NA'; 
    end
end

if strcmp(loadAction,'New') == 1  || strcmp(updateAction,'Retrain network') == 1

    % CNN network architecture, the brains of the operation
    numUnits = 30;
    depth = 5;
    branchDepth = 3;
    numTarget = length(Targets);
    numClass = length(Classes);

    PartDetectorLayers = imageInputLayer([inputSize,3],'Normalization','none','Name','Input');
    outputName = PartDetectorLayers.Name;
    PartDetectorLayers = layerGraph(PartDetectorLayers);   
            
    for i = 1:depth
        
        if i == 1
            fistCon = 7;
        elseif i == 2
            firstCon = 5;
        else
            firstCon = 3;
        end
    
        hiddenLayers = [
            convolution2dLayer(fistCon,numUnits,'Padding','same',Name="conv_"+i+".0")
            batchNormalizationLayer(Name="BN_"+i+".0")
            reluLayer(Name="reLu_"+i+".0")];

        PartDetectorLayers = addLayers(PartDetectorLayers,hiddenLayers);
        PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"conv_"+i+".0");
        outputName = "reLu_"+i+".0";

        for ii = 1:branchDepth

            branchLayers = [
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_"+i+"_"+ii+".0")
                batchNormalizationLayer(Name="BN_BN"+i+"_"+ii+".0")
                reluLayer(Name="reLu_BN_"+i+"_"+ii+".0")
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_"+i+"_"+ii+".1")
                batchNormalizationLayer(Name="BN_BN"+i+"_"+ii+".1")
                reluLayer(Name="reLu_BN_"+i+"_"+ii+".1")
                additionLayer(2,Name="add_BN_"+i+"_"+ii+".0")];

            PartDetectorLayers = addLayers(PartDetectorLayers,branchLayers);
            PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"conv_BN_"+i+"_"+ii+".0");
            PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"add_BN_"+i+"_"+ii+".0/in2");
            outputName = "add_BN_"+i+"_"+ii+".0";
        end
            
        hiddenLayers = [
            convolution2dLayer(3,numUnits,'Padding','same',Name="conv_"+i+".1")
            batchNormalizationLayer(Name="BN_"+i+".1")
            reluLayer(Name="reLu_"+i+".1")
            maxPooling2dLayer(2,'Stride',2,Name="maxpool_"+i+".1")];

        PartDetectorLayers = addLayers(PartDetectorLayers,hiddenLayers);
        PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"conv_"+i+".1");
        outputName = "maxpool_"+i+".1";         
        numUnits = numUnits*2;
    end
    
    bridgeLayers = [
        convolution2dLayer(3,numUnits,'Padding','same','Name','conv_Bridge_1.0')
        batchNormalizationLayer('Name','BN_Bridge_1.0')
        reluLayer('Name','reLu_Bridge_1.0')];

    PartDetectorLayers = addLayers(PartDetectorLayers,bridgeLayers);
    PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,'conv_Bridge_1.0');
    outputName = 'reLu_Bridge_1.0';

        for ii = 1:branchDepth

            branchLayers = [
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_Bridge_"+ii+".0")
                batchNormalizationLayer(Name="BN_BN_Bridge_"+ii+".0")
                reluLayer(Name="reLu_BN_Bridge_"+ii+".0")
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_Bridge_"+ii+".1")
                batchNormalizationLayer(Name="BN_BN_Bridge_"+ii+".1")
                reluLayer(Name="reLu_BN_Bridge_"+ii+".1")
                additionLayer(2,Name="add_BN_Bridge_"+ii+".0")];

            PartDetectorLayers = addLayers(PartDetectorLayers,branchLayers);
            PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"conv_BN_Bridge_"+ii+".0");
            PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"add_BN_Bridge_"+ii+".0/in2");
            outputName = "add_BN_Bridge_"+ii+".0";
        end

    bridgeLayers = [
        convolution2dLayer(3,numUnits,'Padding','same','Name','conv_Bridge_1.1')
        batchNormalizationLayer('Name','BN_Bridge_1.1')
        reluLayer('Name','reLu_Bridge_1.1')];
    
    PartDetectorLayers = addLayers(PartDetectorLayers,bridgeLayers);
    PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,'conv_Bridge_1.1');
    outputName = 'reLu_Bridge_1.1';
    
    Count = 1;
    for i = depth:-1:1
    
        numUnits = numUnits/2;
        hiddenLayers = [
            resize2dLayer('Scale',2,Name="resize_"+i)
            convolution2dLayer(3,numUnits,'Padding','same',Name="conv_up_"+i+".0")
            additionLayer(2,Name="Fuse_"+Count)
            batchNormalizationLayer(Name="BN_up_"+i+".0")
            reluLayer(Name="reLu_up_"+i+".0")];

        PartDetectorLayers = addLayers(PartDetectorLayers,hiddenLayers);
        PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"resize_"+i);
        outputName = "reLu_up_"+i+".0";
        
        for ii = 1:branchDepth

            branchLayers = [
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_up_"+i+"_"+ii+".0")
                batchNormalizationLayer(Name="BN_BN_up_"+i+"_"+ii+".0")
                reluLayer(Name="reLu_BN_up_"+i+"_"+ii+".0")
                convolution2dLayer(3,numUnits,'Padding','same',Name="conv_BN_up_"+i+"_"+ii+".1")
                batchNormalizationLayer(Name="BN_BN_up_"+i+"_"+ii+".1")
                reluLayer(Name="reLu_BN_up_"+i+"_"+ii+".1")
                additionLayer(2,Name="add_BN_up_"+i+"_"+ii+".0")];

            PartDetectorLayers = addLayers(PartDetectorLayers,branchLayers);
            PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"conv_BN_up_"+i+"_"+ii+".0");
            PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"add_BN_up_"+i+"_"+ii+".0/in2");
            outputName = "add_BN_up_"+i+"_"+ii+".0";
        end

        hiddenLayers = [
            convolution2dLayer(3,numUnits,'Padding','same',Name="conv_up_"+i+".1")
            batchNormalizationLayer(Name="BN_up_"+i+".1")
            reluLayer(Name="reLu_up_"+i+".1")];
    
        PartDetectorLayers = addLayers(PartDetectorLayers,hiddenLayers);
        PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,"conv_up_"+i+".1");
        outputName = "reLu_up_"+i+".1";
        Count = Count+1;
    end
    
    Count = depth;
    for i = 1:depth
        PartDetectorLayers = connectLayers(PartDetectorLayers,"reLu_"+i+".1","Fuse_"+Count+"/in2");
        Count = Count-1;
    end
    
    outputLayers = [
        convolution2dLayer(1,numTarget,'Padding','same','Name','convFC_1.0')
        depthConcatenationLayer(2,'Name','catInput')
        convolution2dLayer(1,numTarget,'Padding','same','Name','convFC_1.1')
        regressionLayer('Name','RegOut')];
    PartDetectorLayers = addLayers(PartDetectorLayers,outputLayers);
    PartDetectorLayers = connectLayers(PartDetectorLayers,outputName,'convFC_1.0');
    PartDetectorLayers = connectLayers(PartDetectorLayers,'Input','catInput/in2'); 
    
    RegressionLayers = [   
        depthConcatenationLayer(2,'Name','Infuse')
        
        convolution2dLayer(5,numUnits,'Padding','same','Name','convRegBN_1.0')
        reluLayer('Name','reLuRegBN_1.0')
        convolution2dLayer(3,numUnits,'Padding','same','Name','convRegBN_1.1')
        reluLayer('Name','reLuRegBN_1.1')
        convolution2dLayer(3,numUnits,'Padding','same','Name','convRegBN_1.2')
        reluLayer('Name','reLuRegBN_1.2')
        
        maxPooling2dLayer(2,'Stride',2,'Name','maxpoolReg_1')
        convolution2dLayer(3,numUnits*2,'Padding','same','Name','convRegBN_2.0')
        reluLayer('Name','reLuRegBN_2.0')
        convolution2dLayer(3,numUnits*2,'Padding','same','Name','convRegBN_2.1')
        reluLayer('Name','reLuRegBN_2.1')
            
        resize2dLayer('Scale',2,'Name','resizeReg_4')
        convolution2dLayer(3,numUnits,'Padding','same','Name','convRegUp_3.0')
        reluLayer('Name','reLuRegBN_3.0')
        convolution2dLayer(3,numUnits,'Padding','same','Name','convRegUp_3.1')
        reluLayer('Name','reLuRegBN_3.1')
        additionLayer(2,'Name','addReg_1.0')
        
        convolution2dLayer(1,numTarget,'Stride',1,'Name','conv2dR_final')
        reluLayer('Name','reLuRegBN_Out')
        additionLayer(2,'Name','FuseReg_1.0')
        regressionLayer('Name','RegOut')];
    
    % LSTM network
    numFeatures = (6*actLayer_numFilt)+1;        
    LSTMlayers = [
        sequenceInputLayer(numFeatures,'Name','sequence')
        gruLayer(2500,'OutputMode','sequence','Name','bilstm1')
        dropoutLayer(0.4,'Name','drop1')
        gruLayer(1250,'OutputMode','last','Name','bilstm2')
        dropoutLayer(0.4,'Name','drop1')
        fullyConnectedLayer(200,'Name','fc_1')
        fullyConnectedLayer(numClass,'Name','fc_final')
        softmaxLayer('Name','softmax')
        classificationLayer('Name','classification')]; 
    
    save([ProjectPathName,'\NetworkLayers.mat'],'PartDetectorLayers','RegressionLayers','LSTMlayers');
end
%%
if strcmp(labelAction_CNN,'Yes') == 1 
    while 1        
        if exist('VidFilePath','var') == 1 && ischar(VidFilePath) == 1
            cd(VidFilePath)
        end
    
        [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select training videos','MultiSelect','off');
        if ischar(VidFilePath) == 0
            break
        end
        cd(VidFilePath)   
            
        if strcmp(updateAction,'Refine network') == 1    
            load([ProjectPathName,'\netCNN.mat']);
            blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort',false,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',10);
        end   
    
        FullName = {fullfile(VidFilePath,VidFileName)};
        cd(VidFilePath)
        
        % Create video reader
        disp('Reading frames.')
        vr = VideoReader(FullName{1,1}); %#ok<TNMLP> 
        H = vr.Height;
        W = vr.Width;
        C = 3;
        fps = vr.FrameRate;
        numFrames = vr.NumFrames; 
        
        % Get frames from selected video
        if strcmp(updateAction,'Refine network') == 0

            idx = round(linspace(1,numFrames-1,ClusterPool));    

            % Pull frames from video and cluster        
            frameTrain = uint8(zeros(H,W,3,ClusterPool));
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
                    frameTrain(:,:,:,i) = frame;
                    GS = rgb2gray(frame);
                    GS = imresize(GS,[50,50]);
                    frameVect(:,i) = double(GS(:));
                end
            end
            frameTrain = frameTrain(:,:,:,1:i);
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
            frameTrain = frameTrain(:,:,:,idx);

        elseif strcmp(updateAction,'Refine network') == 1

            definput = {'0','0','1','5'};
            frameTrain = uint8(zeros(H,W,3,0));
            Count = 1;
            figure
            set(gcf, 'Position', get(0, 'Screensize'));
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
                frameR = (double(imresize(frame,inputSize))/255);
                reg = predict(netCNN,frameR);
                
                xy = zeros(fps*L,2,numTarget);
                Con = zeros(numTarget,(fps*L));
                for i = 1:size(reg,4)
                    
                    clf
                    subplot(2,1,1)
                    title('CNN Confidence')
                    hold on
                    for ii = 1:numTarget
                        plot(Con(ii,:),'Color',CB(ii,:)/255)
                    end
                    hold off
                    axis tight
            
                    subplot(2,1,2)
                    f = frame(:,:,:,i);
                    imshow(f)
                    hold on
                    for ii = 1:numTarget
                    
                        P = reg(:,:,ii,i);
                        P(P < 0.3) = 0;
                        P(P > 1) = 1;
                        mask = logical(P);
                        [area,centroids] = blobAnalyser(mask);
                        [~,aMax] = max(area);
                        
                        if isempty(centroids) == 0
            
                            xy(i,:,ii) = round(centroids(aMax,:));
                            xy(i,:,ii) = RemapPoint(xy(i,:,ii),[inputSize(1,1),inputSize(1,2)],[H,W],0,[1,1]);
                            plot(xy(i,1,ii),xy(i,2,ii),'o-','MarkerFaceColor',CB(ii,:)/255,'MarkerEdgeColor',CB(ii,:)/255);
                        else
                            xy(i,:,ii) = [nan,nan];
                        end
                    
                        [counts,centers] = hist(P(P > 0),linspace(0,1,20)); %#ok<HIST> 
                        [~,idx] = max(counts);
                        Con(ii,i) = centers(idx);
                    end
                    drawnow
                end
            
                while 1
            
                    [sel(1,1),sel(1,2),button] = ginput(1);
                    x = round(sel(1,1));
                    
                    subplot(2,1,2)
                    f = frame(:,:,:,x);
                    imshow(f)
                    hold on
                    for ii = 1:numTarget     
                        plot(xy(x,1,ii),xy(x,2,ii),'o-','MarkerFaceColor',CB(ii,:)/255,'MarkerEdgeColor',CB(ii,:)/255);
                    end
                    hold off
            
                    promptMessage = sprintf('Add frame to training data?');
                    titleBarCaption = 'settings';
                    answer = questdlg(promptMessage,titleBarCaption,'Yes','No','Exit','Yes');
                    
                    if strcmp(answer,'Yes') == 1
                        frameTrain(:,:,:,Count) = f;
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
            numTrainingFrames = size(frameTrain,4);
        end

        % Begin UI labeling
        xyTrain = zeros(numTarget,2,numTrainingFrames);
        goBack = 0;
        Stop = 0;
        f = figure('Resize','on');
        f.Position = [0,41,1920,963];
        i = 0;            
        while i < numTrainingFrames

            if goBack == 0
                i = i+1;
            else
                i = i-1;  
            end
            
            frame = frameTrain(:,:,:,i);
            if any(size(frame) ~= [BaseSize,3])
                frame = imresize(frame,BaseSize);
            end
            imshow(frame);
            hold on
            f.Position = [0,41,1920,963];    
    
            if goBack == 1
                for i3 = 1:numTarget-1
                     PlotPoint(1,i3) = plot(xyTrain(i3,1,i),xyTrain(i3,2,i),'o-','MarkerFaceColor',CB(i3,:)/255,'MarkerEdgeColor',CB(i3,:)/255); %#ok<SAGROW> 
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
                        xyTrain(ii,:,i) = sel;
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
        
        % Create vector of all points in image for alphashape query
        QP = zeros(H*W,2);
        Count = 1;
        for i4 = 1:W
            for i5 = 1:H
                QP(Count,:) = [i4,i5]; 
                Count = Count+1;
            end
        end
        
        % If new trainging data is being added to an existing model, set
        % aside for to sort into training and validation before saving
        numTrainingFrames = i;
        if strcmp(loadAction,'Update pose network') == 1
            if inputSize(1,1) == inputSize(1,1)
                m = 8;
            else
                m = 1;
            end
            vX = zeros([inputSize,3,numTrainingFrames*m]);
            vY = single(zeros([inputSize,numTarget,numTrainingFrames*m])); 
        end
        
        % Generate training data
        disp('Writting training data')
        Count = 1;
        for i = 1:numTrainingFrames
            
            frame = frameTrain(:,:,:,i);
            frame = (double(imresize(frame,inputSize))/255);
            xy = xyTrain(:,:,i);
            xy = RemapPoint(xy,BaseSize,inputSize,0,[1,1]);
            
            Y = single(zeros([inputSize,numTarget]));
            for ii = 1:numTarget
            
                x = xy(ii,2);
                y = xy(ii,1);
        
                if x == 0 && y == 0
                    continue
                end
                
                % Label ground truth frames
                p = nsidedpoly(1000,'Center',[x,y], 'Radius',3.5);        
                shp = alphaShape(p.Vertices);
                a = criticalAlpha(shp,'all-points')*2;
                shp.Alpha = a;    
                
                IS = inShape(shp,QP);
                IS = find(IS == 1);
               
                for i3 = 1:length(IS)                         
                    row = QP(IS(i3,1),1); 
                    col = QP(IS(i3,1),2);                     
                    Y(row,col,ii) = 1;
                end               
            end
            
            % Augment trainging frames if frames have equal dimentions
            if size(frame,1) == size(frame,2)
    
                frameO = frame;
                YO = Y;
    
                for ii = 1:4
                   
                    if strcmp(loadAction,'New') == 1
                        totalFrames = totalFrames+1;
                        XY = {frame,Y};
                        save([ProjectPathName,'\Training\Data\image_',num2str(totalFrames),'.mat'],'XY')
                    else
                        vX(:,:,:,Count) = frame;
                        vY(:,:,:,Count) = Y;
                        Count = Count+1;
                    end
                            
                    % reflect
                    frame = flip(frameO,1);
                    Y = flip(YO,1);
                    if size(SymPair,2) > 1
                        for i3 = 1:size(SymPair,1)                                
                            y = Y(:,:,SymPair(1,1));
                            Y(:,:,SymPair(1,1)) = Y(:,:,SymPair(1,2));
                            Y(:,:,SymPair(1,2)) = y;
                        end
                    end
                    
                    if strcmp(loadAction,'New') == 1
                        totalFrames = totalFrames+1;
                        XY = {frame,Y};
                        save([ProjectPathName,'\Training\Data\image_',num2str(totalFrames),'.mat'],'XY')
                    else
                        vX(:,:,:,Count) = frame;
                        vY(:,:,:,Count) = Y;
                        Count = Count+1;
                    end
    
                    %Rotate 90
                    frameO = imrotate(frameO,90);
                    frame = frameO;
                    YO = imrotate(YO,90);
                    Y = YO;
                end 
            else                
                if strcmp(loadAction,'New') == 1
                    totalFrames = totalFrames+1;
                    XY = {frame,Y};
                    save([ProjectPathName,'\Training\Data\image_',num2str(totalFrames),'.mat'],'XY')
                else
                    vX(:,:,:,Count) = frame;
                    vY(:,:,:,Count) = Y;
                    Count = Count+1;
                end
            end
        end
        save([ProjectPathName,'\ProjectData.mat'],'inputSize','Targets','totalFrames');
        
        if strcmp(loadAction,'New') == 1
            promptMessage = sprintf('Add another video?');
            titleBarCaption = 'settings';
            AddVideo = questdlg(promptMessage, titleBarCaption, 'Yes','No','Yes');
    
            if strcmp(AddVideo,'No') == 1
                break
            end
        else
            break
        end
    end
    
    if exist('vX','var') == 1 || strcmp(loadAction,'New') == 1
        promptMessage = sprintf('Would you like to train now?');
        titleBarCaption = 'settings';
        trainAction_CNN = questdlg(promptMessage, titleBarCaption, 'Yes','Later','Yes');    
        
        if exist('vX','var') == 1 && strcmp(trainAction_CNN,'Later') == 1
            for i = 1:size(vX,4)
                    
                XY = {vX(:,:,:,i),vY(:,:,:,i)};
                totalFrames = totalFrames+1;
                save([ProjectPathName,'\Training\Data\image_',num2str(totalFrames),'.mat'],'XY')
            end
            save([ProjectPathName,'\ProjectData.mat'],'inputSize','Targets','totalFrames');
        end
    else
        trainAction_CNN = 'NA'; 
    end
    trainAction_LSTM = 'NA'; 

elseif strcmp(labelAction_LSTM,'Yes') == 1 
    
    load([ProjectPathName,'\Training\SequenceData.mat']);
    load([ProjectPathName,'\netCNN.mat'])
    Classes = sequenceData{2,1}(:,1);
    numClass = size(sequenceData{2,1},1);
    
    if strcmp(reactAction,'No') == 1
        for i = 1:numClass
    
            numFiles = dir(fullfile([ProjectPathName,'\Training\Sequences\',Classes{i,1}],'*.mat'));
            numFiles = size(numFiles,1);
            sequenceData{2,1}(i,3) = {numFiles};
            sequenceData{2,1}(i,2) = {numFiles/8};
        end
    else
        sequenceData{2,1}(:,2:3) = {0};
    end
    
    while 1
        if strcmp(reactAction,'No') == 1

            if exist('VidFilePath','var') == 1 && ischar(VidFilePath) == 1
                cd(VidFilePath)
            end
        
            [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select training videos','MultiSelect','off');
            if ischar(VidFilePath) == 0
                break
            end
            vidFile = fullfile([VidFilePath,VidFileName]);
            cd(VidFilePath) 
    
            if strcmp(loadAction,'Refine network') == 1    
                [CoorFileName,CoorFilePath] = uigetfile('*.mat','Select coordinate data','MultiSelect','off');
                load(fullfile(CoorFilePath,CoorFileName));    
            end   
        else
            if iFile > length(sequenceData{2,2})
                break
            end
            disp(['File ',num2str(iFile),'/',num2str(length(sequenceData{2,2}))]);
            vidFile = sequenceData{2,2}{iFile,1};                               
        end                

        % Read Video
        disp('Reading video file...')
        vr = VideoReader(vidFile); %#ok<TNMLP> 
        H = vr.Height;
        W = vr.Width;
        C = 3;
        fps = round(vr.FrameRate);          
        
        if strcmp(reactAction,'No') == 1
            
            numSeq = floor(vr.NumFrames/vr.FrameRate)-1;
            video = zeros(H,W,C,framesPerObs);
            idx = 1:numSeq;
            idxPast = zeros(numSeq,1);
            seqClass = zeros(numSeq,1);
            
            i = 1;
            idxLoc = 1;
            randNext = 0;
            exitPick = 0;
            skip = 0;
            goBack = 0;
            f = figure('Resize','on');
            f.Position = [553,393,814,573];
            while i <= numSeq                       
                if i > 1
                    idxPast(i,1) = loc;
                end
                if goBack == 1 
                    if seqClass(i,1) > 0
                        sequenceData{2,1}(seqClass(i,1),2) = sequenceData{2,1}(seqClass(i,1),2);
                    end
                    goBack = 0;
                else                       
                    if randNext == 0
                        try
                            loc = idx(1,idxLoc);
                        catch
                            idxLoc = 1;
                            loc = idx(1,idxLoc);
                        end                           
                        idx(:,idxLoc) = [];
                    elseif randNext == 1
                        [loc,idxLoc] = datasample(idx,1);
                        idx(:,idxLoc) = [];
                    elseif randNext == 2

                        timePrompt{1,1} = 'Minutes:';
                        timePrompt{1,2} = 'Seconds:';
                        dlgtitle = ' ';
                        defInput{1,1} = '0';
                        defInput{1,2} = '0';
                        skip2Time = inputdlg(timePrompt,dlgtitle,[1,30],defInput);

                        if isempty(skip2Time) == 0

                            skip2Time = (str2double(skip2Time{1,1})*60)+str2double(skip2Time{2,1});                            
                            idxLoc = find(idx == skip2Time);
                            if isempty(idxLoc) == 0
                                loc = skip2Time;
                            else
                                continue
                            end
                            idx(:,idxLoc) = [];
                        else
                            continue
                        end
                    end
                end           
                
                for ii = 0:framesPerObs-1                 
                    frame = read(vr,(loc*fps)+ii);
                    video(:,:,:,ii+1) = frame;        
                end            
    
                ClassCount = cell(numClass,1);
                for ii = 1:numClass
                    ClassCount{ii,1} = [Classes{ii,1},' (',num2str((sequenceData{2,1}{ii,2})),')'];
                end            
                
                while 1
                    for ii = 1:framesPerObs
                        frame = video(:,:,:,ii)/255;
                        frame = imresize(frame,BaseSize);
                        imshow(frame);
                        drawnow
                        f.Position = [553,393,814,573];                    
                    end
    
                    [pick,tf] = listdlg('PromptString',{['Sequence: ',num2str(i),'/',num2str(numSeq)];'Select a class,';'Press cancel to exit.'}, ... 
                    'ListString',[{'Play Again'};ClassCount;{'Skip Next'};{'Skip Random'};{'Skip to time'};{'Go back'}],'SelectionMode','single','InitialValue',1);
    
                    if isempty(pick)
                        exitPick = 1;
                        break
                    elseif pick > 1 && pick < numClass+2
                        seqClass(i,1) = pick-1;
                        sequenceData{2,1}{pick-1,2} = sequenceData{2,1}{pick-1,2}+1;
                        break
                    elseif pick >= numClass+2 && pick <= numClass+4
                        skip = 1;
                        if pick == numClass+2
                            randNext = 0;
                        elseif pick == numClass+3
                            randNext = 1;
                        elseif pick == numClass+4
                            randNext = 2;
                        end
                        break
                    elseif pick == numClass+5
                        goBack = 1;
                        break
                    end
                end
                if exitPick == 1
                    break
                elseif skip == 1
                    skip = 0;
                    i = i+1; 
                    continue
                elseif goBack == 1
                    if i > 1
                        i = i-1;
                    end
                else           
                    i = i+1; 
                    randNext = 0;
                end
            end 
            close(f)      
            
            labeled = seqClass > 0;
            seqClass(seqClass == 0,:) = [];
            labeled = idxPast(labeled,1);
            [labeled,idx] = sort(labeled);
            seqClass = seqClass(idx,1);
        else
            labeled = sequenceData{2,3}{iFile,1};
            seqClass = sequenceData{2,4}{iFile,1}; 
            iFile = iFile+1;  
        end
        numLabeled = length(labeled);
        
        Cvideo = zeros([inputSize,C,framesPerObs]);
        wb = waitbar(0/numLabeled,{'Writing training data...'});
        for i = 1:numLabeled    
            
            infoAct = zeros((actLayer_numFilt*6)+1,framesPerObs);
            pixelDiff = zeros(1,framesPerObs);
            for ii = 1:framesPerObs  
                frame = read(vr,(labeled(i,1)*fps)+ii);      
                Cvideo(:,:,:,ii) = (double(imresize(frame,inputSize))/255);
                if ii > 1
                    pixelDiff(1,ii) = mean(Cvideo(:,:,:,ii)-Cvideo(:,:,:,ii-1),'all');
                end
            end  
            infoAct(1,:) = pixelDiff;
            
            Y = categorical(seqClass(i,1),1:numClass);

            for ii = 1:4                
                for iii = 1:2

                    sequenceData{2,1}{seqClass(i,1),3} = sequenceData{2,1}{seqClass(i,1),3}+1;
                    Act = activations(netCNN,Cvideo,actLayer);
                    
                    for i3 = 1:size(Act,4)
    
                        info = zeros(5,size(Act,3));
                        for i4 = 1:size(Act,3)
    
                            temp = Act(:,:,i4,i3);
                            
                            info(1,i4) = mean(temp,'all')/10;
                            
                            numZero = length(find(temp == 0));
                            info(2,i4) = numZero/(size(Act,1)*size(Act,2));
    
                            if mean(temp,'all') > 0
                                temp(temp == 0) = nan;
                            end
                    
                            [m,locm] = min(temp,[],'all');
                            info(3,i4) = m/10;
                            info(4,i4) = locm/(size(Act,1)*size(Act,2));
                    
                            [M,locM] = max(Act(:,:,i4,i3),[],'all');
                            info(5,i4) = M/10;
                            info(6,i4) = locM/(size(Act,1)*size(Act,2));
                        end
                        info = info-0.5;
                        infoAct(2:end,i3) = info(:);
                    end  

                    X = infoAct;
                    save([ProjectPathName,'\Training\Sequences\',Classes{seqClass(i,1),1},'\Sequence_',num2str(sequenceData{2,1}{seqClass(i,1),3}),'.mat'],'X','Y');
                    
                    Cvideo2 = flip(Cvideo,1);
                end            
                Cvideo = imrotate(Cvideo,90);
            end            
            waitbar(i/numLabeled,wb,{'Writing training data...'});
        end
        close(wb)
        
        if strcmp(updateAction,'Retrain network') == 0
            if isempty(sequenceData{2,2})
                sequenceData{2,2} = {[VidFilePath,VidFileName]};
                sequenceData{2,3} = {labeled};
                sequenceData{2,4} = {seqClass};
            else
                sequenceData{2,2}(end+1,1) = {[VidFilePath,VidFileName]};
                sequenceData{2,3}(end+1,1) = {labeled};
                sequenceData{2,4}(end+1,1) = {seqClass};
            end                
    
            promptMessage = sprintf('Add another video?');
            titleBarCaption = 'settings';
            AddVideo = questdlg(promptMessage, titleBarCaption, 'Yes','No','Yes');
    
            if strcmp(AddVideo,'No') == 1
                break
            end
        else
            for i = 1:numClass
                sequenceData{2,1}{i,2} = sequenceData{2,1}{i,3}/8;
            end
        end            
    end
    save([ProjectPathName,'\Training\SequenceData.mat'],'sequenceData');
    
    if strcmp(reactAction,'No') == 1
    promptMessage = sprintf('Would you like to train now?');
    titleBarCaption = 'settings';
    trainAction_LSTM = questdlg(promptMessage, titleBarCaption, 'Yes','Later','Yes');
    else
        trainAction_LSTM = 'Yes';
    end
    trainAction_CNN = 'NA';
end
%%
if strcmp(trainAction_CNN,'Yes') == 1
    if strcmp(loadAction,'New') == 0 && strcmp(updateAction,'Retrain network') == 0

        L = size(vX,4);
        idx = randperm(L);
        vX = vX(:,:,:,idx);
        vY = vY(:,:,:,idx);

        numVal = rem(L,miniBatchSizeCNN);
        if numVal < 80
            numVal = numVal+(miniBatchSizeCNN*round(100/miniBatchSizeCNN));
        end

        if strcmp(updateAction,'Refine network') == 0
            for i = 1:L-numVal
                
                frame = vX(:,:,:,i);
                Y = vY(:,:,:,i);
                totalFrames = totalFrames+1;
                XY = {frame,Y};
                save([ProjectPathName,'\Training\Data\image_',num2str(totalFrames),'.mat'],'XY')
            end
            vX(:,:,:,1:i) = [];
            vY(:,:,:,1:i) = [];

        elseif strcmp(updateAction,'Refine network') == 1

            tX = vX(:,:,:,1:end-numVal);
            tY = vY(:,:,:,1:end-numVal);
            vX(:,:,:,1:end-numVal) = [];
            vY(:,:,:,1:end-numVal) = [];
        end
    end

    % Create datastores for training frames and their labels
    %%
    imDirTrain = [ProjectPathName,'\Training\Data'];
    imds = fileDatastore(imDirTrain,'ReadFcn',@load2VarLabeled);
    imds = shuffle(imds);
    numObs = length(imds.Files);
    %%

    if numObs > 200
        
        % Pull observations for validation and make training data
        % divisible by the minibach size    
        if strcmp(updateAction,'Refine network') == 0  

            numVal = rem(numObs,miniBatchSizeCNN);

            if strcmp(loadAction,'New') == 1 || strcmp(updateAction,'Retrain network') == 1
            
                if numVal < 80
                    numVal = numVal+(miniBatchSizeCNN*round(100/miniBatchSizeCNN));
                end
                vX = zeros([inputSize,3,numVal]);
                vY = single(zeros([inputSize,numTarget,numVal]));   
            end                     
            
            for i = 1:numVal                
                if strcmp(loadAction,'New') == 1 || strcmp(updateAction,'Retrain network') == 1
                    readDS = read(imds);
                    vX(:,:,:,i) = readDS{1,1};
                    vY(:,:,:,i) = readDS{1,2};
                else
                    readDS = read(imds);
                    vX = cat(4,vX,readDS{1,1});
                    vY = cat(4,vY,readDS{1,2});
                end
    
                delimg = matches(imds.Files,imds.Files{1,1});
                imds = subset(imds,~delimg);
            end
            numObs = length(imds.Files);
        elseif strcmp(updateAction,'Refine network') == 1

            numVal = 0;
            
            tX2 = zeros(size(tX));
            tY2 = zeros(size(tY));
            for i = 1:size(vX,4)
                readDS = read(imds);
                tX2(:,:,:,i) = readDS{1,1};
                tY2(:,:,:,i) = readDS{1,2};
            end
            tX2 = cat(4,tX,tX2);
            tY2 = cat(4,tY,tY2);                
        end            

        if strcmp(loadAction,'New') == 1 || strcmp(updateAction,'Retrain network') == 1
            S = 1; %[1,2,3];
            loadNet = 0;
        elseif strcmp(updateAction,'Refine network') == 1
            S = [4,5];
            loadNet = 1;
        else
            S = 6;
            loadNet = 1;
        end

        for i = S
        
            % Set options for network training
            if i == 1
                numEpoch = round(numIter/(numObs/miniBatchSizeCNN));
                learnRate = 0.001;
                learnRateDrop = round((numIter*0.33)/(numObs/miniBatchSizeCNN));
                outPut = 'best-validation-loss';
            elseif i == 2
                numObs = length(imds.Files);
                numEpoch = 10;
                learnRate = 0.001;
                learnRateDrop = 8;
                outPut = 'last-iteration';
            elseif i == 3 || i == 6
                numObs = length(imds.Files);
                numEpoch = 10; 
                learnRate = 0.0005;
                learnRateDrop = 8;
                outPut = 'last-iteration';                
            elseif i == 4
                numObs = size(tX2,4);
                numEpoch = 2;
                learnRate = 0.0005;
                learnRateDrop = 2;
                outPut = 'best-validation-loss';
            elseif i == 5
                numObs = length(imds.Files);
                numEpoch = round(numIter/(numObs/miniBatchSizeCNN));
                learnRate = 0.001;
                learnRateDrop = round((numIter*0.95)/(numObs/miniBatchSizeCNN));
                outPut = 'last-iteration';
            end

            options = trainingOptions('sgdm', ...
            'Momentum',0.95, ...
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
            'ValidationPatience',20, ...
            'Plots','training-progress', ...
            'Verbose',true, ...
            'VerboseFrequency',floor((numObs/miniBatchSizeCNN)/4), ...
            'ExecutionEnvironment',execution, ...
            'OutputNetwork',outPut);            
            
            if i == 1
                 layersCNN = PartDetectorLayers;
            elseif i == 2
               
                layersCNN = replaceLayer(layersCNN,layersCNN.Layers(end).Name,reluLayer('Name','reLu_Out'));
                layersCNN = setLearnRate(layersCNN,0);
                layersCNN = addLayers(layersCNN,RegressionLayers);
                layersCNN = connectLayers(layersCNN,'reLuRegBN_1.2','addReg_1.0/in2');
                layersCNN = connectLayers(layersCNN,'reLu_Out','FuseReg_1.0/in2');
                layersCNN = connectLayers(layersCNN,'reLu_Out','Infuse/in1');
                layersCNN = connectLayers(layersCNN,'Input','Infuse/in2');
                
            elseif i > 3
                if loadNet == 1
                    load([ProjectPathName,'\netCNN.mat']);
                    loadNet = 0;
                end
                layersCNN = layerGraph(netCNN);
            end

            % Train network
            if strcmp(execution,'gpu') == 1
                reset(gpuDevice(1))
            end

            if i ~= 4

                netCNN = trainNetwork(imds,layersCNN,options);
                layersCNN = layerGraph(netCNN);

                if i == 2
                    layersCNN = setLearnRate(layersCNN,1);
                end

            elseif i == 4

                netCNN = trainNetwork(tX2,tY2,layersCNN,options);
                for ii = 1:size(tX,4)
                
                    totalFrames = totalFrames+1;
                    XY = {tX(:,:,:,i),tY(:,:,:,i)};
                    save([ProjectPathName,'\Training\Data\image_',num2str(totalFrames),'.mat'],'XY')
                end
                save([ProjectPathName,'\ProjectData.mat'],'inputSize','Targets','totalFrames');

                imDirTrain = [ProjectPathName,'\Training\Data'];
                imds = fileDatastore(imDirTrain,'ReadFcn',@load2VarLabeled);
                imds = shuffle(imds);
            end
            if strcmp(execution,'gpu') == 1
                reset(gpuDevice(1))
            end
        end
        save([ProjectPathName,'\netCNN.mat'],'netCNN');
        disp('Your network is now trained =D')

        if strcmp(loadAction,'New') == 0 && strcmp(updateAction,'Retrain network') == 0              
            for i = 1:size(vX,4)
                
                totalFrames = totalFrames+1;
                XY = {vX(:,:,:,1),vY(:,:,:,1)};
                save([ProjectPathName,'\Training\Data\image_',num2str(totalFrames),'.mat'],'XY')
            end
            save([ProjectPathName,'\ProjectData.mat'],'inputSize','Targets','totalFrames');
        end
    else
        disp('Please add more training frames.')
    end

elseif strcmp(trainAction_LSTM,'Yes') == 1

    trainingData = fileDatastore([ProjectPathName,'\Training\Sequences'],"ReadFcn",@load2VarLabeled,'IncludeSubfolders',true);
    trainingData = shuffle(trainingData);
    numObs = length(trainingData.Files); 
    
    % Pull observations for validation and make training data
    % divisible by 100    
    numVal = rem(numObs,miniBatchSizeLSTM);
    
    if numVal < 80
        numVal = numVal+(miniBatchSizeLSTM*round(100/miniBatchSizeLSTM));
    end
    vX = cell(numVal,1);
    vY = categorical(ones(1,numVal),1:numClass)';
    
    for i = 1:numVal    
        
        T = read(trainingData);
        vX{i,1} = T{1,1};
        vY(i,1) = T{1,2};
    
        delimg = matches(trainingData.Files,trainingData.Files{1,1});
        trainingData = subset(trainingData,~delimg);
    end
    numObs = length(trainingData.Files); 
    
    options = trainingOptions('sgdm', ...
    'MaxEpochs',200, ...
    'MiniBatchSize',miniBatchSizeLSTM, ...
    'InitialLearnRate',0.001, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',90, ...
    'GradientThreshold',2, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{vX,vY}, ... 
    'ValidationFrequency',numObs/miniBatchSizeLSTM, ...
    'ValidationPatience',50, ...
    'Plots','training-progress', ...
    'Verbose',false, ...
    'ExecutionEnvironment','gpu', ...
    'OutputNetwork','last-iteration');
    
    reset(gpuDevice(1))
    netLSTM = trainNetwork(trainingData,LSTMlayers,options);  

    save([ProjectPathName,'\netGRU.mat'],'netLSTM');
    disp('Your network is now trained =D')
end
%%
if strcmp(loadAction,'Analyse video') == 1 
    
    % UI video select
    if exist('VidFilePath','var') == 1 && ischar(VidFilePath) == 1
        cd(VidFilePath)
    end

    vidFiles = cell(1,1);
    dataFiles = cell(1,1);
    i = 1;
    while 1

        [VidFileName,VidFilePath] = uigetfile('*.mp4;*.avi;*.mov','Select video file','MultiSelect','off');
        vidFiles{i,1} = fullfile(VidFilePath,VidFileName);
        cd(VidFilePath);       
        
        if Sync == 1
    
            promptMessage = sprintf('Would you like to sync to neural data?');
            titleBarCaption = 'settings';
            syncNeural = questdlg(promptMessage, titleBarCaption, 'Yes','No','Yes');
        
            if strcmp(syncNeural,'Yes') == 1         
                [DigFileName,DigFilePath] = uigetfile('*.dat','Select digitalin.dat file','MultiSelect','off');
                dataFiles{i,1} = fullfile(DigFilePath,DigFileName);
            else
                dataFiles{i,1} = 'N/A';
            end
        end

        promptMessage = sprintf('Would you like to add another file?');
        titleBarCaption = 'settings';
        addFile = questdlg(promptMessage, titleBarCaption, 'Yes','No','Yes');
    
        if strcmp(addFile,'Yes') == 1 
            i = i+1;
        else
            break
        end
    end

    for I = 1:size(vidFiles,1)

        disp("Working on file: "+I+"/"+size(vidFiles,1))
    
        % create video reader
        disp('Reading video file...')
        vr = VideoReader(vidFiles{I,1}); %#ok<TNMLP> 
        H = vr.Height;
        W = vr.Width;
        C = 3;
        fps = round(vr.FrameRate);
        linSync = 0;
        
        if Sync == 1   
    
            disp('Finding LED')
    
            % Create object detector to detect LED pulses
            blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort',true,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',150,'MaximumBlobArea',800);    
            
            % Gathering frames from first part of video and find bounding
            % box of LED
            numTrainingFrames = (fps*LEDtime)*2;
            vr.CurrentTime = 0;
            lastFrame = nan;
            i = 1;
            BBox = 0;
            Ck = 1;
            while i <= numTrainingFrames
                
                % Read frame and crop to ROI
                frame = rgb2gray(readFrame(vr));
                if mean(frame,'all') > 100
                    
                    frame(1:end,100:end-99) = 0;
                    frame(1:end-99,1:end) = 0;
    
                    if all(isnan(lastFrame),'all')
                        lastFrame = frame;
                    else
                        frameLED = frame-lastFrame;
                        frameLED(frameLED < 50) = 0;
                        frameLED(frameLED > 0) = 1;
                        
                        [A,centroids,BB] = blobAnalyser.step(logical(frameLED));
                        if isempty(A) == 0
                            
                            BB = double(BB);
                            if Ck == 1
                                BBox = BB;
                                Ck = 0;
                            else
                                BBox = round(mean([BBox;BB],1));
                            end 
                        else
                            lastFrame = frame;
                            if Ck == 0
                                break
                            end
                        end
                    end
                    i = i+1;
                end
            end   
            if ~all(BBox == 0,'all')
                BBox = bbox2points(BBox);  
                frameDims = size(lastFrame);
                BBox(BBox(:,1) > frameDims(1,2),1) = frameDims(1,2);
                BBox(BBox(:,2) > frameDims(1,1),2) = frameDims(1,1);
                lastFrame = lastFrame(min(BBox(:,2)):max(BBox(:,2)),min(BBox(:,1)):max(BBox(:,1)),:);
            else
                linSync = 1;
                disp('LED not dectected, switching to linear sync')
            end
    
            if strcmp(dataFiles{I,1},'N/A') == 0   
    
                disp('Analyzing digitalin.dat')
                
                samplesPerFrame = floor(Sfreq/fps);
                s = dir(dataFiles{I,1});
                numSamples = s.bytes/2;

                if linSync == 1

                    theoreticalFrames = floor(numSamples/samplesPerFrame);
                    frameSpacing = round(linspace(1,vr.NumFrame,theoreticalFrames));
                else
                
                    fileID = fopen(dataFiles{I,1});
                    fseek(fileID,0,'bof');
                    samplesPerPulse = 0;
                    Ck = 0;
                    Count = 1;
                    while 1
        
                        temp = fread(fileID,[1,Sfreq*100],'int16');
                        if isempty(temp) == 0
                            for ii = 1:length(temp)
                                if temp(1,ii) ~= LEDchannel
                                    samplesPerPulse(Count,1) = samplesPerPulse(Count,1)+1; %#ok<SAGROW> 
                                    Ck = 0;
                                elseif temp(1,ii) == LEDchannel && Ck == 0
                                    Count = Count+1;
                                    samplesPerPulse(Count,1) = 1; %#ok<SAGROW> 
                                    Ck = 1;
                                elseif temp(1,ii) == LEDchannel && Ck == 1
                                    samplesPerPulse(Count,1) = samplesPerPulse(Count,1)+1; %#ok<SAGROW> 
                                end
                            end
                        else
                            break
                        end
                    end            
                    framesPerPulse = round(samplesPerPulse/samplesPerFrame);      
                end
            end  
            
            % Create video writer objects
            writerObj_Sync = VideoWriter([vidFiles{I,1}(1,1:end-4),'_Sync.mp4'],'MPEG-4'); %#ok<TNMLP> 
            writerObj_Sync.FrameRate = fps;
            open(writerObj_Sync);
        end  
    
        if Pose == 1
    
            load([ProjectPathName,'\netCNN.mat'])
    
            % Create object detector to detect CNN parts
            blobAnalyserCNN = vision.BlobAnalysis('BoundingBoxOutputPort',false,'AreaOutputPort',true,'CentroidOutputPort',true,'MinimumBlobArea',10);
    
            % Create vector of all points in image for alphashape query
            QP = zeros(H*W,2);
            Count = 1;
            for i4 = 1:W
                for i5 = 1:H
                    QP(Count,:) = [i4,i5]; 
                    Count = Count+1;
                end
            end
            
            % Create array to store pose coordinates
            allXY = cell(2,numTarget+1);
            for i = 1:numTarget
                allXY{1,i} = Targets{i,1};
            end
            allXY{1,end} = 'Behavior';
    
            if isempty(gcp('nocreate')) == 1 && numWorkers > 0
                parpool('local',numWorkers);
            end 
            
            % Create video writer objects
            if labelVidPose == 1 || labelVidBehav == 1
    
                writerObj_Pose = VideoWriter([vidFiles{I,1}(1,1:end-4),'_Pose.mp4'],'MPEG-4'); %#ok<TNMLP> 
                writerObj_Pose.FrameRate = fps;
                open(writerObj_Pose);
            end
        end   
        
        % Begin reading frames from video    
        
        if linSync == 0
            numFrame = vr.NumFrames;
        elseif linSync == 1
            numFrame = theoreticalFrames;
        end        

        disp('Analyzing video file')
        vr.CurrentTime = 0;
        iFrame = 1;
        iPulse = 1;
        wb = waitbar(iFrame/numFrame,{'Reading frames...';['Frame ',num2str(1),'/',num2str(numFrame)]});
        Stop = 0;
        while iFrame <= numFrame
            
            % Prealocate frame store 
            if Sync == 0 || linSync == 1
                vid = uint8(zeros(H,W,C,fps*10));
            elseif Sync == 1
                if strcmp(dataFiles{I,1},'N/A') == 1
                    vid = uint8(zeros(H,W,C,(fps*LEDtime)+fps));
                elseif strcmp(dataFiles{I,1},'N/A') == 0
                    if iPulse > size(framesPerPulse,1)
                        break
                    end
                    vid = uint8(zeros(H,W,C,framesPerPulse(iPulse,1)));
                end
            end
           
            % Read frames from video in chunks. Either in 1s intervals if sync
            % is not enabled or all frames between LED pulses
            numFrames = 0;  
            Ck = 0;
            while 1
    
                waitbar(iFrame/numFrame,wb,{'Reading frames...';['Frame ',num2str(iFrame),'/',num2str(numFrame)]});                
                
                try
                    if linSync == 0
                        frame = read(vr,iFrame);
                    elseif linSync == 1
                        frame = read(vr,frameSpacing(1,iFrame));
                    end
                    numFrames = numFrames+1;
                    iFrame = iFrame+1;
                catch
                    Stop = 1;
                    break
                end                            
    
                if Sync == 1 && linSync == 0

                    if numFrames > round(framesPerPulse(iPulse,1)*1.2)
                        numFrames = numFrames-1;
                        iFrame = iFrame-1;
                        break
                    end
    
                    frameLED = rgb2gray(frame(min(BBox(:,2)):max(BBox(:,2)),min(BBox(:,1)):max(BBox(:,1)),:));  
                    frameLED = frameLED-lastFrame;
                    frameLED = imresize(frameLED,[25,25]);   
                    frameLED(frameLED < 50) = 0;
    
                    if sum(frameLED,'all') < 7500

                        vid(:,:,:,numFrames) = frame;
                        if mean(rgb2gray(frame),'all') > 100
                            lastFrame = rgb2gray(frame(min(BBox(:,2)):max(BBox(:,2)),min(BBox(:,1)):max(BBox(:,1)),:));  
                        end
                        Ck = Ck+1;
                    elseif sum(frameLED,'all') >= 5000
                        if Ck > 3
                            numFrames = numFrames-1;
                            iFrame = iFrame-1;
                            break
                        else
                            vid(:,:,:,numFrames) = frame;
                            Ck = 0;
                        end
                    end                
                else
                    if numFrames <= fps*10
                        vid(:,:,:,numFrames) = frame;
                        if numFrames == fps*10
                            break
                        end
                    end                
                end
            end
            vid = vid(:,:,:,1:numFrames);

            % if sync is enabled the number of frames between LED pulses is
            % compared to the expected value and frames are added or subtracted
            % as needed
            if Sync == 1 && Stop == 0
                if linSync == 0                
            
                    waitbar(iFrame/numFrame,wb,{'Syncing video...';['Frame ',num2str(iFrame),'/',num2str(numFrame)]});

                    if strcmp(dataFiles{I,1},'N/A') == 1
                        Expected = fps*LEDtime;
                    else
                        Expected = framesPerPulse(iPulse,1);
                        framesPerPulse(iPulse,2) = numFrames;
                        iPulse = iPulse+1;
                    end
                    Missing = Expected-numFrames;
                
                    vidSync = uint8(zeros(H,W,C,Expected));
                    
                    if Missing > 0
                        if Missing < numFrames
                            Add = sort(datasample(1:numFrames,Missing,'Replace',false));
                        else
                            Add = sort(datasample(1:numFrames,Missing,'Replace',true));
                        end
                        if size(Add,1) > 1
                            Add = Add';
                        end
                        Add = sort([1:numFrames,Add]);
                        for ii = 1:Expected
                            vidSync(:,:,:,ii) = vid(:,:,:,Add(1,ii));
                        end
                    elseif Missing < 0
        
                        Sub = round(linspace(1,numFrames,Expected));
                        for ii = 1:Expected
                            vidSync(:,:,:,ii) = vid(:,:,:,Sub(1,ii));
                        end
                    else
                        for ii = 1:Expected
                            vidSync(:,:,:,ii) = vid(:,:,:,ii);
                        end
                    end
                elseif linSync == 1
                    Expected = numFrames;
                    vidSync = vid;
                end                
            end    
            
            % If pose is enabled each frame is processed by the CNN and the
            % target coordinates are determined and frames are labeled
            if Pose == 1            
                
                waitbar(iFrame/numFrame,wb,{'Detecting targets...';['Frame ',num2str(iFrame),'/',num2str(numFrame)]}); 
                
                if Sync == 1 && Stop == 0 
                    vidPose = vidSync;
                else
                    vidPose = vid;
                    Expected = size(vid,4);
                end
    
                % Make heat maps with CNN
                vidPoseC = (double(imresize(vidPose,inputSize))/255);
                YPred = predict(netCNN,vidPoseC(:,:,:,1:Expected),'ExecutionEnvironment',execution,'Acceleration','auto');      
                
                % Find coordinates from the heatmaps and get confidence
                XY = zeros(Expected,3,numTarget);
                if numWorkers > 0
                    parfor ii = 1:Expected                  
                        for i3 = 1:numTarget
                        
                            P = YPred(:,:,i3,ii);
                            P(P < 0.2) = 0;
                            P(P > 1) = 1;
                            mask = logical(P);
                            [area,centroids] = blobAnalyserCNN(mask); %#ok<PFBNS> 
                            if isempty(centroids) == 0

                                [~,aMax] = max(area);
                                xy = round(centroids(aMax,:)); 
                                xy = RemapPoint(xy,[inputSize(1,1),inputSize(1,2)],[H,W],0,[1,1]); %#ok<PFBNS> 

                                [counts,centers] = hist(P(P > 0),linspace(0.01,1,20)); %#ok<HIST> 
                                [~,idx] = max(counts);
                                Con = centers(idx);
                            else
                                xy = [nan,nan];
                                Con = 0;
                            end
                            XY(ii,:,i3) = [xy(1,1),xy(1,2),Con];                
                        end            
                    end
                else            
                    for ii = 1:Expected                  
                        for i3 = 1:numTarget
                        
                            P = YPred(:,:,i3,ii);
                            P(P < 0.3) = 0;
                            P(P > 1) = 1;
                            mask = logical(P);
                            [area,centroids] = blobAnalyserCNN(mask);
                            if isempty(centroids) == 0

                                [~,aMax] = max(area);
                                xy = round(centroids(aMax,:)); 
                                xy = RemapPoint(xy,[inputSize(1,1),inputSize(1,2)],[H,W],0,[1,1]);

                                [counts,centers] = hist(P(P > 0),linspace(0.01,1,20)); %#ok<HIST> 
                                [~,idx] = max(counts);
                                Con = centers(idx);
                            else
                                xy = [nan,nan];
                                Con = 0;
                            end                            
                            XY(ii,:,i3) = [xy(1,1),xy(1,2),Con];           
                        end            
                    end
                end
                
                % Apply smoothing to coordinates to reduce jitter
                if xSmooth == 1
                   for i = 1:numTarget
    
                        temp = XY(:,1:2,i);          
                        
                        ii = 1;
                        notNAN1 = nan;
                        numNAN = 0;
                        idx = false(size(temp,1),1);
                        ck = 0;
                        while ii <= size(temp,1)
                            if ~isnan(temp(ii,1))
                                if ck == 0
                                    notNAN1 = temp(ii,:);
                                    numNAN = 0;
                                    idx(ii,1) = true;
                                    ii = ii+1;
                                else
                                    notNAN2 = temp(ii,:);
                                    temp(ii-numNAN-1:ii,1) = round(linspace(notNAN1(1,1),notNAN2(1,1),numNAN+2));
                                    temp(ii-numNAN-1:ii,2) = round(linspace(notNAN1(1,2),notNAN2(1,2),numNAN+2));
                                    idx(ii-numNAN-1:ii,1) = true;
                                    ck = 0;
                                end
                            elseif isnan(temp(ii,1))
                                if ~isnan(notNAN1)

                                    numNAN = numNAN+1;
                                    ck = 1;
    
                                    if numNAN > 5
                                        ck = 0;
                                    end
                                end
                                ii = ii+1;
                            end
                        end

                        for ii = 1:2
                            for i3 = 1:2
                                temp(idx,ii) = smooth(temp(idx,ii));
                            end
                            temp(idx,ii) = floor(temp(idx,ii));
                        end
                    end
                end
    
                if Behavior == 1
                    
                    B = cell(Expected,3);                
                    for i = 1:fps:Expected
    
                        infoAct = zeros((actLayer_numFilt*6)+1,fps);
    
                        pixelDiff = zeros(1,fps);
                        for ii = 2:fps
                            pixelDiff(1,ii) = mean(vidPoseC(:,:,:,ii)-vidPoseC(:,:,:,ii-1),'all');
                        end
                        infoAct(1,:) = pixelDiff;
                        
                        Act = activations(netCNN,vidPoseC,actLayer);
                        for ii = 1:size(Act,4)
        
                            info = zeros(5,size(Act,3));
                            for i3 = 1:size(Act,3)
                        
                                temp = Act(:,:,i3,ii);
    
                                info(1,i3) = mean(temp,'all');
    
                                numZero = length(find(temp == 0));
                                info(2,i3) = numZero/(size(Act,1)*size(Act,2));
                                temp(temp == 0) = nan;
                        
                                [m,locm] = min(temp,[],'all');
                                info(3,i3) = m/10;
                                info(4,i3) = locm/(size(Act,1)*size(Act,2));
                        
                                [M,locM] = max(Act(:,:,i3,ii),[],'all');
                                info(5,i3) = M/10;
                                info(6,i3) = locM/(size(Act,1)*size(Act,2));
                            end
                            info = info-0.5;
                            infoAct(2:end,ii) = info(:);
                        end
                        
                        YPred = predict(netLSTM,infoAct,'ExecutionEnvironment','gpu');                   
                        [con,class] = max(YPred);
                        B(i:i+fps-1,1) = {class};
                        B(i:i+fps-1,2) = Classes(class,1);
                        B(i:i+fps-1,3) = {con};
                    end
                    allXY{2,end} = [allXY{2,end};B];
                else
                    B = 0;
                end
                
                % Label video frames
                if numWorkers > 0
                    parfor ii = 1:Expected
        
                        vidFrame = vidPose(:,:,:,ii);
        
                        if labelVidPose == 1
                            for i3 = 1:numTarget 
                                
                                xy = XY(ii,1:2,i3); %#ok<PFBNS> 
                                Con = XY(ii,3,i3); 
            
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
                                        vidFrame(row,col,:) = CB(i3,:); %#ok<PFBNS> 
                                    end
                                end                     
                            end
                        end
                        if Behavior == 1 && labelVidBehav == 1
                            vidFrame = insertText(vidFrame,[0,0],[B{ii,2},' ',num2str(B{ii,3})],'FontSize',18,'BoxColor','blue','BoxOpacity',0.4,'TextColor','white'); %#ok<PFBNS> 
                        end
                        vidPose(:,:,:,ii) = vidFrame;
                    end
                else
                    for ii = 1:Expected
        
                        vidFrame = vidPose(:,:,:,ii);
        
                        if labelVidPose == 1
                            for i3 = 1:numTarget 
                                
                                xy = XY(ii,1:2,i3);
                                Con = XY(ii,3,i3); 
            
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
                                        vidFrame(row,col,:) = CB(i3,:); 
                                    end
                                end                     
                            end
                        end
                        if Behavior == 1 && labelVidBehav == 1
                            vidFrame = insertText(vidFrame,[0,0],[B{ii,2},' ',num2str(B{ii,3})],'FontSize',18,'BoxColor','blue','BoxOpacity',0.4,'TextColor','white');
                        end
                        vidPose(:,:,:,ii) = vidFrame;
                    end
                end
                
                % Add coordinate data to data store
                for ii = 1:numTarget
                    allXY{2,ii} = [allXY{2,ii};XY(:,:,ii)];
                end            
            end
            
            % Write frames to new video file
            waitbar(iFrame/numFrame,wb,{'Writting video...';['Frame ',num2str(iFrame),'/',num2str(numFrame)]});
            for ii = 1:Expected
                if Sync == 1 && Stop == 0 
                    writeVideo(writerObj_Sync,vidSync(:,:,:,ii));
                end
                if Pose == 1 && labelVidPose == 1
                    writeVideo(writerObj_Pose,vidPose(:,:,:,ii));
                end
            end
        end
        if Sync == 1
            close(writerObj_Sync);
        end
        if Pose == 1 && labelVidPose == 1
            close(writerObj_Pose);
        end
        close(wb);
        
        % Compute euclidean distance for all target coordinates
        if Pose == 1
            for i = 1:numTarget
                allXY{2,i}(:,4) = 0;
                for ii = 1:size(allXY{2,i},1)
                    if ii > 1
                        allXY{2,i}(ii,4) = pdist([allXY{2,i}(ii-1,1:2);allXY{2,i}(ii,1:2)],'euclidean');
                    end
                end
            end
            save([vidFiles{I,1}(1:end-4),'_PoseData.mat'],'allXY');
        end
    end
    disp('Finished')
end
%%
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

% Datastore read function
function Out = load2VarLabeled(In)
     
    data = load(In);
    name = fieldnames(data);
    Out = data.(name{1,1});   
end

function lgraph = setLearnRate(lgraph,rate)

    numLayer = length(lgraph.Layers);
    for i = 1:numLayer
    
        layer = lgraph.Layers(i,1);                   
        try
            layer.WeightLearnRateFactor = rate;
            layer.BiasLearnRateFactor = rate;
        catch
            continue
        end
        lgraph = replaceLayer(lgraph,lgraph.Layers(i,1).Name,layer);
    end
end