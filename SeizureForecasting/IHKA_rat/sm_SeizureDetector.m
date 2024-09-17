
%get all 5.1 files
fils = getAllExtFiles('R:\DGregg\NeuralData\EDS_Cohort2','mp4',1);
fils = fils(contains(fils,'4.0'));

for i = 15:length(fils)
    try
        
        f   = getAllExtFiles(fileparts(fils{i}),'dat',1);
        kp = find(contains(f,'amplifier.dat'));
        fname = f{kp};
        if ~exist([fname(1:end-3) 'lfp']);
            fxml = [fname(1:end-3) 'xml'];
            if ~exist(fxml)
                fxml = 'R:\DGregg\NeuralData\EDS_Cohort2\Bipolar\10-24-2023(12.58)\RHS_231024_130000\amplifier.xml';
                
            end
            bz_LFPfromDat(fileparts(fname),'basename','amplifier','fxml',fxml);
            
        end
        i
    end
end


%%

fils = getAllExtFiles('R:\DGregg\NeuralData\EDS_Cohort2','mp4',1);
fils = fils(contains(fils,'4.0'));
ix=1;
clear fname
for i = 1:length(fils)
    
    f   = getAllExtFiles(fileparts(fils{i}),'lfp',1);
    kp = find(contains(f,'amplifier.lfp')& ~contains(f,'DoseResponse'));
    if any(kp)
        fname{ix} = f{kp};
        ix=ix+1;
    else
        f   = getAllExtFiles(fileparts(fils{i}),'dat',1);
        kp = find(contains(f,'amplifier.dat')& ~contains(f,'DoseResponse'));
        
        if length(kp)==1
            fxml = [f{kp}(1:end-3) 'xml'];
            if ~exist(fxml)
                fxml = 'R:\DGregg\NeuralData\EDS_Cohort2\Bipolar\10-24-2023(12.58)\RHS_231024_130000\amplifier.xml';
                
            end
            bz_LFPfromDat(fileparts(f{kp}),'basename','amplifier','fxml',fxml);
        end
    end
    
    
    
end



%%


% get all extreme events
k  = gaussian2Dfilter([10000 1],1250*5);
for i = 33
    cd(fileparts(fname{i}))
    dat  = LoadBinary(fname{i},'frequency',1250,'nchannels',32,'channels',6);
    filtsig = nanconvn(InstAmplitude(BandpassFilter(double(dat(:,1)),1250,[5 20])),k);
    
    tims = find(diff(zscore(filtsig)>3)>0)/1250
    fname{i}
    
    
end
    
%%

%autodetect EDS 5.1


%%
%%

seizfils = getAllExtFiles('R:\DGregg\NeuralData\EDS_Cohort2','szr',1);
ix = 1;clear gdses
for i = 1:length(seizfils)
    ev = LoadEvents(seizfils{i});
    
    subjkp = contains(ev.description(:,1),'4.0');
    on = contains(ev.description(:,1),'on');
    off = contains(ev.description(:,1),'ff');
    tim = [ev.time(subjkp&on) ev.time(subjkp&off)];
    
    kp = diff(tim,[],2)>5;
    
    if any(kp)
        gdses{ix,1} = seizfils{i};
        gdses{ix,2} = tim(kp,:);
        ix = ix+1;
    end
    
    %get all EDS 4.0 seisure
    
end


%%
% list of annotated files
annotations = [...
    {'R:\DGregg\NeuralData\EDS\12-23-2022(14.13)\RHD_221223_141500\manualDetect.evt.szr'          };...
    {'R:\DGregg\NeuralData\EDS\OL\Week1\1-3-2023(12.53)\RHD_230103_130000\manualDetect.evt.szr'   };...
    {'R:\DGregg\NeuralData\EDS\OL\Week1\12-28-2022(13.28)\RHD_221228_132904\manualDetect.evt.szr' };...
    {'R:\DGregg\NeuralData\EDS\OL\Week1\12-30-2022(13.1)\RHD_221230_130149\manualDetect.evt.szr'  };...
    {'R:\DGregg\NeuralData\EDS\OL\Week2\1-4-2023(12.50)\RHD_230104_130000\manualDetect.evt.szr'   };...
    {'R:\DGregg\NeuralData\EDS\OL\Week2\1-5-2023(12.47)\RHD_230105_130000\manualDetect.evt.szr'   };...
    {'R:\DGregg\NeuralData\EDS\OL\Week2\1-6-2023(12.58)\RHD_230106_130000\manualDetect.evt.szr'   };...
    {'R:\DGregg\NeuralData\EDS\Prophylactic3\4-13-2023(12.59)\RHS_230413_130000\amplifier.evt.szr'};...
    {'R:\DGregg\NeuralData\EDS\Prophylactic3\4-14-2023(12.59)\RHS_230414_130000\amplifier.evt.szr'};...
    
    ];

datfil = cellfun(@(a) [fileparts(a) filesep 'amplifier.lfp'],annotations,'uni',0);

training = 1:2:length(datfil);
test = 2:2:length(datfil);


nBin = 16;




Twindow = 2; % Time to analyse in seconds

% (let's use bz_LFPfromDat)
%lowCut = 4; % Lower frequency limit
%hiCut = 200.1; % Upper frequency limit

promptMessage = sprintf('Would you like to make a dataset, or update an exsiting dataset?');
titleBarCaption = 'settings';

ProjectPathName = uigetdir(cd,'Select folder for Datastore');
eval(['mkdir ',ProjectPathName,'\Training']);


TrecFileName = datfil(training);



VrecFileName = datfil(test);



cd(ProjectPathName)


%%


% list of annotated files
annotations = gdses(:,1);

datfil = cellfun(@(a) [fileparts(a) filesep 'amplifier.lfp'],annotations,'uni',0);



% define which subject
subjName = '4.0';
X =[]; Y=[];
%loop over trianing/test






recFileName = datfil;
numFile = length(recFileName);
allSeq = cell(numFile,0);
allTime = cell(numFile,1);




% loop over all files to build training data
numfreq = 10;

for i = 1:numFile
    fil = recFileName{i};
    
    %get meta data in xml
    xmlfil = [fil(1:end-3) 'xml'];
    
    
    if ~exist(xmlfil)
        xmlfil = 'R:\DGregg\NeuralData\EDS_Cohort2\Bipolar\10-24-2023(12.58)\RHS_231024_130000\amplifier.xml';
        
    end
    
    anfil =  annotations(ismember(fileparts(annotations),fileparts(fil)));
    
    % get annotations
    
    ev = LoadEvents(anfil{1});
    disp(['Reading file ', num2str(i),'/',num2str(length(recFileName)),'...']);
    
    
    
    
    
    
    % get channel ID
    channels = 1:8;
    
    
    xmlData = LoadXml(xmlfil);
    
    Fs = xmlData.lfpSampleRate; % Sample frequency
    numChan = xmlData.nChannels;  % Number of channels
    
    % load data
    EEGdata = LoadBinary(fil,'nchannels',numChan,'frequency',Fs,'channels',channels);
    
    
    ts = (1:size(EEGdata,1))/Fs;
    
    % make Twindow cell arrays with and without seizures (this period
    % is taken just before each seizure)
    
    kp_on = contains(ev.description,subjName) & contains(ev.description,'sz_on');
    kp_off = contains(ev.description,subjName) & contains(ev.description,'sz_off');
    
    
    sz_on = ev.time(kp_on);
    sz_off = ev.time(kp_off);
    
    
    
    % put all channels in a cell array
    for ii = 1:size(X1,1)
        Xt{ii,1} = cell2mat(X1(ii,:))';
        
        
    end
    ops.freqs = logspace(log(2),log10(300),numfreq);
    
    
    
    %get spectrogram
    for ii = 1:size(Xt,1)
        nCh = size(Xt{ii},1);
        Xtt{ii} = nan(numfreq*nCh,nBin);
        for j = 1:size(Xt{ii},1)
            tmp = abs(awt_freqlist(Xt{ii}(j,:),Fs,ops.freqs))';
            tmp1 = nan(size(tmp,1),nBin);
            for jj = 1:size(tmp,1)
                t = avghist(1:size(tmp(jj,:),2),tmp(jj,:),linspace(1,nSamp,nBin+1));
                tmp1(jj,:)  = t(1:end-1);
            end
            
            
            %resample
            
            if any(isnan(tmp1(:)))
                error('here')
            end
            Xtt{ii}(j:nCh:end,:) = tmp1;
        end
        
    end
    X = [X;Xtt'];
    Y = [Y;Yt];
    i
end




%% set up ANN

% this code shows an error with unconnected layers


useLSTM = false; % is it?
numBlocks = 3;
numFilters = 16;
filterSize = 16;
dilationFactor = 1;
dropoutFactor = 0.005;

layer = [
    sequenceInputLayer(nCh*numfreq,Normalization="none",MinLength=nBin,Name="input")];
lgraphCNN = layerGraph(layer);
outputName = layer(end).Name;

for i = 1:numBlocks
    for ii = 1:length(filterSize)
        
        % convolutional block
        layers = [
            convolution1dLayer(filterSize(1,ii),numFilters,DilationFactor=dilationFactor,Padding="causal",Name="conv1_"+"f"+filterSize(1,ii)+"_"+i)
            layerNormalizationLayer(Name="lNorm1_"+"f"+filterSize(1,ii)+"_"+i)
            dropoutLayer(dropoutFactor,Name="dropOut1_"+"f"+filterSize(1,ii)+"_"+i)
            convolution1dLayer(filterSize(1,ii),numFilters,DilationFactor=dilationFactor,Padding="causal",Name="conv2_"+"f"+filterSize(1,ii)+"_"+i)
            layerNormalizationLayer(Name="lNorm2_"+"f"+filterSize(1,ii)+"_"+i)
            reluLayer(Name="relu_"+"f"+filterSize(1,ii)+"_"+i)
            dropoutLayer(dropoutFactor,Name="dropOut2_"+"f"+filterSize(1,ii)+"_"+i)];
        
        % Add and connect layers.
        lgraphCNN = addLayers(lgraphCNN,layers);
        lgraphCNN = connectLayers(lgraphCNN,outputName,"conv1_"+"f"+filterSize(1,ii)+"_"+i);
    end
    
    % Add addtion and attention layers and connect convovultional blocks
    layers = [
        additionLayer(length(filterSize)+1+useLSTM,Name="add_"+i)
        maxPooling1dLayer(2,"Stride",2,Name="maxPool_"+i)];
    lgraphCNN = addLayers(lgraphCNN,layers);
    for ii = 1:length(filterSize)
        lgraphCNN = connectLayers(lgraphCNN,"dropOut2_"+"f"+filterSize(1,ii)+"_"+i,"add_"+i+"/in"+ii);
    end
    
    % Skip connection.
    layer = convolution1dLayer(1,numFilters,Name="convSkip"+i);
    
    lgraphCNN = addLayers(lgraphCNN,layer);
    lgraphCNN = connectLayers(lgraphCNN,outputName,"convSkip"+i);
    lgraphCNN = connectLayers(lgraphCNN,"convSkip"+i,"add_"+i+"/in"+(length(filterSize)+1));
    
    % Update layer output name.
    outputName = "maxPool_"+i;
    numFilters = numFilters*2;
    %dilationFactor = dilationFactor*2;
end

% Add output layers
layers = [
    globalAveragePooling1dLayer('Name','globalPool')
    fullyConnectedLayer(numFilters,'Name','fc')
    fullyConnectedLayer(2,'Name','fc_class')
    softmaxLayer('Name','softMax')
    classificationLayer('Name','Classify')];
lgraphCNN = addLayers(lgraphCNN,layers);
lgraphCNN = connectLayers(lgraphCNN,outputName,'globalPool');

analyzeNetwork(lgraphCNN)

%%

idx = randperm(length(X),length(X))';
X1 = X(idx,1);
Y1 = Y(idx,1);

numVal = floor(length(X1)*.05);
Xval = cell(numVal,1);

clear YVal
Xval(:,1) = X1(1:numVal,1);
Yval = Y1(1:numVal,1);

X1(1:numVal,:) = [];
Y1(1:numVal,:) = [];

MBS = 10;
iLR = 0.001;
numEpoch = 4;
numObs = length(X1);

options = trainingOptions('sgdm', ...
    'MiniBatchSize',MBS, ...
    'MaxEpochs',numEpoch, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',iLR, ...
    'Shuffle','every-epoch', ...
    'ValidationData',{Xval,Yval}, ...
    'ValidationFrequency', 1, ...
    'Plots','training-progress', ...
    'Verbose',true, ...
    'VerboseFrequency',1, ...
    'ExecutionEnvironment','gpu');

netCNN = trainNetwork(X1,Y1,lgraphCNN,options);
lgraphCNN = layerGraph(netCNN);%
reset(gpuDevice(1))

YPred = predict(netCNN,Xval,'ExecutionEnvironment','gpu');


save('R:\DGregg\SeizureForecast\Seizuredetect_demo\Training\netCNN2.mat','netCNN','-v7.3')

clear X Y

%%




%define the groups labels
group = Y;
training = cell2mat(cellfun(@(a) nanmean(a,2)',X,'uni',0));

% randomly sort for cross validation
ops.rix  = randsample(1:length(group),length(group));
training = training(ops.rix,:);
group = group(ops.rix);

%define subset of data to train (odd sessions)
kp = mod(1:length(group),2)==0;


%set up classifer
ops.N = round(size(training,1));         % Number of observations in the training sample
ops.t = templateTree('MaxNumSplits',ops.N/10);
ops.NumLearningCycles = 500;
ops.Learners = ops.t;
ops.LearnRate = 0.1;
ops.Method = 'RUSBoost';


% train model
rusTree = fitcensemble(training,group,'Method',ops.Method, ...
    'NumLearningCycles',ops.NumLearningCycles,'Learners',ops.Learners,'LearnRate',ops.LearnRate,'ScoreTransform','logit');

%get prediction on held out data

[pred,score] = predict(rusTree,training(~kp,:));
save('R:\DGregg\SeizureForecast\Seizuredetect_demo\Training\randomForest4.0_wstim.mat','rusTree','-v7.3')
% SAM STOPPED HERE 9/5/2023


%%

fils = getAllExtFiles('R:\DGregg\NeuralData\EDS_Cohort2','mp4',1);
fils = fils(contains(fils,'4.0'));
ix=1;
for i = 1:length(fils)
    try
        
        f   = getAllExtFiles(fileparts(fils{i}),'lfp',1);
        sz_detect{ix} = f{1};
        ix = ix+1;
    end
end

%%
ops.numfreq  = [4.9334 7.7869 12.2910 19.4002 30.6214 48.3330 76.2892 120.4155 190.0649 300.0000]; % frequencies for wavelet decomp.
ops.nchanFil= 32; % number of channel in dat file
ops.nChanSubj = 8; % # channels for subj
ops.channels = 1:8; % channels for subject
ops.Fs =  1250; % sampling rate
ops.nTempBin =  32; % ?
ops.data =  15; % ?
ops.Twin =  2; % window size
ops.feature_fun = @sm_GetDataFeature2; % function for featurizing the data
model_fname  = 'R:\DGregg\SeizureForecast\Seizuredetect_demo\Training\randomForest4.0_wstim.mat';
v = load(model_fname);

for i = 1:length(sz_detect)
    
    
    
    %for retrospective
    fil = sz_detect{i};
    tmp = dir(fil);
    
    filedur = tmp.bytes/2/ops.nchanFil/ops.Fs;
    
    
    
    %get sz prob in 1s increments
    clear tim
    
    Pr = nan(ceil(filedur-ops.Twin/.1),1);
    ix = 1;
    Fstest = 1;
    for j = 1:(1/Fstest):filedur-ops.Twin
        tic
        data = LoadBinary(fil,'nchannels',ops.nchanFil,'frequency',ops.Fs,'channels',ops.channels,'duration',2 ,'start',j);
        Prob = sm_getSeizProb(data,v.rusTree,ops);
        Pr(ix) = Prob(1);
        tim(ix) = toc;
        ix = ix+1
        
    end
    
    
    % save auto-detected seizures
    sz_ep = [find(diff([0;Pr]==2)>0) find(diff([0;Pr]==2)<0)]/Fstest;
    ev.time  = [sz_ep(:,1);sz_ep(:,2)];
    ev.description = [repmat({'sz_on'},size(sz_ep,1),1);repmat({'sz_off'},size(sz_ep,1),1)];
    
    [~,b] = sort(ev.time(:,1));
    ev.time = ev.time(b,:);
    ev.description = ev.description(b);
    filename = '\autodetect.evt.szr';
    SaveEvents(filename,ev);
end



%%

% test on held out data



% define which subject
subjName = '2.3';
X =[]; Y=[];







recFileName = VrecFileName;
numFile = length(VrecFileName);
allSeq = cell(numFile,0);
allTime = cell(numFile,1);




% loop over all files to build training data
numfreq = 10;

for i = 1:numFile
    fil = recFileName{i};
    
    %get meta data in xml
    xmlfil = [fil(1:end-3) 'xml'];
    
    anfil =  annotations(ismember(fileparts(annotations),fileparts(fil)));
    
    % get annotations
    
    ev = LoadEvents(anfil{1});
    disp(['Reading file ', num2str(i),'/',num2str(length(recFileName)),'...']);
    
    
    
    
    
    
    % get channel ID
    matfil = getAllExtFiles(fileparts(fil),'mat',0);
    kp = contains(matfil,subjName);
    gdfil = matfil(kp);
    if ~isempty(gdfil)
        
        v = load(gdfil{1});
        channels = v.sessiondata.channelID;
    else
        error('missing session info')
    end
    
    
    xmlData = LoadXml(xmlfil);
    
    Fs = xmlData.lfpSampleRate; % Sample frequency
    numChan = xmlData.nChannels;  % Number of channels
    
    % load data
    EEGdata = LoadBinary(fil,'nchannels',numChan,'frequency',Fs,'channels',channels);
    
    
    ts = (1:size(EEGdata,1))/Fs;
    
    % make Twindow cell arrays with and without seizures (this period
    % is taken just before each seizure)
    
    kp_on = contains(ev.description,subjName) & contains(ev.description,'sz_on');
    kp_off = contains(ev.description,subjName) & contains(ev.description,'sz_off');
    
    
    sz_on = ev.time(kp_on);
    sz_off = ev.time(kp_off);
    
    
    
    newEp = [sz_on-10 sz_off+10];
    
    newEp = cell2mat(cellfun(@(a) [(a(1):Twindow:a(2)-Twindow)' (a(1)+Twindow:Twindow:a(2))'],num2cell(newEp,2),'uni',0));
    newEp(:,1) = newEp(:,1)+.001;
    newEp = MergeEpochs2(newEp);
    
    kp = diff(newEp,[],2)<=Twindow;
    newEp = newEp(kp,:);
    seizureSeq = InIntervals(ts,[sz_on sz_off]);
    
    %this builds the cell array
    [~,Yt] = Epoch2Mat(ts,newEp,seizureSeq);
    Yt = cellfun(@any,Yt);
    Yt =  categorical(Yt);
    clear X1
    for ii = 1:size(EEGdata,2)
        [~,tmp] = Epoch2Mat(ts,newEp,EEGdata(:,ii));
        
        % resample to be 1024
        
        ts1 = (1:length(tmp{1}))'/Fs;
        Fs2 = range(ts1)/1024;
        ts2= ts1(1):Fs2:(ts1(end) - Fs2);
        nSamp = length(ts1);
        % tmp = cellfun(@(a) nanPad(a',nSamp)',tmp,'uni',0);
        
        %X1(:,ii) = cellfun(@(a) interp1(ts1,double(a),ts2),tmp,'UniformOutput',0);
        X1(:,ii) = tmp;
    end
    
    
    
    clear Xt Xtt
    
    % put all channels in a cell array
    for ii = 1:size(X1,1)
        Xt{ii,1} = cell2mat(X1(ii,:))';
        
        
    end
    
    %get spectrogram
    for ii = 1:size(Xt,1)
        nCh = size(Xt{ii},1);
        Xtt{ii} = nan(numfreq*nCh,nBin);
        for j = 1:size(Xt{ii},1)
            tmp = abs(awt_freqlist(Xt{ii}(j,:),Fs,logspace(log(2),log10(300),numfreq)))';
            tmp1 = nan(size(tmp,1),nBin);
            for jj = 1:size(tmp,1)
                t = avghist(1:size(tmp(jj,:),2),tmp(jj,:),linspace(1,nSamp,nBin+1));
                tmp1(jj,:)  = t(1:end-1);
            end
            
            
            %resample
            
            if any(isnan(tmp1(:)))
                error('here')
            end
            Xtt{ii}(j:nCh:end,:) = tmp1;
        end
        
    end
    X = [X;Xtt'];
    Y = [Y;Yt];
    i
end



%%






% idx = randperm(length(X),300000)';
%
% if I == 1
%     allX_T = [allX_T;X(idx,1)]; %#ok<AGROW>
%     allY_T = [allY_T;Y(idx,1)]; %#ok<AGROW>
% else
%     allX_V = [allX_V;X(idx,1)]; %#ok<AGROW>
%     allY_V = [allY_V;Y(idx,1)]; %#ok<AGROW>
% end



%save('trainingData.mat','allX_T','allY_T','allX_V','allY_V-v7.3')

%%
%
%
% i = 1;
%
% clf
% subplot(2,1,1)
% plot(Xval{i,1})
% subplot(2,1,2)
% plot(YPred{i,1})
% hold on
% plot(Yval{i,1})
% ylim([-0.2,1.2])
% i = i+1;
%
% %%
%
% L = 5000;
% seizureSeq(:,3) = 0;
% for i = 1:L:size(EEGdata,1)
%
%     YPred = clip(predict(netCNN,EEGdata(i:i+L-1,1)','ExecutionEnvironment','gpu'),0,1);
%     seizureSeq(i:i+L-1,3) = YPred;
% end
%
% clf
% plot(seizureSeq(1:500:end,3))
% hold on
% plot(seizureSeq(1:500:end,1),'Color','k')
%
% idx = find(seizureSeq(:,3) >= 0.95);
% idxDiff = diff(idx);
% maxSkip = 1000;
% startP = [0,0];
% count = 1;
% track = 0;
% for i = 1:length(idx)-1
%
%     if track == 0
%         startP(count,1) = idx(i,1);
%     end
%
%     if idxDiff(i) <= maxSkip
%         track = track+1;
%     elseif idxDiff(i) > maxSkip
%         if track >= 500
%             startP(count,2) = idx(i,1);
%             count = count+1;
%         end
%         track = 0;
%     end
% end
% startP(startP(:,2) == 0,:) = [];
%
% seizureSeq(:,4) = 0;
% for i = 1:size(startP,1)
%     seizureSeq(startP(i,1):startP(i,2),4) = 1;
% end
%
% clf
% plot(seizureSeq(1:500:end,4))
% hold on
% plot(seizureSeq(1:500:end,1)*.8,'Color','k')
%
% i = 1;
% L = 3000;
%
% subplot(2,1,1)
% plot(seizureSeq(startP(i,1)-L:startP(i,2)+L,3)')
% axis tight
% subplot(2,1,2)
% plot(EEGdata(startP(i,1):startP(i,2)+L,1)')
% axis tight
% ylim([-0.5,0.5])
% i = i+1;
%
%
%
% %%
%
% Fs = 20000; % Sample frequency
% DsF = 500; % Downsample frequency
%
% numChan = 32;
% useChan = 1;
%
% [file,path] = uigetfile('.dat');
% cd(path)
% s = dir(path);
% numSamples = s(3).bytes/(2*numChan);
%
% fileID = fopen([path,'\amplifier.dat']);
% fseek(fileID,0,'bof');
%
% X = nan(floor(numSamples/(Fs/DsF)),1);
%
% sT = 10;
% count = 1;
% wb = waitbar(0,'Reading file...');
% for i = 1:Fs*sT:numSamples
%
%     data = fread(fileID,[numChan,(sT*Fs)],'int16');
%     data = convLPFilt(data(useChan,:)',Fs,DsF,[]);
%     try
%         X(count:count+(sT*DsF)-1) = data;
%     catch
%         X = [X;data]; %#ok<AGROW>
%         break
%     end
%     count = count+(sT*DsF);
%     waitbar(i/numSamples,wb,'Reading file...');
% end
% close(wb);
% fclose(fileID);
% X(isnan(X),:) = [];
%
% M = round(std(X)*30);
% X2 = X/M;
%
% L = 5000;
% X2(:,2) = 0;
% wb = waitbar(0,'Making predictions');
% for i = 1:L:length(X2)
%
%     try
%         YPred = clip(predict(netCNN,X2(i:i+L-1,1)','ExecutionEnvironment','gpu'),0,1);
%         X2(i:i+L-1,2) = YPred;
%     catch
%         YPred = clip(predict(netCNN,X2(i:end,1)','ExecutionEnvironment','gpu'),0,1);
%         X2(i:i+length(YPred)-1,2) = YPred;
%     end
%     waitbar(i/length(X2),wb,'Making predictions');
% end
% close(wb);
%
% idx = find(X2(:,2) >= 0.95);
% idxDiff = diff(idx);
% maxSkip = 1000;
% startP = [0,0];
% count = 1;
% track = 0;
% for i = 1:length(idx)-1
%
%     if track == 0
%         startP(count,1) = idx(i,1);
%     end
%
%     if idxDiff(i) <= maxSkip
%         track = track+1;
%     elseif idxDiff(i) > maxSkip
%         if track >= 1000
%             if idx(i,1)-startP(count,1) >= maxSkip
%                 startP(count,2) = idx(i,1);
%                 count = count+1;
%             else
%                 startP(count,:) = [];
%             end
%         end
%         track = 0;
%     end
% end
% startP(startP(:,2) == 0,:) = [];
%
% for i = 1:size(startP,1)
%
%     x = X2(startP(i,1):startP(i,2),2);
%     startP(i,3) = AOC(x);
% end
%
% X2(:,3) = 0;
% for i = 1:size(startP,1)
%     if startP(i,2)-startP(i,1) >= 1000 && startP(i,3) > 10000
%         X2(startP(i,1):startP(i,2),3) = 1;
%     end
% end
%
% clf
% plot(X2(1:500:end,3))
% hold on
% plot(X2(1:500:end,1),'Color','k')
%
% if size(startP,1) > 0
%     i = 1;
%     L = 4000;
%
%     figure
%     subplot(2,1,1)
%     plot(X2(startP(i,1)-L:startP(i,2)+L,2)')
%     axis tight
%     subplot(2,1,2)
%     plot(X2(startP(i,1)-L:startP(i,2)+L,1)')
%     axis tight
%     ylim([-0.5,0.5])
%     i = i+1;
%
%     TS4NS = cell(size(startP(startP(:,3) >= 10000),1)*2,2);
%     count = 1;
%     for i = 1:size(startP,1)
%
%         if startP(i,3) >= 10000
%
%             TS4NS{count,1} = startP(i,1)*(1000/DsF);
%             TS4NS{count,2} = seizureStart;
%             TS4NS{count+1,1} = startP(i,2)*(1000/DsF);
%             TS4NS{count+1,2} = seizureStop;
%             count = count+2;
%         end
%     end
% end
%
%





