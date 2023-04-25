clear all %#ok<CLALL> 

StartTime = [13,0]; % Time to begin recording in 24hr time, format = [hr,min]

% Block schedule 
firstBL = 120; % Pre-stim baseline duration in minutes
stimDur = 60; % Duration of stim period in minutes
finalBL = 120; % Post-stim baseline duration in minutes

% Stim parameters
stimI = 100; % Stim 1 intensities in uA, randomly selected
thresholdTD = 1; % Theda/delta threshold for stim
ISI = 8:12; % Inter stim interval in seconds, randomly selected.

% Intan channel ID
chanID(1,1:2) = [{'000'},{'PrL(L)'}];
chanID(2,1:2) = [{'001'},{'PrL(R)'}];
chanID(3,1:2) = [{'002'},{'AVT(L)'}];
chanID(4,1:2) = [{'003'},{'BLA(R)'}];
chanID(5,1:2) = [{'004'},{'CA1(L)'}];
chanID(6,1:2) = [{'005'},{'CA1(R)'}];
chanID(7,1:2) = [{'006'},{'GRS(L)'}];
chanID(8,1:2) = [{'007'},{'LTD(L)'}];

% Intan TCP
intanIP = '127.0.0.1'; % IP for RHS.
intanPort1 = 5000; % Port for RHS commands.
intanPort2 = 5001; % Port for RHS data output.

dataChan(1,1) = 0; % Port A measure channel
dataChan(2,1) = 6; % Port B measure channel
dataChan(3,1) = 0; % Port C measure channel
dataChan(4,1) = 0; % Port D measure channel

stimChan(1,1) = NaN; % Port A stim channel
stimChan(2,1) = 8; % Port B stim channel
stimChan(3,1) = NaN; % Port C stim channel
stimChan(4,1) = NaN; % Port D stim channel

framesPerBlock = 128;
blocksPerRead = 125;
numChan = 1;
numBandsPerChan = 1; % Amplifier channels can be displayed as LOW, WIDE, or HIGH

% Arduino Settings
ArdPin{1,1} = 'D13'; % Arduino pin for LED.
ArdPin{2,1} = 'D8'; % Arduino pin for Timestamp.
LEDtime = 10; % Time between LED blinks is seconds

% Camera Settings
fps = 25;

% Analysis Settings
sF = 20000; % Recording sample frequency
dSF = 1000; % Downsample frequency
stimTrig = 1;
Pre = 0;
Post = ISI(1,1);

activeChan(1,1:3) = [{'000'},{'PrL(L)'},{1}];
activeChan(2,1:3) = [{'001'},{'PrL(R)'},{1}];
activeChan(3,1:3) = [{'002'},{'AVT(L)'},{1}];
activeChan(4,1:3) = [{'003'},{'BLA(R)'},{1}];
activeChan(5,1:3) = [{'004'},{'CA1(L)'},{1}];
activeChan(6,1:3) = [{'005'},{'CA1(R)'},{1}];
activeChan(7,1:3) = [{'006'},{'GRS(L)'},{1}];
activeChan(8,1:3) = [{'007'},{'LTD(L)'},{1}];

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
    prompt = 'Label WebCam';
    dlgtitle = 'Input';
    WL(i,2) = inputdlg(prompt,dlgtitle,[1,40],WL(i,1));
end

[camChoice,~] = listdlg('PromptString','Select cameras.','ListString',WL(:,2),'SelectionMode','multiple');
numCam = length(camChoice);

Res = cam.AvailableResolutions;
[resChoice,~] = listdlg('PromptString','Select resolution.','ListString',Res,'SelectionMode','single');
Dim = symsepchar(Res{1,resChoice},'x');
Dim = [str2double(Dim{1,2}),str2double(Dim{1,1})];
clear cam

PathName = uigetdir(cd,'Select folder for Datastore');
cd(PathName);
%%
global dataFig videoFig %#ok<NUSED,GVMIS> 
matObj = matfile([tempdir,'Go.mat'],'Writable',true);

Now = clock;
Time = [num2str(Now(1,2)),'-',num2str(Now(1,3)),'-',num2str(Now(1,1)),'(',num2str(Now(1,4)),'.',num2str(Now(1,5)),')'];

PathName = [cd,'\',Time,'_',chanID{stimChan(stimChan > 0),2},'_RHS'];
eval(['mkdir ',PathName]);

vStartTime = (StartTime(1,1)*60)+StartTime(1,2);

hardwareTimer = 1;
intanWriter = 2;
intanReader = 3;
dataHandler = 4;
camReader = 5;
camSaver = 6;
        
% Calculations for accurate parsing
numDsFrames = (framesPerBlock*blocksPerRead)/(sF/dSF);
readTime = numDsFrames/dSF;
numAmplifierBands = numBandsPerChan*numChan;
waveformBytesPerFrame = 4+2*numAmplifierBands;
waveformBytesPerBlock = framesPerBlock*waveformBytesPerFrame+4;
framesPerRead = framesPerBlock*blocksPerRead;
waveformBytesBlocks = blocksPerRead*waveformBytesPerBlock;

disp(['Will start at ',num2str(StartTime(1,1)),':',num2str(StartTime(1,2))]); 
while 1    
    if getTime >= vStartTime-1
       break
    end
end

disp('Connecting to hardware...');

textOut = parallel.pool.DataQueue;
afterEach(textOut,@disp);

videoOut = parallel.pool.DataQueue;
afterEach(videoOut,@showVideo);

dataOut = parallel.pool.DataQueue;
afterEach(dataOut,@plotTD);

if isempty(gcp('nocreate')) == 1
    parpool('local',6);
end 

spmd (6)
    if labindex == hardwareTimer        
        while 1
            try
                Ard = arduino('com3','uno');                          
                break
            catch
                disp('Issue connecting to Arduino')
                pause(3)
            end   
        end                  

        matObj.Go = 1; 

        % sync with system clock
        send(textOut,'Syncing hardware');
        while 1
            if getTime >= vStartTime
                break
            end
        end   

        labBarrier
        send(textOut,'Pre-stime baseline');
        
        Count = 0;
        waitT1 = 1/60;
        holdT1 = 0;
        ck = [0,0];
        seshStart = getTime;
        labSend(seshStart,intanWriter);
        labSend(seshStart,intanReader);
        labSend(seshStart,dataHandler);
        while 1 
            while 1
                t = getTime;
                if rem(t,waitT1) == 0 && t > holdT1
                    holdT1 = t;
                    break
                end
            end
            Count = Count+1;             

            if Count == 1                  
                writeDigitalPin(Ard,ArdPin{1,1},1);                
            elseif Count == LEDtime
                Count = 0;
            end       
            
            if t >= seshStart+firstBL && ck(1,1) == 0
                writeDigitalPin(Ard,ArdPin{2,1},1);
                send(textOut,'Begin stim');
                ck(1,1) = 1;
            elseif t >= seshStart+firstBL+stimDur && ck(1,2) == 0
                writeDigitalPin(Ard,ArdPin{2,1},1);                
                send(textOut,'Post-stime baseline');
                ck(1,2) = 1;
            elseif t >= seshStart+firstBL+stimDur+finalBL
                send(textOut,'Finished');
                break
            end

            pause(0.2)
            writeDigitalPin(Ard,ArdPin{1,1},0);
            writeDigitalPin(Ard,ArdPin{2,1},0);
        end
        
        matObj.Go = 0;

    elseif labindex == intanWriter
        while 1
            try
                intanSend = tcpclient(intanIP,intanPort1);
                write(intanSend,uint8(['set Filename.Path ',PathName]));
                write(intanSend,uint8('set FileFormat OneFilePerSignalType'));
                break
            catch
                disp('Issue connecting to Intan commands')
                pause(3)
            end
        end        

        stim1 = ['(s1-',num2str(stimI),'-'];    
        write(intanSend,uint8(['set Filename.BaseFilename ',stim1(1,1:end-1),')uA_']));        
       
        for ii = 1:size(stimChan,1)
            if ii == 1
                port = 'a';
            elseif ii == 2
                port = 'b';
            elseif ii == 3
                port = 'c';
            elseif ii == 4
                port = 'd';
            end   

            if ~isnan(stimChan(ii,1))    
                if dataChan(ii,1) > 0
                    write(intanSend,uint8(['set ',port,'-',chanID{dataChan(ii,1),1},'.tcpdataoutputenabled true;']));
                end
                
                if stimChan(ii,1) > 0

                    write(intanSend,uint8(['set ',port,'-',chanID{stimChan(ii,1),1},'.stimenabled true']));
                    write(intanSend,uint8(['set ',port,'-',chanID{stimChan(ii,1),1},'.source KeyPressF1']));
                    write(intanSend,uint8(['set ',port,'-',chanID{stimChan(ii,1),1},'.FirstPhaseAmplitudeMicroAmps ',num2str(stimI)]));
                    write(intanSend,uint8(['set ',port,'-',chanID{stimChan(ii,1),1},'.SecondPhaseAmplitudeMicroAmps ',num2str(stimI)]));
                else
                    write(intanSend,uint8(['set ',port,'-',chanID{stimChan(ii,1),1},'.stimenabled false']));
                end
                write(intanSend,uint8(['execute uploadstimparameters ',port,'-',chanID{stimChan(ii,1),1}]));
                send(textOut,['Stim amp set to ',num2str(stimI),'uA.']) 
                pause(5)
            end                    
        end        

        labBarrier
        seshStart = labReceive(hardwareTimer);

        write(intanSend,uint8('set runMode record'));

        while getTime <= seshStart+firstBL
            pause(0.001)
        end

        while 1
            stim = labReceive(dataHandler);
            if stim == 1
                write(intanSend,uint8('execute ManualStimTriggerPulse F1'));
            elseif stim == -1
                break
            end
        end

        while getTime <= seshStart+firstBL+stimDur+finalBL
            pause(0.001)
        end        

        write(intanSend,uint8('set runMode stop')); 
    elseif labindex == intanReader        
        while 1
            try
                intanRead = tcpclient(intanIP,intanPort2);
                break
            catch
                disp('Issue connecting to Intan data output')
                pause(3)
            end   
        end

        flush(intanRead);

        % Pre-allocate memory for blocks of waveform data (the amount that's
        % plotted at once)
        amplifierData = 32768*ones(numChan,framesPerRead);
        ampTimewindow = zeros(1,numDsFrames*5);
    
        % Initialize counters
        chunkCount = 0;
        blockCount = 0;
        amplifierTimestampsIndex = 1;
        Stopper = 0;

        labBarrier
        seshStart = labReceive(hardwareTimer);

        while getTime <= seshStart+firstBL
            flush(intanRead);
            pause(0.001)
        end
        
        flush(intanRead);
        while getTime <= seshStart+firstBL+stimDur

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

                %Scale these data blocks and downsample
                amplifierData = downsample(0.195*(amplifierData-32768),sF/dSF);

                %Shift time window to next step
                ampTimewindow(:,1:numDsFrames) = [];
                ampTimewindow = [ampTimewindow,amplifierData]; %#ok<AGROW>

                labSend(ampTimewindow,dataHandler);          
                
                %Reset index
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

        while getTime <= seshStart+firstBL+stimDur+finalBL
            flush(intanRead);
            pause(0.001)
        end   
    elseif labindex == dataHandler
        
        bankISI = (round(ISI(1,1)/readTime)):round((ISI(1,end)/readTime));
        nextISI = bankISI(1,1);
        readCount = 0;

        labBarrier
        seshStart = labReceive(hardwareTimer);

        while getTime <= seshStart+firstBL
            pause(0.001)
        end

        while getTime <= seshStart+firstBL+stimDur

            ampTimewindow = labReceive(intanReader);
            send(dataOut,[{ampTimewindow};{dSF}]); 
            readCount = readCount+1;

            if readCount >= nextISI

                % Calculate theta/delta ratio
                theta = bandpass(ampTimewindow,[5,12],dSF);
                delta = bandpass(ampTimewindow,[1,4],dSF);
                thetaAmp = abs(hilbert(theta));
                deltaAmp = abs(hilbert(delta));
                theta_delta = mean(thetaAmp(1,end-dSF+1:end))/mean(deltaAmp(1,end-dSF+1:end));

                if theta_delta <= thresholdTD
                    labSend(1,intanWriter);
                    send(textOut,['Stim @ ',num2str(theta_delta),' ISI: ',num2str(readCount)]);
                    nextISI = datasample(bankISI,1);
                    readCount = 0;
                else
                    labSend(0,intanWriter);
                end
            end
        end
        labSend(-1,intanWriter);

        while getTime <= seshStart+firstBL+stimDur+finalBL
            pause(0.001)
        end     

    elseif labindex == camReader 
        
        % Create all webcam objects           
        cam = cell(numCam,1); 
        for ii = 1:numCam                
            cam{ii,1} = webcam(camChoice(1,ii));
        end           
        frames = uint8(zeros(Dim(1,1),Dim(1,2),3,numCam));     

        labBarrier
        
        skipFrame = 0;
        waitT2 = (1/60)/fps;
        holdT2 = 0;
        while matObj.Go == 1
            while 1
                t2 = getTime;
                if rem(t2,waitT2) == 0 && t2 > holdT2
                    holdT2 = t2;
                    break
                end
            end

            for ii = 1:numCam
                frames(:,:,:,ii) = snapshot(cam{ii,1});
            end

            labSend(frames,camSaver);

            % Send video to client for on-screen display
            if skipFrame == 0
                send(videoOut,frames);
                skipFrame = 1;
            else
                skipFrame = 0;
            end    
        end 
        labSend([],camSaver);
    elseif labindex == camSaver       
        
        VidFileName = cell(numCam,1);
        writerObj = cell(numCam,1);

        % Create all video writer objects
        for ii = 1:numCam   
         
            stim1 = ['(s1-',num2str(stimI),'-'];           
            VidFileName{ii,1} = [WL{camChoice(1,ii),2},'_',stim1(1,1:end-1),')uA.mp4'];

            writerObj{ii,1} = VideoWriter([PathName,'\',VidFileName{ii,1}],'MPEG-4'); %#ok<TNMLP>
            writerObj{ii,1}.FrameRate = fps;
            open(writerObj{ii,1});
        end

        videos = uint8(zeros(Dim(1,1),Dim(1,2),3,fps,numCam));

        labBarrier
        while 1
            for ii = 1:fps  
                frames = labReceive(camReader);
                if isempty(frames)
                    break
                end
                for iii = 1:numCam
                    videos(:,:,:,ii,iii) = frames(:,:,:,iii);
                end
            end
            if isempty(frames)
                break
            end           
            for ii = 1:numCam
                for iii = 1:fps
                    writeVideo(writerObj{ii,1},videos(:,:,:,iii,ii));
                end
            end
        end
        for ii = 1:numCam
             close(writerObj{ii,1});
        end
    end
end
close(gcf)

delete(gcp('nocreate'))
disp('Finsished');

%%

disp('Starting analysis');

listing = dir(PathName);
dirRemove = false(length({listing.name}),1);
for i = 1:length({listing.name})
    if strcmp(listing(i).name,'.') == 1 || strcmp(listing(i).name,'..') == 1
        dirRemove(i,1) = true;
    end
end
listing(dirRemove) = [];

recInfo = [listing.isdir]';
recInfo = {listing(recInfo).name}';
numFolders = length(recInfo);
dirRemove = false(numFolders,1);
for i = 1:numFolders
        
    if contains(recInfo{i,1},'s1') == 1 && contains(recInfo{i,1},'s2') == 0

        loc1 = strfind(recInfo{i,1},'s1-');
        loc2 = strfind(recInfo{i,1},')uA');
        recInfo{i,2} = str2double(recInfo{i,1}(1,loc1+3:loc2-1));
        recInfo{i,3} = NaN;
    elseif contains(recInfo{i,1},'s1') == 1 && contains(recInfo{i,1},'s2') == 1
        
        loc1 = strfind(recInfo{i,1},'s1-');
        loc2 = strfind(recInfo{i,1},'-s2');
        recInfo{i,2} = str2double(recInfo{i,1}(1,loc1+3:loc2-1));

        loc1 = strfind(recInfo{i,1},'s2-');
        loc2 = strfind(recInfo{i,1},')uA');
        recInfo{i,3} = str2double(recInfo{i,1}(1,loc1+3:loc2-1));
    else
        recInfo(i,2:3) = {NaN};
        dirRemove(i,1) = true;
    end    
    recInfo{i,1} = [listing(1).folder,'\',recInfo{i,1}];
end
recInfo(dirRemove,:) = [];
recInfo = sortrows(recInfo,2);
numFolders = size(recInfo,1);

recInfo = [{'File'},{'Stim1 amp (uI)'},{'Stim2 amp (uI)'},{'Stim start and stop'},{'Theta/delta'};recInfo,{[]},{[]}];

for i = 2:numFolders+1

    if isnan(recInfo{i,2})
        continue
    end
    
    % Digital in
    s = dir([recInfo{i,1},'\digitalin.dat']);
    numSamples = s.bytes/2;   
    
    fileID = fopen([recInfo{i,1},'\digitalin.dat']);
    fseek(fileID,0,'bof');
    stim = fread(fileID,[1,numSamples],'int16');
    stim = downsample(stim,sF/dSF);
    stim(stim ~= stimTrig) = 0;
    
    % Amplifier in
    numChan = find(cell2mat(activeChan(:,3)) == 1);
    for ii = 1:length(numChan)
        iChan = numChan(ii,1);
        recInfo{i,5}{ii,1} = activeChan{iChan,2};
    end

    numChan = length(numChan);
    s = dir([recInfo{i,1},'\amplifier.dat']);
    numSamples = s.bytes/(2*numChan);   
    
    fileID = fopen([recInfo{i,1},'\amplifier.dat']);
    fseek(fileID,0,'bof');
    data = fread(fileID,[numChan,numSamples],'int16');
    data = downsample(data',sF/dSF)';
    
    timeStamps = 0;
    stimTS = 0;
    ck = 0;
    count = 1;
    for ii = 1:length(stim)
        if stim(1,ii) == stimTrig && ck == 0            

            timeStamps(count,1) = ii; %#ok<SAGROW> 
            count = count+1;
            ck = 1;
        elseif stim(1,ii) == 0 && ck == 1
            ck = 0;
        end
    end
    recInfo{i,4} = timeStamps/dSF;
    clear stim      

    for ii = 1:size(data,1)
        
        theta = bandpass(data(ii,:),[5,12],dSF);
        delta = bandpass(data(ii,:),[1,4],dSF);

        thetaAmp = abs(hilbert(theta));
        deltaAmp = abs(hilbert(delta));

        theta_deltaMean = zeros(1,floor(length(thetaAmp)/dSF));
        for i3 = 1:length(theta_deltaMean)
            try
                tWin_thetaAmp = thetaAmp(1,i3+(dSF*(i3-1)):i3+(dSF*(i3)-1));
                tWin_deltaAmp = deltaAmp(1,i3+(dSF*(i3-1)):i3+(dSF*(i3)-1));
                theta_deltaMean(1,i3) = mean(tWin_thetaAmp)/mean(tWin_deltaAmp);
            catch
                tWin_thetaAmp = thetaAmp(1,i3+(dSF*(i3-1)):end);
                tWin_deltaAmp = deltaAmp(1,i3+(dSF*(i3-1)):end);
                theta_deltaMean(1,i3) = mean(tWin_thetaAmp)/mean(tWin_deltaAmp);
            end
        end
        recInfo{i,5}{ii,2} = theta_deltaMean;
    end   

%     for ii = 1:numChan
%         [wave,PC,frex] = waveletize(tWin(:,:,ii),dSF,50,0.1,100,0,0);
%         wave = zscore(wave);
%         recInfo{i,6}{ii,2} = wave;
%         recInfo{i,7}{ii,2} = PC;
%     end
    disp([num2str(i),'/',num2str(numFolders)])
end

% count = 1;
% for i = 1:numFolders
%     if isnan(recInfo{i,2})
%         continue
%     end
% 
%     for ii = 1:numChan 
% 
%         subplot(numFolders,numChan,count)
%         contourf(recInfo{i,6}{ii,2},100,'linecolor','none');
%         set(gca,'xlim',[1,Post*dSF],'XTick',1:(Post*dSF)/10:Post*dSF,'XTickLabel',(1:(Post*dSF)/10:Post*dSF)/dSF)
%         set(gca,'ylim',[1,50],'YTick',1:10:50,'YTickLabel',frex(1:10:end))
%         title([recInfo{i,6}{ii,1},' ',num2str(recInfo{i,2}),'uA'])
%         count = count+1;
%     end
% end    
save([PathName,'\recInfo.mat'],'recInfo')
%%

subplot(3,1,1)
hold on
plot(recInfo{2,5}{1,2})
axis tight
xline(recInfo{2,4}(1,1),'LineWidth',3)
xline(recInfo{2,4}(2,1),'LineWidth',3)
ylim([0,15])
title(recInfo{2,5}{3,1})
hold off

subplot(3,1,2)
hold on
plot(recInfo{2,5}{5,2})
axis tight
xline(recInfo{2,4}(1,1),'LineWidth',3)
xline(recInfo{2,4}(2,1),'LineWidth',3)
ylim([0,15])
title(recInfo{2,5}{5,1})
hold off

subplot(3,1,3)
hold on
plot(recInfo{2,5}{6,2})
axis tight
xline(recInfo{2,4}(1,1),'LineWidth',3)
xline(recInfo{2,4}(2,1),'LineWidth',3)
ylim([0,15])
title(recInfo{2,5}{6,1})
hold off

%%
d = recInfo{2,5}{6,2}(1:round(recInfo{2,4}(1,1)));
numD9 = round(length(d)*0.9);
h1 = histogram(d,50);
bins = h1.BinEdges;
h1 = h1.BinCounts;

d = recInfo{2,5}{6,2}(1,round(recInfo{2,4}(1,1)):round(recInfo{2,4}(2,1)));
numD9 = round(length(d)*0.9);
h2 = histogram(d,'BinEdges',bins);
h2 = h2.BinCounts;

d = recInfo{2,5}{6,2}(1,round(recInfo{2,4}(2,1)):end);
numD9 = round(length(d)*0.9);
h3 = histogram(d,'BinEdges',bins);
h3 = h3.BinCounts;

medianBin = zeros(1,50);
count = 1;
for i = 1:length(bins)-1

    medianBin(1,count) = mean(bins(1,i:i+1));
    count = count+1;
end

figure
subplot(2,1,1)
bar(medianBin,h2-h1)
title('stim - pre-baseline hist')

subplot(2,1,2)
bar(medianBin,h3-h1)
title('post-baseline - pre-baseline hist')


%%

function timeVec = getTime()
    Now = clock;
    timeVec = (Now(1,4)*60)+Now(1,5)+(Now(1,6)/60);
end

function showVideo(frame)

    global videoFig
    if isempty(videoFig) == 1 || ishandle(videoFig) == 0
        videoFig = figure;
        videoFig = axes(videoFig);
    else
        cla(videoFig)
    end

    matObj = matfile([tempdir,'Go.mat'],'Writable',true);
    numCam = size(frame,4);
    if matObj.Go == 1
        for i = 1:numCam
            subplot(1,numCam,i,videoFig)
            imshow(frame(:,:,:,i),'Parent',videoFig);
        end
        drawnow
    end
end

function plotTD(data)

    global dataFig
    if isempty(dataFig) == 1 || ishandle(dataFig) == 0
        dataFig = figure;
        dataFig = axes(dataFig);
    end

    theta = bandpass(data{1,1},[5,12],data{2,1});
    delta = bandpass(data{1,1},[1,4],data{2,1});
    thetaAmp = abs(hilbert(theta));
    deltaAmp = abs(hilbert(delta));
    theta_delta = thetaAmp./deltaAmp;
    plot(dataFig,theta_delta)
    title(num2str(mean(theta_delta(end-data{2,1}+1:end))))
    drawnow
end

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