% talk to RHX
useRHS = 1;
intanIP_RHS = '127.0.0.1'; % IP for RHS.
intanPort1_RHS = 5000; % Port for RHS commands.
intanPort2_RHS = 5001; % Port for RHS data output.

%%
%model info
model_fname  = 'R:\McKenzieLab\DGregg\SeizureForecast\Seizuredetect_demo\Training\randomForest.mat';
v = load(model_fname);
ops.numfreq  = [4.9334 7.7869 12.2910 19.4002 30.6214 48.3330 76.2892 120.4155 190.0649 300.0000];
ops.nchanFil= 16;
ops.nChanSubj = 8;
ops.channels = 8:15;
ops.Fs =  20000;
ops.nTempBin =  32;
ops.data =  15;
ops.Twin =  2;
ops.feature_fun = @sm_GetDataFeature2;

%%

sF = 20000; % Recording sample frequency
dSF = 1000; % Downsample frequency

framesPerBlock = 128;
blocksPerRead = 300;
numBandsPerChan = 1; % Amplifier channels can be displayed as LOW, WIDE, or HIGH

% Calculations for accurate parsing
numChan = 8;
numDsFrames = (framesPerBlock*blocksPerRead);
readTime = numDsFrames/sF;
numAmplifierBands = numBandsPerChan*numChan;
waveformBytesPerFrame = 4+2*numAmplifierBands;
waveformBytesPerBlock = framesPerBlock*waveformBytesPerFrame+4;
framesPerRead = framesPerBlock*blocksPerRead;
waveformBytesBlocks = blocksPerRead*waveformBytesPerBlock;
%%
intanRead = tcpclient(intanIP_RHS,intanPort2_RHS);

%%
flush(intanRead);
        
% Pre-allocate memory for blocks of waveform data (the amount that's
% plotted at once)
amplifierData = 32768*ones(numChan,framesPerRead);
ampTimewindow = zeros(numChan,2*sF);

% Initialize counters
chunkCount = 0;
blockCount = 0;
amplifierTimestampsIndex = 1;
Stopper = 0;



while 1

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
        
        % Scale these data blocks and downsample
        %amplifierData = downsample(0.195*(amplifierData'-32768),sF/dSF)';

        % Shift time window to next step
        ampTimewindow(:,1:numDsFrames) = [];
        ampTimewindow = [ampTimewindow,amplifierData]; %#ok<AGROW>

        Prob = sm_getSeizProb(amplifierData',v.rusTree,ops);
        clf
        imagesc(amplifierData(:,1:100:end))

        title(Prob)
        % Reset index
        drawnow
        blockCount = 0;
        amplifierTimestampsIndex = 1;
        Stopper = 0;
    end
end