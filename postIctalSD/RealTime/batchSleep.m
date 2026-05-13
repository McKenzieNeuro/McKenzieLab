
[file,path] = uigetfile('*.xlsx');
t = readtable([path,file],'VariableNamingRule','preserve');
t = convertvars(t,{'Sleep Score'},'cellstr');
t_cell = table2cell(t);

%%

% Channels to look at, in base 1
useChan_delta = [1,2];
useChan_theta = [5,6];
for i = 1:size(t_cell,1)

    disp("Working on file: "+i+"/"+size(t_cell,1))

    if ~isempty(t_cell{i,3}) && ~isnan(t_cell{i,4})
        
        subjectName = t_cell{i,3}(1,1:end-4);
        pFile = [t_cell{i,2},'\',subjectName,'_poseData.mat'];
        if isfile(pFile)

            poseData = matfile(pFile,'Writable',true);
            if isprop(poseData,'neuralTS') && sum(poseData.neuralTS,'omitnan') > 0                

                if ~ischar(poseData.tsPulse_video)
                    t(i,11) = {'Yes'};
                else
                    t(i,11) = {'No'};
                end

                t(i,12) = {'Yes'};

                if ~ischar(poseData.tsPulse_video)  
                    t.('numPulse neural/video')(i) = {[num2str(length(poseData.tsPulse_neural)),'/',num2str(length(poseData.tsPulse_video))]};
                    t.('neural/video similarity')(i) = {round(similarity(poseData),3)};
                else
                    t.('numPulse neural/video')(i) = {[num2str(length(poseData.tsPulse_neural)),'/',num2str(0)]};
                    t.('neural/video similarity')(i) = {nan};
                end

                if isprop(poseData,'velocity') && sum(poseData.velocity) > 0

                    t.('Pose')(i) = {'Yes'};
                    
                    basePath = [t_cell{i,2},'\',t_cell{i,5}];
                    if isfile([basePath,'\amplifier.lfp'])

                       t(i,16) = {'Yes'};

                        if ~isfile([t_cell{i,2},'\',t_cell{i,3}(1,1:end-4),'_sleepData2.mat'])

                            port = t_cell{i,8};
                            numChan = t_cell{i,9};
                            switch port
                                case 'A'
                                    ckChan = 0:7;
                                case 'B'
                                    ckChan = 8:15;
                                case 'C'
                                    ckChan = 16:23;
                                case 'D'
                                    ckChan = 24:31;
                            end
                            fileChan_delta = ckChan(useChan_delta);
                            fileChan_theta = ckChan(useChan_theta);
    
                            try
                                savefolder = [fileparts(basePath),'\sleepData'];
                                if ~exist(savefolder,'dir')
                                    mkdir(savefolder)
                                end
                                SleepState = SleepScoreMaster(basePath,'savedir',[savefolder,'\',subjectName],'SWChannels',fileChan_delta,'ThetaChannels',fileChan_theta,'MotionSource','velocity','velocityFile',pFile,'overwrite',true);
                                save([t_cell{i,2},'\',t_cell{i,3}(1,1:end-4),'_sleepData2.mat'],'SleepState')
                                t(i,17) = {'Yes'};
                            catch
                                t{i,17} = {'No'};
                                disp('Problem with data')
                                
                            end
                        else
                            t(i,17) = {'Yes'};
                            disp('Sleep score done')
                        end
                    else
                        t(i,16) = {'No'};
                        t(i,17) = {'No'};
                        disp('No .lfp file')
                    end
                else
                    t(i,13) = {'No'};
                    disp('Pose not done')
                end
            else
                t(i,12) = {'No'};
                disp('Sync not done')
            end
        else
            t(i,11) = {'No'};
            t(i,12) = {'No'};
            t(i,13) = {'No'};
            t(i,17) = {'No'};
            t(i,18) = {'NaN'};
            t.('neural/video similarity')(i) = nan;
            disp('File not processed')
        end
    else
        t(i,11) = {'No'};
        t(i,12) = {'No'};
        t(i,13) = {'No'};
        t(i,17) = {'No'};
        t(i,18) = {'NaN'};
        t.('neural/video similarity')(i) = {nan};
        disp('File not processed')
    end
end

function out = similarity(poseData)

    digData = poseData.digData1KHz;
    conLED = poseData.conLED;
    neuralTS = poseData.neuralTS;

    if ischar(conLED)
        conLED = zeros(size(poseData.diffPix));
    else
        conLED(conLED > 0.9) = 1;
        conLED(conLED < 1) = 0;
    end

    neuralTS_ms = round(neuralTS*1000)+1;
    neuralTS_ms(isnan(neuralTS_ms)) = [];
    
    if length(digData) < neuralTS_ms(end,1)
        digData(end+1:neuralTS_ms(end,1)) = 0;
    end

    digData = digData(neuralTS_ms,1);
   
    out = sqrt(sum((digData-conLED).^2));
end