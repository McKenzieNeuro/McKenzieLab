
% Use this to check the alighntment of the video sync and sleep score.
% Select the excel sheet and enter the row number for i

[file,path] = uigetfile('*.xlsx');
t = readtable([path,file],'VariableNamingRule','preserve');
t_cell = table2cell(t);
numFile = size(t_cell,1);
%%
i = 1;

i = i-1;
disp(['Loading File ',num2str(i),'/',num2str(numFile)])

digInFile = [t_cell{i,2},'\',t_cell{i,5},'\',t_cell{i,6}];
pFile = [t_cell{i,2},'\',t_cell{i,3}(1,1:end-4),'_poseData.mat'];
sFile = [t_cell{i,2},'\',t_cell{i,3}(1,1:end-4),'_sleepData2.mat'];
disp(pFile)

poseData = matfile(pFile,'Writable',true);

digData = poseData.digData1KHz;
conLED = poseData.conLED;
neuralTS = poseData.neuralTS;
velocity = poseData.velocity;

if ischar(conLED)
    conLED = zeros(size(poseData.diffPix));
end

neuralTS_ms = round(neuralTS*1000)+1;
neuralTS_ms(isnan(neuralTS_ms)) = [];

try
    sleepData = load(sFile);
    sleepStates = sleepData.SleepState.idx.states;
    sleepTS = sleepData.SleepState.idx.timestamps;

    idxTS = zeros(length(sleepTS),1);
    for ii = 1:length(sleepTS)
        [~,idxTS(ii,1)] = min(abs(neuralTS-ii));
    end
catch
    sleepStates = zeros(length(neuralTS),1);
    idxTS = 1:length(neuralTS);
end

if length(digData) < neuralTS_ms(end,1)
    digData(end+1:neuralTS_ms(end,1)) = 0;
end

figure
clf
subplot(2,1,1)
    hold on
    plot(digData(neuralTS_ms,1)*2,'k')
    plot(conLED,'r')
    axis tight
    xlabel('Frame')
    yticks([1,2])
    yticklabels({'LED','Intan'})
subplot(2,1,2)
    hold on
    plot(smooth(velocity(idxTS,1),20),'r')
    plot(sleepStates,'k')
    axis tight
    xlabel('Time (s)')
    ylim([0,5])
    yticks([1,2,3,5])
    yticklabels({'Active','Awake','SWS','REM'})
