function recInfo = sm_recfileAnalysis(ProjectPathName)
%%

load([ProjectPathName,'/config.mat'])


if ~exist(PathName)
    PathName = ProjectPathName;
end

stimTrig = 1;

disp('Starting analysis');

recInfo = [{'File'},{'System'},{'Data'}];

listing = dir(PathName);
dirRemove = false(length({listing.name}),1);
for i = 1:length({listing.name})
    if strcmp(listing(i).name,'.') == 1 || strcmp(listing(i).name,'..') == 1 || listing(i).isdir == 0
        dirRemove(i,1) = true;
    end
end
listing(dirRemove) = [];

numFolders = size(listing,1);
for i = 1:size(listing,1)
    
    recInfo{i+1,3} = cell(5,7);
    recInfo{i+1,3}(1,:) = [{'Port'},{'Subject'},{'Stim amp (uI)'},{'Stim Chan'},{'Stim start/stop'},{'Stim TS'},{'Theta/delta'}];
    recInfo{i+1,3}{2,1} = 'A';
    recInfo{i+1,3}{3,1} = 'B';
    recInfo{i+1,3}{4,1} = 'C';
    recInfo{i+1,3}{5,1} = 'D';
    
    if contains(listing(i).name,'RHS') == 1
        
        recInfo{i+1,1} = [listing(i).folder,'\',listing(i).name];
        recInfo{i+1,2} = 'RHS';
        
        for ii = 1:4
            if boxRHS(ii,1) > 0
                recInfo{i+1,3}{ii+1,2} = boxSubject(boxRHS(ii,1),1);
                recInfo{i+1,3}{ii+1,3} = stimI(ii,1);
                recInfo{i+1,3}{ii+1,4} = stimChan(ii,:);
            end
        end
    elseif contains(listing(i).name,'RHD') == 1
        
        recInfo{i+1,1} = [listing(i).folder,'\',listing(i).name];
        recInfo{i+1,2} = 'RHD';
        
        for ii = 1:4
            if boxRHD(ii,1) > 0
                recInfo{i+1,3}{ii+1,2} = boxSubject(boxRHD(ii,1),1);
            end
        end
    end
end

for i = 2:numFolders+1
    
    if strcmp(recInfo{i,2},'RHD') == 1
        
        
        continue
    end
    
    % Digital in
    disp('Working on digitalin.dat')
    s = dir([recInfo{i,1},'\digitalin.dat']);
    numSamples = s.bytes/2;
    
    fileID = fopen([recInfo{i,1},'\digitalin.dat']);
    fseek(fileID,0,'bof');
    L = floor(numSamples/(sF*600));
    digIN = nan(1,L*600*dSF);
    
    Count = 1;
    for ii = 1:L
        
        temp = fread(fileID,[1,sF*600],'int16');
        temp = downsample(temp,sF/dSF);
        temp(temp ~= stimTrig) = 0;
        digIN(:,Count:Count+size(temp,2)-1) = temp;
        Count = Count+size(temp,2);
    end
    
    timeStamps = 0;
    ck = 0;
    count = 1;
    for ii = 1:length(digIN)
        if digIN(1,ii) == stimTrig && ck == 0
            
            timeStamps(1,count) = ii; %#ok<SAGROW>
            count = count+1;
            ck = 1;
        elseif digIN(1,ii) == 0 && ck == 1
            ck = 0;
        end
    end
    recInfo{i,3}([false;stimI > 0],5) = {round(timeStamps/dSF)};
    clear digIN
    
    % Find number of active channels
    numChan = length(find(cell2mat(activeChan(:,3,:)) == 1));
    chan = cell(numChan,2);
    count = 1;
    for ii = 1:size(activeChan,3)
        for i3 = 1:size(activeChan,1)
            if activeChan{i3,3,ii} == 1
                chan(count,:) = activeChan(i3,1:2,ii);
                count = count+1;
            end
        end
    end
    
    % Stim TS
    disp('Working on stim.dat')
    stimfil = [recInfo{i,1},'\stim.dat'];
    
    if exist(stimfil)
        s = dir(stimfil);
        numSamples = s.bytes/(2*numChan);
        
        fileID = fopen([recInfo{i,1},'\stim.dat']);
        fseek(fileID,0,'bof');
        L = floor(numSamples/(sF*600));
        stimTS = nan(numChan,L*600*dSF);
        
        Count = 1;
        for ii = 1:L
            
            temp = fread(fileID,[numChan,sF*600],'int16');
            temp = downsample(temp',sF/dSF)';
            temp(temp > 0) = 1;
            stimTS(:,Count:Count+size(temp,2)-1) = temp;
            Count = Count+size(temp,2);
        end
        
        
        
        count = 1;
        for ii = 1:size(activeChan,3)
            for i3 = 1:size(activeChan,1)
                if ~isempty(recInfo{i,3}{ii+1,3}) && activeChan{i3,3,ii} == 1
                    if any(recInfo{i,3}{ii+1,4} == i3)
                        
                        temp = stimTS(count,:);
                        temp = find(diff(temp) > 0);
                        recInfo{i,3}{ii+1,6} = temp';
                    end
                    count = count+1;
                end
            end
        end
        clear stimTS
    end
    count = zeros(4,1);
    for ii = 1:size(chan,1)
        
        if contains(chan{ii,1},'a')
            idx = 1;
        elseif contains(chan{ii,1},'b')
            idx = 2;
        elseif contains(chan{ii,1},'c')
            idx = 3;
        elseif contains(chan{ii,1},'d')
            idx = 4;
        end
        count(idx,1) = count(idx,1)+1;
        recInfo{i,3}{idx+1,7}(count(idx,1),1:2) = chan(ii,:);
        
        
    end
end
save([PathName,'\recInfo.mat'],'recInfo')

beep
disp('Finsished');