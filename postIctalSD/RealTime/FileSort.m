
sortDir = uigetdir;
d = dir(sortDir);
d(1:2) = [];

%%

tS = {'Animal ID','Root Folder','Video File','FPS','Data Folder','Data File','Sample Rate','Port','Channel Number'};
tX = tS;

countS = 1;
countX = 1;
for i = 1:size(d,1)

    disp(['Folder ',num2str(i),'/',num2str(size(d,1))])
    if d(i).isdir == 1

        d2 = dir([d(i).folder,'\',d(i).name]);
        d2(1:2) = [];
    
        if ~isempty(d2) && sum([d2.isdir]) > 0
    
            cd(d2(1).folder)
            
            try
                load('recInfo.mat')
            catch
                continue
            end
    
            for ii = 1:size(recInfo.fileData,2)
                if ~isempty(recInfo.fileData(ii).subject)
        
                    t = cell(1,9);
                    t{1,1} = recInfo.fileData(ii).subject;
                    t{1,2} = d2(1).folder;
                    t{1,3} = [t{1,1},'.mp4'];
        
                    vr = VideoReader([t{1,2},'\',t{1,3}]); %#ok<TNMLP>
                    t{1,4} = vr.FrameRate;
                    clear vr
        
                    t{1,5} = d2([d2.isdir] == 1).name;
        
                    d3 = dir([d2(1).folder,'\',d2([d2.isdir] == 1).name]);
                    d3(1:2) = [];
        
                    for i3 = 1:size(d3,1)
                        if contains(d3(i3).name,'digitalin') && contains(d3(i3).name,'.dat')
                            t{1,6} = d3(i3).name;
                        end
                    end
        
                    t{1,8} = recInfo.fileData(ii).port;
                    t{1,9} = sum(~cellfun(@isempty,{recInfo.fileData.subject}))*8;
        
                    blocks = recInfo.fileData(ii).blockDuration(:,1);
                    S = 1;
                    if length(blocks) > 1
                        
                        for i3 = 1:length(blocks)
                            if contains(blocks{i3,1},'X')
                                S = 0;
                            end
                        end
                    else
                        if strcmp(blocks{1,1},'BL')
                            S = 0;
                        end
                    end
        
                    if S == 1
                        countS = countS+1;
                        tS = [tS;t]; %#ok<AGROW>
                    else
                        countX = countX+1;
                        tX = [tX;t]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

            




            



        



