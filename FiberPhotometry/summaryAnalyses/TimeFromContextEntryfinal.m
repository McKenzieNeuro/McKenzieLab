fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'Novel Env'));

fils = fils(kp);

%%For 2022 use :[dirs] = fileparts(fils);
[dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);


for ii = 1:length(dirs)
    ii
    cd(dirs{ii})

    if exist([dirs{ii} filesep 'sessiondataR.mat'])

        load('sessiondataR.mat')
        %%
       
        timeFrmCChange = nan(size(sessiondata.behavior.ts_video(:,1)));
        sessiondata.behavior.timeFrmEntry = nan(size(sessiondata.behavior.ts_video(:,1)));
               
        % loop over each environment
    
        for i = 1:size(sessiondata.contextEntry,1)
         
            edges = sessiondata.contextEntry{i,4};

            [ptStartts,ptStartInd] = bestmatch(sessiondata.contextEntry{i,2},sessiondata.behavior.ts_video);    

            if i<size(sessiondata.contextEntry,1)   
                [ptEndts, ptEndInd] = bestmatch(sessiondata.contextEntry{i+1,2},sessiondata.behavior.ts_video);
                ptEndInd = ptEndInd- 1;
            else
                ptEndts = sessiondata.behavior.ts_video(end);
                ptEndInd = length(sessiondata.behavior.ts_video);
            end
            kp = ptStartInd:ptEndInd;        


            for i3 = kp
                if i3 == 1
                    timeFrmCChange(i3,1) = 0;
                else
                    timeFrmCChange(i3,1) = (sessiondata.behavior.ts_video(i3) - cell2mat(sessiondata.contextEntry(i,2)));
                end
            end

            sessiondata.behavior.timeFrmEntry = timeFrmCChange;
            save('sessiondataR.mat','sessiondata','-v7.3')
        end
       
    end

end
