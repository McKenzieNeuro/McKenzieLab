function sm_MergeTimeFromTransition(dirName)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here



    cd(dirName)

    if exist([dirName filesep 'sessiondata.mat'])

        load('sessiondata.mat')
        %%
       
        timeFrmCChange = nan(size(sessiondata.behavior.ts_video(:,1)));
        sessiondata.behavior.timeFrmEntry = nan(size(sessiondata.behavior.ts_video(:,1)));
               
        % loop over each environment
    
        for i = 1:size(sessiondata.contextEntry,1)
         
            %set the boundaries/edges of the environment
            edges = sessiondata.contextEntry{i,4};

            %Compare the timestamp from the context transition to the array
            %of timestamps per frame to find the best fitting frame index
            %at the point of transition
            [ptStartts,ptStartInd] = bestmatch(sessiondata.contextEntry{i,2},sessiondata.behavior.ts_video);    

            %Using similar logic as above, end index should be near the
            %beginning of the following context transition OR the last
            %frame
            if i<size(sessiondata.contextEntry,1)   
                [ptEndts, ptEndInd] = bestmatch(sessiondata.contextEntry{i+1,2},sessiondata.behavior.ts_video);
                ptEndInd = ptEndInd- 1;
            else
                ptEndts = sessiondata.behavior.ts_video(end);
                ptEndInd = length(sessiondata.behavior.ts_video);
            end
            kp = ptStartInd:ptEndInd;        


            %Time post context transition is using simple subtraction
            for i3 = kp
                if i3 == 1
                    timeFrmCChange(i3,1) = 0;
                else
                    timeFrmCChange(i3,1) = (sessiondata.behavior.ts_video(i3) - cell2mat(sessiondata.contextEntry(i,2)));
                end
            end

        end
        sessiondata.behavior.timeFrmEntry = timeFrmCChange;
        save('sessiondata.mat','sessiondata','-v7.3')

    end

end




