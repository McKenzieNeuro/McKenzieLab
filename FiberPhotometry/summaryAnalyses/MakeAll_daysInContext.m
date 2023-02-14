
fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'Novel Env'));


fils = fils(kp);

%%For 2022 use :
[dirs, name, ext] = fileparts(fils);
% [dirs] = cellfun(@fileparts,fils,'UniformOutput',false);
dirs =  unique(dirs);


%%
ii = fullfile(dirs);
[ii, nameWDate, ext] = fileparts(ii);
nameWDate

AllContexts = [ ...
    {'home'}; ...
    {'context1'} ; ...
    {'context2'} ; ...
    {'BB'} ; ...
    {'CB'} ; ...
    {'RB'} ; ...
    {'WB'} ; ...
    {'CC'} ; ...
    {'OTT'};...
    {'homecage'};...
    {'homeContext'};...
    ];

dayHomeCount = 0;
daytext1Count = 0;
daytext2Count = 0;
dayBBCount = 0;
dayCBCount = 0;
dayRBCount = 0;
dayWBCount = 0;
dayCCCount = 0;
dayOTTCount = 0;

%Loop through all directories
for i = 1:length(dirs)
    cd(dirs{i})

    if exist('sessiondata1.mat')

        %load tracking data
        load('sessiondata1.mat');
    
        totalDaysperContext = nan(size(sessiondata.behavior.ts_video(:,1)));
        sessiondata.behavior.dayNum = nan(size(sessiondata.behavior.ts_video(:,1)));
        
        for i3 = 1:size(sessiondata.contextEntry,1)
            [ptStartts,ptStartInd] = bestmatch(sessiondata.contextEntry{i3,2},sessiondata.behavior.ts_video);    

            if i3<size(sessiondata.contextEntry,1)   
                [ptEndts, ptEndInd] = bestmatch(sessiondata.contextEntry{i3+1,2},sessiondata.behavior.ts_video);
                ptEndInd = ptEndInd- 1;
            else
                ptEndts = sessiondata.behavior.ts_video(end);
                ptEndInd = length(sessiondata.behavior.ts_video);
            end
            kp = ptStartInd:ptEndInd; 

            %Find if the subject is the same as previous
            %First directory exception
            if i == 1
                if contains(sessiondata.contextEntry{i3,1}, 'ome')
                    dayHomeCount = 1;
                    totalDaysperContext(kp,1) = 100;
                elseif contains(sessiondata.contextEntry{i3,1}, 'text1')
                    daytext1Count = 1;
                    totalDaysperContext(kp,1) = daytext1Count;
                elseif contains(sessiondata.contextEntry{i3,1}, 'text2')
                    daytext2Count = 1;
                    totalDaysperContext(kp,1) = daytext2Count;
                elseif contains(sessiondata.contextEntry{i3,1}, 'BB')
                    dayBBCount = 1;
                    totalDaysperContext(kp,1) = dayBBCount;
                elseif contains(sessiondata.contextEntry{i3,1}, 'CB')
                    dayCBCount = 1;
                    totalDaysperContext(kp,1) = dayCBCount;
                elseif contains(sessiondata.contextEntry{i3,1}, 'RB')
                    dayRBCount = 1;
                    totalDaysperContext(kp,1) = dayRBCount;
                elseif contains(sessiondata.contextEntry{i3,1}, 'WB')
                    dayWBCount = 1;
                    totalDaysperContext(kp,1) = dayWBCount;
                elseif contains(sessiondata.contextEntry{i3,1}, 'CC')
                    dayCCCount = 1;
                    totalDaysperContext(kp,1) = dayCCCount;
                elseif contains(sessiondata.contextEntry{i3,1}, 'OTT')
                    dayOTTCount = 1;
                    totalDaysperContext(kp,1) = dayOTTCount;
                end

            %Find if the subject is the same as previous
            else
                %Checks if parent folder are of the same Sub
                if strcmp(nameWDate{i}(1:5), nameWDate{i-1}(1:5))        
                    dateNTime1 = extractAfter(nameWDate{i},'-');
                    dateNTime2 = extractAfter(nameWDate{i-1},'-');
                    %Check if the dates are different, if they are add one
                    %for the new day in the context
                    if ~strcmp(dateNTime1(1:6), dateNTime2 (1:6))
        
                        if contains(sessiondata.contextEntry{i3,1}, 'ome')
                            dayHomeCount = dayHomeCount+1;
                            totalDaysperContext(kp,1) = 100;
                        elseif contains(sessiondata.contextEntry{i3,1}, 'text1')
                            daytext1Count = daytext1Count+1;
                            totalDaysperContext(kp,1) = daytext1Count;
                        elseif contains(sessiondata.contextEntry{i3,1}, 'text2')
                            daytext2Count = daytext2Count+1;
                            totalDaysperContext(kp,1) = daytext2Count;
                        elseif contains(sessiondata.contextEntry{i3,1}, 'BB')
                            dayBBCount = dayBBCount+1;
                            totalDaysperContext(kp,1) = dayBBCount;
                        elseif contains(sessiondata.contextEntry{i3,1}, 'CB')
                            dayCBCount = dayCBCount+1;
                            totalDaysperContext(kp,1) = dayCBCount;
                        elseif contains(sessiondata.contextEntry{i3,1}, 'RB')
                            dayRBCount = dayRBCount+1;
                            totalDaysperContext(kp,1) = dayRBCount;
                        elseif contains(sessiondata.contextEntry{i3,1}, 'WB')
                            dayWBCount = dayWBCount+1;
                            totalDaysperContext(kp,1) = dayWBCount;
                        elseif contains(sessiondata.contextEntry{i3,1}, 'CC')
                            dayCCCount = dayCCCount+1;
                            totalDaysperContext(kp,1) = dayCCCount;
                        elseif contains(sessiondata.contextEntry{i3,1}, 'OTT')
                            dayOTTCount = dayOTTCount+1;
                            totalDaysperContext(kp,1) = dayOTTCount;
                        end

                    end

                else
                    if i3 == 1
                        dayHomeCount = 0;
                        daytext1Count = 0;
                        daytext2Count = 0;
                        dayBBCount = 0;
                        dayCBCount = 0;
                        dayRBCount = 0;
                        dayWBCount = 0;
                        dayCCCount = 0;
                        dayOTTCount = 0;
                    end
                    if contains(sessiondata.contextEntry{i3,1}, 'ome')
                        dayHomeCount = 1;
                        totalDaysperContext(kp,1) = 100;
                    elseif contains(sessiondata.contextEntry{i3,1}, 'text1')
                        daytext1Count = 1;
                        totalDaysperContext(kp,1) = daytext1Count;
                    elseif contains(sessiondata.contextEntry{i3,1}, 'text2')
                        daytext2Count = 1;
                        totalDaysperContext(kp,1) = daytext2Count;
                    elseif contains(sessiondata.contextEntry{i3,1}, 'BB')
                        dayBBCount = 1;
                        totalDaysperContext(kp,1) = dayBBCount;
                    elseif contains(sessiondata.contextEntry{i3,1}, 'CB')
                        dayCBCount = 1;
                        totalDaysperContext(kp,1) = dayCBCount;
                    elseif contains(sessiondata.contextEntry{i3,1}, 'RB')
                        dayRBCount = 1;
                        totalDaysperContext(kp,1) = dayRBCount;
                    elseif contains(sessiondata.contextEntry{i3,1}, 'WB')
                        dayWBCount = 1;
                        totalDaysperContext(kp,1) = dayWBCount;
                    elseif contains(sessiondata.contextEntry{i3,1}, 'CC')
                        dayCCCount = 1;
                        totalDaysperContext(kp,1) = dayCCCount;
                    elseif contains(sessiondata.contextEntry{i3,1}, 'OTT')
                        dayOTTCount = 1;
                        totalDaysperContext(kp,1) = dayOTTCount;
                    end
                end
            end
        end
    end
    sessiondata.behavior.dayNum = totalDaysperContext;
    save('sessiondata1.mat','sessiondata','-v7.3')
end
