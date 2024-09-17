function sm_MergeDist2Edge_as(dirName)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here


% Defining the shape of each environment
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
    {' large'};...
    {' pink'};...
    {' rat'};...
    {'rat'};...
    {'BT'};...
    {'RR'};...
    ];

circles = AllContexts([4,7,12]);
ellipse = AllContexts([6]);
rectangle = AllContexts([1,2,3,5,8,9,10,11,13,14,15,16,17]);


%% Loop through all directories and find the shortest distance from the
% right ear, left ear, and tail base to the environments edge

% for ii = 1
    cd(dirName)
    if exist([dirName filesep 'sessiondata.mat'])

            load('sessiondata.mat');
            totalEdgDist = nan(size(sessiondata.behavior.ts_video(:,1)));

            % loop over each environment
            for i = 1:size(sessiondata.contextEntry,1)
%             for i = 2

                % Define edges for context i
                edges = sessiondata.contextEntry{i,4};

                %Defining kp as the range of indices belonging to the
                %subjects position while in context i by matching the
                %indices to the timestamp where the subject was placed and
                %removed in its respective environment
               
                [ptStartts,ptStartInd] = bestmatch(sessiondata.contextEntry{i,2},sessiondata.behavior.ts_video);

                if i<size(sessiondata.contextEntry,1)
                    [ptEndts, ptEndInd] = bestmatch(sessiondata.contextEntry{i+1,2},sessiondata.behavior.ts_video);
                    ptEndInd = ptEndInd- 1;
%                                         ptEndInd = ptEndInd - 1500

                else
                    ptEndts = sessiondata.behavior.ts_video(end);
                    ptEndInd = length(sessiondata.behavior.ts_video);

                    %CHECK
%                     ptEndInd = ptEndInd - 1500
                end
                kp = ptStartInd:ptEndInd;

                %Defining the position of the left ear (LE), right ear (RE),
                % and tailbase (TB) while in environment i by looking at the
                %kp indices only
                pointsLE = sessiondata.behavior.position.left_ear_cor(kp,:);
                pointsRE = sessiondata.behavior.position.right_ear_cor(kp, :);
                pointsTB = sessiondata.behavior.position.tail_base_cor(kp, :);
                
                
                %test if environment i is a circle
                %totalEdgDist is an array that will hold the distances from 
                % the subjects LE, RE, TB position to the edge of its
                % environment in cm. The first row is distance from LE,
                % second row is distance from RE, and third row is distance
                % from TB
                
                %Find distance while in circular environment
                if ismember(sessiondata.contextEntry{i,1},circles)

                    for j = 1:3
                        if j == 1
                            points = pointsLE;
                            totalEdgDist(kp,j)  =  ptoc_distance_sm(edges, points);
                        end
                        if j == 2
                            points = pointsRE;
                            totalEdgDist(kp,j)  =  ptoc_distance_sm(edges, points);
                        end
                        if j == 3
                            points = pointsTB;
                            totalEdgDist(kp,j)  =  ptoc_distance_sm(edges, points);
                        end

                    end


                %Find distance while in rectangular environment
                elseif ismember(sessiondata.contextEntry{i,1},rectangle)

                    for j = 1:3
                        if j == 1
                            points = pointsLE;
                            totalEdgDist(kp,j)  = ptor_distance_sm(edges, points);
                        end
                        if j == 2
                            points = pointsRE;
                            totalEdgDist(kp,j)  = ptor_distance_sm(edges, points);
                        end
                        if j == 3
                            points = pointsTB;
                            totalEdgDist(kp,j)  = ptor_distance_sm(edges, points);
                        end
                    end

                %Find distance while in ellipse environment
                elseif ismember(sessiondata.contextEntry{i,1},ellipse)

                    for j = 1:3
                        if j == 1
                            points = pointsLE;
                            totalEdgDist(kp,j)  = ptoe_distance_sm_progress(edges, points);
                        end
                        if j == 2
                            points = pointsRE;
                            totalEdgDist(kp,j)  = ptoe_distance_sm_progress(edges, points);
                        end
                        if j == 3
                            points = pointsTB;
                            totalEdgDist(kp,j)  = ptoe_distance_sm_progress(edges, points);
                        end
                    end

                end
        %check if point is in polygon
        totalEdgDistT = totalEdgDist(kp,:);
        
        in = inpolygon(pointsLE(:,1),pointsLE(:,2),edges(:,1),edges(:,2));
        totalEdgDistT(~in,1) =  -totalEdgDistT(~in,1);
      
        
        in = inpolygon(pointsRE(:,1),pointsRE(:,2),edges(:,1),edges(:,2));
       totalEdgDistT(~in,2)  =  -totalEdgDistT(~in,2);
        
        in = inpolygon(pointsTB(:,1),pointsTB(:,2),edges(:,1),edges(:,2));
      totalEdgDistT(~in,3) =  -totalEdgDistT(~in,3);
        totalEdgDist(kp,:) = totalEdgDistT;
             
            end
        sessiondata.behavior.totalEdgDist = totalEdgDist;
        save('sessiondata.mat','sessiondata','-v7.3');
    end
end


