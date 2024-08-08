fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
% fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'Novel Env'));

fils = fils(kp);

%%For 2022 use :[dirs] = fileparts(fils);
[dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);

%%

dirs = [...
        {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210819-123801'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210820-140801'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210821-081259'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210822-094611'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210823-165050'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210824-124626'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210825-101604'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210826-105416'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210827-104741'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210830-130702'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210831-135215'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210901-103622'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210902-114718'};...
    {'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210907-135954'}];
%%
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

circles = AllContexts([4,7]);
ellipse = AllContexts([6]);
rectangle = AllContexts([1,2,3,5,8,9,10,11]);

for ii = 1:length(dirs)
    ii
    cd(dirs{ii})
    if exist([dirs{ii} filesep 'sessiondataR.mat'])

            load('sessiondataR.mat')
%             load('arena_edges.mat')
            %%


            %             totalEdgDist = nan(size(sessiondata.behavior.position.left_ear(:,1)));
            totalEdgDist = nan(size(sessiondata.behavior.ts_video(:,1)));

            % loop over each environment

            for i = 1:size(sessiondata.contextEntry,1)

                edges = sessiondata.contextEntry{i,4};

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

                pointsLE = sessiondata.behavior.position.left_ear_cor(kp,:);
                pointsRE = sessiondata.behavior.position.right_ear_cor(kp, :);
                pointsTB = sessiondata.behavior.position.tail_base_cor(kp, :);
                %test if environment i is a circle

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


                    % FIND TIMES WHEN MOUSE IS IN CONTEXT
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

                for jjj = kp
                
                        if sessiondata.behavior.position.left_ear_cor(jjj,1) > max(edges(:,1)) || sessiondata.behavior.position.left_ear_cor(jjj,2) > max(edges(:,2))
                            totalEdgDist(jjj,1) = -totalEdgDist(jjj,1);
                        elseif sessiondata.behavior.position.left_ear_cor(jjj,1) < min(edges(:,1)) || sessiondata.behavior.position.left_ear_cor(jjj,2) < min(edges(:,2))
                            totalEdgDist(jjj,1) = -totalEdgDist(jjj,1);
                        end
                        if sessiondata.behavior.position.right_ear_cor(jjj,1) > max(edges(:,1)) || sessiondata.behavior.position.right_ear_cor(jjj,2) > max(edges(:,2))
                            totalEdgDist(jjj,2) = -totalEdgDist(jjj,2);
                        elseif sessiondata.behavior.position.right_ear_cor(jjj,1) < min(edges(:,1)) || sessiondata.behavior.position.right_ear_cor(jjj,2) < min(edges(:,2))
                            totalEdgDist(jjj,2) = -totalEdgDist(jjj,2);
                        end

                        if sessiondata.behavior.position.tail_base_cor(jjj,1) > max(edges(:,1)) || sessiondata.behavior.position.tail_base_cor(jjj,2) > max(edges(:,2))
                            totalEdgDist(jjj,3) = -totalEdgDist(jjj,3);
                        elseif sessiondata.behavior.position.tail_base_cor(jjj,1) < min(edges(:,1)) || sessiondata.behavior.position.tail_base_cor(jjj,2) < min(edges(:,2))
                            totalEdgDist(jjj,3) = -totalEdgDist(jjj,3);
                        end

                end
            end
        save('totalDistanceFrmEdgR', 'totalEdgDist');
    end
end