fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'inear'));

fils = fils(kp);

%%For 2022 use :[dirs] = fileparts(fils);
[dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);
%%

%%Load the edge coordinates for the linear track
load('R:\McKenzieLab\DANEHippocampalResponse\linearTrack_edges.mat')


for ii = 1:length(dirs)
    ii
    cd(dirs{ii})
    if exist([dirs{ii} filesep 'sessiondata2.mat'])
        load('sessiondata2.mat')
        %%

        %             totalEdgDist = nan(size(sessiondata.behavior.position.left_ear(:,1)));
        totalEdgDist = nan(length(sessiondata.behavior.ts_video(:,1)),3);

        % loop over each environment

        for i = 1:size(sessiondata.contextEntry,1)

            if contains(sessiondata.contextEntry{i,1}, 'near' )
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

                pointsLE = sessiondata.behavior.position.left_ear_cor(kp,:);
                pointsRE = sessiondata.behavior.position.right_ear_cor(kp, :);
                pointsTB = sessiondata.behavior.position.tail_base_cor(kp, :);
                %test if environment i is a circle


                    % FIND TIMES WHEN MOUSE IS IN CONTEXT
                if contains(sessiondata.contextEntry{i,1}, 'inear')

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

        end

        save('totalDistanceFrmEdg', 'totalEdgDist');
    end
end