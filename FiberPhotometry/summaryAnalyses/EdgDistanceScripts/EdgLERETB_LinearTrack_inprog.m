
fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'inear'));

fils = fils(kp);

[dirs] = fileparts(fils);
%%For 2020 use :[dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);
%%


for ii = 1:length(dirs)
    ii
    cd(dirs{ii})
    if exist([dirs{ii} filesep 'sessiondata3.mat'])

        %load tracking data
        load([dirs{ii} filesep 'sessiondata3.mat']);


        %         totalEdgDist = nan(length(sessiondata.behavior.ts_video(:,1)),3);

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


                positionName = fieldnames(sessiondata.behavior.position);

                for j = 1:length(positionName)
                    if contains(positionName{j}, [strrep(sessiondata.contextEntry{i,1},' ', ''), '_cor'])
                        %Assign LE, RE, and TB
                        if contains(positionName{j}, 'LE')
                            pointsLE = sessiondata.behavior.position.(positionName{j});
                        elseif contains(positionName{j}, 'RE')
                            pointsRE = sessiondata.behavior.position.(positionName{j});
                        elseif contains(positionName{j}, 'TB')
                            pointsTB = sessiondata.behavior.position.(positionName{j});
                        end
                    end
                end

                %                 pointsLE = sessiondata.behavior.position.left_ear_cor(kp,:);
                %                 pointsRE = sessiondata.behavior.position.right_ear_cor(kp, :);
                %                 pointsTB = sessiondata.behavior.position.tail_base_cor(kp, :);
                %test if environment i is a circle


                %%CORRECT FROM HERE UP
                %%EDIT BELOW

                % FIND TIMES WHEN MOUSE IS IN CONTEXT
                if (contains(sessiondata.contextEntry{i,1}, 'near' ) || contains(sessiondata.contextEntry{i,1}, 'rack' ))

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

                %%%REPLACE THE left_ear_cor names below

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

                edgDist_name =  ['edgDist',  strrep(sessiondata.contextEntry{ii,1},' ','')];

               edgDistCor



            end

        end

        save('totalDistanceFrmEdg', 'totalEdgDist');
    end
end