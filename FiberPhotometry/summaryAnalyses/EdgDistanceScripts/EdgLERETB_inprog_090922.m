fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'Novel Env'));

fils = fils(kp);

%%For 2022 use :[dirs] = fileparts(fils);
[dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);
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
    if exist([dirs{ii} filesep 'sessiondata.mat'])
        if exist('arena_edges.mat')

            load('sessiondata1.mat')
            load('arena_edges.mat')
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
                else
                    ptEndts = sessiondata.behavior.ts_video(end);
                    ptEndInd = length(sessiondata.behavior.ts_video);
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
                            totalEdgDist(kp,j)  = ptoe_distance_sm_progress_inprog_090922(edges, points);
                       end
                       if j == 2
                            points = pointsRE;
                            totalEdgDist(kp,j)  = ptoe_distance_sm_progress_inprog_090922(edges, points);
                       end
                       if j == 3
                            points = pointsTB;
                            totalEdgDist(kp,j)  = ptoe_distance_sm_progress_inprog_090922(edges, points);
                       end
                       
                    end
                end
        
            end
             
        end
    end
    %save('totalDistanceFrmEdg', 'totalEdgDist');
end
