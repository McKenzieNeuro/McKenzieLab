% %function MakeAll_AlignCoordinatesScratch(ok, LE, RE, TB)
% fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
% % kp = cellfun(@any,regexp(fils,'Novel Env'));
% % kp = cellfun(@any,regexp(fils,));
% % fils = fils(kp);
% 
% [dirs] = fileparts(fils);
% 
% dirs =  unique(dirs);
% kp = ~contains(dirs,'Linear');
% dirs=  dirs(kp);
% 

fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'Novel Env'));

fils = fils(kp);

[dirs] = fileparts(fils);
%%For 2020 use :[dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);
%

%%
for i = 1:length(dirs)
    cd(dirs{i})
    i

    %%Checking for sessiondata
    if exist([dirs{i} filesep 'sessiondata.mat'])
        if exist('arena_edges.mat')

    
            %load tracking data
            load([dirs{i} filesep 'sessiondata.mat']);
    
            %load context edges
    
            if exist([dirs{i} filesep 'contextTransition1.mat'])
                load([dirs{i} filesep 'contextTransition1.mat']);
            else
                load([dirs{i} filesep 'arena_edges.mat']);
                load([dirs{i} filesep 'contextTransition.mat']);
    
                if size(data,1)==2 && any(contains(data(:,1),'home'))
                    hm_edge = contextEdges.edges{contains(contextEdges.contextName,'home')};
                    other_edge = contextEdges.edges{~contains(contextEdges.contextName,'home')};
    
                    data{contains(data(:,1),'home'),3} = hm_edge;
                    data{~contains(data(:,1),'home'),3} = other_edge;
    
                elseif size(data,1)==1
    
    
                    data{1,3} = contextEdges.edges{1};
                else
                    error('fix me')
                    % [a,b] = ismember(contextEdges.contextName,data(:,1));
    
                end
    
            end
    
    
            sessiondata.contextEntry = data;
            %get times for each context entry
            context_entry = cell2mat(data(:,2));
    
            %sort and get which index goes in which order (b)
            [context_entry,b] = sort(context_entry);
            %sort the edge/transition by time
            data = data(b,:);
    
            %define entry times
            epochs_on  = context_entry;
    
            %define exit times
            epochs_off = [context_entry(2:end); sessiondata.behavior.ts_video(end)];
    
            %build matrix of onsets and offsets
            epochs = [epochs_on epochs_off];
    
    
            sessiondata.contextEntry =  sortrows(sessiondata.contextEntry, 2);
            
            
            nSamples = length(sessiondata.behavior.position.left_ear);
            sessiondata.behavior.position.left_ear_cor = nan(nSamples,2);
            sessiondata.behavior.position.right_ear_cor = nan(nSamples,2);
            sessiondata.behavior.position.tail_base_cor = nan(nSamples,2);
            sessiondata.behavior.position.context = cell(nSamples,1);
            %loop over contexts
            for j = 1:size(epochs,1)
    
   
                [ptStartts,ptStartInd] = bestmatch(sessiondata.contextEntry{j,2},sessiondata.behavior.ts_video);    

                if j<size(sessiondata.contextEntry,1)   
                    [ptEndts, ptEndInd] = bestmatch(sessiondata.contextEntry{j+1,2},sessiondata.behavior.ts_video);
                    ptEndInd = ptEndInd- 1;
                else
                    ptEndts = sessiondata.behavior.ts_video(end);
                    ptEndInd = length(sessiondata.behavior.ts_video);
                end
                kp = ptStartInd:ptEndInd;  
               
                %keep data for those subsets of time (kp)
                LE =  sessiondata.behavior.position.left_ear(kp,:);
               
                RE =  sessiondata.behavior.position.right_ear(kp,:);
              
                TB = sessiondata.behavior.position.tail_base(kp,:);
                
    
    
                %define the edges of the context for that subst of time (note
                %we need to flip by 480 since the edges are defined in flipped
                %coordinates  (for some reason) and the height of the video is 480px


                edg = [data{j,3}(:,1) 480-data{j,3}(:,2)];

                if max(pdist(edg(:,2))) > max(pdist(edg(:,1)))
                    tempedg = edg(:,2);
                    edg(:,2)= edg(:,1);
                    edg(:,1) = tempedg;
    
                    tempLE = LE(:,2);
                    LE(:,2)= LE(:,1);
                    LE(:,1) = tempLE;
    
                    tempRE = RE(:,2);
                    RE(:,2)= RE(:,1);
                    RE(:,1) = tempRE;                
    
                    tempTB = TB(:,2);
                    TB(:,2)= TB(:,1);
                    TB(:,1) = tempTB;
    
                end
               
                % rescale by the long edge (y-axis)
                LE = LE/range(edg(:,2));
                RE = RE/range(edg(:,2));
                TB = TB/range(edg(:,2));
                edg = edg/range(edg(:,2));



                %Projective Transformation
                movingPoints = [edg(:,1), edg(:,2)];
                fixedPoints = [0 0; max(pdist(edg(:,1))) 0; max(pdist(edg(:,1))) max(pdist(edg(:,2))); 0 max(pdist(edg(:,2)))];
                tform = fitgeotrans(movingPoints , fixedPoints , "projective");
                
                [edg(:,1),edg(:,2)] = transformPointsForward(tform, movingPoints (:,1)',movingPoints (:,2)');
                [LE(:,1),LE(:,2)] = transformPointsForward(tform, LE(:,1)', LE(:,2)');
                [RE(:,1),RE(:,2)] = transformPointsForward(tform, RE(:,1)', RE(:,2)');
                [TB(:,1),TB(:,2)] = transformPointsForward(tform, TB(:,1)', TB(:,2)');

                   
                %check for negative values:
                if min(edg(:,1)) < 0 
                    edg(:,1) = edg(:,1) + abs(min(edg(:,1)));
                    LE(:,1) = LE(:,1) + abs(min(edg(:,1)));
                    RE(:,1) = RE(:,1) + abs(min(edg(:,1)));
                    TB(:,1) = TB(:,1) + abs(min(edg(:,1)));
    
    
                end
    
                if min(edg(:,2)) < 0 
                    edg(:,2) = edg(:,2) + abs(min(edg(:,2)));
                    LE(:,2) = LE(:,2) + abs(min(edg(:,2)));
                    RE(:,2) = RE(:,2) + abs(min(edg(:,2)));
                    TB(:,2) = TB(:,2) + abs(min(edg(:,2)));
                end

                %%order the X and Ys for the edges
                sortEdgX = sort(edg(:,1));
                sortEdgY = sort(edg(:,2));
    
    %%
                sessiondata.behavior.position.left_ear_cor(kp,:) = LE;
                sessiondata.behavior.position.right_ear_cor(kp,:) = RE;
                sessiondata.behavior.position.tail_base_cor(kp,:) = TB;
                sessiondata.behavior.position.context(kp) = data(j,1);
                sessiondata.contextEntry{j,4} = edg;
    
                
    
               
            end
            
            
            save('sessiondata1.mat','sessiondata','-v7.3')
            i
        end
    end
end

