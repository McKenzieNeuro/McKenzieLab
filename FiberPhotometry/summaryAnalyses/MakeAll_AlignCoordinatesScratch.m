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
               
                
                %find the left edge
                [~,b] = sort(edg(:,1));
                left = edg(b(1:2),:);
                
                %find the bottom left point
                [~,b] = min(left(:,2));
                BL = left(b,:);
    
                % translate the position and edges so the bottom left is at the
                % origin
    %%%%%%%%%%%%Where negatives appear %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
                LE(:,1) = LE(:,1) - BL(1);
                LE(:,2) = LE(:,2) - BL(2);
    
                RE(:,1) = RE(:,1) - BL(1);
                RE(:,2) = RE(:,2) - BL(2);
    
                TB(:,1) = TB(:,1) - BL(1);
                TB(:,2) = TB(:,2) - BL(2);
    
                edg(:,1) = edg(:,1) - BL(1);
                edg(:,2) = edg(:,2) - BL(2);
    
                
                %now find top left
                [~,b] = sort(edg(:,1));
                left = edg(b(1:2),:);
                [~,b] = max(left(:,2));
                TL = left(b,:);
                
                % get angle of offset
                deg_offset = (atan(TL(1)/TL(2)));
                
                %find the rotation matrix to offset that angle
                R = rot2d( -deg_offset );
             
                %rotate the position data
                LE = LE*R;
                RE = RE*R;
                TB = TB*R;
    
    
                %rotate the edges
                ok =(edg*R);
    
                
                
                %keep long edge up/down
                xyrange = range(ok);
    
                
                %check if max x value is bigger than max y value
                %if it is, lets make the x values the y values
    
                if max(ok(:,2)) > max(ok(:,1))
                    tempok = ok(:,2);
                    ok(:,2)= ok(:,1);
                    ok(:,1) = tempok;
    
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
    
    
                %check for negative values:
                if min(ok(:,1)) < 0 
                    ok(:,1) = ok(:,1) + abs(min(ok(:,1)));
                    LE(:,1) = LE(:,1) + abs(min(ok(:,1)));
                    RE(:,1) = RE(:,1) + abs(min(ok(:,1)));
                    TB(:,1) = TB(:,1) + abs(min(ok(:,1)));
    
    
                end
    
                if min(ok(:,2)) < 0 
                    ok(:,2) = ok(:,2) + abs(min(ok(:,2)));
                    LE(:,2) = LE(:,2) + abs(min(ok(:,2)));
                    RE(:,2) = RE(:,2) + abs(min(ok(:,2)));
                    TB(:,2) = TB(:,2) + abs(min(ok(:,2)));
                end




                 
    %             if xyrange(1)<xyrange(2)
    %                 R = rot2d( 90 );
    %                
    %                 LE = LE*R;
    %                 RE = RE*R;
    %                 TB = TB*R;
    %                 ok =(ok*R);
    %                 
    %                 %shift over by new width
    %                 w = range(ok(:,1));
    %                 
    %                 LE(:,1) = LE(:,1) + w;
    %                 RE(:,1) = RE(:,1) + w;
    %                 TB(:,1) = TB(:,1) + w;
    %                 ok(:,1) = ok(:,1) +w;
    %                 
    %             end
                
                % rescale by the long edge (y-axis)
                LE = LE/range(ok(:,2));
                RE = RE/range(ok(:,2));
                TB = TB/range(ok(:,2));
                ok = ok/range(ok(:,2));
    
    
                % rescale
    
                %%order the X and Ys for the edges
                sortEdgX = sort(ok(:,1));
                sortEdgY = sort(ok(:,2));
    
    %             %%find min X values to find left most coordinates
    % 
    %             for i = 1:(size(sortEdgX)-1)
    % 
    %                 AB = 
    % 
    % 
    %             end
    %             %%for edges, find the shortest distance to determine which are 
    %             %%sharing a vertex distance between X and Ys
    %             
    %             distanceXs
    % 
    % 
    %             %%Find X and Y distance between node  to determine the longest 
    %             %%edge
    % 
    % 
    %             %%Divide all edges by the longest edge distance (height-wise)
    %             % to scale the rectangle to a height of one
    % 
    %             %%Move the rectangle
    % 
    % 
    % 
    % 
    % 
    %             %%%%%%
    
    
    %%




                sessiondata.behavior.position.left_ear_cor(kp,:) = LE;
                sessiondata.behavior.position.right_ear_cor(kp,:) = RE;
                sessiondata.behavior.position.tail_base_cor(kp,:) = TB;
                sessiondata.behavior.position.context(kp) = data(j,1);
                sessiondata.contextEntry{j,4} = ok;
    
                
    
               
            end
            
            
            save('sessiondata1.mat','sessiondata','-v7.3')
            i
        end
    end
end



%%

%sanity check (do this per context)

% binX = ,-150:5:150; % this is the resolution to bin the rows
% binY = ,-150:5:150; % this is the resolution to bin the columns

% [n,~,~,ix] = histcn([X Y],binX,binY); % bin your 2D data
% kp  = all(ix>0,2); % exclude all points that are outside of your bins
% mean_dist = accumarray(ix(kp,:),distanceCircle(kp),[],@nanmean,nan); % calculate the mean distance from the edge at every binned position
%imagesc(mean_dist)