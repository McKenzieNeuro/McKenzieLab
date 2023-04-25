masterDir = 'R:\McKenzieLab\DANEHippocampalResponse';
fils = getAllExtFiles(masterDir,'mat',1);

%kp_edges = cellfun(@any,regexp(fils,'arena_edges'));
%kp_trans = cellfun(@any,regexp(fils,'transition'));
kp_Novel = cellfun(@any,regexp(fils,'Novel Env'));

fils_N = fils(kp_Novel);
[dirs] = fileparts(fils_N);

dirs = unique(dirs);
%%
for i = 1:length(dirs)
    i
    cd(dirs{i})
    
    %%
    % NOTE TO INFANIA: TAKE EVERYTHING BELOW AND TURN INTO ITS OWN SCRIPT
    if exist('arena_edges.mat') && exist('contextTransition.mat')
        v  = load('arena_edges.mat');
        v1 = load('contextTransition.mat');
        
        if length(v.contextEdges.contextName) > length(v.contextEdges.edges)
            %ASSUME THAT EXTRA NAMES ARE LEFT BLANK
            v.contextEdges.contextName = v.contextEdges.contextName(1:length(v.contextEdges.edges));
        end
        
        
        
        % do something
        if size(v.contextEdges.contextName,2) == size(v1.data,1)
            
            
            % deal with home cage for edges
            % WE ASSUME THAT THE FIRST HOMECAGE HAPPENED FIRST IN THE
            % RECORDING
            kp_home_edge = cellfun(@any,regexp(v.contextEdges.contextName,'home'));
            home_edges = v.contextEdges.edges(kp_home_edge);
            
            %match with transition
            kp_home = find(cellfun(@any,regexp(v1.data(:,1),'home')));
            [~,b] = sort(cell2mat(v1.data(kp_home,2)));
            
            if length(home_edges) == length(kp_home)
                
                %match each homecage with each edge
                for j = 1:length(kp_home)
                    v1.data{kp_home(b(j)),3} = home_edges{j};
                    
                end
                
            else
                error('wrong number of homecages')
                
            end
            
            
            % deal with other context
            
            kp_context_edge = ~kp_home_edge;
            context_edges = v.contextEdges.edges(kp_context_edge);
            context_names = v.contextEdges.contextName(kp_context_edge);
            contextNums = regexp(context_names,['[0-9]'],'match');
            contextNums = cellfun(@(a) str2num(cell2mat(a)),contextNums);
            [~,b] = sort(contextNums);
            context_edges = context_edges(b); % ensures that the context edges are sorted by context number
            
            
            kp_context = setdiff(1:size(v1.data,1),kp_home);
            
            
            
            [~,b] = sort(cell2mat(v1.data(kp_context,2))); % sorts context by time
            
            
            if length(context_edges) == length(kp_context)
                
                %match each homecage with each edge
                for j = 1:length(kp_context)
                    v1.data{kp_context(b(j)),3} = context_edges{j};
                    
                end
                
            else
                error('wrong number of contexts')
                
            end
            data = v1.data;
            save('contextTransition1.mat','data')
            
        else
            
            error('maybe missing first homecage?')
            
            
        end % end test for equal number of edges transitions
        
        
        
    else
        
        error([dirs{i} ' missing transitions or edges'])
    end % end test for edges/transitions
    
    %%
    
end % end directories


