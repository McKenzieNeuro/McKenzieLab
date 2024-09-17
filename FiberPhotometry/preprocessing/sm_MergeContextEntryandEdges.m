function sm_MergeContextEntryandEdges(dirName)
%%Merges output of context transition and context edge and saves in sessiondata.contextEntry

%-This script will need to load neural data (sessiondata), context transition and context edges
%-then match which edges goes with which context
%-then save this as a field in sessiondata: sessiondata.contextEntry which has N rows (each row is a context) and three columns:
%Col1 = context name
%Col2 = time of entrys
%Col3 = context corners in pixels

cd(dirName)

if exist(['arena_edges.mat']) && exist(['contextTransitionRevised.mat']) && exist(['sessiondata.mat'])
    
    arenaEdg  = load('arena_edges.mat'); 
    contextTrans = load('contextTransitionRevised.mat');
    load('sessiondata.mat');
    
    %%This removes the extra context names that were never used
    if length(arenaEdg.contextEdges.contextName) > length(arenaEdg.contextEdges.edges)
        %ASSUME THAT EXTRA NAMES ARE LEFT BLANK
        arenaEdg.contextEdges.contextName = arenaEdg.contextEdges.contextName(1:length(arenaEdg.contextEdges.edges));
    end
    
    
    % this will output the sessiondata.contextEntry array
    if size(arenaEdg.contextEdges.contextName,2) == size(contextTrans.data,1)
        
        
        %% deal with home cage for edges
        % WE ASSUME THAT THE FIRST HOMECAGE HAPPENED FIRST IN THE RECORDING
        kp_home_edge = cellfun(@any,regexp(arenaEdg.contextEdges.contextName,'home')); %logical array, 1 = in home context, 0 means in other context
        home_edges = arenaEdg.contextEdges.edges(kp_home_edge);
        
        %match the context data with transition
        kp_home = find(cellfun(@any,regexp(contextTrans.data(:,1),'home')));
        [~,b] = sort(cell2mat(contextTrans.data(kp_home,2)));
        
        if length(home_edges) == length(kp_home)
            
            %match each homecage with each respective edge
            for jj = 1:length(kp_home)
                contextTrans.data{kp_home(b(jj)),3} = home_edges{jj};
                
            end
            
        else
            error('wrong number of homecages')
            
        end
        
        
        %% deal with other context (not homecage)
        
        kp_context_edge = ~kp_home_edge; %logical array
        context_edges = arenaEdg.contextEdges.edges(kp_context_edge);
        context_names = arenaEdg.contextEdges.contextName(kp_context_edge);
        
        %contexts were labeled with a context name and ended with the
        %number which indicated the order of the context transition
        contextNums = regexp(context_names,['[0-9]'],'match');
        contextNums = cellfun(@(a) str2num(cell2mat(a)),contextNums);
        [~,b] = sort(contextNums);
        context_edges = context_edges(b); % ensures that the context edges are sorted by context number
        
        
        kp_context = setdiff(1:size(contextTrans.data,1),kp_home);
        
        
        [~,b] = sort(cell2mat(contextTrans.data(kp_context,2))); % sorts context by time
        
        
        if length(context_edges) == length(kp_context)
            
            %match each homecage with each edge
            for jj = 1:length(kp_context)
                contextTrans.data{kp_context(b(jj)),3} = context_edges{jj};
                
            end
            
        else
            error('wrong number of contexts')
            
        end
        data = contextTrans.data;
        sortedContextTrans = sortrows(data, 2);
        sessiondata.contextEntry(:,1:3) = sortedContextTrans;
        
        save('sessiondata.mat', 'sessiondata','-v7.3')
        
    else
        
        error('maybe missing first homecage?')
        
    end % end test for equal number of edges transitions
    
    
    
else
    
    error([dirName ' missing transitions or edges'])
end % end test for edges/transitions



end
