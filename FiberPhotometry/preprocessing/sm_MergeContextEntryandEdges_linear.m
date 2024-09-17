function sm_MergeContextEntryandEdges_linear(dirName)
%%Merges output of context transition and context edge and saves in sessiondata.contextEntry

%-This script will need to load neural data (sessiondata), context transition and context edges
%-then match which edges goes with which context
%-then save this as a field in sessiondata: sessiondata.contextEntry which has N rows (each row is a context) and three columns:
%Col1 = context name
%Col2 = time of entrys
%Col3 = context corners in pixels

cd(dirName)

if (exist(['contextTransition.mat'])|| exist(['trackON.mat'])) && exist(['sessiondata.mat'])
    
    if exist(['contextTransition.mat'])
        contextTrans = load('contextTransition.mat');
        
    else
        contextTrans = load('trackON.mat');
    end
    load('sessiondata.mat');
    
    if size(contextTrans.data,1)==2
        [~,b]= sort(cell2mat(contextTrans.data(:,2)));
        
        sortedContextTrans = contextTrans.data(b,:);
        
        % add first home cage
        if ~contains(sortedContextTrans{1,1},'home')
            sortedContextTrans = [{'home_cage'} ,{0};sortedContextTrans];
            
        end
        sessiondata.contextEntry(:,1:2) = sortedContextTrans;
    elseif size(contextTrans.data,1)>3
         [~,b]= sort(cell2mat(contextTrans.data(:,2)));
        
        sortedContextTrans = contextTrans.data(b,:);
        
        % add first home cage
        if ~contains(sortedContextTrans{1,1},'home')
            sortedContextTrans = [{'home_cage'} ,{0};sortedContextTrans];
            
        end
        sessiondata.contextEntry(:,1:2) = sortedContextTrans;
        
        
    end
    
    
    % get track coordinates
    
    kp_lin = find(contains( sessiondata.contextEntry(:,1),'track','IgnoreCase',true));
    
    if length(kp_lin)==1
    linEpoch = cell2mat(sessiondata.contextEntry(kp_lin:kp_lin+1,2))';
    else
            linEpoch = [cell2mat(sessiondata.contextEntry(kp_lin(1):kp_lin(1)+1,2))';... 
                cell2mat(sessiondata.contextEntry(kp_lin(2):kp_lin(2)+1,2))'];
    end
        
    kp = InIntervals(sessiondata.behavior.ts_video,linEpoch);
    k  = gaussian2Dfilter([100 100],2);
    sessiondata.behavior.position.left_ear(sessiondata.behavior.position.left_ear(:,1)<15,1) = nan;
    ok = histcn([sessiondata.behavior.position.left_ear(kp,1),sessiondata.behavior.position.left_ear(kp,2)],1:500,1:500);
    ok = nanconvn(ok,k)'>.5;
    [ix1,ix2] = find(ok);
    
    %left
    idx1 = min(ix2);
    idx2 = max(ix2);
    NW = find(sum(ok(:,idx1:idx1+5),2),1,'last');
    SW = find(sum(ok(:,idx1:idx1+5),2),1,'first');
    NE = find(sum(ok(:,idx2-5:idx2),2),1,'last');
    SE = find(sum(ok(:,idx2-5:idx2),2),1,'first');
    edges = [idx1 SW;idx2 SE;idx2 NE;idx1 NW];
    
    %flip it to match DLC
    
    edges(:,2) = 480- edges(:,2);
    
    
    sessiondata.contextEntry(kp_lin,3) = {edges};
    %bin and find bound
    save('sessiondata.mat', 'sessiondata','-v7.3')
    
    
    
else
    
    error([dirName ' missing transitions or edges'])
end % end test for edges/transitions



end
