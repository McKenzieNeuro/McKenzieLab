function sm_alignEdgPosition_pix2cm_as(dirName)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


cd(dirName)
%%Checking for sessiondata
if exist([dirName filesep 'sessiondata.mat'])
    
    
    %load tracking data
    load([dirName filesep 'sessiondata.mat']);

    %define entry times
    context_entry  = cell2mat(sessiondata.contextEntry(:,2));
    
    %define exit times
    epochs_off = [context_entry(2:end); sessiondata.behavior.ts_video(end)];
    
    %build matrix of onsets and offsets
    epochs = [context_entry epochs_off];
    
    
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
            ptEndInd = ptEndInd- 1501; % backwards from next context entry
        else
            ptEndts = sessiondata.behavior.ts_video(end);
            ptEndInd = length(sessiondata.behavior.ts_video);
            ptEndInd = ptEndInd - 1500;
        end
        kp = ptStartInd:ptEndInd;
        sessiondata.behavior.position.context(kp) = sessiondata.contextEntry(j,1);
        %keep data for those subsets of time (kp)
        LE =  sessiondata.behavior.position.left_ear(kp,:);
        
        RE =  sessiondata.behavior.position.right_ear(kp,:);
        
        TB = sessiondata.behavior.position.tail_base(kp,:);
        
        %define the edges of the context for that subst of time (note
        %we need to flip by 480 since the edges are defined in flipped
        %coordinates  (for some reason) and the height of the video is 480px
        if ~isempty(sessiondata.contextEntry{j,3})
            edg = [sessiondata.contextEntry{j,3}(:,1) 480-sessiondata.contextEntry{j,3}(:,2)];
            
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
            
            %Projective Transformation
            fixedPoints = nan;
            movingPoints = [edg(:,1), edg(:,2)];
            if contains(sessiondata.contextEntry{j,1}, 'ome' )
                fixedPoints = [0 0; 24.88 0; 24.88 14.85; 0 14.85];
            elseif contains(sessiondata.contextEntry{j,1}, 'ontext1')
                fixedPoints = [0 0; 24.88 0; 24.88 14.85; 0 14.85];
            elseif contains(sessiondata.contextEntry{j,1}, 'ontext2')
                fixedPoints = [0 0; 39.43 0; 39.43 18.1; 0 18.1];
            elseif contains(sessiondata.contextEntry{j,1}, 'BB')
                fixedPoints = [0 0; 16.47 0; 16.47 16.47; 0 16.47];
            elseif contains(sessiondata.contextEntry{j,1}, 'CB')
                fixedPoints = [0 0; 9.97 0; 9.97 9.97; 0 9.97];
            elseif contains(sessiondata.contextEntry{j,1}, 'RB')
                fixedPoints = [0 0; 14.67 0; 14.67 9.90; 0 9.90];
            elseif contains(sessiondata.contextEntry{j,1}, 'CC')
                fixedPoints = [0 0; 13.77 0; 13.77 13.77; 0 13.77];
            elseif contains(sessiondata.contextEntry{j,1}, 'OTT')
                fixedPoints = [0 0; 25.60 0; 25.60 15.60; 0 15.60];
            elseif contains(sessiondata.contextEntry{j,1}, 'WB')
                fixedPoints = [0 0; 10.00 0; 10.00 10.00; 0 10.00];
            elseif contains(sessiondata.contextEntry{j,1}, 'pink')
                fixedPoints = [0 0; 25.60 0; 25.60 15.60; 0 15.60];
            elseif contains(sessiondata.contextEntry{j,1}, 'rat')
                fixedPoints = [0 0; 39.43 0; 39.43 18.1; 0 18.1];
            elseif contains(sessiondata.contextEntry{j,1}, 'large')
                fixedPoints = [0 0; 49.0 0; 49.0 49.0;0 49.0 ];
            elseif contains(sessiondata.contextEntry{j,1}, 'RR')
                fixedPoints = [0 0; 21 0; 21 13; 0 13];
            elseif contains(sessiondata.contextEntry{j,1}, 'BT')
                fixedPoints = [0 0; 45 0; 45 32.5; 0 32.5];
            elseif contains(sessiondata.contextEntry{j,1}, 'track','IgnoreCase',true)
                fixedPoints = [0 0; 121.92 0;121.92 7.62;0 7.62];
                
            end
            tform = fitgeotrans(movingPoints , fixedPoints , "projective");
            
            [edg(:,1),edg(:,2)] = transformPointsForward(tform, edg(:,1)',edg(:,2)');
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
            
            
            
            for ixx = 1:length(LE)
                
                
                if LE(ixx,1) >= (max(edg(:,1))+20) || LE(ixx,1) <= (min(edg(:,1))-20) || LE(ixx,2) >= (max(edg(:,2))+20) || LE(ixx,2) <= (min(edg(:,2))-20)
                    LE(ixx,1) = nan;
                elseif RE(ixx,1) >= (max(edg(:,1))+20) || RE(ixx,1) <= (min(edg(:,1))-20) || RE(ixx,2) >= (max(edg(:,2))+20) || RE(ixx,2) <= (min(edg(:,2))-20)
                    RE(ixx,1) = nan;
                elseif TB(ixx,1) >= (max(edg(:,1))+20) || TB(ixx,1) <= (min(edg(:,1))-20) || TB(ixx,2) >= (max(edg(:,2))+20) || TB(ixx,2) <= (min(edg(:,2))-20)
                    TB(ixx,1) = nan;
                end
                
            end
            
            
            
            %%
            sessiondata.behavior.position.left_ear_cor(kp,:) = LE;
            sessiondata.behavior.position.right_ear_cor(kp,:) = RE;
            sessiondata.behavior.position.tail_base_cor(kp,:) = TB;
            
            sessiondata.contextEntry{j,4} = edg;
            
        else
            sessiondata.behavior.position.left_ear_cor(kp,:) = nan;
            sessiondata.behavior.position.right_ear_cor(kp,:) = nan;
            sessiondata.behavior.position.tail_base_cor(kp,:) = nan;
            
            sessiondata.contextEntry{j,4} = nan;
        end
        
    end
    
    % clean up
    kp = find(cellfun(@isempty,sessiondata.behavior.position.context));
    
    
    
    sessiondata.behavior.position.context(kp) = {'None'};
    
    save('sessiondata.mat','sessiondata','-v7.3')
    
end

end



