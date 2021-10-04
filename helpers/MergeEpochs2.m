function epochs = MergeEpochs2(epochs, mergetouching)
%
% merge epochs of arbirary precision without creating binary vector
%
% returned sorted
%

if ~exist('mergetouching', 'var'), mergetouching = 1; end

if size(epochs,1)<=1, return; end

if any(epochs(:,2)-epochs(:,1)<0), error('epochs must be zero or positive windows of time like [tstart tstop]'); end

[~, inds] = sort(epochs(:,1));



epochs = epochs(inds,:); % sort epochs by start time

i = 1;

while i < size(epochs,1)
    if mergetouching
        if epochs(i,2) >= epochs(i+1,1)
            
            epochs(i, 2) = epochs(i+1, 2);
            
            epochs(i+1, :) = [];
            
        else
            
            i = i+1;
            
        end
    else
        if epochs(i,2) > epochs(i+1,1)
            
            epochs(i, 2) = epochs(i+1, 2);
            
            epochs(i+1, :) = [];
            
        else
            
            i = i+1;
            
        end
        
    end
end
    
    % size_limit = 10000;
    %
    % if size(epochs, 1) > size_limit, epochs = ForBigEpochsLists(epochs, mergetouching); return; end
    %
    % [~, inds] = sort(epochs(:,1));
    %
    % epochs = epochs(inds,:); % sort epochs by start time
    %
    % S = repmat(epochs(:,1), 1, length(epochs));
    %
    % E = repmat(epochs(:,2), 1, length(epochs));
    %
    % eS = repmat(epochs(:,1)', length(epochs), 1);
    %
    % eE = repmat(epochs(:,2)', length(epochs), 1);
    %
    % if mergetouching
    %
    %     sTF = eS >= S & eS <= E; % start times which exist within another epoch
    %
    %     sTF(1:length(sTF)+1:end) = 0; % elimate diagonal
    %
    %     eTF = eE <= E & eE >=S;
    %
    %     eTF(1:length(sTF)+1:end) = 0; % end times which exist within another epoch
    %
    % else
    %
    %     sTF = eS > S & eS < E; % start times which exist within another epoch
    %
    %     sTF(1:length(sTF)+1:end) = 0; % elimate diagonal
    %
    %     eTF = eE < E & eE > S;
    %
    %     eTF(1:length(sTF)+1:end) = 0; % end times which exist within another epoch
    %
    % end
    %
    % eDEL = sum(eTF,1)>=1; % which epochs ends need to be deleted
    %
    % sDEL = sum(sTF,1)>=1; % which epochs ends need to be deleted
    %
    % epochs = [epochs(~sDEL,1), epochs(~eDEL,2)];
    %
    % function epochs = ForBigEpochsLists(epochs, mergetouching)
    %
    % [~, inds] = sort(epochs(:,1));
    %
    % if mergetouching
    %
    %     op = {'>=', '<='};
    %
    % else
    %
    %     op = {'>', '<'};
    %
    % end
    %
    % epochs = epochs(inds,:); % sort epochs by start time
    %
    % last_i = 1;
    %
    % i = 1;
    %
    % while i < length(epochs)
    %
    %     if eval(['epochs(i,2)' op{1} 'epochs(i+1,1)'])
    %
    %         epochs(i, 2) = epochs(i+1, 2);
    %
    %         epochs(i+1, :) = [];
    %
    %     else
    %
    %         i = i+1;
    %
    %     end
    %
    % end
    %
    % end