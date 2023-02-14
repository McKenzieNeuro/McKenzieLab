function [mat,v2] = Epoch2Mat(ts,epochs,varargin)
v2 =[];
epochs1 = MergeEpochs_local(epochs);

if size(epochs1,1)~=size(epochs,1) || ~all(epochs(:) ==epochs1(:))
    error('overlapping or non sorted epochs, epochs must be in order and non-overlapping')
end


[status, interval] = InIntervals(ts,epochs);
ts1 = ts(status);

n1 = histoc(interval(interval>0),1:size(epochs,1));
mat = mat2cell( ts1(:),n1);
if ~isempty(varargin)
    v2 = varargin{1};
    v2 = v2(status);
    v2 = mat2cell( v2(:),n1);
end
end

function epochs = MergeEpochs_local(epochs, mergetouching)
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

end