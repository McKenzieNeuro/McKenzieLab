function Y = getXPctTim(X,prc,Tim)
% function gets random time bins within the bounds of X

% Inputs
% X = [M x 2], M(:,1) = start, M(:,2) = end
% prc = percentage of time point to take in this window
% Tim = resolutation of time bins

% Returns
% Y = [M x N] matrix,  M is size(X,2) if X is not nan, N = number of time bins

%%
% check inputs are valid

% return empty matrix if no input
if all(isnan(X(:)))
    Y = [];
    return
end

% only consider inputs with onset and offset
X(any(isnan(X),2),:) = [];

if isempty(X)
    Y = [];
    return
end

%%

%calculate the number of events to return
n = max(ceil(floor(diff(X,[],2)/Tim)*prc));
Y = nan(size(X,1),n);

% loop over event
for i = 1:size(X,1)
    
    %potential time bins
    tims = X(i,1):Tim:X(i,2)-Tim;
    
    %numer of events to return
    n = ceil(length(tims)*prc);
    
    %randomly chosen events
    ix = randsample(1:length(tims),n);
    
    %save events to return
    Y(i,1:length(ix)) = tims(ix)';
end

end

