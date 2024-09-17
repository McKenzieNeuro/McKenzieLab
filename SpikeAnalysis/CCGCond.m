

function [ccg,n,t] = CCGCond(times,groups,conditions,varargin)

%CCG - Compute multiple cross- and auto-correlograms
%
%  USAGE
%
%    [ccg,n,t] = CCGCond(times,groups,conditions, <options>)
%    ccg = NGroups x Ngroups x NConditions x nBins - coincidence counts
%    n = NGroups x Ngroups x NConditions x nBins - # reference observations
%    t = nBins - temporal lag

%    times          times of all events (sorted)
%    groups         group IDs for each event in time list
%    conditions     categorical label for each time point
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'        bin size in s (default = 0.01)
%     'duration'       duration in s of each xcorrelogram (default = 2)
%     'Fs'             sampling rate of time series (defualt = 1/20000)
%     'across_groups'  logical flag for computing CCG across conditions
%    =========================================================================
%
%  

% Default values
duration = 2;
binSize = 0.01;
Fs = 1/20000;
across_groups = true;
% Check parameters
if nargin < 3,
    error('Incorrect number of parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end


if ~isdscalar(groups) && ~isdvector(groups),
    error('Parameter ''groups'' is not a real-valued scalar or vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(groups) && length(times) ~= length(groups),
    error('Parameters ''times'' and ''groups'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end


groups = groups(:);
times = times(:);
conditions = conditions(:);


% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'binsize'
            binSize = varargin{i+1};
           
        case 'across_groups'
            across_groups = varargin{i+1};
    
        case 'duration'
            duration = varargin{i+1};
            if ~isdscalar(duration,'>0')
                error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
            end
            
        case 'Fs'
            Fs = varargin{i+1};
            if ~isdscalar(Fs,'>0')
                error('Incorrect value for property ''Fs'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
            end
            
    end
end



% Number of groups, number of bins, etc.
if length(groups) == 1,
    groups = ones(length(times),1);
    nGroups = 1;
else
    nGroups = max(unique(groups));
end


halfBins = round(duration/binSize/2);
nBins = 2*halfBins+1;

times = round(times/Fs);
binSize = round(binSize/Fs);

if length(times) <= 1,
    % ---- MODIFIED BY EWS, 1/2/2014 ----
    % *** Use unsigned integer format to save memory ***
    ccg = uint16(zeros(nBins,nGroups,nGroups));
    % -----------------------------------
    return
end

% Compute CCGs
nEvents = length(times);
%
if across_groups
    %the CCG will be calculated separately for the condition of each
    %reference/target combination
    
    [counts,cnt] = (CCGHeartCrossGroup(times,uint32(groups),uint32(conditions),binSize,uint32(halfBins)));
    nGroup = max(conditions)*max(conditions);
else
    %the CCG will be calculated accorgint to just the reference condition
    [counts,cnt] = (CCGHeartWithinGroup(times,uint32(groups),uint32(conditions),binSize,uint32(halfBins)));
    nGroup = max(conditions);
end
% -----------------------------------
%
% Reshape the results
counts = double(counts);
n = max(groups);


t = Fs*((-halfBins:halfBins)'*binSize);


cnt = reshape(cnt,[n nGroup]);
counts = reshape(counts,[nBins n n nGroup]);


ccg = permute(counts,[ 1 3 2 4]);
%counts = flipud(counts);



n = nan(size(ccg));
for i = 1:size(cnt,2)
    [~,nn2] = meshgrid(cnt(:,i));
    
    n(:,:,:,i) = permute(repmat(nn2,1,1,size(ccg(:,:,:,i),1)),[3 1 2]);
    
    
end


