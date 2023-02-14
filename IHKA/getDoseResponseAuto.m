DirectoryName = 'R:\McKenzieLab\DGregg\NeuralData\EDS\OL\Week2';
files = getAllExtFiles(DirectoryName,'mat',1);
keepFiles = contains(files,'recInfo');
files = files(keepFiles);
for ii = 1:length(files)
    files(ii) = erase(files(ii),'\recInfo.mat');
end

xCategories = categorical(["pre-stim","stim incl.","stim excl.","post-stim"]); % create categorical values for the x axis (note: automatically alphabetizes)
xCategories = reordercats(xCategories, string(xCategories)); % re-order the categories so they are no longer in alphabetical order

for kk = 1:length(files)
    dir = files(kk);
    cd(dir{1});
    clear recInfo
    load('recInfo.mat')
    ISI = load('config.mat','ISI');
    ISI = mean(ISI.ISI);
    clear TDdata
    TDdata = struct;
    %% Loop through animals (only two)
    for ii = 1:2
        %% Loop through channels (8)
        if ~isempty(recInfo{3,3}{ii+1,6})
            for jj = 1:length(recInfo{3,3}{1+ii,end})
                channelName = recInfo{3,3}{1+ii,end}{jj,2};
                subjectName = recInfo{3,3}{1+ii,2}{1};
                TD = recInfo{3,3}{1+ii, end}{jj,3}; % take data from rec info. Cycle through the animals and the channels
                TD = TD(~isnan(TD));
                ts = 1:length(TD);
                stimON = recInfo{3,3}{1+ii, 6}/1000; % stim data encoded in ms. convert to seconds

                %% find the theta delta ratio while ignoring the stim and 2 seconds following stim
                epoch = [stimON(1:end-1)+2 stimON(2:end)-0.1; stimON(end)+2 stimON(end)+10]; % select the times 2 seconds after stim until 100 ms before the next stim
                keep = InIntervals(ts, epoch); % keep only the times that overlap with the select inter-stim window
                u_TD_woStim = mean(TD(keep)); % average the TD ratio for the selected times

                %% find the theta delta ratio including the stim window
                epoch_stim = [stimON(1:end-1) stimON(2:end)]; % select times that are within the stim period. do not reject the stim window
                keep_stim = InIntervals(ts, epoch_stim); % keep only the times that overlap with the stimulation period
                u_TD_wStim = mean(TD(keep_stim)); % average the TD ratio for the stim window

                %% pre-baseline theta delta ratio
                preBLend = stimON(1)-.1; % select times that occur up to 100 ms before the first stimulation
                pre_keep = InIntervals(ts,[ts(1) preBLend]); % keep only the times that overlap with the pre-stimulation baseline
                pre_TD = mean(TD(pre_keep)); % average the TD ratio for the pre-stim baseline

                %% post-baseline theta delta ratio
                postBLstart = stimON(end)+10; % select times that occur 10 seconds or more after the final stim
                post_keep = InIntervals(ts, [postBLstart ts(end)]); % keep only the times that overlap with the post-stimulation baseline
                post_TD = mean(TD(post_keep)); % average the TD ratio for the post-stim baseline

                %% get the stim information
                stimChannels = recInfo{3,3}{1+ii,4}; % load the active stimulation channels for comparison
                if stimChannels(2) == 0
                    stimType = 'mono';
                elseif stimChannels(2) == 8
                    stimType = 'bi';
                else
                    stimType = 'error in code';
                end
                stimAmp = recInfo{3,3}{1+ii,3}; % load the stimulation amplitude

                %% save the data
                TDdata(ii).subject = subjectName;
                TDdata(ii).data(jj).channel = channelName;
                TDdata(ii).data(jj).preTD = pre_TD;
                TDdata(ii).data(jj).TDwithoutStim = u_TD_woStim;
                TDdata(ii).data(jj).TDwithStim = u_TD_wStim;
                TDdata(ii).data(jj).postTD = post_TD;
                TDdata(ii).StimType = stimType;
                TDdata(ii).StimAmp = stimAmp;

            end
            
        end
    end
    save('adjustedThetaDelta.mat',"TDdata");
    close all
    for ii = 1:length(TDdata)
        if ~isempty(TDdata(ii).data)
            figure(ii)
            hold on
            titlestr = sprintf('%s at %spolar stim intensity %duA with %d sec ISI',TDdata(ii).subject,TDdata(ii).StimType,TDdata(ii).StimAmp, ISI);
            title(titlestr)
            %% Loop through channels (8)
            for jj = 1:length(TDdata(ii).data)
                plot(xCategories, [TDdata(ii).data(jj).preTD/TDdata(ii).data(jj).preTD, TDdata(ii).data(jj).TDwithStim/TDdata(ii).data(jj).preTD, TDdata(ii).data(jj).TDwithoutStim/TDdata(ii).data(jj).preTD, TDdata(ii).data(jj).postTD/TDdata(ii).data(jj).preTD],'DisplayName',TDdata(ii).data(jj).channel)
                ylabel('Theta/Delta')
                legend('Location','northeastoutside') % don't block the data with the legend!!!!

            end
            set(gca,'linestyleorder',{'-','--','-.'},'colororder',[1 0 0; 0 0 1; 0 0 0]) % specify the line colors and line style. The matlab presets suck.
            savename = sprintf('%s normalized.png',titlestr); % create a file name to save the figure.
            saveas(gcf,savename)
        end
    end
end
















%% function InIntervals

function [status,interval,index] = InIntervals(values,intervals,varargin)

%InIntervals - Test which values fall in a list of intervals.
%
%  USAGE
%
%    [status,interval,index] = InIntervals(values,intervals,<options>)
%
%    values         values to test (these need not be ordered)
%    intervals      list of (start,stop) pairs
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'verbose'     display information about ongoing processing
%                   (default = 'off')
%    =========================================================================
%
%  OUTPUT
%
%    status         logical indices (1 = belongs to one of the intervals,
%                   0 = belongs to none)
%    interval       for each value, the index of the interval to which
%                   it belongs (0 = none)
%    index          for each value, its index in the interval to which
%                   it belongs (0 = none)
%
%  NOTE
%
%    If the intervals overlap, the outputs 'interval' and 'index' refer to the
%    last overlapping interval (i.e. if one value belongs to intervals #7 and #8,
%    it will be listed as belonging to interval #8).
%
%  SEE
%
%    See also ConsolidateIntervals, SubtractIntervals, ExcludeIntervals,
%    Restrict, FindInInterval, CountInIntervals, PlotIntervals.

% Copyright (C) 2004-2011 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
verbose = false;

% Check number of parameters
if nargin < 2 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).');
end

% Check parameters
if isempty(intervals)
    status = logical(zeros(size(values)));
    interval = zeros(size(values));
    index = zeros(size(values));
    return
end


intervals = double(intervals);
values = double(values);
if ~isdmatrix(intervals) || size(intervals,2) ~= 2,
  error('Incorrect intervals (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).');
end

if isempty(values)
    warning('values is an empty vector, returning nothing..')
    status=[];
    interval=[];
    index=[];
    return
end

if size(values,1) == 1,
	values = values(:);
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'verbose',
			verbose = varargin{i+1};
			if ~isstring_FMAT(verbose,'on','off'),
				error('Incorrect value for property ''verbose'' (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).');
			end
			verbose = strcmp(verbose,'on');

		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help InIntervals">InIntervals</a>'' for details).']);
	end
end

[values,order] = sortrows(values(:,1));

% Max nb of digits for display
l = int2str(floor(log10(max(max(intervals*100))))+2);

% Determine if intervals overlap (in which case we must use a 'slow' algorithm)
[intervals,intorder] = sortrows(intervals,1);
di = intervals(2:end,1)-intervals(1:end-1,2);
overlap = any(di<0);
if ~overlap,
	% Fast algorithm: for the next interval, start from the end of the previous interval
	k = 2;
else
	% Slow algorithm: for the next interval, start from the beginning of the previous interval
	k = 1;
end

% Retrieve values in intervals
previous = 1;
n = size(intervals,1);
status = logical(zeros(size(values)));
interval = zeros(size(values));
index = zeros(size(values));
times = values;
for i = 1:n,
	from = intervals(i,1);
	to = intervals(i,2);
	timeString = sprintf(['%' l '.2f %' l '.2f (%' l '.2f)'],from,to,to-from);
	% Get values
	more = FindInInterval(values,[from to],previous);
	if ~isempty(more),
		previous = more(k); % See note above about algorithm
		nMore = more(2)-more(1)+1;
		interval(more(1):more(2)) = intorder(i);
		status(more(1):more(2)) = 1;
		index(more(1):more(2)) = (1:nMore);
	end
	if verbose, disp([timeString ' - ' int2str(nMore) ' values']); end
end

status(order) = status;
interval(order) = interval;
index(order) = index;
end

function fileList = getAllExtFiles(dirName,ext,varargin)

if ~isempty(varargin)
    rec = varargin{1};
else
    rec = 1;
end

dirData = dir(dirName);      %# Get the data for the current directory
dirIndex = [dirData.isdir];  %# Find the index for directories
fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
fileList = fileList(cellfun(@length,fileList)>3);
fileList=fileList(cellfun(@(a) all(a(end-2:end)==ext),fileList)); %only take matlab files
if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
        fileList,'UniformOutput',false);
end
subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
if rec                                          %#   that are not '.' or '..'
    for iDir = find(validIndex)                  %# Loop over valid subdirectories
        nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
        fileList = [fileList; getAllExtFiles(nextDir,ext,1)];  %# Recursively call getAllFiles
    end
end
end
