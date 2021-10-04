%% Online Raster PSTH Example
%
% <html>
% Import snippet and epoc data into Matlab using SynapseLive while the
% experiment is running <br>
% Generate peri-event raster and histogram plots over all trials <br>
% Good for stim-response experiments, such as optogenetic or electrical
% stimulation, where you need immediate visual feedback <br>
% </html>

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
% path.
close all; clc;
[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Setup
% Setup the variables for the data you want to extract. We will extract
% channel 1 from the eNe1 snippet data store, created by the PCA Sorting
% gizmo, and use the Tick as our stimulus onset.
REF_EPOC = 'Tick';
EVENT = 'eNe1';
CHANNEL = 1;
SORTCODE = 0; % set to 0 to use all sorts
TRANGE = [-0.3, 0.8]; % start time, duration
DO_RASTER = 1; % set to 0 to only see histogram

%% Setup SynapseLive
t = SynapseLive('MODE', 'Preview', 'EXPERIMENT', 'RasterPSTHdemo'); % we will default to 'Preview' mode
t.NEWONLY = 0;  % read all available events in the block on every iteration
t.TIMESTAMPSONLY = 1;  % we don't care what the snippets look like, just their timestamps
t.TYPE = {'snips', 'epocs', 'scalars'}; % we only care about these types of events
t.VERBOSE = false;

%% The Main Loop

% Set figure size base on number of plots
if DO_RASTER
    h = figure('Position',[100, 100, 500, 800]);
else
    h = figure('Position',[100, 100, 500, 400]);
end
while 1
    
    % slow it down a little
    pause(1)
    
    % get the most recent data, exit loop if the block has stopped.
    if isempty(t.update)
        break
    end

    % read the snippet and event timestamp data.
    r = t.get_data(EVENT);
    if isstruct(r)
        if ~isnan(r.ts)
            
            % do our timestamp filtering
            if DO_RASTER
                data = TDTfilter(t.data, REF_EPOC, 'TIME', TRANGE);
            else
                data = TDTfilter(t.data, REF_EPOC, 'TIME', TRANGE, 'TIMEREF', 1);
            end
            
            % do our channel and sort code filtering
            if SORTCODE ~= 0
                i = find(data.snips.(EVENT).chan == CHANNEL & data.snips.(EVENT).sortcode == SORTCODE);
            else
                i = find(data.snips.(EVENT).chan == CHANNEL);
            end
            
            % grab only the relevant timestamps
            try
                TS = data.snips.(EVENT).ts(i);
                if isempty(TS)
                    continue
                end
            catch
                continue
            end
            
            % that's it for the data extraction, now we plot
            num_trials = size(data.time_ranges, 2);
            if DO_RASTER
                % match timestamp to its trial
                all_TS = cell(num_trials, 1);
                all_Y = cell(num_trials, 1);
                for trial = 1:num_trials
                    trial_TS = TS(TS >= data.time_ranges(1, trial) & TS < data.time_ranges(2, trial));
                    all_TS{trial} = trial_TS - data.time_ranges(1, trial) + TRANGE(1);
                    all_Y{trial} = trial * ones(numel(trial_TS), 1);
                end
                all_X = cat(1, all_TS{:});
                all_Y = cat(1, all_Y{:});

                % plot raster
                subplot(2,1,1)
                hold on;
                plot(all_X, all_Y, '.', 'MarkerEdgeColor','k', 'MarkerSize',10)
                line([0 0], [1, trial], 'Color','r', 'LineStyle','-', 'LineWidth',3)
                axis tight; axis square;
                set(gca, 'XLim', [TRANGE(1), TRANGE(1)+TRANGE(2)]);
                ylabel('trial number')
                xlabel('time, s')
                title(sprintf('Raster ch=%d sort=%d, %d trials', CHANNEL, SORTCODE, num_trials))
                hold off;
                TS = all_X;
                subplot(2,1,2)
            end
            
            % plot PSTH
            NBINS = floor(numel(TS)/10);
            if NBINS > 0
                hist(TS, NBINS);
                hold on;
                N = hist(TS, NBINS);
                line([0 0], [0, max(N)*1.1], 'Color','r', 'LineStyle','-', 'LineWidth',3);
                hold off;
                axis tight; axis square;
                set(gca, 'XLim', [TRANGE(1), TRANGE(1)+TRANGE(2)]);
                ylabel('number of occurrences')
                xlabel('time, s')
                title(sprintf('Histogram ch=%d sort=%d, %d trials', CHANNEL, SORTCODE, num_trials))
            end
            
            % force the plots to update
            try
                snapnow
            catch
                drawnow
            end
            
            % for publishing, end early
            if exist('quitEarly','var') && num_trials > 30
                break
            end
        end
    end
end