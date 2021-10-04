%% Online Signal Averaging Example
%
% <html>
% Import strobe store gizmo data into Matlab using SynapseLive during the
% experiment <br>
% Plot the average waveform <br>
% Good for Evoked Potential visualization <br>
% </html>

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
% path.
close all; clc;
[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Variable Setup
% Set up the variables for the data you want to extract. We will extract
% a single channel from a fixed duration strobed storage gizmo.
EVENT = 'StS1';
CHANNEL = 1;

%%
% show the last N waveforms in the plot.
N = 5; 

%%
% Set KEEPALL to 0 to only show the running average of the last N waveforms.
% Otherwise, all waveforms in the block are included in the average.
KEEPALL = 0; 

%%
% Setup SynapseLive
t = SynapseLive('MODE', 'Preview', 'EXPERIMENT', 'OnlineAveragingDemo'); % we will default to 'Preview' mode
t.TYPE = {'snips'}; % we only care about these types of events
t.VERBOSE = false;
first_pass = true;

%% The Main Loop
prevWaves = cell(1,N);
nsweeps = 0;
while 1
    
    % slow it down a little
    pause(1)
    
    % get the most recent data, exit loop if the block has stopped
    if isempty(t.update)
        break
    end

    % read the snippet events
    r = t.get_data(EVENT);
    if isstruct(r)
        if ~isnan(r.data)
            % get our channel of data
            chan_data = r.data(r.chan == CHANNEL,:);
            nsize = size(chan_data,1);
            
            % cache the waveforms in our circular buffer
            prevWaves = circshift(prevWaves, -nsize);
            for i = 1:(min(nsize, N))
                prevWaves{i} = chan_data(end-(i-1),:);
            end
            
            % find average signal
            cache_ind = ~cellfun('isempty', prevWaves);
            if KEEPALL == 0
                % if we are only keeping the previous N, do average on just those
                avg_data = mean(cell2mat(prevWaves(cache_ind)'), 1);
            else
                if first_pass
                    first_pass = false;
                    nsweeps = nsize;
                    avg_data = new_mean;
                else
                    new_mean = mean(chan_data, 1);
                    % add new average into the old average
                    avg_data = (avg_data .* nsweeps + new_mean * nsize) / (nsweeps + nsize);
                end
            end
            
            nsweeps = nsweeps + nsize;
            
            % plot the preview N waves in gray
            t_ms = 1000*(1:numel(avg_data)) / r.fs;
            plot(t_ms, cell2mat(prevWaves(cache_ind)')','Color', [.85 .85 .85]); hold on;
            
            % plot the average signal in thick blue
            plot(t_ms, avg_data, 'b', 'LineWidth', 3); hold off;
            
            % finish up plot
            title(sprintf('nsweeps = %d, last %d shown', nsweeps, N));
            xlabel('Time, ms','FontSize',12)
            ylabel('V', 'FontSize', 12)
            temp_axis = axis;
            temp_axis(1) = t_ms(1);
            temp_axis(2) = t_ms(end);
            axis(temp_axis);
            
            % force the plots to update
            try
                snapnow
            catch
                drawnow
            end
            
            % for publishing, end early
            if exist('quitEarly','var') && nsweeps > 30
                break
            end
        end
    end
end