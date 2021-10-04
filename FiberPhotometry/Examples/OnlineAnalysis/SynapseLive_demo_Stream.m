%% Stream Example
%
% <html>
% Read a streaming data store into Matlab from Tank server during a recording <br>
% </html>

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
% path.
close all; clc;
[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Setup
EVENT = 'EEG1';
t = SynapseLive('MODE', 'Preview', 'EXPERIMENT', 'OnlineStreamDemo');
t.VERBOSE = false;

%% Main Loop
first_pass = true;
while 1
    
    % slow it down
    pause(1)
    
    % get the most recent data, exit loop if the block has stopped.
    if isempty(t.update)
        break
    end

    % grab the latest events
    r = t.get_data(EVENT);
    if isstruct(r)
        
        % plot them
        ts = linspace(t.T1, t.T2, max(size(r.data))-1);
        plot(ts, r.data(:,1:end-1)')
        title(EVENT)
        xlabel('Time, s')
        ylabel('V')
        axis tight
        
        % force the plots to update
        try
            snapnow
        catch
            drawnow
        end
    end
    
    % for publishing, end early
    if exist('quitEarly','var') && t.T2 > 30
        break
    end
end