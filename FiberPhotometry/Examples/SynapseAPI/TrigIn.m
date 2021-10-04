%% Trigger Input
%
% <html>
% This examples shows how to start a Synapse recording using an external
% trigger into your processor <br>
% Matlab monitors the bit input and controls Synapse state <br>
% </html>

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
close all; clear all; clc;
[APIEXAMPLESPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
[SDKPATH,name,ext] = fileparts(APIEXAMPLESPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Setup
% This example uses the 'TrigIn' experiment included in the download.
syn = SynapseAPI();
EXPERIMENT = 'TrigIn';
syn.setCurrentExperiment(EXPERIMENT);

%% Runtime
% Move into Standby Mode to wait for external trigger.
% Note: You must have Standby Mode enabled in Menu > Preferences.
syn.setMode(1);

%%
% Wait for external trigger to move to Preview mode
prevTime = 0; tic;
while syn.getParameterValue('Trigger1','TrigIn') == 0
    currTime = floor(toc);
    if currTime > prevTime
        fprintf('Waiting %d seconds for trigger\n', currTime);
    end
    prevTime = currTime;
end
fprintf('%.1f seconds elapsed.\nGo!\n', toc);
syn.setModeStr('Preview');

% Experiment Loop
pause(5);

% Move from Preview to Idle Mode when done
syn.setModeStr('Idle');

