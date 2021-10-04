%% Run Experiment For Set Duration
%
% <html>
% This examples shows how to use Matlab to control the Synapse mode and
% monitor recording status <br>
% This uses the 'ExperimentSetTime' example experiment, but it could run with any
% experiment <br>
% </html>

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
close all; clear all; clc;
[APIEXAMPLESPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
[SDKPATH,name,ext] = fileparts(APIEXAMPLESPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Setup
% Choose which experiment to run and the duration. It could be anything.
% This example uses a simple experiment with just a Tick store.
% To see full list of available experiments use syn.getKnownExperiments()
EXPERIMENT = 'ExperimentSetTime';
TOTAL_TIME = 30;

% Connect to Synapse
syn = SynapseAPI();

% Set your experiment
syn.setCurrentExperiment(EXPERIMENT);
 
%% Runtime
% Set the system to 'Preview' mode
syn.setModeStr('Preview');

% Wait five seconds to give 'getSystemStatus' time to update internally
pause(5);

%% Main Loop
currTime = 0; prevTime = 0;
% Poll the system status until it reaches the desired state
while currTime < TOTAL_TIME
    % Add any additional API controls here
    currTime = syn.getSystemStatus.iRecordSecs;
    if prevTime ~= currTime
        fprintf('Current elapsed time: %ds\n', currTime);
    end
    prevTime = currTime;
end

%%
% Our desired elapsed time has passed, switch to Idle mode
syn.setModeStr('Idle');
disp('done');
