%% pulse_control.m
% MATLAB control for Arduino pulse generator

clear;
clc;

%% Serial Configuration

port = "COM3";      % CHECK THIS
baud = 115200;

arduinoObj = serialport(port, baud);

configureTerminator(arduinoObj, "LF");

pause(2); % allow Arduino reset

%% Pulse Parameters

train_duration = 10; % complete duration (seconds)
frequency_Hz = 10;   % pulse repetition frequency
trainLength   = round(frequency_Hz*train_duration);    % number of pulses
duty = 90; % percent
pulseWidth_ms = round((1/frequency_Hz)*(duty/100),2)*1000;   % pulse width in microseconds
pulseWidth_s = pulseWidth_ms/1000;

%% Construct Command String

cmd = sprintf('%d,%.3f,%d\n', ...
    pulseWidth_ms, ...
    frequency_Hz, ...
    trainLength);

fprintf('Sending command:\n%s\n', cmd);
%% Start background flashes
write(arduinoObj, "START_BG\n", "string");

response = readline(arduinoObj);

disp(response)

pause(60)

%% Stop Background flashes 
write(arduinoObj, "STOP_BG\n", "string");

%% Send Command 
write(arduinoObj, cmd, "string"); 

%% Wait for Arduino Response 
response = readline(arduinoObj); fprintf('Arduino response: %s\n', response); %% Cleanup clear arduinoObj;