%% Read Continuous Multi-Channel Data Directly Into Matlab
%
% <html>
% Use SynapseAPI to read multi-channel data buffer into Matlab array <br>
% Online process the data in Matlab <br>
% Requires the StreamMC user gizmo <br>
% </html>

%% Configuration
close all; clear all; clc;
GIZMO_NAME = 'StreamMC1';
HISTORY_SECONDS = 5; % history to keep, in seconds

%% Setup SynapseAPI
syn = SynapseAPI();
fprintf('\n\n')
if syn.getMode < 2
    fprintf('Waiting for Synapse to enter run mode before continuing.\n\n');
    while syn.getMode < 2
        pause(.5)
    end
end

%% The big loop
TIME_ARRAY = zeros(1,50); % holds time delay information
FIRST_PASS = true; % initializes variables when recording starts
ind = 1;
while true
    tic
    if FIRST_PASS
        % Set up variables, determine sampling rate
        SAMP_RATES = syn.getSamplingRates();
        PARENT     = syn.getGizmoParent(GIZMO_NAME);
        BUFF_SIZE  = syn.getParameterValue(GIZMO_NAME, 'size');
        NCHAN      = syn.getParameterValue(GIZMO_NAME, 'numchan');
        DECIMATION = syn.getParameterValue(GIZMO_NAME, 'decimation');
        wav.fs = SAMP_RATES.(PARENT) / DECIMATION;
        fprintf('%d channels in %d sample buffer at %.2f Hz\n', NCHAN, BUFF_SIZE, wav.fs);
        
        SAMPLE_LIMIT = floor(wav.fs / 4);
        SAMPLE_LIMIT = max(1000, SAMPLE_LIMIT - mod(SAMPLE_LIMIT, 1000));
        
        % Fetch the first data points and set up Matlab memory buffer
        curr_index = syn.getParameterValue(GIZMO_NAME, 'index');
        curr_cycle = syn.getParameterValue(GIZMO_NAME, 'loops');
        prev_index = curr_index;
        
        wav.data = zeros(1, NCHAN * round(HISTORY_SECONDS * wav.fs));
        FIRST_PASS = false;
    end
    
    % Look for new data
    curr_index = syn.getParameterValue(GIZMO_NAME, 'index');
    
    %fprintf('curr:%d prev:%d\n', curr_index, prev_index);
    if curr_index ~= prev_index
        if curr_index > prev_index
            npts = curr_index - prev_index;
        elseif prev_index > curr_index
            % buffer wrapped back to the beginning
            % just read up until the end of the buffer this time around
            npts = BUFF_SIZE - prev_index;
            curr_cycle = curr_cycle + 1;
        end
        
        % Read the new data and rotate the Matlab memory buffer
        npts = npts - mod(npts, NCHAN); % make sure we read a multiple of NCHAN
        npts = min(npts, NCHAN*SAMPLE_LIMIT); % read no more than NCHAN*SAMPLE_LIMIT points
        
        temp = syn.getParameterValues(GIZMO_NAME, 'data', npts, prev_index)';
        wav.data = circshift(wav.data, -npts);
        wav.data(end-npts+1:end) = temp;

        % DO PROCESSING HERE
        ddd = reshape(wav.data, NCHAN, []);
        curr_time = (curr_index + curr_cycle * BUFF_SIZE) / wav.fs / NCHAN;
        size_time = numel(wav.data) / wav.fs / NCHAN;
        ts = linspace(curr_time-size_time, curr_time, numel(wav.data) / NCHAN);
        
        % Update TDT buffer index variable for next loop
        prev_index = prev_index + npts;
        if prev_index >= BUFF_SIZE
            prev_index = prev_index - BUFF_SIZE;
            fprintf('buffer loop\n');
        end
    else
        if syn.getMode < 2
            fprintf('Waiting for Synapse to enter run mode before continuing.\n\n');
            while syn.getMode < 2
                pause(.5)
            end
            FIRST_PASS = true;
        end
    end
    
    % Get/show delay stats
    TIME_ARRAY = circshift(TIME_ARRAY, -1);
    TIME_ARRAY(end) = toc;
    ind = ind + 1;
    if mod(ind, numel(TIME_ARRAY)) == 0
        fprintf('max delay: %.2fms min delay: %.2f, mean: %.2fms\n', max(TIME_ARRAY)*1000, min(TIME_ARRAY)*1000, mean(TIME_ARRAY)*1000)
    end
end