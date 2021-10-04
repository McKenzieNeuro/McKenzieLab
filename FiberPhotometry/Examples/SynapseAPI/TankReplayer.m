%% Export Continuous Data To Corpus Emulator
%
% <html>
% Import continuous data into Matlab using TDTbin2mat <br>
% Export to binary files that are read into Corpus as PZ5 data <br>
% Use SynapseAPI to generate new recordings with this data <br>
% Concatenate the resulting data sets into one master data structure <br>
% </html>

%% Housekeeping
% Clear workspace and close existing figures. Add SDK directories to Matlab
% path.
close all; clear all; clc;
[MAINEXAMPLEPATH,name,ext] = fileparts(cd); % \TDTMatlabSDK\Examples
DATAPATH = fullfile(MAINEXAMPLEPATH, 'ExampleData'); % \TDTMatlabSDK\Examples\ExampleData
[SDKPATH,name,ext] = fileparts(MAINEXAMPLEPATH); % \TDTMatlabSDK
addpath(genpath(SDKPATH));

%% Importing the Data
% This example assumes you downloaded our example data sets
% (<https://www.tdt.com/files/examples/TDTExampleData.zip link>) and extracted
% it into the \TDTMatlabSDK\Examples\ directory. To import your own data, replace
% |BLOCKPATH| with the path to your own data block.
%
% In Synapse, you can find the block path in the database. Go to Menu > History. 
% Find your block, then Right-Click > Copy path to clipboard.
BLOCKPATH = fullfile(DATAPATH,'Algernon-180308-130351');

%% Setup the variables for the data you want to extract
% We will extract a stream store and output it to a Mat file
STORE = 'Wav1';
% the circuit loads and begins running before Synapse switches to record
% mode, so add 5 seconds of zeros to the front of it to make sure we are
% recording all of our data.
PADDING = 5; 

%% Experiment setup
% This example uses the '_Demo_TankReplayer' experiment included in the download.
% In the RZ HAL, PortC0 is enabled as an input and its epoc store is set to
% 'Full'. This is critical to capture when the PZ data is valid. This will
% create the PC0_ epoc event that we can time filter on during extraction.
EXPERIMENT = '_Demo_TankReplayer';
TARGET_DEVICE = 'PZ5'; % PZ5 or PZ2
if ~strcmp(EXPERIMENT(1:6), '_Demo_')
    error('Experiment name must start with _Demo_ to automatically import Corpus data')
end

% Connect to Synapse
syn = SynapseAPI();

% Set appropriate scale factor depending on our device type
if strcmpi(TARGET_DEVICE, 'PZ5')
    SCALE = 1e9 / 4;
elseif strcmpi(TARGET_DEVICE, 'PZ2') || strcmpi(TARGET_DEVICE, 'PZ3') || strcmpi(TARGET_DEVICE, 'PZ4')
    SCALE = 1e9;
end

%% 
% Read the headers to see how much data we have
heads = TDTbin2mat(BLOCKPATH, 'TYPE', {'streams'}, 'STORE', STORE, 'HEADERS', 1);

% get the total data size first
nchan = double(max(heads.stores.(STORE).chan));
npts = double(heads.stores.(STORE).size - 10);

MAX_SAMPLES = 256e6 / 2; % Corpus can load up to 256M samples total. Cut it in half to be safe.

% leave room for our digital trigger as well
nsamples = npts * double(numel(heads.stores.(STORE).chan)) * (nchan+1) / nchan;
FS = heads.stores.(STORE).fs;
ZEROPAD = zeros(nchan, ceil(PADDING * FS));

% if it's too big, we need to split it up
iter = ceil((nsamples + numel(ZEROPAD)) / MAX_SAMPLES);

% figure out what T1 and T2 should be on each iteration
T1 = 0;
% subtract our pad off of it, plus a little extra to be safe
DELTA = MAX_SAMPLES / FS / (nchan + 1) - (PADDING + 0.5);

FOLDERNAME = EXPERIMENT(7:end);
ROOT = sprintf('C:\\TDT\\Corpus\\Simulation Files\\Synapse Demos\\%s', FOLDERNAME);

blockPaths = cell(1,iter);
for ii = 1:iter
    % first set it to the default experiment
    syn.setCurrentExperiment('Experiment');
    refreshExperiment = tic;
    
    % read our data into Matlab and create MAT file for PZ5 emulated data
    data = TDTbin2mat(BLOCKPATH, 'TYPE', {'streams'}, 'STORE', STORE, 'HEADERS', heads, 'T1', T1, 'T2', T1 + DELTA);
    folder = fullfile(ROOT, 'PZ1');
    mkdir(folder);
    OUTFILE = fullfile(folder, 'chan.mat');
    arr = [ZEROPAD data.streams.(STORE).data] * SCALE;
    
    % log it
    DURATION = size(arr, 2) / FS;
    fprintf('Wrote %s to output file %s\n', STORE, OUTFILE);
    fprintf('Sampling Rate:\t%.6f Hz\n', FS);
    fprintf('Num Channels:\t%d\n', size(data.streams.(STORE).data, 1));
    fprintf('Duration:\t\t%.1fs (includes padding)\n', DURATION);

    % make the digital input trigger that tells when the data is valid
    bits = [zeros(1, size(ZEROPAD, 2)) ones(1, size(data.streams.(STORE).data,2))];
    
    clear data;
    save(OUTFILE, 'arr');
    clear arr;
    
    % save the digital input trigger data into Corpus directory
    folder = fullfile(ROOT, 'RZ1');
    mkdir(folder);
    OUTFILE = fullfile(folder, 'bits.i16');
    fileID = fopen(OUTFILE,'w');
    fwrite(fileID, bits, 'int16');
    fclose(fileID);
    fprintf('Wrote Sync Bit0 to output file %s\n', OUTFILE);
    clear bits;
    
    % Reload your experiment to force Corpus to reload the new data
    % then to our experiment
    
    % make sure enough time has elapsed that Corpus detected the change
    while tic - refreshExperiment < 3
        pause(.1);
    end
    
    % load our experiment so Corpus catches it
    syn.setCurrentExperiment(EXPERIMENT);
    
    disp('Wait until data has loaded into Corpus (estimated)')
    pause(max(DURATION/10, 5))
    syn.setModeStr('Record');
    
    fprintf('Wait until data has completely cycled through (%d seconds)\n', ceil(DURATION))
    pause(PADDING);
    while syn.getSystemStatus.iRecordSecs < ceil(DURATION)
        pause(1)
    end
    
    % store this block name
    thisBlockPath = fullfile(syn.getCurrentTank, syn.getCurrentBlock);
    fprintf('Done recording into: %s\n', thisBlockPath)
    blockPaths{ii} = thisBlockPath;
    
    % get ready for next iteration
    syn.setModeStr('Idle');
    T1 = T1 + DELTA;
end

%% Import the new data back into Matlab
%
disp('Reading the data back into Matlab')
for ii = 1:numel(blockPaths)
    tempdata = TDTbin2mat(blockPaths{ii}, 'TYPE', {'epocs'}); % get just the epoc events
    
    % we just want the first valid epoc event
    tr = [tempdata.epocs.PC0_.onset(1); tempdata.epocs.PC0_.offset(1)];
    
    % extract our data, just on that time range
    tempdata = TDTbin2mat(blockPaths{ii}, 'RANGES', tr);
    
    % your custom analysis code here
end