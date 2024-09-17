%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: preprocess to select/rearrange channels, then resample to 200Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;
addpath('R:\Analysis\McKenzieLab\detectors\seizure\A_RapidIIIC_pipeline_B');

%% Get the list of files
dataDir = 'R:\Analysis\McKenzieLab\detectors\seizure\nedc_eeg_smile\edf\MAT\';
files = struct2cell(dir([dataDir, '*.mat']))';

%% Define output folder 
trtDir = '.\Data\processed\';
if exist(trtDir, 'dir')~=7
    mkdir(trtDir)
end

%% Load pre-defined mapping table
tmp = load('channel_mappings.mat');
ch_map = lower(tmp.LUT);

%% Main loop for each MAT file 
for i = 1:size(files, 1)
    file = files{i, 1};

    % Load data
    tmp = load([dataDir, file]);
    data0 = tmp.data;
    channels0 = strrep(strrep(tmp.channels, 'EEG ', ''), '-Ref', '');
    channels0 = strrep(lower(channels0), ' ', '');
    Fs0 = round(tmp.Fs);
    startTime = datetime(tmp.startTime);
    
    % Pre-process: channel-select and resample, and denoise
    Fs = 200;
    Fc = 50; % powerline interference frequency, note: it is 60Hz for US data
    [data, channels, ch_miss] = fcn_preprocess(data0, channels0, Fs0, Fs, ch_map, Fc);  
    
    % Export
    if ~isempty(ch_miss)
        disp('Missing channels!')
        keyboard;
    else
        save([trtDir, file], 'data', 'channels', 'Fs', 'startTime', '-v7.3')
    end
end
