%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: Compute spectrograms                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear

addpath('.\Tools\qEEG\')

%% Get the list of files
dataDir = '.\Data\processed\' ;
files = struct2cell(dir([dataDir , '*.mat']))';

%% Define output folder 
trtDir = '.\Data\Spectrograms\';
if exist(trtDir, 'dir')~=7
    mkdir(trtDir)
end

%% Spectrogram parameters
Fs = 200;
params.movingwin = [4, 2];       
params.tapers    = [2, 3];
params.fpass     = [.5, 20];
params.Fs        = Fs;       

%% Main loop
for i = 1:size(files, 1)
    
    fileName = strrep(files{i, 1}, '.mat', '');
    specFile = [fileName, '_spect.mat'];
  
    %% Load data
    tmp = load([dataDir, fileName]);
    data = tmp.data;
    data(isnan(data)) = 0;

    %% convert to L-bipolar montage
    eeg_bi = fcn_Bipolar(data(1:19, :));

    %% Get regional average spectrograms
    [Sdata, stimes, sfreqs] = fcn_computeSpec(eeg_bi, params);
    stimes = round(stimes);

    %% Export
    save([trtDir, specFile], 'Sdata', 'stimes', 'sfreqs', 'params', '-v7.3');
   
end