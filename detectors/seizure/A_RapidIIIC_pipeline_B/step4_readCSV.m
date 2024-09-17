%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: Read CSV to MAT and pad for two ends                                                       
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear

%% Get list of files
dataDir = '.\Data\processed\';
files = struct2cell(dir([dataDir, '*.mat']))';

modelDir = '.\Data\iiic\';
trtDir = '.\Data\iiic\model_prediction\';
if exist(trtDir, 'dir')~=7
    mkdir(trtDir)
end

Fs = 200;

%% main loop 
for i = 1:size(files, 1)
    file = strrep(files{i, 1}, '.mat', '');

    %% Get EEG duration 
    matObj = matfile([dataDir, file, '.mat']);
    [~, M] = size(matObj, 'data');
    mm = ceil(M/(2*Fs));

    %% Read CSV
    Y = cell2mat(table2cell(readtable([modelDir, file, '_score.csv'])));
    
    % pad for the 1st 4sec in time 
    Y = [repmat(Y(1, :), 2, 1) ; Y]; 
    
    % pad the rest in the end
    nn = size(Y, 1); 
    if nn>mm
        keyboard
    else
        dd = mm-nn;
        Y = [Y; repmat(Y(end, :), dd, 1)];
    end

    %% Export 
    Y_model = Y;
    save([trtDir, file, '_score.mat'], 'Y_model');
 
end
