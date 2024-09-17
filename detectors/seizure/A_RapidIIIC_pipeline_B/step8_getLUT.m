%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: Get the look-up-table for labelling GUI                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear

%% Get the list of files
dataDir = '.\Data\CP_centers\';
files = struct2cell(dir([dataDir, '*.mat']))';
files = files(:, 1);

%% Initialization 
pp = {'Seizure', 'LPD', 'GPD', 'LRDA', 'GRDA', 'Other'};
LUT_all = cell(size(files, 1), 5);

%% Main loop 
for i = 1:size(LUT_all, 1)
    
    file = files{i, 1};
    disp(['add ', num2str(i), ' ', file])
 
    tmp = load([dataDir, file]);
    scr = tmp.scores;
    cpr = tmp.idx_CPrange;
    fileKey = tmp.fileKey;
    segindx = tmp.idx_CPcenter;
    
    [~, idx] = max(scr([2:6, 1]));
    prediction = pp{idx};
    
    LUT_all(i, :) = {fileKey, segindx, cpr, prediction, scr};
end
save('.\Data\GUI_LUT.mat', 'LUT_all')