%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Goal: Get PaCMAP                                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear

embedDir_trt = './Data/PaCMAP/';
if exist(embedDir_trt, 'dir')~=7
    mkdir(embedDir_trt)
end
pacmapFile_in  = 'pacmap_input.mat';
pacmapFile_out = 'pacmap_output.mat';

%% Load look-up-table
tmp = load('./Data/GUI_LUT.mat');
LUT = tmp.LUT_all;
Y_model = cell2mat(LUT(:, end));
X = [Y_model, sum(Y_model(:, 2:5), 2)];
save([embedDir_trt, pacmapFile_in], 'X', 'Y_model')  
 
%% Compute PaCMAP
str = strrep(strrep([embedDir_trt,  pacmapFile_in, ' ', embedDir_trt, pacmapFile_out], '.mat', ''),'\','/');   % no whitespace plz!!%
commandStr = ['python ./Tools/computePaCMAP.py ', str];
disp(commandStr)
system(commandStr);

%% Export 
tmp = load([embedDir_trt, pacmapFile_out]);
Vxy = double(tmp.X_pacmap);
save([embedDir_trt, pacmapFile_out], 'Vxy', 'Y_model')  


