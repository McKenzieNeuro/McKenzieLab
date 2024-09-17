%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: EEG segmentation via CP detection                                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear
addpath('.\Tools\')

%% Get the list of files
specDir = '.\Data\spectrograms\' ;
dataDir = '.\Data\processed\' ;
files = struct2cell(dir([dataDir , '*.mat']))';

%% Define output folder
trtDir = '.\Data\CPDs\';
if exist(trtDir, 'dir')~=7
    mkdir(trtDir)
end

%% CPD parameters
Fs = 200;
alpha_cpd = .1; 

%% Main loop
for i = 1:size(files, 1)
    
    fileName = files{i, 1};
    
    %% find total number of 2-sec intervals in each file 
    matObj = matfile([dataDir, fileName]);
    [~, N] = size(matObj,'data'); 
    nn = ceil(N/(2*Fs));
        
    %% find the spectrogram length
    tmp = load([specDir, strrep(fileName, '.mat', '_spect.mat')]);
    Sdata = tmp.Sdata(:, [2 1]);
    mm = size(Sdata{1}, 2);
    
    %% adjust the length to align in time
    if mm>=nn
        for kk = 1:size(Sdata, 1)
            Sdata{kk, 1} = Sdata{kk, 1}(:, 1:nn);
        end
    else
        dd1 = size(Sdata{1, 1}, 1); dd2 = nn-mm;
        for kk = 1:4
            Sdata{kk, 1} = [Sdata{kk, 1}, eps+zeros(dd1, dd2)];
        end
    end
    
    %% CP deteciton 
    [isCPs, isCPcenters] = fcn_cpd(Sdata, alpha_cpd);

    % CP range + center infor %
    idx_rise = find(isCPs==1);
    idx_fall = unique([idx_rise(2:end)-1; length(isCPs)]);
    idx_cpc  = find(isCPcenters==1);  
    
    %% Export look-up-table
    lut_cpd = [idx_cpc idx_rise idx_fall];
    
    disp(['-add ', fileName, ': found ', num2str(size(lut_cpd, 1)), ' cp segments'])
    save([trtDir, strrep(fileName, '.mat', '_cpc.mat')], 'isCPs', 'isCPcenters', 'lut_cpd');
end
