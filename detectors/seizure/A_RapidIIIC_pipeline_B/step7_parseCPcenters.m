%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: Parse data - CP centers only (14sec EEG, 10min spectrogram etc.)                                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear
addpath('.\Tools\')

%% Get the list of files
cpcDir  = '.\Data\CPDs\';
speDir  = '.\Data\spectrograms\';
scrDir  = '.\Data\iiic\model_prediction\';
dataDir = '.\Data\processed\' ;
files = struct2cell(dir([dataDir, '*.mat']))';

%% Define output folder
trtDir = '.\Data\CP_centers\';
if exist(trtDir, 'dir')~=7
    mkdir(trtDir)
end

%% Parameters
Fs = 200;
channels = {'fp1';'f3';'c3';'p3';'f7';'t3';'t5';'o1';'fz';'cz';'pz';'fp2';'f4';'c4';'p4';'f8';'t4';'t6';'o2';'EKG'};   
ww_eeg = 14;        % sec 
ww_spe = (10*60)/2; % 2sec sgements

%% Main loop
for k = 1:size(files, 1)
  
    fileKey = strrep(files{k, 1}, '.mat', '');
    
    %% load cp look-up-table
    tmp = load([cpcDir, fileKey, '_cpc.mat']);
    lut = tmp.lut_cpd;  % [1]idx_center [2]idx_start [3]idx_end
    
    %% load spectrogram
    tmp = load([speDir, fileKey, '_spect.mat']); 
    SDATA = tmp.Sdata(:, [2 1]); sfreqs = tmp.sfreqs;
    M = size(SDATA{1, 1}, 2);
    nFreqs = size(SDATA{1, 1}, 1);
    
    %% load model predicted probabilities
    tmp = load([scrDir, fileKey, '_score.mat']); 
    Y = tmp.Y_model;
    
    %% load EEG data
    matObj = matfile([dataDir, fileKey, '.mat']); matObj.Properties.Writable = true;
    size_data = size(matObj, 'data'); 
    N = size_data(2);
    nCh = size_data(1);
 
    %% loop for each CP center - 14sec EEG and 10min Spectrograms 
    for kk = 1:size(lut, 1)
        
        idx_CPcenter = lut(kk, 1);     % in 2sec
        idx_CPrange  = lut(kk, [2 3]);
        scores = Y(lut(kk, 1), :);
        
        % naming convention for each sample
        sampleKey = [fileKey, '_', num2str(idx_CPcenter)];       
        disp(['- parse ', sampleKey])
            
        % parse 14sec EEG
        t_center = idx_CPcenter*2-1;                % in sec 
        t_left   = (t_center*Fs)-(ww_eeg/2)*Fs+1;   % in data points
        t_right  = (t_center*Fs)+(ww_eeg/2)*Fs; 
        seg  = matObj.data(:, max(t_left, 1):min(t_right, N));     
        SEG = fcn_parseData([nCh, N], seg, Fs, t_left, t_right, ww_eeg);
        
        % parse 10min spectrograms 
        st_center = idx_CPcenter;       
        st_left  = (st_center)-(ww_spe/2)+1;  
        st_right = (st_center)+(ww_spe/2); 

        Sdata = SDATA;
        for ii = 1:size(Sdata, 1)
            ss = Sdata{ii, 1}(:, max(st_left, 1):min(st_right, M));
            Sdata{ii, 1} = fcn_parseData([nFreqs, M], ss, 1, st_left, st_right, ww_spe);
        end
        
        % export 
        save([trtDir, sampleKey, '.mat'], 'fileKey', 'SEG', 'Sdata', 'sfreqs', 'scores', 'idx_CPrange', 'idx_CPcenter')   
   
    end
end

