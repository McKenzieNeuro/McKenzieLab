%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: Get the BoW model for template matching in GUI                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear

%% Load LUT
dataDir = '.\Data\CP_centers\';
tmp = load('.\Data\GUI_LUT.mat');
LUT = tmp.LUT_all;

%% Get BoW input 
SS =  cell(size(LUT, 1), 1);
event_idx = [];
for i = 1:size(LUT, 1)
    
    disp(['-add ', num2str(i)])
    file = [LUT{i, 1}, '_', num2str(LUT{i, 2})];
    dd = LUT{i, 3}(2)-LUT{i, 3}(1)+1;
    
    % load sprectrogram 
    tmp = load([dataDir, file]);
    sdata = cell2mat(tmp.Sdata(:, 1));
    idx1 = max(1, size(sdata, 2)/2-round(dd/2)+1);
    idx2 = min(size(sdata, 2), idx1+dd-1);
    SS{i} = sdata(:, idx1:idx2)';
    
    dd = idx2-idx1+1;
    event_idx = [event_idx; repmat(i, dd, 1)];
end
 
%% BoW model - learn 500 words via K-means
K_bow = 500;
X = cell2mat(SS);
s = pow2db(X+eps);
s(s<-10) = -10; s(s>25) = 25;

rng('default')
bow = kmeans(s, K_bow);

%% BoW model - get word distribution in each segment
bow_vec = NaN(size(LUT, 1), K_bow);
for i = 1:size(LUT, 1)
    disp(['-add ', num2str(i)])

    idx = find(event_idx==i);
    y = bow(idx);
    
    bow_vec(i, :) = hist(y, 1:K_bow);
end
bow_vec = bow_vec./repmat(sum(bow_vec, 2), 1, size(bow_vec, 2));
 
%% Export
save('.\Data\BoW.mat', 'bow_vec')
