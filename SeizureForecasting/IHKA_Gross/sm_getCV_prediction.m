%load all features from all sessions
v=load('G:\data\IHKA_gross\features.mat','sessions');

%define group ID for each observation
%group_all = cell2mat(cellfun(@(a,b) b*ones(size(a,1),1),v.dat,num2cell(1:length(v.dat)),'uni',0)');
%sesID_all = cell2mat(v.sesID');
%dat_all = cell2mat(v.dat');
[~,mouse] = fileparts(fileparts(v.sessions'));
[umouse,~,c] = unique(mouse);
%%
cd('G:\data\IHKA_gross\classification')

%fils = getAllExtFiles(pwd,'mat',0);
C_norm = [];ix=1;
%%
C_norm1 = [];ix=1;
for i = 1:110
    clear C_norm
    fname = ['classification_' num2str(i) '.mat'];
    if exist(fname)
        try
            load(fname,'C_norm','C')
            if exist('C_norm')
                C_norm1(:,:,i) =  C_norm;
            else
                C_norm1(:,:,i) = nan(6);
            end
        end
    else
        C_norm1(:,:,i) = nan(6);
    end
    i
end

%%
lab = [{'-inf to -3600s'},{'-3600s to -100s'},{'-100 to -10'},{'-10 to seiz.'},{'seiz'},{'post'}];
close all
figure
for i = 1:max(c)
    subplot(3,3,i)
imagesc(1:6,1:6,100*(nanmedian(C_norm1(:,:,i==c),3)-(1/6))./(1/6))
set(gca,'xtick',1:6,'xticklabel',lab,'ytick',1:6,'yticklabel',lab,'ydir','normal')
xlabel('Predicted class')
ylabel('Real class')
colormap('bluewhitered')
colorbar
title(umouse(i))
end
%%

close all
figure

imagesc(1:6,1:6,100*(nanmean(C_norm1,3)-(1/6))./(1/6))
set(gca,'xtick',1:6,'xticklabel',lab,'ytick',1:6,'yticklabel',lab,'ydir','normal')
xlabel('Predicted class')
ylabel('Real class')
colormap('bluewhitered')
colorbar
title('All mice')

%%
figure
ok=[];
for i = 1:7
tmp  =100*(C_norm1(:,:,i)-(1/6))./(1/6);
ok(i,:)  =tmp(eye(6)==1);
end
plotMeanSEM(1:6,ok,'k')
set(gca,'xtick',1:6,'xticklabel',lab)
ylabel('% improvement over chance')