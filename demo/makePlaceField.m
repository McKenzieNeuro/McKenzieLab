% load spikes
cd('R:\WSun\BiconBehav\RAM11\T1')

dirN = 'R:\WSun\BiconBehav\RAM11\T1';
[spikes] = sm_LoadPhy('basepath',dirN,'forceReload',true,'getWaveforms',false,'basename','amplifier_analogin_auxiliary_int16');

%%

% load position
csvF = 'T1_AMDLC_resnet50_2nd trialOct4shuffle1_100000.csv';
DLC_pos = readtable(csvF);
LE = cell2mat([table2cell(DLC_pos(:,2)) table2cell(DLC_pos(:,3))]);

conf_LE = cell2mat([table2cell(DLC_pos(:,4))]);
LE(conf_LE<.5,:) = nan;

RE = cell2mat([table2cell(DLC_pos(:,5)) table2cell(DLC_pos(:,6))]);

conf_RE = cell2mat([table2cell(DLC_pos(:,7))]);
RE(conf_RE<.5,:) = nan;


hat = cell2mat([table2cell(DLC_pos(:,8)) table2cell(DLC_pos(:,9))]);
conf_hat = cell2mat([table2cell(DLC_pos(:,10))]);
hat(conf_hat<.5,:) = nan;


X  = nanmean([LE(:,1) RE(:,1) hat(:,1)],2);
Y  = nanmean([LE(:,2) RE(:,2) hat(:,2)],2);

%%
%get sync info
% ts_syn = [Intan_time video_frame estimate_video_frame pulse#]

load('R:\WSun\BiconBehav\RAM11\T1\OptoRAM11_240813_094944\ts_syn.mat')

% grab up until the last pulse

good_ix = 1:max(ts_syn(:,2));
X = X(good_ix);
Y = Y(good_ix);
pos = [X Y];

ts_vid = interp1(ts_syn(:,2),ts_syn(:,1),good_ix);
kp = find(~isnan(ts_vid),1,'first'):length(ts_vid);
ts_vid = ts_vid(kp);
pos = pos(kp,:);

%%
k = gaussian2Dfilter([1000 1000],1);
dt = mode(diff(ts_vid));


bin_X = min(X):range(X)/100:max(X);
bin_Y = min(Y):range(Y)/100:max(Y);
% get spikes in each time bin


%get spatial bin at each time
[occ,~,~,b] = histcn(pos,bin_X,bin_Y);

occ = nanconvn(occ,k);
occ = occ*dt; % convert to seconds


%get mean spike count
kp = all(b>0,2);
clear ratemap
figure
for i = 1:spikes.numcells
    n = histc(spikes.times{i},ts_vid);
    tmp = accumarray(b(kp,:),n(kp),[101 101],@sum,nan);
    tmp = nanconvn(tmp,k);
    tmp = tmp./occ;
    tmp(occ<.1) = nan;
    ratemap(:,:,i) = tmp;
    
    imagesc(ratemap(:,:,i))
    waitforbuttonpress
    close all
end