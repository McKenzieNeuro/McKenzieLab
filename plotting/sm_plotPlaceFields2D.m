%get positio
%clear all

warning off
spikes = bz_GetSpikes;


basepath  =pwd;
basename = bz_BasenameFromBasepath(basepath);



load('position_info.mat')



%%
kp = cellfun(@range,pos_inf.len_ep)>50;

pos_inf.in_eps = pos_inf.in_eps(kp);

pos_inf.out_eps = pos_inf.out_eps(kp);

pos_inf.len_ep = pos_inf.len_ep(kp);

pos_inf.ts_ep = pos_inf.ts_ep(kp);


ineps = cell2mat(cellfun(@(a) [a(1) a(end) ] ,pos_inf.ts_ep(pos_inf.in_eps),'uni',0));
outeps = cell2mat(cellfun(@(a) [a(1) a(end) ] ,pos_inf.ts_ep(pos_inf.out_eps),'uni',0));






t= pos_inf.ts;
%kp = InIntervals(sti{1}(:,1),[min(t) max(t)]);
len_in_st =[];
len_in_fin =[];
len_out_st =[];
len_out_fin =[];

interval_out =[];
len = pos_inf.lin_pos;



t= pos_inf.ts;

kp = cellfun(@range,pos_inf.ts_ep)<30;

pos_inf.len_ep = pos_inf.len_ep(kp);
pos_inf.ts_ep = pos_inf.ts_ep(kp);
pos_inf.in_eps = pos_inf.in_eps(kp);
pos_inf.out_eps = pos_inf.out_eps(kp);
ineps = cell2mat(cellfun(@(a) [a(1) a(end) ] ,pos_inf.ts_ep(pos_inf.in_eps),'uni',0));
outeps = cell2mat(cellfun(@(a) [a(1) a(end) ] ,pos_inf.ts_ep(pos_inf.out_eps),'uni',0));





%%

k  = gaussian2Dfilter([1000 1000],[1 1]);
dt = mode(diff(pos_inf.ts));

%plot data for each neu1ron
binsX = min(pos_inf.x):range(pos_inf.x)/100:max(pos_inf.x);
binsY = min(pos_inf.y):range(pos_inf.y)/100:max(pos_inf.y);
% get occupancy

[n,~,~,b] = histcn([pos_inf.x pos_inf.y],binsX,binsY);
n(n<10) = nan;
occ = nanconvn(n,k,'nanout',true)*dt;
 [t,x,y,vx,vy,ax,ay] = KalmanVel(pos_inf.x,pos_inf.y,pos_inf.ts,2);
 
 kp = sqrt(vy.^2+vx.^2)>50;
 
 [n,~,~,b] = histcn([pos_inf.x pos_inf.y],binsX,binsY);
n(n==0) = nan;
b(~kp,:) = -1;

%%
for i = 	1:length(spikes.times)
    h = figure('visible','on');
    
    %get binned spike times
    n_spk  =  histc(spikes.times{i},pos_inf.ts);
    binned_spikeCount = accumarray(b(kp,:),n_spk(kp),[101 101],@sum);
    binned_spikeCount(isnan(occ)) = nan;
    binned_spikeCount = nanconvn(binned_spikeCount,k,'nanout',true);
    FR = binned_spikeCount./occ;
    imagesc(FR)%,[0 20])
    
    set(h,'Units','Inches');
    
    pos = get(h,'Position');
    
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    set(groot,'DefaultFigureRenderer','painters')
    colorbar
    % print(h, '-dpsc2', [basename '.ps'], '-append');
    
    uiwait
  %  close all
end

