
function [outputs] = sm_plot_fiber_photometry_neural(baseDirectory,varargin)



% parse inputs


p = inputParser;


addParameter(p,'ephysDirectory',[],@ischar)
addParameter(p,'TDTDirectory',[],@ischar)

addParameter(p,'photoBleachCorrection','quadratic',@ischar)
addParameter(p,'savefig',false,@islogical)
addParameter(p,'skiptime',0,@isnumeric)
addParameter(p,'streams',{'x465A','x405A'},@iscell)
addParameter(p,'isosbestic','x405A',@ischar)
addParameter(p,'endTime',[],@isnumeric)
addParameter(p,'plotIntervals',[50 50],@isvector);
addParameter(p,'figureName','myPETH.fig',@ischar);
addParameter(p,'returnedDataType','corrected',@ischar);


parse(p,varargin{:})


photoBleachCorrection = p.Results.photoBleachCorrection;
savefig = p.Results.savefig;
skiptime = p.Results.skiptime;
streams = p.Results.streams;
isosbestic = p.Results.isosbestic;
endTime = p.Results.endTime;
plotIntervals = p.Results.plotIntervals;
figureName = p.Results.figureName;
returnedDataType = p.Results.returnedDataType;

% directory with both Intan and TDT
dirN_TDT = 'R:\McKenzieLab\DACSD\DACSD1\DACSD1_220527\Subject1-220527-160113';
cd(dirN_TDT)

% load TDT for pulses
tdt_data =  TDTbin2mat(pwd);
tdt_ups = tdt_data.epocs.Pu1_.onset;

clear tdt_data


% load TDT GEFI data, corrected
[signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(pwd,'streams',{'Dv4D','Dv5D'},'isosbestic','Dv4D');

%define smoothing kernel
k  = gaussian2Dfilter([10000 1],500);
signal_DFoF = nanconvn(signal_DFoF,k');

%load video scoring (should be one file)
v1 = load('objs.mat');
v2 = load('objs2.mat');
data = [v1.data;v2.data];

%%

dirN = 'R:\DACSD\DACSD1\DACSD1_220527\DACSD1_220527_160109';
[ups,dwns]  = sm_getDigitalin(dirN,'digitalin',1,30000);




aux  =LoadBinary('amplifier_analogin_auxiliary_int16_median.lfp','nchannels',35,'frequency',1250,'channels',35);
%%
freqs = logspace(log10(.5),log10(200),50);
d  =LoadBinary('amplifier_analogin_auxiliary_int16_median.lfp','nchannels',35,'frequency',1250,'channels',7);
ts_lfp = (1:length(d))/1250;
ts_fixed = interp1(ups,tdt_ups,ts_lfp,'linear','extrap');


%%

%
d_ts = ts_data;
obj_ts = cell2mat(data(:,2));
d_ts(d_ts<obj_ts(1)) = nan;
clear dat
for i = 1:length(obj_ts )-1
    kp = d_ts>obj_ts(i) & d_ts<obj_ts(i+1);
    d_ts(kp) = d_ts(kp) - obj_ts(i);
    dat(:,i) = avghist(d_ts(kp),signal_DFoF(kp),0:10:300);
    hold on
end
kp=d_ts>obj_ts(end);
    d_ts(kp) =  d_ts(kp) - obj_ts(end);
    
      dat = [dat avghist(d_ts(kp),signal_DFoF(kp),0:10:300)'];
      %%
      close all
      col  = linspecer(length(b));
      figure
      
      
     b= bar(dat');
      
     for i = 1:length(b)
         b(i).FaceColor = col{i};
     end
    %  legend({'Obj1','Obj2','Obj3','Obj4','Obj5','Obj6','Obj7','Obj8','Obj9'})
      set(gca,'xtick',1:9,'xticklabel',1:9)
      ylabel('DA response')
      xlabel('Object ID')
      colormap(cell2mat(col))
  
    
    h = colorbar;
set(get(h,'label'),'string','time from obj. (min)');
h.Ticks = 0:.2:1;
h.TickLabels= 0:5;

%%


%%

load('amplifier_analogin_auxiliary_int16_median.ripples.events.mat')

[~,b] = histc(ripples.peaks,ts_lfp);
ts_ripples = ts_fixed(b);
%%
d = awt_freqlist(double(d),1250,freqs);
d = zscore(abs(d));

%%
bints = -10:5:max(ts_fixed)+60;
 [~,b] = histc(ts_fixed,bints);
 clear d1
 for i = 1:size(d,2)
     
    d1(:,i) =  accumarray(b',d(:,i),[],@nanmean,nan);
     
 end
 aux1 =  accumarray(b',double(abs(aux)),[],@nanmean,nan);
 
 
  [~,b] = histc(ts_data,bints);
  
   DA =  accumarray(b',double(signal_DFoF),[],@nanmean,nan);

%%
bins = prctile(DA,0:5:100);
[~,b1] = histc(DA,bins);
clear d2
for i = 1:50
d2(:,i) = accumarray(b1(b1~=0),d1(b1~=0,i),[],@nanmean,nan);

end

aux2  = accumarray(b1(b1~=0),aux1(b1~=0),[],@nanmean,nan);
%%
figure
col  = linspecer(size(d2,1),'jet');
for i = 1:size(d2,1)-1
    semilogx(freqs,d2(i,:),'color',col{i})

hold on
end
semilogx(freqs,d2(end,:),'color','k','linewidth',3)

%legend({'Lowest','Low','Mid','High','Highest'})
title('St. PYR')
xlabel('Frequency')
ylabel('zscore power')
%%
close all
k  = gaussian2Dfilter([10000 1],50);
signal_DFoF = nanconvn(signal_DFoF,k');
[~,b] = histc(ts_ripples,ts_data);
ts_b  =-20*1250:1250*20;
ix = repmat(b',1,length(ts_b)) + repmat(ts_b,length(b),1);

ix = ix(all(ix>0,2) & ~any(ix>length(signal_DFoF),2),:);
plot(ts_b/1250,nanmean(signal_DFoF(ix)))



%%
