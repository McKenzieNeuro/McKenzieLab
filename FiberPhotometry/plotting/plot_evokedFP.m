FP_dir = 'R:\DANEHippocampalResponse\lc_cHR4\LC_Chr4_240417_085356\DACSD1-240417-085738';
Intan_dir = 'R:\DANEHippocampalResponse\lc_cHR4\LC_Chr4_240417_085356';


[ups,dwns]  = sm_getDigitalin(Intan_dir,'digitalin.dat',30000,16);
[signal_DFoF1,ts_data1,fs] = sm_getSignal_DFoF(FP_dir,'streams',{'x500D','x450D'},'isosbestic','x450D','baseline',[5 500]);
data = TDTbin2mat(FP_dir);
k = gaussian2Dfilter([ 1 fs*100],fs);
signal_DFoF2 = nanconvn(signal_DFoF1,k);
%%


fs = data.streams.x500D.fs;


stim_ts_FP = data.epocs.PC0_.onset;
%stim_ts_Intan = ups{1};

kp = diff([0;stim_ts_FP])>1 ;

[ix,early,late,ts_ev] = sm_getIndicesAroundEvent(stim_ts_FP(kp),30,60,fs,length(signal_DFoF1));

ix = ix(~early &~late,:);
figure

plotMeanSEM(ts_ev,zscore(signal_DFoF2(ix),[],2),'k')
%%
cd(Intan_dir)
kp = diff([0;ups{1}])>1 ;
ts = ups{1}(kp);
clear d dd
k = gaussian2Dfilter([ 10000 1],1000);
for j = 1:8
    for i = 1:sum(kp)
        tmp = LoadBinary('amplifier_analogin_auxiliary_int16_median.dat','nchannels',11,'channels',j,'frequency',30000,'start',ts(i)-2,'duration',5);
        %d(i,:) = nanconvn(InstAmplitude(BandpassFilter(double(tmp),30000,[500 5000])),k');
        d(i,:) = double(tmp);
    end
    dd(j,:) = nanmean(d);
    j
end

%%
cd('R:\DANEHippocampalResponse\LC_Chr2\LC_Chr2_240425_090551')

xml =  LoadXml('amplifier_analogin_auxiliary_int16.xml');

sh1 = 2;
sh = cell2mat({xml.AnatGrps(sh1).Channels})+1;
ev = LoadEvents('amplifier_analogin_auxiliary_int16.evt.ctx');

up_bad = (abs(tt'-bestmatch(tt,ups{1}))<.002);
dwn_bad = (abs(tt'-bestmatch(tt,dwns{1}))<.002);
%%

clear spect d csd volt
ix = 1;
for i = 1:length(ts)
  
   
        tmp = LoadBinary('amplifier_analogin_auxiliary_int16.dat','nchannels',22,'channels',5,'frequency',30000,'start',ts(i)-.1,'duration',.2);
        tmp = double(tmp);
       % tmp(up_bad|dwn_bad) = nan;
       % tt = (1:length(tmp))/30000 +ts(i)-1;
        d(ix,:) = tmp;
        
        %d1 = awt_freqlist(tmp,1250,logspace(log10(1),log10(200),100));
        
        %spect(:,:,i) = abs(d1);
   
   
      ix = ix+1;
   % [ csd(:,:,i), CSDelecinds] = csd5pt(d,.025,0);
   %volt(:,:,i) = d;
    i
end

%%

plot(((1:size(csd,2))/30000)-1,nanmedian(csd(3,:,ts<ev.time(1)),3))
hold on
plot(((1:size(csd,2))/30000)-1,nanmedian(csd(6,:,ts<ev.time(1)),3))
plot([0 2],[100 100],'k','linewidth',6)
%%


%eStim
cd('R:\DANEHippocampalResponse\lc_cHR4\LC_Chr4_240426_153810')
d = LoadBinary('stim.dat','nchannels',16,'channels',2,'frequency',20000);
stim = find(diff([0;d]>.5e4)>0)/20000;
%%
k = gaussian2Dfilter([ 1 10000],200);
for i = 1:length(stim)
    
   LC(i,:) = LoadBinary('amplifier.dat','nchannels',16,'channels',13,'frequency',20000,'start',stim(i)-2,'duration',4); 
   fil(i,:) = nanconvn(InstAmplitude(BandpassFilter(double(LC(i,:)),20000,[100 5000])),k);
   %fil(i,:) = BandpassFilter(double(LC(i,:)),20000,[100 5000]);
   %fil_notch(i,:) = sm_notch_filter(double(LC(i,:)),20000,60);
  
   i
end


