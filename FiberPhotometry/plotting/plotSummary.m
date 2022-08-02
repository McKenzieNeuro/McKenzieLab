function plotSummary(fname_GEFI,evFile)

k  = gaussian2Dfilter([1 10000],1000);

load(fname_GEFI)

if ~isempty(evFile)
e = load(evFile);
end


[a,b] = fileparts(fname_GEFI);


[ix1, early, late, ts_PETH] = sm_getIndicesAroundEvent(GEFI.reward1,10,10,GEFI.fs,length(GEFI.corrected));
[ix2, early, late, ts_PETH] = sm_getIndicesAroundEvent(GEFI.reward2,10,10,GEFI.fs,length(GEFI.corrected));
ix = [ix1;ix2];

signal_smoothed = nanconvn(GEFI.corrected, k);


subplot(2,2,1)
text(0,1,['Subj: ' GEFI.subjName])
text(0,.9,['Date: ' GEFI.date])
text(0,.8,['Power: ' num2str(GEFI.power_signal) ' / ' num2str(GEFI.power_isos)])
text(0,.7,['wavelength: ' GEFI.signal ' / ' GEFI.isosbestic])


subplot(2,2,2)
hold on



plot(GEFI.ts_data,signal_smoothed)
plot([cell2mat(e.data(:,2)) cell2mat(e.data(:,2))]',[-2 2],'k')
text(cell2mat(e.data(:,2)),2.5*ones(size(e.data,1),1),e.data(:,1)) 
ylim([-2 3])
title('Iso corrected')



subplot(2,2,3)
plotMeanSEM(ts_PETH,signal_smoothed(ix),'k')
title('Mean response to reward')

subplot(2,2,4)
imagesc(ts_PETH,[],signal_smoothed(ix))
title('Reward PETH')


print(f1, '-dpsc2',filenameps ,'-append')




  %  ps2pdf('psfile',filenameps , 'pdffile',filenamepdf, 'gspapersize', 'a4','deletepsfile',1)
    
    