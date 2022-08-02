function [data,ix,signal_smoothed] = sm_QC_DANE(dirName,sensor)
k = gaussian2Dfilter([100000 1], 1000);
data = TDTbin2mat(dirName);

switch sensor
    case 'DA'
        stream1 = 'x450D';
        stream2 = 'x500D';
        [signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(dirName,'streams',{'x450D','x500D'},'isosbestic','x450D');
    case 'NE'
        stream1 = 'x405A';
        stream2 = 'x465A';
        [signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(dirName,'streams',{'x405A','x465A'},'isosbestic','x405A');
end

[ix1, early, late, ts_PETH] = sm_getIndicesAroundEvent(data.epocs.U11_.onset,10,10,fs,length(signal_DFoF));
[ix2, early, late, ts_PETH] = sm_getIndicesAroundEvent(data.epocs.U12_.onset,10,10,fs,length(signal_DFoF));
ix = [ix1;ix2];

signal_smoothed = nanconvn(signal_DFoF, k');
subplot(2,2,1)
plot(ts_data,data.streams.(stream1).data,'k')
hold on
plot(ts_data,data.streams.(stream2).data,'r')
title('Raw data')

subplot(2,2,2)
plot(ts_data,signal_smoothed)
title('Iso corrected')

subplot(2,2,3)
plotMeanSEM(ts_PETH,signal_smoothed(ix),'k')
title('Mean response to reward')

subplot(2,2,4)
imagesc(ts_PETH,[],signal_smoothed(ix))
title('Reward PETH')

end