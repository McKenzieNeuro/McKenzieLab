function feat = sm_GetDataFeature1(data,ops)



nsample = size(data,1);
nfreq = length(ops.numfreq);
feat= nan(nfreq*ops.nChanSubj,ops.nTempBin);
for j = 1:ops.nChanSubj
    tmp = abs(awt_freqlist(double(data(:,j)),ops.Fs,ops.numfreq))';


    for jj = 1:nfreq
        t = avghist(1:nsample,tmp(jj,:),linspace(1,nsample,ops.nTempBin+1));
        feat((jj-1)*ops.nChanSubj+j,:)  = t(1:end-1);
    end


end
end