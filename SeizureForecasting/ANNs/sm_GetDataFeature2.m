function feat = sm_GetDataFeature2(data,ops)
% this function takes the wavelet spectra of data according to the options in ops
% output will be a vector that is the number of channels x number of 

nsample = size(data,1); % 
nfreq = length(ops.numfreq);

%initialize 
feat= nan(1,nfreq*ops.nChanSubj);

%loop
for j = 1:ops.nChanSubj
    
    % wavelet decomposition
    tmp = abs(awt_freqlist(double(data(:,j)),ops.Fs,ops.numfreq))';

    %loop over the frequencies
    for jj = 1:nfreq
        
        % for the sampling window, take the mean power at each frequency.
        tmp = mean(tmp,2);
        feat((jj-1)*ops.nChanSubj+j)  = tmp(jj);
    end


end
end