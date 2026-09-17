function feat = sm_GetDataFeature_rat_normalized_sk(data,tim,ops)
% data is either the full path to a binary file or an Nxch int16 matrix
% where N = number of samples. If data is a matrix, tim is not used
%
% tim = time point to read from file
% ops = options describing how to calculate features (wavelet)
%
% ops.art_thres = threshold for considering data sample artifact
% ops.freqs = frequencies with which to calculate wavelet spectra
% ops.nCh_raw = number of channels in data matrix
% ops.nCh_featureFile = number of channels in binary file
% ops.Fs = sampling rate
% ops.ch_subj = which channels to load from binary file
% ops.durFeat = time window to calculate features

if isstr(data)
    
    
    fxml = strrep(data,'lfp','xml');
    xml = LoadXml(fxml);
    data = LoadBinary(data,'nchannels',xml.nChannels,'frequency',ops.Fs,'channels',ops.ch_subj,'duration', ops.durFeat,'start',tim);
end
stim_art = any(any(abs(data)>ops.art_thres));
nfreq = length(ops.freqs);
nCh_subj = length(ops.ch_subj);
%initialize
feat= nan(1,(nfreq*nCh_subj));
if ~stim_art
    
    
    
    %loop
    for j = 1:size(data,2)
        
        % wavelet decomposition
        tmp = abs(awt_freqlist(double(data(:,j)),ops.Fs,ops.freqs))';
        
        %loop over the frequencies
        for jj = 1:nfreq+2
            
            % for the sampling window, take the mean power at each frequency.
            if jj<=nfreq
                tmp = mean(tmp,2);
                feat((jj-1)*nCh_subj+j)  = tmp(jj);
            end
            if jj==(nfreq+1)
                feat((jj-1)*nCh_subj+j)  = skewness(double(data(:,j)));
            end
            if jj==(nfreq+2)
                feat((jj-1)*nCh_subj+j)  = kurtosis(double(data(:,j)));
            end
        end


    end
    feat=[feat,sum(feat(1:2:nfreq*2)),sum(feat(2:2:nfreq*2))];
    feat(1:2:nfreq*2)=feat(1:2:nfreq*2)/sum(feat(1:2:nfreq*2));
    feat(2:2:nfreq*2)=feat(2:2:nfreq*2)/sum(feat(2:2:nfreq*2));
else
    %error('here')
end
end