function feat = sm_GetDataFeature_rat(data,tim,ops)
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
feat= nan(1,nfreq*nCh_subj);
if ~stim_art
    
    
    
    %loop
    for j = 1:size(data,2)
        
        % wavelet decomposition
        tmp = abs(awt_freqlist(double(data(:,j)),ops.Fs,ops.freqs))';
        
        %loop over the frequencies
        for jj = 1:nfreq
            
            % for the sampling window, take the mean power at each frequency.
            tmp = mean(tmp,2);
            feat((jj-1)*nCh_subj+j)  = tmp(jj);
        end
        
        
    end
else
    %error('here')
end
end