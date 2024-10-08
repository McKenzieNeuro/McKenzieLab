function sm_dat_from_nk(fname,fnameOut,scaleFactor)
%INPUT
%fname = full path of *.EEG file
chunk = 1e6;
[sFile, ChannelMat] = in_fopen_nk(fname);

nSample = sFile.epochs.times(2)*sFile.header.sample_rate;
nLoop = floor(nSample/chunk);
fidO = fopen(fnameOut,'w');


idx = [-1 -1];
for ch = 1:nLoop
    
    idx = [0 chunk-1]+((ch-1)*chunk);
    [F, TimeVector] = in_fread(sFile, ChannelMat,1,idx);
    EEG_kp = contains({ChannelMat.Channel.Type},'EEG');
    F(EEG_kp,:) = F(EEG_kp,:) - repmat(median(F(EEG_kp,:)),sum(EEG_kp),1);
    F = F*scaleFactor;
    fwrite(fidO,F(:),'int16');
end

idx = [idx(2)+1 nSample];
[F, TimeVector] = in_fread(sFile, ChannelMat,1,idx);
EEG_kp = contains({ChannelMat.Channel.Type},'EEG');
F(EEG_kp,:) = F(EEG_kp,:) - repmat(median(F(EEG_kp,:)),sum(EEG_kp),1);
F = F*scaleFactor;

fwrite(fidO,F(:),'int16');

% deal with remainder


fclose(fidO);

