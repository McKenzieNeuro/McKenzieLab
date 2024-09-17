function data=fcn_preprocess(file_path)
    tmp=load(file_path);
    record=tmp.data;Fs=tmp.Fs;
    channels_org=tmp.channels;
    
    % select channels
    channels={'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};  
    idx_ch=NaN(size(channels,1),1);
    for j=1:size(channels,1)
        jj=find(ismember(lower(channels_org),lower(channels{j})));
        if isempty(jj)
            disp(['warning: missing channel: ',channels{j},'!']);
            keyboard;
        else
            idx_ch(j)=jj;
        end
    end
    data=record(idx_ch,:); 
    
    % resample to 200Hz
    if round(Fs)~=200
        data=resample(data',200,Fs)';
        Fs=200;
    end
    
    % denoise
    [B1,A1]=butter(3,[.5,70]/(Fs/2));
    [B2,A2]=butter(3,[60-2.5,60+2.5]/(Fs/2),'stop');
    data(isnan(data))=0;data(isinf(data))=0;
    data=filtfilt(B1,A1,data')';
    data=filtfilt(B2,A2,data')';
end