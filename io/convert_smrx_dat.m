function convert_smrx_dat(fname,fnameOut)
% inputs
% fname = matlab files with structs for each channel
% fnameout = *.dat file with int16 for each channel merged




v = load(fname);

field_names = fields(v);
kp = contains(field_names,'Ch');



field_names = field_names(kp);

nCh = sum(kp);
ch = nan(nCh,1);


% put all on same clock







for i = 1:length(field_names)
    
    fs(i) = 1./v.(field_names{i}).interval;
    
end

[fs_M,b] = max(fs);
nSamples = length(v.(field_names{b}).values);
ts = (1:nSamples)/fs_M;

data(nCh, nSamples) = int16(0);
for i = 1:length(field_names)
    
    
    ix = regexp(field_names{i},'_Ch');
    
    ch(i) = str2num(field_names{i}(ix+3:end));
    
    if fs(i)~=fs_M
        dat = v.(field_names{i}).values;
        ts1 = (1:length(dat))/fs(i);
        dat  = int16(interp1(ts1,double(dat),ts));
         data(i,:)  =dat';
        
    else
    data(i,:) = v.(field_names{i}).values';
    end
end

[ch,b] = sort(ch);

data = data(b,:);


sm_Mat2Dat(data,fnameOut)

