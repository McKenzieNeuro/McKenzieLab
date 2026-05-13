function fileDur = sm_getFileDur(fname,varargin)

p = inputParser;


addParameter(p,'nChannels',[],@isnumeric) %method to correct for decaying signal
addParameter(p,'fs',[],@isnumeric) % how much time (s) to skip before calculating mean/std/decay



parse(p,varargin{:})


nChannels = p.Results.nChannels;
fs = p.Results.fs;

%fname is the dat file with the raw neural data
%fileDur is time in seconds (assumed 16bit rez)
[~,~,c] = fileparts(fname);

if strcmp(c,'.dat')
    xml_fname = strrep(fname,'dat','xml');
    
elseif  strcmp(c,'.lfp')
    xml_fname = strrep(fname,'lfp','xml');
    
else
    error('needs dat or lfp')
    
end


if exist(xml_fname)
    xml = LoadXml(xml_fname);
        
    
    if strcmp(c,'.dat')
    fs = xml.SampleRate;
    
    elseif  strcmp(c,'.lfp')
        fs = xml.lfpSampleRate;
    end
    
    nChannels = xml.nChannels;
    
end


ok = dir(fname);
siz = ok.bytes;


fileDur = siz/2/fs/nChannels;



end

