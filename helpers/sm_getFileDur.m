function fileDur = sm_getFileDur(fname)

%fname is the dat file with the raw neural data
%fileDur is time in seconds (assumed 16bit rez)
xml_fname = strrep(fname,'dat','xml');

xml = LoadXml(xml_fname);
ok = dir(fname);
siz = ok.bytes;


fileDur = siz/2/xml.SampleRate/xml.nChannels;



end

