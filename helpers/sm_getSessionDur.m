function dur = sm_getSessionDur(fnameIn)

 fxml = [fnameIn(1:end-4) '.xml'];
    if ~exist(fnameIn,'file') && ~exist(fxml,'file')
        error('Dat file and/or Xml file does not exist')
    end
    
    syst = LoadXml(fxml);
    
    fInfo = dir(fnameIn);
    dur  = fInfo.bytes/2/syst.nChannels/syst.SampleRate;
  
end

