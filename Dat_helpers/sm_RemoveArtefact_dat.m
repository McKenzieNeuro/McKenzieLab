function sm_RemoveArtefact_dat(fnameIn,artefact)

%
% INPUT
%     - fnameIn     : source file
%     - artefact   : Nx2 matrix of start end times to interpolate over


fxml = [fnameIn(1:end-4) '.xml'];
if ~exist(fnameIn,'file') && ~exist(fxml,'file')
    error('Dat file and/or Xml file does not exist')
end
sizeInBytes = 2; % change it one day...

syst = LoadXml(fxml);


nbChan = syst.nChannels;


fidI = fopen(fnameIn,'r+');
fidO = fidI;

for i = 1:size(artefact,1)
    duration = artefact(i,2)-artefact(i,1);
    nBytes = round(syst.SampleRate*duration);
    
    start = floor(artefact(i,1)*syst.SampleRate)*syst.nChannels*sizeInBytes;
    status = fseek(fidI,start,'bof');
    if status ~= 0,
        error('Could not start reading (possible reasons include trying to read a closed file or past the end of the file).');
    end
    
    
    
    dat = fread(fidI,nbChan*nBytes,'int16');
    dat = reshape(dat,[nbChan nBytes]);
    n = size(dat,2);
    
    for ii = 1:size(dat,1)
        dat(ii,:) =  interp1([1 n],[dat(ii,1) dat(ii,end)],1:n);
    end
   
    dat = int16(dat(:));
 
    status = fseek(fidO,start,'bof');
    fwrite(fidO,dat,'int16');
end


fclose(fidI);





