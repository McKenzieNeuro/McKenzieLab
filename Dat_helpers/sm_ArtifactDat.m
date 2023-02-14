function sm_ArtifactDat(fnameIn,artifact,ch,start,xml)

% USAGE
%     CopyDat(fnameIn,fnameOut,'optionNames',optionalues)
% copy a given number of seconds of data from one dat file into another one
%
% INPUT
%     - fnameIn     : source file
%     - fnameOIut   : target file
% options:
%     - 'duration'  : duration in seconds
%     - 'start'     : start time in seconds
%     - 'channelIx' : subset of channels to be copied

% Adrien Peyrache, 2012

channelIx = [];



fxml = [fnameIn(1:end-4) '.xml'];
if ~exist(fnameIn,'file') && ~exist(fxml,'file')
    error('Dat file and/or Xml file does not exist')
end
sizeInBytes = 2; % change it one day...
chunk = 1e5; % depends on the system... could be bigger I guess

syst = LoadXml(fxml);
duration = length(artifact)/syst.SampleRate;
fInfo = dir(fnameIn);

nbChan = syst.nChannels;
if ~isempty(duration)
    nBytes = syst.SampleRate*duration*sizeInBytes*nbChan;
    if fInfo.bytes<nBytes
        error('Duration longer than total duration of file')
    end
else
    nBytes = fInfo.bytes;
end
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunk));

%fidO = fopen(fnameOut,'r+');
fidI = fopen(fnameIn,'r+');

fidO = fidI;

nBytes = round(syst.SampleRate*duration);

startbit = floor(start*syst.SampleRate)*syst.nChannels*sizeInBytes;
status = fseek(fidI,startbit,'bof');
if status ~= 0,
    error('Could not start reading (possible reasons include trying to read a closed file or past the end of the file).');
end



dat = fread(fidI,nbChan*nBytes,'int16');
dat = reshape(dat,[nbChan nBytes]);



dat1 = dat;
correct = dat(ch,:);
correct = double(correct) - repmat(artifact,length(ch),1);
dat1(ch,:)  = correct;
dat1 = int16(dat1);

status = fseek(fidO,startbit,'bof');
fwrite(fidO,dat1(:),'int16');


fclose(fidI);




