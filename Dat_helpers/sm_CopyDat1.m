function sm_CopyDat(fnameIn,fnameOut,varargin)

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



p = inputParser;


addParameter(p,'duration',[],@ischar)
addParameter(p,'start',0,@isnumeric)
addParameter(p,'channelIx',[],@isvector);
addParameter(p,'nbChan',1,@isnumeric);
addParameter(p,'SampleRate',[],@isnumeric);



parse(p,varargin{:})


duration = p.Results.duration;
start = p.Results.start;
channelIx = p.Results.channelIx;
nbChan = p.Results.nbChan;
SampleRate = p.Results.SampleRate;

sizeInBytes = 2; % int16 file
chunk = 1e5; % depends on the RAM



if ~isempty(fnameIn)
    
    
    fxml = [fnameIn(1:end-4) '.xml'];
    if ~exist(fnameIn,'file') && ~exist(fxml,'file')
        error('Dat file and/or Xml file does not exist')
    end
    
    syst = LoadXml(fxml);
    
    fInfo = dir(fnameIn);
    
    
    if isempty(nbChan)
    nbChan = syst.nChannels;
    end
    
    if isempty(SampleRate)
        SampleRate = syst.SampleRate;
    end
    
    if ~isempty(duration)
        nBytes = SampleRate*duration*sizeInBytes*nbChan;
        if fInfo.bytes<nBytes
            error('Duration longer than total duration of file')
        end
    else
        nBytes = fInfo.bytes;
    end
    
    
else
    
     if ~isempty(nbChan) && ~isempty(SampleRate)  && ~isempty(duration)  
         nBytes = SampleRate*duration*sizeInBytes*nbChan;
     end
     
end


nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunk));

fidO = fopen(fnameOut,'w');

if ~isempy(fnameIn)
    fidI = fopen(fnameIn,'r');
else
    fidI = nan;
end
start = floor(start*SampleRate)*nbChan*sizeInBytes;
status = fseek(fidI,start,'bof');
if status ~= 0,
    error('Could not start reading (possible reasons include trying to read a closed file or past the end of the file).');
end

for ii=1:nbChunks
    h=waitbar(ii/(nbChunks+1));
  if ~isnan(fidI)
    dat = fread(fidI,nbChan*chunk,'int16');
  else
      
      dat = zeros(nbChan*chunk,1);
  end
  
  
    if ~isempty(channelIx)
        dat = reshape(dat,[nbChan chunk]);
        dat = dat(channelIx,:);
        dat = dat(:);
    end
    fwrite(fidO,dat,'int16');
end

remainder = nBytes/(sizeInBytes*nbChan) - nbChunks*chunk;
if ~isempty(remainder)
    
     if ~isnan(fidI)
    dat = fread(fidI,nbChan*remainder,'int16');
     else
          dat = zeros(nbChan*remainder,1);
     end
    if ~isempty(channelIx)
        dat = reshape(dat,[nbChan remainder]);
        dat = dat(channelIx,:);
        dat = dat(:);
    end
    fwrite(fidO,dat,'int16');
end

close(h);

fclose(fidI);
fclose(fidO);




