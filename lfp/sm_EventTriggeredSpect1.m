function sm_EventTriggeredSpect(lfpfile,events,channel,varargin)




p = inputParser;


addParameter(p,'freqs',[],@isvector)
addParameter(p,'savefig',false,@islogical)
addParameter(p,'plotIntervals',[50 50],@isvector);
addParameter(p,'figureName','myPETH.fig',@ischar);
addParameter(p,'edge',[],@isnumeric);
addParameter(p,'fs',1250,@isnumeric);
addParameter(p,'channel',1,@isnumeric);

parse(p,varargin{:})


freqs = p.Results.freqs;
savefig = p.Results.savefig;
plotIntervals = p.Results.plotIntervals;
figureName = p.Results.figureName;
edge = p.Results.edge;
channel = p.Results.channel;





%%
if ~exist(lfpfile)
    
    error('LFP file does not exist')
    
    
end

xmlF = [lfpfile(1:end-3) 'xml'];

if ~exist(xmlF)
    
    error('XML file does not exist')
    
    
end
xml = LoadXml(xmlF);

if isempty(fs)
    
    fs = xml.lfpSampleRate;
    
end

if isempty(nchannels)
    nchannels = xml.nChannels;
end
ok = dir(lfpfile);
fileDur = (ok.bytes/2/nchannels/fs);


if isempty(freqs)
    
    freqs = logspace(log10(1),log10(1000),100);
    
end

edge = freqs(1)*3;

numpts = round(2*edge*fs);


wavspec = nan(length(freqs),numpts,length(events));

for i = 1:length(events)
    
    if (events(i)-edge >0) && (events(i)+edge <= fileDur)
        dat = LoadBinary(lfpfile,'nchannels',nchannels,'channels',channel+1,'frequency',fs,'onset',events(i)-edge,'duration',2*edge);
        d = awt_freqlist(double(dat),fs,freqs);
        wavspec(:,:,i) = abs(d);
        
    end
    
end