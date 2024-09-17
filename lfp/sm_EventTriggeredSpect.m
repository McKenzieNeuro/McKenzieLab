function h = sm_EventTriggeredSpect(lfpfile,events,varargin)
% This function calculated the event triggered spectrogram from a single
%channels of a flat binary file and returns that spectra and makes a plot
%
%inputs
%lfpfile = filename of a *.lfp binary file
%events = Nx1 matrix of event time stamps (s)
%
%optional inputs
%freqs = list of frequencies to calculate the wavelet spectra
%savefig = true/(false) to save a file
%plotIntervals = 1x2 matrix of [pre post] times to get spectrogram (s)
%figureName = filename of figure
%fs = sampling rate of file (1250)
%channel = channel to read from binary, base 1 (default = 1)
%
% example
% h = sm_EventTriggeredSpect(mydat.dat,events,'fs',20000,'channel',6);
%%

p = inputParser;


addParameter(p,'freqs',[],@isvector)
addParameter(p,'savefig',false,@islogical)
addParameter(p,'plotIntervals',[50 50],@isvector);
addParameter(p,'figureName','myPETH.fig',@ischar);
addParameter(p,'fs',1250,@isnumeric);
addParameter(p,'channel',1,@isnumeric);
addParameter(p,'nchannels',[],@isnumeric);

parse(p,varargin{:})


freqs = p.Results.freqs;
savefig = p.Results.savefig;
plotIntervals = p.Results.plotIntervals;
figureName = p.Results.figureName;
channel = p.Results.channel;
nchannels = p.Results.nchannels;

fs = p.Results.fs;





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


dur = plotIntervals(2) + plotIntervals(1);
numpts = round(dur*fs);

% row = freq, column = time, 3rd dim = # events
wavspec = nan(length(freqs),numpts,length(events));

%loop over events
for i = 1:length(events)
    
    % silently skip events that fall too close to the edge of the file
    % FUTURE - DON'T MAKE THIS SILENT
    
    
    if (events(i)-plotIntervals(1) >0) && (events(i)+plotIntervals(2) <= fileDur)
        
        %load the data at the right time
        dat = LoadBinary(lfpfile,'nchannels',nchannels,'channels',channel,'frequency',fs,'start',events(i)-plotIntervals(1),'duration',dur);
        
        %take wavelet decomposition
        d = awt_freqlist(double(dat),fs,freqs);
        
        %take the power
        wavspec(:,:,i) = abs(d)';
        
    end
    
end


nfreq = length(freqs);
ts = ((1:numpts)/fs)  - plotIntervals(1);
h = figure;
imagesc(ts,[],nanmean(wavspec,3))
set(gca,'ytick',1:10:nfreq,...
    'yticklabel',round(freqs(1:10:end)*10)/10,'ydir','normal')
%
if savefig
    
    saveas(h,figureName)
end
end
