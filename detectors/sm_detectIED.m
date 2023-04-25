function st = sm_detectIED(dirN,varargin)


events.description =[];
events.time =[];
p = inputParser;
addParameter(p,'basename','amplifier',@isstr);
addParameter(p,'baselineDur',3600,@isnumeric);
addParameter(p,'st',[],@isnumeric);


parse(p,varargin{:});

basename = p.Results.basename;
baselineDur = p.Results.baselineDur;
st = p.Results.st;


lfpfil = [dirN filesep basename '.lfp'];
xmlfil = [dirN filesep basename '.xml'];


if ~exist(xmlfil)
    error('make xml file')
end
xml = LoadXml(xmlfil);


if ~exist(lfpfil)
    bz_LFPfromDat(dirN,'basename',basename)
end

fs = xml.lfpSampleRate;
inf = dir(lfpfil);
durFil =inf.bytes/2/fs/xml.nChannels;

if isempty(st) || (length(st) ~=xml.nChannels)
    get_std = true;
else
     get_std = false;
end



%indiced for baseline calc



k  = gaussian2Dfilter([ 2*fs 1],fs/500);
nSTD = 10;
% loop throug channels
for i = 1:xml.nChannels
    
    
    
    
    tmp = LoadBinary(lfpfil,'nchannels',xml.nChannels,'channels',i,'frequency',fs);
    tmp = BandpassFilter(double(tmp),fs,[20 200]);
    tmp = InstAmplitude(tmp);
    if i ==1
        ts = (1:length(tmp))/fs;
    end
    
    if get_std
        
        if isinf(baselineDur)
            baselineDur = ts(end);
        end
        
        ixx = floor(1:baselineDur*fs);
        st(i) = std((double(tmp(ixx))));
    end
    
    IED = ts((tmp)>nSTD*st(i));
    
    
    
    
    %get width
    width = nan(length(IED),1);
    for k1 = 1:length(IED)
        tmp =  LoadBinary(lfpfil,'nchannels',xml.nChannels,'channels',i,'frequency',fs,'start',IED(k1)-.25,'duration',.5);
        if max(InstAmplitude(BandpassFilter(double(tmp),1250,[200 500])))<1000
            
            tmp = BandpassFilter(double(tmp),fs,[20 200]);
            tmp = InstAmplitude(tmp);
            ii = round(length(tmp)/2);
            %determine if its a peak or valley
            
            [height,loc,w] = findpeaks(double(tmp));
            ix = bestmatch(ii,loc);
            if any(loc)
                width(k1) =  w(loc==ix)/fs;
            end
            
        end
        
    end
    
    IED1 = IED(width<.075 & width>.01);
    % only keep IEDs
    kp = diff([IED1])>.2 ;
    IED1 = IED1(kp);
    
    
    
    
    
    
    % determine if it is noise, IED, or spike train
    
    
    
    
    events.description = [  events.description;repmat({['IED: Ch' num2str(i)]},length(IED1),1)];
    events.time = [ events.time;IED1(:)];
    
    
    
    
end
[events.time,b] = sort(events.time);
events.description =  events.description(b);
SaveEvents([dirN filesep 'autoDetect.evt.IED'],events)
end

