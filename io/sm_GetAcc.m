function [acc,frequency] = sm_GetAcc(basepath,varargin)

p = inputParser;


addParameter(p,'channels',[],@isnumeric) % base 0
addParameter(p,'st',.1,@isnumeric) 

parse(p,varargin{:})


channels = p.Results.channels;
st = p.Results.st;




basename = bz_BasenameFromBasepath(basepath);

aux = [basename '_auxiliary.dat'];


if exist(aux)
    channels = 1:3;
    xmlf = [aux(1:end-3) 'xml'];
    if ~exist(xmlf)
        
        system(['neuroscope ' aux])
    end
    xml = LoadXml(xmlf);
    
else
    
    if isempty(channels)
        
        error('user must enter channels (base 0)')
    end
    
    
    aux = [basename '.lfp'];
    
    if exist(aux)
        xmlf = [aux(1:end-3) 'xml'];
        xml = LoadXml(xmlf);
        if ~exist(xmlf)
            
            system(['neuroscope ' aux])
        end
        frequency  = xml.lfpSampleRate;
        
    else
        aux = [basename '.dat'];
        
        
        xmlf = [aux(1:end-3) 'xml'];
        
        if ~exist(xmlf)
            
            system(['neuroscope ' aux])
        end
        xml = LoadXml(xmlf);
        
        frequency  = xml.SampleRate;
    end
    
    
    
    
    
    
    
    
    
    
end

acc = LoadBinary(aux,'nchannels',xml.nChannels,'channels',channels+1,'frequency',frequency);

acc = sum(double(acc).^2,2).^.5;

st = st*frequency;
k  = gaussian2Dfilter([100*st 1],[st 1]);
acc = nanconvn(acc,k);
end