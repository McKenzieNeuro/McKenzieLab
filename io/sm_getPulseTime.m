function [onset,offset,pulse,center] = sm_getPulseTime(fname,ch,type,varargin)

%fname = name of the analogin dat file
%ch = which channel to detect spikes
%type  = pulse or sine (smooth input)
%varargin -
%   thres - absolute threshold for detection
%   events - start stop events for the experiment, so you don't have to
%            look over entire session


onset = [];
offset = [];
pulse = [];
center = [];
thres = [];
thres1 = [];
evs = zeros(1,2);
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'thres'
            thres = varargin{i+1};
        case 'thres1'
            thres1 = varargin{i+1};
        case 'events'
            evs = varargin{i+1};
    end
    
end
if exist([fname(1:end-3) 'xml'],'file')
    [xml, rxml] = LoadXml([fname(1:end-3) 'xml']);
    
else
    
    system(['neuroscope ' fname]);
    [xml, rxml] = LoadXml([fname(1:end-3) 'xml']);
end

nChannels = xml.nChannels;
Fs = xml.SampleRate;

for ii = 1:size(evs,1)
    
    if ~any(evs)
        d = LoadBinary(fname,'nchannels',nChannels,'channels',ch,'precision','uint16');
        t1 = 0;
    else
        
        
        
        d = LoadBinary(fname,'nchannels',nChannels,'channels',ch,'start',evs(ii,1),'duration',diff(evs(ii,:)),'frequency',Fs,'precision','uint16');
        t1 = evs(ii,1);
    end
    %
%      if Fs > 3000
%         Fs = Fs/10;
%          d = d(1:10:end);
%      end
%     


   d = Smooth(double(d - median(d)),10);
    %d = double(d - median(d(1:10)));
 %   subtract baseline (may not be at zero)
    
    
    
    
    
    if strmatch(type,'pulse')
        if isempty(thres)
            thres = .10*max(d);
        end
        onset{ii} = find(diff([d]>thres)>0)-1;
        offset{ii} = find(diff([d]>thres)<0)-1;
        
        
        if any(onset{ii})
            for i = 1:length(onset{ii})
                pulse{ii}(i) = mean(d(onset{ii}(i):offset{ii}(i)));
                
            end
            onset{ii} = onset{ii}/Fs +t1;
            offset{ii} = offset{ii}/Fs +t1;
            center{ii} = mean([offset{ii}  onset{ii}],2);
            
        end
    elseif strmatch(type,'sine')
        
        k = gaussian2Dfilter([Fs 1],Fs*.01);
        
        d = nanconvn(d,k);
        
        %sometimes the spike2 software was switched on part way
        %so to find baseline, first see fi there are two stable values
        if isempty(thres)
            thres = .10*max(d);
        end
        
        
        
        [pulse{ii},pks] = findpeaks(d,'MinPeakDistance',Fs*.05,'MinPeakHeight',thres);
        
        
        
        
        %get noise around peaks
        dt  = -round(Fs*2):round(Fs*2);
        pk_id = repmat(pks,1,length(dt)) + repmat(dt,length(pks),1);
        pk_id(pk_id < 1) = 1;
        pk_id(pk_id > length(d)) =  length(d);
        idx = setdiff(1:length(d),pk_id(:));
        if isempty(thres1)
            thres1 = prctile(d(idx),99.99);
        end
        
        onset{ii} = [];
        offset{ii} = [];
        for i = 1:length(pulse{ii})
            onset_t = find(diff(d(pk_id(i,:))>thres1)>0);
            offset_t = find(diff(d(pk_id(i,:))>thres1)<0);
            
            if length(onset_t) >1 || length(offset_t) >1
                
                if length(offset_t) == length(onset_t)
                    [~,kp] = max(offset_t - onset_t);
                elseif offset_t(end)>onset_t(end) && onset_t(1) < offset_t(1)
                    [~,kp] = max(offset_t(end-1) - onset_t);
                elseif onset_t(end)>offset_t(end) && onset_t(1) < offset_t(1)
                    [~,kp] = max(offset_t - onset_t(1:end-1));
                elseif onset_t(1) > offset_t(1) && length(offset_t) > length(onset_t)
                    [~,kp] = max(offset_t(2:end) - onset_t);
                end
                onset_t = onset_t(kp);
                offset_t = offset_t(kp);
            end
            
            if any(onset_t)
                try
            onset{ii}(i) = [pks(i)/Fs + dt(onset_t)/Fs]';
                catch
                    disp('here')
                end
            offset{ii}(i) = [pks(i)/Fs + dt(offset_t)/Fs]';
            end
        end
        
        if ~isempty(offset{ii}) &&  ~isempty(onset{ii})
            kp = [offset{ii} - onset{ii}]>0;
            pks= pks(kp);
            pulse{ii} = pulse{ii}(kp);
            
            
            onset{ii}= onset{ii}(kp) +t1;
            offset{ii}= offset{ii}(kp)+t1;
            center{ii} =   (onset{ii} + offset{ii})/2;
            
        else
            
            onset{ii} = [];
            offset{ii}= [];
            center{ii} = [];
        end
        
        
    end
    
end

end