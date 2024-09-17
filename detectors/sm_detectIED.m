function st = sm_detectIED(dirN,varargin)
%% this function will detect IEDs from a binary file
% Input
% dirN = directory name where a *.lfp file is stored (must have
% accompanying xml to parse file with the same basename
% Optional inputs
% basename = the name of the lfp and xml file (e.g. basename.xml)
% baselineDur = duration (s) in which to calculate mean and std over which IEDs
%               are detected
% st = standard deviation of the baseline over which an IED is detected
% nSTD = # of standard deviations over which an IED should be detected
% ch = [1 N] list of channels (base 1) to read from basename.lfp
% ampThres = maximum size of event (uV) over which it is noise
% start_time = time (s) to begin reading data
% analyze_duration = length (s) to analyze relative to start_time
% pass_band = bandpass limits to detect IED
% conversionConstant = constant to get binary file into mV
% high_pass = bandpass to find sharp noise events
% width_thres = [minW maxW] min and max duration of IED
% IED_ISI_thres = minimum inter IED interval
%
% OUTPUT
% st = standard deviation of the baseline
% saves and FMA event file 'autoDetect.evt.IED'
%     to load these ev = LoadEvents('autoDetect.evt.IED');
%
% example call
% st = sm_detectIED(dirN,'basename','continuous','ch',60);
%%

% parse inputs


p = inputParser;
addParameter(p,'basename','amplifier',@isstr);
addParameter(p,'baselineDur',3600,@isnumeric);
addParameter(p,'start_time',0,@isnumeric);
addParameter(p,'analyze_duration',[],@isnumeric);
addParameter(p,'st',[],@isnumeric);
addParameter(p,'nSTD',8,@isnumeric);
addParameter(p,'ch',[],@isnumeric);
addParameter(p,'ampThres',195,@isnumeric);
addParameter(p,'pass_band',[20 200],@isnumeric);
addParameter(p,'conversionConstant',.195,@isnumeric);
addParameter(p,'high_pass',[200 500],@isnumeric);
addParameter(p,'width_thres',[.01 .075],@isnumeric);
addParameter(p,'IED_ISI_thres',.2,@isnumeric);


parse(p,varargin{:});

basename = p.Results.basename;
baselineDur = p.Results.baselineDur;
st = p.Results.st;
ch = p.Results.ch;
nSTD = p.Results.nSTD;
ampThres = p.Results.ampThres;
start_time = p.Results.start_time;
analyze_duration = p.Results.analyze_duration;
pass_band = p.Results.pass_band;
conversionConstant = p.Results.conversionConstant;
high_pass = p.Results.high_pass;
width_thres = p.Results.width_thres;
IED_ISI_thres = p.Results.IED_ISI_thres;
%%

% initialize output struct
events.description =[];
events.time =[];

% define lfp and xml files
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

if isempty(analyze_duration)
   analyze_duration = durFil - start_time;
end
    
if isempty(ch)
    
    ch = 1:xml.nChannels;
end
ch = ch(:)';

if isempty(st) || (length(st) ~= length(ch))
    get_std = true;
else
     get_std = false;
end



%indiced for baseline calc



k  = gaussian2Dfilter([ 2*fs 1],fs/500);

% loop through channels
idx = 1;
for i = ch
    
    
    
    % load data
    tmp = LoadBinary(lfpfil,'nchannels',xml.nChannels,'channels',i,'frequency',fs,'start',start_time,'duration',analyze_duration);
    
    %bandpass filter data
    tmp = BandpassFilter(double(tmp),fs,pass_band);
    
    %get instantaneous amplitude of filtered data
    tmp = InstAmplitude(tmp);
    
    %define time stamps
    if idx ==1
        ts = (1:length(tmp))/fs+start_time;
    end
    
    %get standard deviation
    if get_std
        
        if isinf(baselineDur) || baselineDur > analyze_duration
            baselineDur = ts(end) - start_time;
        end
        
        ixx = floor(1:baselineDur*fs);
        st(idx) = std((double(tmp(ixx))));
    end
    
    
    % get the time points where the data is nSTD times over baseline s.t.d.
    kp = (tmp)>nSTD*st(idx);
    kp = diff([0;kp])>0;
    IED = ts(kp);
    
    
    
    % filter by width and amplidtude 
    %get width
    width = nan(length(IED),1);
    for k1 = 1:length(IED)
        tmp =  LoadBinary(lfpfil,'nchannels',xml.nChannels,'channels',i,'frequency',fs,'start',IED(k1)-.25,'duration',.5);
        tmp = double(tmp)*conversionConstant; % mulitply by intan conversion rate
        if max(InstAmplitude(BandpassFilter(tmp,fs,high_pass)))<ampThres
            
            tmp = BandpassFilter(double(tmp),fs,pass_band);
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
    
    IED1 = IED(width<width_thres(2) & width>width_thres(1));
    % only keep IEDs
    kp = diff([IED1])>IED_ISI_thres ;
    IED1 = IED1(kp);
    
    
    
    
    
    
    % determine if it is noise, IED, or spike train
    
    
    
    
    events.description = [  events.description;repmat({['IED: Ch' num2str(i)]},length(IED1),1)];
    events.time = [ events.time;IED1(:)];
    
    
    
    idx = idx+1;
end
[events.time,b] = sort(events.time);
events.description =  events.description(b);
SaveEvents([dirN filesep 'autoDetect.evt.IED'],events)
end

