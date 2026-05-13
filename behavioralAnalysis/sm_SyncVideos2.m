 %run from directory that has all of yoru data
clear all
dirN = pwd;

cd(dirN)

s = 3;

switch s
    
    case 1
        
        % one video/one intan
        cond = '1Intan1Vid';
    case 2
        % two videos/two merged intans
        cond = '2Intan2Vid';
    case 3
        % two videos/one intan
        cond = '1Intan2Vid';
end

%%

% merge digitalins
fout = [dirN filesep 'digitalin.dat'];

fils = getAllExtFiles(dirN,'dat',1);
kp = contains(fils,'digitalin.dat');

fils = fils(kp);

if length(fils)>1 & ~exist(fout)
    ts_dat = nan(size(fils));
    
    for i = 1:length(fils)
        
        tmp = dir(fils{i});
        
        ts_dat(i)= datenum(tmp.date);
    end
    
    [~,b] = sort(ts_dat);
    
    fils = fils(b);
    
    sm_ConcatDats(fils,fout)
end
%%

dirN = pwd;
% get the Intan sync date
LED_sync_ch = 1; % change this if different

if ~exist('TTL_pulse.mat')
    [ups,dwns]  = sm_getDigitalin(dirN,'digitalin.dat',30000,16);
else
    load('TTL_pulse.mat')
end


if length(ups{LED_sync_ch}) ~= length(dwns{LED_sync_ch})
    
    warning('up/dwn mismatch, I guess you left the LEDs on at the beginning or end of the recording')
    
    
    u = ups{LED_sync_ch};
    d = nan(size(u));
    for i = 1:length(u)
        d(i) = dwns{LED_sync_ch}(find(dwns{1}>u(i),1,'first'));
    end
    dur = d-u;
    m_dur = mode(d-u);
    kp  = dur>.8*m_dur & dur<1.2*m_dur;
    
    ups{LED_sync_ch} = u(kp);
    dwns{LED_sync_ch} = d(kp);
    
    save([dirN filesep 'TTL_pulse.mat'],'ups','dwns')
end

%%
clear whl
% get the video sync data by labeling the ROI for the LED

fils = getAllExtFiles(dirN,'mp4',1);
kp = ~contains(fils,'snapshot'); % exclude DLC
fil = fils(kp);


ts_vid = nan(size(fils));

for i = 1:length(fils)
    
    tmp = dir(fils{i});
    
    ts_vid(i)= datenum(tmp.date);
end

[~,b] = sort(ts_vid);

fils = fils(b);

for i = 1:length(fils)
    
    fout = [fils{i}(1:end-4) '_LED.mat'];
    if ~exist(fout)
        [whl{i},in,threshF,fs] = sm_ExtractLed3(fils{i});
    else
        
        tmp = load(fout,'whl','fs');
        whl{i} = tmp.whl;
        
        if isfield(tmp,'fs')
            fs = tmp.fs;
        else
            
            fs = 27.2093;
        end
    end
    
end
vidFils = fils;

%%

%get durations of all files that need to be merged
    fs_intan = 30000;
    nCh = 1;
%amplifier basename
keyword = 'digitalin.dat';

%get amplifier data

fils = getAllExtFiles(dirN,'dat',1);
kp = contains(fils,keyword);
fils = fils(kp);
if length(fils)>1
    
    
    %exclude the merged in the top directory
    subdirs = fileparts(fils);
    
    
    kp = ~strcmp(dirN,subdirs);
    fils = fils(kp);
    
    subdirs = fileparts(fils);
    
    
    
    %sort by time
    tims = cellfun(@(a) str2num(cell2mat(regexp(a,'[0-9]','match'))),subdirs);
    [~,b] = sort(tims);
    fils = fils(b);
    


    
end


for i = 1:length(fils)
    
    
    ok = dir(fils{i});
    siz = ok.bytes;
    dur(i) = siz/fs_intan/nCh/2;
    
end

%%
if ~exist('fs')
    fs = 27.2093;
end
threshF = 10; % change if it is not lining up
sm_match_LED_digitalin(dirN,ups{LED_sync_ch},dwns{LED_sync_ch},whl,dur,fs,vidFils,'merge_cond',cond,'threshF',threshF)