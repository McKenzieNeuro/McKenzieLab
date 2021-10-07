
%%rename intan files (cd into directory and ctr -enter)

warning off

dirName  = pwd;

basename = bz_BasenameFromBasepath(dirName);

sm_Process_Intan_kilosort(basename,dirName);

%%
% RUN KILOSORT
dirName  = pwd;

savepath = sm_KiloSortWrapper('basepath',dirName,'config_version','McKenzie','intan_version','combined');

%%
%make an event file  with basename.evt.sti

dirName  = pwd;
basename = bz_BasenameFromBasepath(dirName);
analogfils_dat = [basename '_analogin.dat'];
%define experiment start and stop based off of pulses -  pulse_on pulse_off

if ~exist(analogfils_dat)
    
    analogfils_dat = [basename '.dat'];
    
end

system(['neuroscope ' analogfils_dat])



%%
%%get pulse times

dirName  = pwd;
basename = bz_BasenameFromBasepath(dirName);
analogfils_dat = [basename '_analogin.dat'];
filename = getAllExtFiles(pwd,'sti',0);
events = LoadEvents(filename{1});

if ~exist(analogfils_dat)
    
    analogfils_dat = [basename '.dat'];
    xml = LoadXml(fullfile(dirName,[basename '.xml']));
    ch = xml.AnatGrps(end-1).Channels+1;
    sm_MakeEventTime([dirName filesep analogfils_dat],events,'ch',ch(1:4),'thres1',500,'thres2',300); %(1:# is channels 1-#, or [2 4] if only 2 and 4)
    
else
    
    sm_MakeEventTime([dirName filesep analogfils_dat],events,'ch',[ 1 ],'thres1',3000,'thres2',300); %(1:# is channels 1-#, or [2 4] if only 2 and 4)
    
    
end



%%

%save spikes into matlab
spikes = bz_GetSpikes;

%%

ch1 = 25;%ripple channel (base 0)
chN = 8;
ripples = sm_FindRipples(pwd,ch1,'noise', chN,'restrict',eps);

%%

sm_TheStateEditor
%%

%line up video



basepath = pwd;
if basepath(end) == filesep
    basepath =  basepath(1:end-1);
end
% find all avis
avifils = getAllExtFiles(basepath,'avi',1);

%get blinking light
LED = [];
for i = 1:length(avifils)
    
    if i==1
        [temp,in,threshF] = sm_GetLEDfromAVI(avifils{i},['blink_ON-' sprintf('%02.0f',i) '.mat'],[],[],[],false);
    else
        [temp,in,threshF] = sm_GetLEDfromAVI(avifils{i},['blink_ON-' sprintf('%02.0f',i) '.mat'],in,threshF,false);
    end
    
    LED = [LED;temp];
end

%%

dirName = pwd;
sl = regexp(dirName,filesep);
basename = dirName(sl(end)+1:end);

% kp = (dwnLED-upLED)>1;
% upLED= upLED(kp);
% dwnLED= dwnLED(kp);
%  kp = diff([0;upLED])>100;
%  upLED= upLED(kp);
% dwnLED= dwnLED(kp);
digitalin_ch =1;

[ups,dwns]  =  sm_getBaslerPos(dirName,basename,digitalin_ch,30000);



%%


%get LED files

fils =  getAllExtFiles(pwd,'mat',0) ;
kp = cellfun(@any,regexp(fils,'LED'));
pos_fils = fils(kp);

kp = cellfun(@any,regexp(fils,'blink'));
LED_fils = fils(kp);


kp = cellfun(@any,regexp(fils,'TTL_pulse'));
TTL = fils(kp);
v = load(TTL{1});
ups = v.ups;
dwns = v.dwns;


LED = [];
for i = 1:length(LED_fils)
    v = load(LED_fils{i});
    LED = [LED;v.whl(:,1)];
end


whl = [];
for i = 1:length(pos_fils)
    v = load(pos_fils{i});
    whl = [whl;v.whl];
end

%%

%%
upLED = find(diff(smooth(LED(:,1),10) > 10)>0);
dwnLED = find(diff(smooth(LED(:,1),10) > 10)<0);
%kp = (dwnLED-upLED)>10;
%kp = (dwns-ups)  >1;
kp = [true;diff(ups)>5];
ups = ups(kp);
dwns = dwns(kp);

%upLED = upLED(3:end-1);
%ups = ups(1:end-1);
figure
plot(diff(upLED)/28,'r') % plot from video
hold on
plot(diff(ups(1:end)),'k') % plot fro Intan
%% load DLC
fils = getAllExtFiles(pwd,'.h5',1);
whl =[];
for i = 1:length(fils)
    data = h5read(fils{i},'/df_with_missing/table');
    
    pt1 = [data.values_block_0(1,:);data.values_block_0(2,:)]';
    pt2 = [data.values_block_0(4,:);data.values_block_0(5,:)]';
    whl = [whl;pt1 pt2];
end
%%

[LEDts,b] = sort([upLED]);
pulsets = [ups];


if length(pulsets)~=length(LEDts)
    error('LED pulse does not match TTL pulse')
end


pulsets = pulsets(b);
if ~all(diff(pulsets)>0)
    disp('temporal mismatch between LED up/LEDdown')
end



%%
basename = bz_BasenameFromBasepath(pwd);
filename = [basename '.dat'];
dat = LoadBinary(filename,'nchannels',39,'channels',37,'frequency',30000);

%%


dt = mode(diff(upLED));
upLED1 = upLED(3)+(0:dt:dt*93-dt);
upLED2 = upLED(97) + (0:dt:dt*89-dt);
upLED = [upLED1(:);upLED2(:)];

[LEDts,b] = sort([upLED]);
pulsets = [ups];


xx = nanmean(whl(LEDts(1):LEDts(end),[1 3 ]),2);
yy = nanmean(whl(LEDts(1):LEDts(end),[2 4 ]),2);

[x,y] = sm_FixPosition(xx,yy,LEDts(1):LEDts(end));

[t,x,y,vx,vy,ax,ay] = KalmanVel(x,y,1:length(x),2);

if length(pulsets)~=length(LEDts)
    error('LED pulse does not match TTL pulse')
end


pulsets = pulsets(b);
if ~all(diff(pulsets)>0)
    disp('temporal mismatch between LED up/LEDdown')
end

ts1 = interp1(LEDts,pulsets,LEDts(1):LEDts(end));
[~,b] = histc(ups,(1:length(dat))/30000);
[~,b1] = histc(upLED,LEDts(1):LEDts(end));
%%

for i = 1:10
figure

plot(dat(b(i):b(i)+3000))
hold on
plot(vx(b1(i):b1(i)+3000))
uiwait
end
%%


%
% ts1 = ups';
%  x = nanmean(whl(:,[1 3 ]),2);
%   y = nanmean(whl(:,[2 4 ]),2);
whl(whl<=0) = nan;


kp = ~isnan(y);

t = ts1;
x1 = x(kp);
y1 = y(kp);
t1 = t(kp);
x = interp1(t1,x1,t);
y = interp1(t1,y1,t);
%%

[len,pos,xp,yp,del] =sm_linearizePositionTrack(x(:),y(:));

%len = interp1(t1,len1,t)';
idx = [find(diff(isnan([nan;len]))<0) find(diff(isnan([len;nan]))>0)];


[status, interval] = InIntervals(ts1,ts1(idx));


n1 = histoc(interval(interval>0),1:size(idx,1));

len_ep = mat2cell(len(status),n1);
ts_ep = mat2cell(ts1(status)',n1);
kp = cellfun(@(a) range(a,1),len_ep)>10 & cellfun(@range,ts_ep)<30;

len_ep = len_ep(kp);
pos_inf.len_ep = len_ep;
pos_inf.ts_ep = ts_ep(kp);
pos_inf.in_eps =  cellfun(@(a) mean(diff(a))>0,len_ep);
pos_inf.out_eps =  cellfun(@(a) mean(diff(a))<0,len_ep);
pos_inf.x = x(:);
pos_inf.y = y(:);
pos_inf.lin_pos = len;
pos_inf.ts = ts1;
%   ed


save('position_info.mat','pos_inf')

%%

sm_plotPlaceFieldsLinear
%%
sm_plotPlaceFields2D

