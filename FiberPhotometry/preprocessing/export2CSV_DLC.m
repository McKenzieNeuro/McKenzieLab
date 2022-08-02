% get position from DLC
dirName = 'R:\DANEHippocampalResponse\NE2h2 (Named NE2m3)\Novel Environment\NE2m3\NE2m3-210819-123801';
dirName = 'R:\DANEHippocampalResponse\DA2h1\Food Reward\DA2h1-211112-133652';

TDTdata = TDTbin2mat(dirName);

[signal_DFoF,ts_neural,fs] = sm_getSignal_DFoF(pwd,'streams',{'x465A','x405A'},'isosbestic','x405A');
%%
nFrames = length(TDTdata.epocs.Cam1.onset);



ts_video = TDTdata.epocs.Cam1.onset;

good_ix = 1:length(ts_video); % assume extra frames are at end and drop

%upSample position



TDTdata = TDTbin2mat(dirName);

%%

% find DLC CSV file (unfiltered)
fils = dir(pwd);
kp = cellfun(@any,regexp({fils.name}','DLC')) & cellfun(@any,regexp({fils.name}','.csv')) & ~cellfun(@any,regexp({fils.name}','filtered'));

if isempty(kp) || sum(kp)>1
    error('DLC csv file not found or multiple found')
end

csvF = fils(kp);
DLC_pos = readtable(csvF.name);
LE = cell2mat([table2cell(DLC_pos(:,2)) table2cell(DLC_pos(:,3))]);
RE = cell2mat([table2cell(DLC_pos(:,5)) table2cell(DLC_pos(:,6))]);
TB = cell2mat([table2cell(DLC_pos(:,8)) table2cell(DLC_pos(:,9))]);

%% get velocity for each body part
[t,x1,y1,vx1,vy1,ax1,ay1] = KalmanVel(LE(good_ix,1),LE(good_ix,2),ts_video,2);
[t,x2,y2,vx2,vy2,ax2,ay2] = KalmanVel(RE(good_ix,1),RE(good_ix,2),ts_video,2);
[t,x3,y3,vx3,vy3,ax3,ay3] = KalmanVel(TB(good_ix,1),TB(good_ix,2),ts_video,2);
%% get velocity of head 
LE_v = abs([vx1+(vy1*j)]);
RE_v = abs([vx2+(vy2*j)]);

vel = mean([LE_v RE_v],2);

LE_a = abs([ax1+(ay1*j)]);
RE_a = abs([ax2+(ay2*j)]);

acc = mean([LE_a RE_a],2);

vel_upsample = interp1(ts_video,vel(good_ix,:),ts_neural);
acc_upsample = interp1(ts_video,acc(good_ix,:),ts_neural);

%%

load('newContext.mat')

ts_homecage = [cell2mat(data(cellfun(@any,regexp(data(:,1),'home')),2))];
ts_new = [cell2mat(data(cellfun(@any,regexp(data(:,1),'new')),2))];


HC = [0 ts_new(1);ts_homecage(1) ts_new(2); ts_homecage(2) ts_new(3);ts_homecage(3) ts_neural(end)];
ConA = [ ts_new(1) ts_homecage(1)];
ConB = [ ts_new(2) ts_homecage(2)];
ConC = [ ts_new(3) ts_homecage(3)];


inHC = InIntervals(ts_neural,HC);
inConA = InIntervals(ts_neural,ConA);
inConB= InIntervals(ts_neural,ConB);
inConC = InIntervals(ts_neural,ConC);

%%
data_mat = [ts_neural(:) double(signal_DFoF(:)) vel_upsample(:) acc_upsample(:) inHC(:) inConA(:) inConB(:) inConC(:)];
%%

% fname = 'DA_animals_NEE-210818-135257_NE2m3-210819-123801_Cam1DLC_resnet50_OLMJan10shuffle1_100000.csv';

% data = readmatrix(fname);
% 
% LE = [data(:,2) data(:,3)];
% RE = [data(:,5) data(:,6)];
% TB = [data(:,8) data(:,9)];
% 
% 
% 
% %%% NOTE THE nframes must be the same as the number of frames analyzed by
% %%% DLC - where did the extra frames come from??




