function h = sm_plotVelAcc(dirName)

% plot Da/NE as function of acceleration and velocity


cd(dirName)
%load data
TDTdata = TDTbin2mat(dirName);

[signal_DFoF,ts_neural,fs] = sm_getSignal_DFoF(pwd,'streams',{'x465A','x405A'},'isosbestic','x405A');
%%

% get video times
nFrames = length(TDTdata.epocs.Cam1.onset);



ts_video = TDTdata.epocs.Cam1.onset;

good_ix = 1:length(ts_video); % assume extra frames are at end and drop

%%

%load CSV data from DLC



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

h = figure;
subplot(2,1,1)
plot(0:5:100,avghist(vel_upsample,double(signal_DFoF),0:5:100))
subplot(2,1,2)
plot(0:5:100,avghist(acc_upsample,double(signal_DFoF),0:5:100))
end

