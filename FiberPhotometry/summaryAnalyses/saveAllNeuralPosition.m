function saveAllNeuralPosition(dirName)

% this script takes the output from TDT and the output from DLC and saves
% into a standard sessiondata struct. The neural data will be saved with the photobleach correction
% and isosbestic correction (see sm_getSignal_DFoF)
% the saved struct (sessiondata) has the following organization:
%
%
% sessiondata.stream = list of stream names (should be LED wavelength) for the neural data
%                      (as set in TDT);
% sessiondata.path = directory path for raw data;
% sessiondata.subject = subject name;
%
% sessiondata.task = what the animal was doing;
%
% sessiondata.neural.fs_neural = sampling rate of the neural acquisition (from TDT);
% sessiondata.neural.signal_DFoF = corrected neural data;
% sessiondata.behavior.ts_video = time stamps (s) of the video (according to TDT clock);
% sessiondata.behavior.position.left_ear = position of the left ear in pixels;
% sessiondata.behavior.position.right_ear = position of the right ear in pixels;
% sessiondata.behavior.position.tail_base = position of the tail base in pixels;
% sessiondata.behavior.vel = velocity (in pixels/s);
%
% sessiondata.behavior.acc = acceleration (in pixels/s);
% sessiondata.behavior.pixel2cm = conversion factor from pixels to cm;





cd(dirName)

[a,subj] = fileparts(dirName);
sl = regexp(dirName,'\');

% on Sam's computer the task name is in between the 3rd and 4th slash -
% check how the R drive is mapped!

task = dirName(sl(3)+1:sl(4)-1);
% check which stream
TDTdata = TDTbin2mat(dirName);

% get the photobleached/ isosbestic corrected neural signal (signal_DFoF)
if isfield(TDTdata.streams,'x465A')
    stream = {'x465A','x405A'};
    [signal_DFoF,ts_neural,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x405A','baseline',[5 500]);
else
    stream = {'x450D','x500D'};
    [signal_DFoF,ts_neural,fs] = sm_getSignal_DFoF(dirName,'streams',stream,'isosbestic','x450D');
    
end
%%

% get the time stamps for each video frame according to the TDT clock
nFrames = length(TDTdata.epocs.Cam1.onset);



ts_video = TDTdata.epocs.Cam1.onset;

good_ix = 1:length(ts_video); % assume extra frames are at end and drop

%upSample position




%%

% find DLC CSV file (unfiltered)
fils = dir(dirName);
kp = cellfun(@any,regexp({fils.name}','DLC')) & cellfun(@any,regexp({fils.name}','.csv')) & ~cellfun(@any,regexp({fils.name}','filtered'));

if ~any(kp) || sum(kp)>1
    warning('DLC csv file not found or multiple found')
    
elseif any(kp)
    % read the DLC CSV file and get the left ear (col 2/3), right ear (col 5/6)
    % and tailbail (col 8/9)
    
    csvF = fils(kp);
    DLC_pos = readtable(csvF.name);
    LE = cell2mat([table2cell(DLC_pos(:,2)) table2cell(DLC_pos(:,3))]);
    RE = cell2mat([table2cell(DLC_pos(:,5)) table2cell(DLC_pos(:,6))]);
    TB = cell2mat([table2cell(DLC_pos(:,8)) table2cell(DLC_pos(:,9))]);
    
    %% get velocity for each body part
    [t,x1,y1,vx1,vy1,ax1,ay1] = KalmanVel(LE(good_ix,1),LE(good_ix,2),ts_video,2);
    [t,x2,y2,vx2,vy2,ax2,ay2] = KalmanVel(RE(good_ix,1),RE(good_ix,2),ts_video,2);
    %[t,x3,y3,vx3,vy3,ax3,ay3] = KalmanVel(TB(good_ix,1),TB(good_ix,2),ts_video,2);
    %% get velocity of head
    LE_v = abs([vx1+(vy1*j)]);
    RE_v = abs([vx2+(vy2*j)]);
    
    vel = mean([LE_v RE_v],2);
    
    LE_a = abs([ax1+(ay1*j)]);
    RE_a = abs([ax2+(ay2*j)]);
    
    acc = mean([LE_a RE_a],2);
    
    %%
    
    % save all of the data into the output struct
    sessiondata.behavior.ts_video = ts_video;
    sessiondata.behavior.position.left_ear = LE;
    sessiondata.behavior.position.right_ear = RE;
    sessiondata.behavior.position.tail_base = TB;
    sessiondata.behavior.vel = vel;
    
    sessiondata.behavior.acc = acc;
    sessiondata.behavior.pixel2cm = nan;
    
    
end
sessiondata.stream = stream;
sessiondata.path = dirName;
sessiondata.subject = subj;

sessiondata.task = task;

sessiondata.neural.fs_neural = fs;
sessiondata.neural.signal_DFoF = signal_DFoF;

save('sessiondata','sessiondata','-v7.3')


%%


%vel_upsample = interp1(ts_video,vel(good_ix,:),ts_neural);
%acc_upsample = interp1(ts_video,acc(good_ix,:),ts_neural);
%
% load('newContext.mat')
%
% ts_homecage = [cell2mat(data(cellfun(@any,regexp(data(:,1),'home')),2))];
% ts_new = [cell2mat(data(cellfun(@any,regexp(data(:,1),'new')),2))];
%
%
% HC = [0 ts_new(1);ts_homecage(1) ts_new(2); ts_homecage(2) ts_new(3);ts_homecage(3) ts_neural(end)];
% ConA = [ ts_new(1) ts_homecage(1)];
% ConB = [ ts_new(2) ts_homecage(2)];
% ConC = [ ts_new(3) ts_homecage(3)];
%
%
% inHC = InIntervals(ts_neural,HC);
% inConA = InIntervals(ts_neural,ConA);
% inConB= InIntervals(ts_neural,ConB);
% inConC = InIntervals(ts_neural,ConC);


%data_mat = [ts_neural(:) double(signal_DFoF(:)) vel_upsample(:) acc_upsample(:) inHC(:) inConA(:) inConB(:) inConC(:)];




