function sm_MergeVelocityAcceleration(dirName)
%UNTITLED7 Summary of this function goes here
%   % Will output an array, 'sessiondata.behavior.vel_cor'  and 'sessiondata.behavior.acc_cor'
% in the sessiondata.m file, that will be the length of the  number of frames in
% the videos submitted to DLC. The row number in the array will correspond with
% the frame number in the video and the contents will be the acceleration
% or velocity in cm distance and sec time

cd(dirName)
%%Checking for sessiondata
if exist([dirName filesep 'sessiondata.mat'])
    
    %load tracking data
    load([dirName filesep 'sessiondata.mat']);
    
    
    ts_video = sessiondata.behavior.ts_video;
    good_ix = 1:length(ts_video); % assume extra frames are at end and drop
    
    LE = sessiondata.behavior.position.left_ear_cor(:,:);
    RE = sessiondata.behavior.position.right_ear_cor(:,:);
    TB = sessiondata.behavior.position.tail_base_cor(:,:);
    
    %% get velocity for each body part
    [t1,x1,y1,vx1,vy1,ax1,ay1] = KalmanVel(LE(good_ix,1),LE(good_ix,2),ts_video,2);
    [t2,x2,y2,vx2,vy2,ax2,ay2] = KalmanVel(RE(good_ix,1),RE(good_ix,2),ts_video,2);
    
    %% get average velocity head
    LE_v= nan(size(sessiondata.behavior.acc));
    RE_v  = nan(size(sessiondata.behavior.acc));
    
    
    LE_v( ismember(sessiondata.behavior.ts_video,t1)) = abs([vx1+(vy1*i)]);
    RE_v( ismember(sessiondata.behavior.ts_video,t2)) = abs([vx2+(vy2*i)]);
    
    vel = mean([LE_v RE_v],2);
    
    LE_a= nan(size(sessiondata.behavior.acc));
    RE_a  = nan(size(sessiondata.behavior.acc));
    
    LE_a( ismember(sessiondata.behavior.ts_video,t1)) = abs([ax1+(ay1*i)]);
    RE_a( ismember(sessiondata.behavior.ts_video,t2))  = abs([ax2+(ay2*i)]);
    
    
    acc = mean([LE_a RE_a],2);
    
    
    
    sessiondata.behavior.vel_cor = vel;
    sessiondata.behavior.acc_cor = acc;
    
    
    
    save('sessiondata.mat','sessiondata','-v7.3')
    
end
end




