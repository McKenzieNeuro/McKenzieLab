fils = getAllExtFiles('R:\McKenzieLab\DANEHippocampalResponse','mat',1);
kp = cellfun(@any,regexp(fils,'Novel Env'));
% kp = cellfun(@any,regexp(fils,'inear'));

fils = fils(kp);

[dirs] = fileparts(fils);
%%For 2020 use :[dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);
%


for i = 1:length(dirs)
    cd(dirs{i})
    i

    %%Checking for sessiondata
    if exist([dirs{i} filesep 'sessiondataR.mat'])
         
                %load tracking data
                load([dirs{i} filesep 'sessiondataR.mat']);
                

                ts_video = sessiondata.behavior.ts_video;
                good_ix = 1:length(ts_video); % assume extra frames are at end and drop        

                LE = sessiondata.behavior.position.left_ear_cor(:,:);
                RE = sessiondata.behavior.position.right_ear_cor(:,:);
                TB = sessiondata.behavior.position.tail_base_cor(:,:);

               %% get velocity for each body part
                [t,x1,y1,vx1,vy1,ax1,ay1] = KalmanVel(LE(good_ix,1),LE(good_ix,2),ts_video,2);
                [t,x2,y2,vx2,vy2,ax2,ay2] = KalmanVel(RE(good_ix,1),RE(good_ix,2),ts_video,2);
                [t,x3,y3,vx3,vy3,ax3,ay3] = KalmanVel(TB(good_ix,1),TB(good_ix,2),ts_video,2);
                %% get average velocity head
                LE_v = abs([vx1+(vy1*i)]);
                RE_v = abs([vx2+(vy2*i)]);
           
                vel = mean([LE_v RE_v],2);
                
                LE_a = abs([ax1+(ay1*i)]);
                RE_a = abs([ax2+(ay2*i)]);

                
                acc = mean([LE_a RE_a],2);

                sessiondata.behavior.vel_cor = vel;
                sessiondata.behavior.acc_cor = acc;


                
                save('sessiondataR.mat','sessiondata','-v7.3')
                i


    end
end
