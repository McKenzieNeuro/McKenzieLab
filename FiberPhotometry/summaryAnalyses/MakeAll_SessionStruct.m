%get all Session structs for novel context

topDir = 'R:\DANEHippocampalResponse\';
fils = getAllExtFiles(topDir,'txt',1);
kp = contains(fils,'Novel Env') & contains(fils,'NEmut');% & ~contains(fils,'NE2h10');
dirs = unique(fileparts(fils(kp)));



%%
for i = 1:length(dirs)
    cd(dirs{i})
    disp(dirs{i})
    if ~isfile([dirs{i},'\ContextTransitionRevised.mat'])
        vidpath=getAllExtFiles(dirs{i},'avi',1);
        videoContactSheet(vidpath{:});
        disp('Label contexts in format: homecage,BB')
        sm_labelTDTvideo();
        pause()
    end
    if ~isfile([dirs{i},'\arena_edges.mat'])
        disp('Label contexts in format: homecage1,BB1,homecage2')
        sm_getArenaEdges(dirs{i});
        pause()
    end
    

end

%%
 makeField = {'MergeDLCwNeural','MergeContextEntryandEdges','alignEdgPosition_pix2cm','MergeVelocityAcceleration',...
        'MergeTimeFromTransition','MergeDist2Edge'};
    
for i = 1:length(dirs)
    cd(dirs{i})
    sm_MakeSessionStruct(dirs{i},'makeField',makeField)
    i
end
%%

for i = 1:length(dirs)
    cd(dirs{i})
    sm_MakeSessionStruct(dirs{i},'makeField',{'newBaseline'})
    i
end

%%

%once all sessions are processed then define # of days per session
cd(dirs{1})
sm_MakeSessionStruct(dirs{1},'makeField',{'MergeDaysInContext'})

%%

topDir = 'R:\DANEHippocampalResponse';
fils = getAllExtFiles(topDir,'mat',1);
kp = contains(fils,'SOR') & contains(fils,'NE2');
dirs = unique(fileparts(fils(kp)));


%%

 makeField = {'MergeDLCwNeural','MergeContextEntryandEdges','alignEdgPosition_pix2cm','MergeVelocityAcceleration',...
        'MergeTimeFromTransition','MergeDist2Edge','MergeTimeFromObjectIntro'};
    
    
for i = 1:length(dirs)
    sm_MakeSessionStruct(dirs{i},'makeField',makeField)
    i
end



%%



%%

topDir = 'R:\DANEHippocampalResponse';
fils = getAllExtFiles(topDir,'mat',1);
kp = contains(fils,'Linear') & contains(fils,'NE2');
dirs = unique(fileparts(fils(kp)));


%%

 makeField = {'MergeDLCwNeural','MergeContextEntryandEdges_linear','alignEdgPosition_pix2cm','MergeVelocityAcceleration',...
        'MergeTimeFromTransition','MergeRewardTime'};
    
    for i = 1:length(dirs)
    sm_MakeSessionStruct(dirs{i},'makeField',makeField)
    i
    end

    %%
    
for i = 1:length(dirs)
    cd(dirs{i})
    
    if exist('sessiondata.mat')
   load('sessiondata.mat')
   
   [~,b1] = histc(sessiondata.behavior.rewardTimeLeft,sessiondata.behavior.ts_video);
   [~,b2] = histc(sessiondata.behavior.rewardTimeRight,sessiondata.behavior.ts_video);
   hold on
   plot(sessiondata.behavior.position.left_ear_cor(:,1),sessiondata.behavior.position.left_ear_cor(:,2),'.')
   plot(sessiondata.behavior.position.left_ear_cor(b1,1),sessiondata.behavior.position.left_ear_cor(b1,2),'x')
    plot(sessiondata.behavior.position.left_ear_cor(b2,1),sessiondata.behavior.position.left_ear_cor(b2,2),'o')
   waitforbuttonpress
    close all
    end
end

%%

% do injections session
topDir = 'R:\ASommer\FP experiments DA-NE\Desipramine injections';
fils = getAllExtFiles(topDir,'mat',1);
kp = cellfun(@any,regexp(fils,'inject')) ;
fils = fils(kp);
[dirs,bi] = fileparts(fils);

   
    
for i = 1:length(dirs)
       sm_MakeSessionStruct(dirs{i},'makeField',{'newBaseline'})

end


%%
%needs DLC

%R:\DANEHippocampalResponse\NE2h7\Linear Track\NE2h7-220720-093937
%R:\DANEHippocampalResponse\NE2h8\LinearTrack\NE2h8-220623-121635

