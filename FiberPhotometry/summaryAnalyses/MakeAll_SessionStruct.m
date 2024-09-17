%get all Session structs for novel context

topDir = 'R:\DANEHippocampalResponse\DA3h7\';
fils = getAllExtFiles(topDir,'mat',1);
kp = contains(fils,'Novel Env') & contains(fils,'DA3');% & ~contains(fils,'NE2h10');
dirs = unique(fileparts(fils(kp)));


%%

 makeField = {'MergeDLCwNeural','MergeContextEntryandEdges','alignEdgPosition_pix2cm','MergeVelocityAcceleration',...
        'MergeTimeFromTransition','MergeDist2Edge'};
    
    
for i = 1:length(dirs)
    sm_MakeSessionStruct(dirs{i},'makeField',makeField)
    i
end
%%

for i = 1:length(dirs)
    sm_MakeSessionStruct(dirs{i},'makeField',{'newBaseline'})
    i
end

%%

%once all sessions are processed then define # of days per session

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

