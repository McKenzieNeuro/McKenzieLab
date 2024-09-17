function sm_LoadSessionStruct(sessionDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% BEFORE RUNNING THIS FUNCTION, THE VIDEOS SHOULD FIRST HAVE THE ARENA
% EDGES AND TRANSITION TIMES ENCODED!!!!!!!!!!!
%
%   (Also, run Deeplabcut first)
%
%
% encode arena edges (arena_edges.mat)
%       "R:\Analysis\McKenzieLab\FiberPhotometry\Video\sm_getArenaEdges.m"
%       arenas should be named home_cage1, context1, home_cage2, ...
%       ... context2, home_cage3, context3, home_cage4, ...etc
% 
% encode transition times (ContextTransitionRevised.mat)
%       "R:\Analysis\McKenzieLab\FiberPhotometry\Video\sm_labelTDTvideo.m"
%       Contexts can be named with any of the names in the file "R:\McKenzieLab\IPimentel\Key\ContextTransitionInstructions.docx"
%
% This function will save the following data: 
%       sessiondata.mat
%       totalDistanceFrmEdgR.mat
%     
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% merge neural data and behavior data from Deeplabcut
sm_MergeDLCwNeural(sessionDir);


%% merge context entry times and arena edges
% R:\McKenzieLab\IPimentel\Matlab\Sessiondata Scripts\forAll\Sessiondata starters\makeAll_mergeContextEntryandEdges.m

fprintf('Merging context entry times and edges...\n\n')
sm_MergeContextEntryandEdges(dirName)

 
 
 
%% align the areda edges, and convert to cm (also does a tranformation to get the skewed arena into a non-skewed form
% R:\McKenzieLab\IPimentel\Matlab\Sessiondata Scripts\forAll\Sessiondata starters\alignEdgPositionFinal_pix2cm_ip.m
fprintf('Converting from pixels to cm...\n\n')
sm_alignEdgPosition_pix2cm(dirName)



%% get velocity and acceleration in cm distance
% R:\McKenzieLab\IPimentel\Matlab\Sessiondata Scripts\Novel Environment\Make Sessiondata Variables\MakeAll_VelAccCalcCorrected_novEnv.m
fprintf('Getting velocity and acceleration...\n\n')
MakeAll_VelAccCalcCorrected_novEnv_as(dirs)


%% get time from context entry
% R:\McKenzieLab\IPimentel\Matlab\Sessiondata Scripts\Novel Environment\Make Sessiondata Variables\TimeFromContextEntryfinal_novEnv.m
fprintf('Getting time from entry...\n\n')
TimeFromContextEntryfinal_novEnv_as(dirs)


%% get number of dats that the animal has experienced the context (homecage is set to 100 days)
% R:\McKenzieLab\IPimentel\Matlab\Sessiondata Scripts\Novel Environment\Make Sessiondata Variables\makeAll_daysInContext_novEnv.m
fprintf('Getting number of days in context...\n\n')
makeAll_daysInContext_novEnv_as(dirs)


%% get how far the animal is from the arena edges in cm
% R:\McKenzieLab\IPimentel\Matlab\Sessiondata Scripts\forAll\EdgDistanceScripts\EdgLERETB.m
fprintf('Finding distance from the arena edges...\n\n')
EdgLERETB_as(dirs)


%% get GMLE models (first best)

% need to separate DA and NE animals from here on.
% fprintf('Calculating initial GLME models and finding best fit...\n\n')
% [GLME, best] = getGLME_novelAndHome_best_as(dirs);


%% compare the GMLE models to find the best fit
% fprintf('Comparing initial GLME models to find best fit...\n\n')
% best = compareGLME_as(GLME);
 
%% find variable with greatest effect 
% fprintf('Finding variable with greatest effect...\n\n')
% GLMEPlots_HCandNovel_allSubs_final_as(dirs,best)

% 
% %% get GMLE models (second best)
% 


% 
% 
% %% compare the GMLE models to find the second best fit
% 
% 
% 

end