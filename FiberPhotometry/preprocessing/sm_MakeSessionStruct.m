function sm_MakeSessionStruct(sessionDir,varargin)
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
%
%           with fields:
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




p = inputParser;
addParameter(p,'makeField',[],@iscell);



parse(p,varargin{:});

makeField = p.Results.makeField;


if isempty(makeField)
    
    makeField = {'MergeDLCwNeural','MergeContextEntryandEdges','alignEdgPosition_pix2cm','MergeVelocityAcceleration',...
        'MergeTimeFromTransition','MergeDist2Edge'};
end

%% merge neural data and behavior data from Deeplabcut

if contains('MergeDLCwNeural',makeField)
    sm_MergeDLCwNeural(sessionDir);
end


%% redefine baseline

if contains('newBaseline',makeField)
    sm_newBaseline(sessionDir,[10 550]);
end


%% merge context entry times and arena edges
if contains('MergeContextEntryandEdges',makeField)
    fprintf('Merging context entry times and edges...\n\n')
    sm_MergeContextEntryandEdges(sessionDir)
end
%%

%% merge context entry times and arena edges for linear track
if any(contains(makeField,'MergeContextEntryandEdges_linear'))
    fprintf('Merging context entry times and edges...\n\n')
    sm_MergeContextEntryandEdges_linear(sessionDir)
end



%% align the areda edges, and convert to cm (also does a tranformation to get the skewed arena into a non-skewed form


if contains('alignEdgPosition_pix2cm',makeField)
    fprintf('Converting from pixels to cm...\n\n')
    sm_alignEdgPosition_pix2cm(sessionDir)
end


%% get velocity and acceleration in cm distance
if contains('MergeVelocityAcceleration',makeField)
    fprintf('Getting velocity and acceleration...\n\n')
    sm_MergeVelocityAcceleration(sessionDir)
end

%% get time from context entry
%
if contains('MergeTimeFromTransition',makeField)
    fprintf('Getting time from entry...\n\n')
    sm_MergeTimeFromTransition(sessionDir)
end


%% get how far the animal is from the arena edges in cm
if contains('MergeDist2Edge',makeField)
    fprintf('Finding distance from the arena edges...\n\n')
    sm_MergeDist2Edge(sessionDir)
end

%% get time from Object intro
if contains('MergeTimeFromObjectIntro',makeField)
    fprintf('Getting time from object intro..\n\n')
  sm_MergeTimeFromObjectIntro(sessionDir)
end


%% get time of reward
if contains('MergeRewardTime',makeField)
    fprintf('Get reward timing .\n\n')
  sm_MergeRewardTime(sessionDir)
end

%% get number of days that the animal has experienced the context (homecage is set to 100 days)
if contains('MergeDaysInContext',makeField)
    fprintf('Getting number of days in context...\n\n')
    topDir = fileparts(sessionDir);
    sm_MergeDaysInContext(topDir)
end


end