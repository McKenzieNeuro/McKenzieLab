%calculate feature files to save to disk as *.dat files for use later
sm_MakeAll_getPowerPerChannel_Gross

%calculate features for the subset of relevant times and save as *.mat for
%more immediate use
sm_PredictIHKA_getAllFeatures

%train model and quantify performance
sm_PredictIHKA

%run model for all times for all files
sm_MakeAllSeizurePred



