function [estimateLabel,trueLabel,inTrainingSet,time2seizure,seizure_start] = sm_getSeizurePred_gross(fname,seizFil,rusTree,trainingTime,ops)

% this function takes the classifier in rusTree and applease the feature
% space specified in ops and classified ever moment in time for the file
% (fname)
%%

Fs =ops.Fs;
bins = ops.bins;
nCh_featureFil = ops.nCh_featureFile;

%%



%these are the subject IDs for all the animals on file
%subject_IDs = {'KA11'};  %for testing, only read in one animals data at a time


sz_files = true;                   %whether to load only the files with seizures (true) or all files (false)



cd(ops.RawDataPath)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add paths to necessary folders

% Here's a url to the CED website where you can download the file (after
% giving them your name/email): https://ced.co.uk/upgrades/spike2matson
% I also included the files on the hard drive (CEDMATLAB folder), but you
% may have to download them separately to get them to work - I followed the
% instructions on the s64mat file and didn't have much issue. Let me know
% if you have trouble.

%load CED MATSON library for reading in .smrx files
cedpath = getenv('R:\Analysis\CEDMATLAB\CEDS64ML');        %this part may be where you need to edit;
cedpath = 'R:\Analysis\CEDMATLAB\CEDS64ML';
%you might be able to switch to the manual path location rather than use
%the 'getenv' function; see pg. 1-6 of the s64mat pdf for more info.
addpath(cedpath)
CEDS64LoadLib(cedpath)

%other libraries - again, may have to adjust when you start
%addpath(genpath('D:\file_for_Sam\'))
%%

dirN = strrep(seizFil,fileparts(fileparts(seizFil)),'');
dirN = dirN(2:end);

 [~,subject_IDs] = fileparts(fileparts(seizFil));

    ix = 1;
    
    [rows_to_extract, stim_metatable] = read_metatable_KAsz(subject_IDs, sz_files);
    
    sz_o = cellfun(@str2num,stim_metatable.Seizures_sec_(rows_to_extract),'UniformOutput',false);
    sz_off = cellfun(@str2num,stim_metatable.SeizuresEnd_sec_(rows_to_extract),'UniformOutput',false);
    ID = [stim_metatable(rows_to_extract,:).AnimalIdentification stim_metatable(rows_to_extract,:).TrialNumber];
    for k2 = 1:size(ID,1)
        
      
        
        str{ix} = [ID{k2,1} filesep ID{k2,1} '_' ID{k2,2}];

        sz{ix} = [sz_o{k2} sz_off{k2}];
        ix   = ix+1;
    end
    
    kp = contains(str,dirN);
    
sz = sz{kp};

if ~isempty(sz)
seizure_start = sz(:,1);
seizure_end = sz(:,2);
end


%%

% get true time to seizure
[~,basename ] = fileparts(fname);
fname = [fname filesep basename];
powerFil = [fname '_1.dat'];
s = dir(powerFil);
disp(powerFil)
disp(s.bytes);
disp(nCh_featureFil);
disp(Fs)
dur = s.bytes/nCh_featureFil/Fs/2;



ts = 0: (dur-ops.durFeat);
time2seizure = ts;
kp = true(size(ts));
for i = 1:length(seizure_start)
    time2seizure(kp&ts<seizure_start(i)) = ts(kp&ts<seizure_start(i)) - seizure_start(i);
    time2seizure(kp& ts>seizure_start(i) & ts<seizure_end(i)) = .5;
    time2seizure(kp& ts>seizure_end(i) & ts<seizure_end(i)+600) = 1.5;
    
    
    kp(ts<seizure_start(i) | ...
        (ts>seizure_start(i) & ts<seizure_end(i)) | ...
        (ts>seizure_end(i) & ts<seizure_end(i)+600)) = false;
end
if ~any(isinf(bins))
[~,trueLabel] = histc(time2seizure,[-inf -(bins) 0 1 2]);
else
    [~,trueLabel] = histc(time2seizure,[ -(bins) 0 1 2]);
end

if isempty(trainingTime)
    inTrainingSet = false(size(ts));
else
    inTrainingSet = histc(trainingTime,ts)>0;
end
%%


%get prediction
estimateLabel =[];
dat1 =[];
for i = ts
    
    
    
    
    tim = i;
    features = ops.features(fname,tim,ops);
    
    
    dat1 = [dat1;features];
    
    if mod(i,100)== 0
        [outpred,conf] = predict(rusTree,dat1);
        estimateLabel = [estimateLabel;outpred conf];
        dat1 =[];
    elseif i > (dur- mod(dur,100))
       [outpred,conf] = predict(rusTree,dat1);
        estimateLabel = [estimateLabel;outpred conf];
         dat1 =[];
    end
    
end




end





