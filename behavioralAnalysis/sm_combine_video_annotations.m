function sm_combine_video_annotations(dirN,keyword)
% SM_COMBINE_VIDEO_ANNOTATIONS
% Combines event annotations from multiple synchronized videos into
% a single file containing Intan timestamps.
%
% INPUTS:
%   dirN    - directory containing video annotation and sync files
%   keyword - string keyword to match specific annotation files
%
% REQUIREMENTS:
%   - You must first run `clickVideo1` and `syncVideo2` to generate the
%     synchronization file `ts_syn.mat`
%   - Annotation files must be named like: [videoBaseName '_event.xls']
%
% OUTPUT:
%   Saves a MAT file containing a table `all_events` with columns:
%       event, Intan_ts, video_frame, video
%



cd(dirN)

fout = [dirN filesep 'annotations_' keyword '.mat'];
%load sync file
if exist('ts_syn.mat')
    
    
    v = load('ts_syn.mat');
else
    error('must sync videos')
end


% get all annotation (must have videoBaseName_event.xls)
fils = getAllExtFiles(dirN,'xls',1);
kp = contains(fils,['_' keyword]);
fils = fils(kp);


if ~isfield(v,'fils')
    
    ts_dat =[];
    for i = 1:length(fils)
        %get videoname
        video_file = strrep(fils{i},['_' keyword '.xls'],'.mp4');
        
        tmp = dir(video_file);
        
        ts_dat(i)= datenum(tmp.date);
    end
    
    [~,b] = sort(ts_dat);
    v.fils = v.fils(b);
else
    [~,tmp] = fileparts(v.fils);
    
    if ~iscell(tmp)
        basename{1} = tmp;
    else
        basename = tmp;
    end
    for i = 1:length(basename)
        idx(i) = find(cellfun(@any,regexp(v.fils,basename{i})));
    end
    
    v.fils = v.fils(idx);
end


if any(contains(v.outTable.Properties.VariableNames,'VideoNumber'))
    numVids = max(v.outTable.VideoNumber);
else
    numVids = 1;
    v.outTable.VideoNumber = ones(length(v.outTable.IntanTime),1);
    
end
all_events =[];
for  i = 1:length(fils)
    
    %get frame number
    tab = readtable(fils{i});
    
    Intan_ts = interp1(v.outTable.VideoFrameNumber(v.outTable.VideoNumber==i),v.outTable.IntanTime((v.outTable.VideoNumber==i)),tab.Var2,'linear','extrap');
    event = tab.Var1;
    video_frame = tab.Var2;
    video = num2cell(repmat(fils{i},length(Intan_ts),1),2);
    all_events = [all_events;table(event,Intan_ts,video_frame,video)];
end

save(fout,'all_events')
end