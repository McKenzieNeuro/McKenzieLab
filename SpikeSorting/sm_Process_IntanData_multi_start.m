function recList = sm_Process_IntanData_multi_start(fbasename,varargin)

% recList = Process_IntanData(fbasename)
% Process Intan data folders starting with 'fbasename'
% Returns a cell array with thee names of the new folders
%
% recList = Process_IntanData(fbasename,newfbasename)
% Change the file basename to 'newfbasename'


eraseDir = 0;
cpVideo = 1;
%nbDigIn = 1;
%camSyncCh = 1;



% parse args
p = inputParser;
addParameter(p,'cpVideo',[],@isstr)
addParameter(p,'intan_version','Combined',@isstr)
addParameter(p,'newfbasename',[],@isstr)
addParameter(p,'dirName',[],@isstr)



parse(p,varargin{:})
cpVideo = p.Results.cpVideo;
intan_version = p.Results.intan_version;
newfbasename = p.Results.newfbasename;
dirName = p.Results.dirName;

if isempty(dirName)
    dirName = pwd;
end


if isempty(newfbasename)
    
    newfbasename = fbasename;
end



if strcmpi(intan_version,'OpenEphys')
   %reorganize directory
   
   %change name of mast directory
  [a,b]  =fileparts(dirName); 
  
    newDirName = [a filesep fbasename];
    cd ..
    
    if ~strcmp(dirName,newDirName)
        movefile(dirName,newDirName);
    end
    
    cd(newDirName)
    %get all sub-directories and rename
    
    ok = dir(newDirName);
    kp = cell2mat({ok.isdir});
    ok = ok(kp);
    dirs =  {ok.name}';
    kp = cellfun(@any,regexp({ok.name}','Record Node'));
    ok = ok(kp);
    
    %rename all subdirectories
    
    sub = [newDirName filesep ok.name];
    
    ok1 = dir(sub);
     
    kp = cell2mat({ok1.isdir});
    ok1 = ok1(kp);
    dirs =  {ok1.name}';
    kp = cellfun(@any,regexp({ok1.name}','experiment'));
    ok1 = ok1(kp);
    
    
    %rename all subdirectories
    
    for ii = 1:length(ok1)
        sub_sub = [sub filesep ok1(ii).name];
  
        %move continuous data into sub directory
        fname = [sub_sub filesep 'recording1' filesep 'continuous' filesep 'Rhythm_FPGA-100.0' filesep 'continuous.dat'];
        new_file =  [newDirName filesep fbasename '_' num2str(ii) filesep 'amplifier_analogin_auxiliary_int16.dat'];
        if ~exist( [newDirName filesep fbasename '_' num2str(ii)])
            mkdir([newDirName filesep fbasename '_' num2str(ii)])
        end
        if exist(fname)
            movefile(fname,new_file)
        end
    end
    
    
end
folders = dir([fbasename '_*']);
token = '_';

if isempty(folders)
    folders = dir([fbasename '-*']);
    token = '-';
end
date = [];
startTime = [];
recName = {};
nRec = length(folders);
for ii=1:nRec
    
    if folders(ii).isdir;
        fname = folders(ii).name;
        recName = [recName;{fname}];
        k = strfind(fname,token);
        if length(k)>1
            date = [date;str2num(fname(k(end-1)+1:k(end)-1))];
        else
            date = [date;str2num(fname(1:k(end)-1))];
        end
        startTime = [startTime;str2num(fname(k(end)+1:end))];
    end
    
end

[startTime,ix] = sort(startTime);
recName = recName(ix);
nRec = length(recName);
recList = cell(nRec,1);






for ii=1:nRec
    
    nber = num2str(ii);
    if ii<10
        nber = ['0' nber];
    end
    
    newFbase = [newfbasename '-' nber];
    recList{ii} = newFbase;
    
    
    if strcmpi(intan_version,'combined') || strcmpi(intan_version,'openephys')
        
        %% combined file
        
        newFbase_dat = [ newFbase '.dat'];
        
        fname = [folders(ii).folder filesep recName{ii} filesep 'amplifier_analogin_auxiliary_int16.dat'];
        if exist(fname,'file')
            movefile(fname, newFbase_dat)
        else
            warning(['Dat file ' fname ' does not exist'])
        end
        
        
        %% amplifier xml
        fname = [folders(ii).folder filesep recName{ii} filesep 'amplifier_analogin_auxiliary_int16.xml'];
        newFbase_dat = [ newFbase '.xml'];
        if exist(fname,'file')
            movefile(fname, newFbase_dat)
        else
            warning(['XML file ' fname ' does not exist. Try animal root folder'])
            fname = fullfile('..','amplifier_analogin_auxiliary_int16.xml');
            if exist(fname,'file')
                movefile(fname, newFbase_dat)
            else
                warning(['XML file ' fname ' does not exist. Skipping it'])
            end
        end
        
        
        
        
    else
        
        
        %% amplifier
        
        newFbase_dat = [ newFbase '.dat'];
        
        fname = [folders(ii).folder filesep recName{ii} filesep 'amplifier.dat'];
        if exist(fname,'file')
            movefile(fname, newFbase_dat)
        else
            warning(['Dat file ' fname ' does not exist'])
        end
        
        %% amplifier xml
        fname = [folders(ii).folder filesep recName{ii} filesep 'amplifier.xml'];
        newFbase_dat = [ newFbase '.xml'];
        if exist(fname,'file')
            movefile(fname, newFbase_dat)
        else
            warning(['XML file ' fname ' does not exist. Try animal root folder'])
            fname = fullfile('..','amplifier.xml');
            if exist(fname,'file')
                movefile(fname, newFbase_dat)
            else
                warning(['XML file ' fname ' does not exist. Skipping it'])
            end
        end
        
        %%  analogin
        
        fname = [folders(ii).folder filesep recName{ii} filesep 'analogin.dat'];
        newFbase_dat = [ newFbase '_analogin.dat'];
        
        
        if exist(fname,'file')
            movefile(fname, newFbase_dat)
        else
            warning(['Analogin ' fname ' does not exist. Try animal root folder'])
            
        end
        %% auxiliary
        
        fname = [folders(ii).folder filesep recName{ii} filesep 'auxiliary.dat'];
        newFbase_dat = [ newFbase '_auxiliary.dat'];
        
        
        
        
        if exist(fname,'file')
            movefile(fname, newFbase_dat)
        else
            warning(['Auxiliary file ' fname ' does not exist. Try animal root folder'])
            
        end
        
        
    end
    
    
    %move digital file
    fname = [folders(ii).folder filesep recName{ii} filesep 'digitalin.dat'];
    newFbase_dat = [ newFbase '_digitalin.dat'];
    if exist(fname,'file')
        movefile(fname, newFbase_dat)
    else
        warning(['digitalin file ' fname ' does not exist'])
    end
    
    
    
    
    % rename folder
    if ~strcmp(recName{ii},newFbase)
        movefile(recName{ii},newFbase)
    end
    
    
end

end
