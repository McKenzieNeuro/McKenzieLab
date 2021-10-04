function recList = Process_IntanData_multi_start(fbasename,varargin)

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

if ~isempty(varargin)
    newfbasename = varargin{1};
    if length(varargin)>1
        cpVideo = varargin{2};
    end
else
    newfbasename = fbasename;
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
    
    newFbase_dat = [ newFbase '.dat'];
    
    fname = [folders(ii).folder filesep recName{ii} filesep 'amplifier.dat'];
    if exist(fname,'file')
        movefile(fname, newFbase_dat)
    else
        warning(['Dat file ' fname ' does not exist'])
    end
    fname = [folders(ii).folder filesep recName{ii} filesep 'digitalin.dat'];
    newFbase_dat = [ newFbase '_digitalin.dat'];
    if exist(fname,'file')
        movefile(fname, newFbase_dat)
    else
        warning(['digitalin file ' fname ' does not exist'])
    end
    
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
    
    
    
    fname = [folders(ii).folder filesep recName{ii} filesep 'analogin.dat'];
    newFbase_dat = [ newFbase '_analogin.dat'];
    
    
    if exist(fname,'file')
        movefile(fname, newFbase_dat)
    else
        warning(['Analogin ' fname ' does not exist. Try animal root folder'])
        
    end
    
    
    fname = [folders(ii).folder filesep recName{ii} filesep 'auxiliary.dat'];
    newFbase_dat = [ newFbase '_auxiliary.dat'];
    
    
    
    
    if exist(fname,'file')
        movefile(fname, newFbase_dat)
    else
        warning(['Auxiliary file ' fname ' does not exist. Try animal root folder'])
        
    end
    
    
    if ~strcmp(recName{ii},newFbase)
        movefile(recName{ii},newFbase)
    end
    
    % fname = fullfile(recName{ii},'digitalin.dat');
    %Convert Baesler video to whl file, warning when coding, we assume
    %fname is the digitalin file here
    %
    %     nbDigIn = 1;
    % camSyncCh = 1;
    %     ConvertBaesler2Whl(videoName,fname)
    %
    
    
end