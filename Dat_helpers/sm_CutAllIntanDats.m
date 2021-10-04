function sm_CutAllIntanDats(basepath,varargin)


% cuts intan files at start/stop location AND DELETES ORIGINAL FILES. HAVE
% A BACK UP

% must have only one amplifier/analogin/auxiliary/digitalin file per
% directory

% should enter a desired duration and start, e.g.

%  CutAllIntanDats(pwd,'start',10,'duration',1600)
%  Will make a new set of dats with the first 10-1610s

%  CutAllIntanDats(pwd,duration',1600)
%  Will make a new set of dats with the first 1600s




datfil = getAllExtFiles(basepath,'dat',0);

% get all files
digitalfils = cellfun(@any,regexp(datfil,'digitalin'));
analogfils = cellfun(@any,regexp(datfil,'analogin'));
auxfils = cellfun(@any,regexp(datfil,'auxiliary'));
datafil =  cellfun(@any,regexp(datfil,'amplifier'));
timefil =  cellfun(@any,regexp(datfil,'time'));
%check to be sure XMLs are made

if ~exist([basepath filesep 'amplifier.xml'],'file')
    h = msgbox('Enter XML info for amplifier file');
    uiwait(h)
    cmd = [' neuroscope ' basepath filesep 'amplifier.dat'];
    system(cmd)
end

if ~exist([basepath filesep 'auxiliary.xml'],'file')
    h = msgbox('Enter XML info for auxiliary file (usually 3 channels)');
    uiwait(h)
    cmd = [' neuroscope ' basepath filesep 'auxiliary.dat'];
    system(cmd)
end

if ~exist([basepath filesep 'digitalin.xml'],'file')
    h=  msgbox('Enter XML info for digitalin file (1 channel)');
    uiwait(h)
    cmd = [' neuroscope ' basepath filesep 'digitalin.dat'];
    system(cmd)
end

if ~exist([basepath filesep 'analogin.xml'],'file')
    h =  msgbox('Enter XML info for analogin file');
    uiwait(h)
    cmd = [' neuroscope ' basepath filesep 'analogin.dat'];
    system(cmd)
end

if ~exist([basepath filesep 'time.xml'],'file')
    h =  msgbox('Enter XML info for time file, two channels');
    uiwait(h)
    cmd = [' neuroscope ' basepath filesep 'time.dat'];
    system(cmd)
end



datafil = datfil(datafil);
auxfils = datfil(auxfils);

analogfils = datfil(analogfils);


outfil_dat = [datafil{1}(1:end-4) '_cut.dat'];
analogfils_dat =  [analogfils{1}(1:end-4) '_cut.dat'];
auxfils_dat =  [auxfils{1}(1:end-4) '_cut.dat'];






CopyDat(datafil{1},outfil_dat,varargin{:})
CopyDat(analogfils{1},analogfils_dat,varargin{:})
CopyDat(auxfils{1},auxfils_dat,varargin{:})


delete datafil{1} analogfils{1} auxfils{1} digfils{1} timefils{1}


movefile(outfil_dat,datafil{1})
movefile(analogfils_dat,analogfils{1})
movefile(auxfils_dat ,auxfils{1})



timefils = datfil(timefil);
digfils = datfil(digitalfils);
digitalfils_dat =  [digfils{1}(1:end-4) '_cut.dat'];
timefils_dat =  [timefils{1}(1:end-4) '_cut.dat'];


CopyDat(digfils{1},digitalfils_dat,varargin{:})
CopyDat(timefils{1},timefils_dat,varargin{:})


delete digfils{1} timefils{1}


movefile(digitalfils_dat, digfils{1})
movefile(timefils_dat,timefils{1})

end