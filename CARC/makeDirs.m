% first collect all directory names and save in text file
% then move that file to CARC
% then fun xargs mkdir -p < "/carc/scratch/projects/mckenzie2016183/code/BASH/directories.txt"




topDir = 'R:\STDP';
fils = getAllExtFiles(topDir,'dat',1);
fil = fils(contains(fils,'_int16.dat'));

[a,b] = fileparts(fils);

fid = fopen('R:\Analysis\McKenzieLab\CARC\directories.txt','wt');


CARC_DIR = '/carc/scratch/projects/mckenzie2016183/data/spikeSorting/';
for j = 1:length(fil)
    
    [a1,b] = fileparts(fil{j});
    tmp_remote = a1(9:end);
    tmp_remote = strrep(tmp_remote,'\','/');
    outDir = [CARC_DIR tmp_remote];
    fprintf(fid, [outDir '\n']);
end

fclose(fid);


%%


topDir = 'R:\STDP';
fils = getAllExtFiles(topDir,'dat',1);
fil = fils(contains(fils,'_int16.dat'));




CARC_DIR = '/carc/scratch/projects/mckenzie2016183/data/spikeSorting/';
%%
fid = fopen('R:\Analysis\McKenzieLab\CARC\copyFiles.txt','wt');
for j = 1:length(fil)
    
    [a1,b] = fileparts(fil{j});
    tmp_remote = a1(9:end);
    tmp_remote = strrep(tmp_remote,'\','/');
   
    
    d = regexp(tmp_remote,'/');
    
    if length(d)>1
        tmp_remote = tmp_remote(1:d(2));
    else
        tmp_remote = [tmp_remote '/'];
    end
    
     outDir = [CARC_DIR tmp_remote];
    tmp_local = strrep(fil{j},'\','\\');
    command = ['scp ' tmp_local ' mckenzie@hopper.alliance.unm.edu:' outDir b c];
    fprintf(fid, [command '\n']);
end

fclose(fid);


