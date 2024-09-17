%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Goal: Run SPaRCNet to get IIIC scores for each 10sec window                                                         
% Todo: Edit line#10 pyExec to local conda env (iiic $where pyhton)     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear

mkdir('R:\Analysis\McKenzieLab\detectors\seizure\nedc_eeg_smile\edf\')

pyExec = 'C:\Users\samckenzie\Anaconda3\envs\iiic\python.exe';

pyRoot = fileparts(pyExec);
p = getenv('PATH');
p = strsplit(p, ';');
addToPath = {
    pyRoot
    fullfile(pyRoot, 'Library', 'mingw-w64', 'bin')
    fullfile(pyRoot, 'Library', 'usr', 'bin')
    fullfile(pyRoot, 'Library', 'bin')
    fullfile(pyRoot, 'Scripts')
    fullfile(pyRoot, 'bin')};
p = [addToPath(:); p(:)];
p = unique(p, 'stable');
p = strjoin(p, ';');
setenv('PATH', p);
 
cmd1 = ['activate iiic';system(cmd1);
cmd2 = 'python R:\\Analysis\\McKenzieLab\\detectors\\seizure\\A_RapidIIIC_pipeline_B\\tools\\runSPaRCNet.py';system(cmd2);
cmd3 = 'conda deactivate';system(cmd3);

