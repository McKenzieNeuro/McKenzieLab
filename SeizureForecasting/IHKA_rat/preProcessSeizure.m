topDir = 'R:\DGregg\NeuralData\PCP\Recordings';
fils = getAllExtFiles(topDir,'evt',1);
kp = contains(fils,'BayesOpt');
fils = fils(~kp);
%%
for  i =1:length(fils)
    dirN = fileparts(fils{i});
    cd(dirN)
    
    if ~exist('amplifier.xml')
        
        copyfile('R:\DGregg\NeuralData\PCP\Recordings\amplifier.xml',[dirN filesep 'amplifier.xml'])
    end
    ev = LoadEvents(fils{i});
    
    if ~isempty(ev.time) && ~exist('amplifier.lfp')
        
       
        bz_LFPfromDat(pwd,'basename','amplifier')
        
    end
    
end