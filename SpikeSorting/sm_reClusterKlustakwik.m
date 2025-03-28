function sm_reClusterKlustakwik(datfil,varargin)
%this function takes output from spyking circus or kilosort and with units
%with more than thres ISI violations runs klustakwik
%
% input
% datafil = pull path of the data file. the samedirectory name that must have:
% basename.dat, basename.xml and
% subdirectory with output of phy




warning off
p = inputParser;


addParameter(p,'recluserThres',.01,@isnumeric) %method to correct for decaying signal
addParameter(p,'reclusterISI',.002,@isnumeric) % how much time (s) to skip before calculating mean/std/decay
addParameter(p,'NPYfil',[],@isstr)

parse(p,varargin{:})


recluserThres = p.Results.recluserThres;
reclusterISI = p.Results.reclusterISI;
NPYfil = p.Results.NPYfil;

%path for program must be in your environment path
program = 'klustakwik';


xml = LoadXml([datfil(1:end-3) 'xml']);
fs = xml.SampleRate;

[dirN] = fileparts(datfil);

%find subdirectory


if isempty(NPYfil)
    fils= getAllExtFiles(dirN,'npy',1);
    
    kp = contains(fils,'pc_features') & contains(fils,'Kilosort_2024');
    fils = fils(kp);
    
    
    
    if ~isempty(fils) && length(fils)==1
        
        NPYfil = fils{1};
        [NPYdir,] = fileparts(fils{1});
    else
        error('speficiy NPY file')
    end
    
else
    NPYdir = fileparts(NPYfil);
end



fet = readNPY(NPYfil);

clufil = [NPYdir filesep 'spike_clusters.npy'];

tempfil = [NPYdir filesep 'spike_templates.npy'];

copyfile(tempfil,clufil)

clu = readNPY(tempfil);

tsfil = [NPYdir filesep 'spike_times.npy'];
ts = readNPY(tsfil);
ts = double(ts)/fs;
uclu = unique(clu);
maxClu = max(clu);


for j = uclu'
    
    tst =  double(ts(clu==j));
    
    
    if mean(diff(tst)<reclusterISI)>recluserThres
        
        %recluster
        kp = any(squeeze(fet(clu==j,1,:))',2);
        
        fett = fet(clu==j,:,kp);
        fett = reshape(fett,size(fett,1),size(fett,2)*size(fett,3));
        
        %rescale for int64
        factor = 2^60;
        factor = int64(factor/max(abs(fett(:))));
        fet2 = int64(fett) * factor;
        
        
        
        fetname = fullfile(NPYdir, ['tmp.fet.1']);
        cluname = fullfile(NPYdir, ['tmp.clu.1']);
        
        
        
        % write temp feature (fet) file
        SaveFetIn(fetname,fet2);
        
        
        
        %recluster
        cmd = [program  ' ' fullfile(NPYdir,'tmp') ' 1'];
        cmd = [cmd ' -UseDistributional 0 -MaxPossibleClusters 20 -MinClusters 20 -UseFeatures ' sprintf('%d',ones(1,size(fet2,2)))];
        
        status = system(cmd);
        
        % relabel
        clut = load(cluname);
        clut = clut(2:end); % toss the first sample to match res/spk files
        clut = int32(maxClu)+int32(clut);
        
        clu(clu==j) = clut;
        writeNPY((clu), clufil)
        maxClu = max(clut);
    end
    
    
end


% exclude noisy units
[good,wf,shank] = sm_assign_noise(datfil);
sm_MergeCluster(datfil,'good',good,'shank',shank,'wf',wf);
%

end


