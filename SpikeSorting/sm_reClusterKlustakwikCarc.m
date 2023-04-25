function sm_reClusterKlustakwikCarc(datfil)
% this function takes output from spyking circus or kilosort and with units
% with more than thres ISI violations runs klustakwik



% edit here for your path
%program = 'C:\Users\alsommer\Documents\MATLAB\McKenzieLab\SpikeSorting\phy2Plugins\klustakwik ';
program = '/carc/scratch/projects/mckenzie2016183/code/software/KlustaKwik-1.5/KlustaKwik';

% free params
recluserThres = .01;
reclusterISI = .002;


dirN = fileparts(datfil);
% find subdirectory
fils= getAllExtFiles(dirN,'npy',1);

kp = contains(fils,'pc_features') & ~contains(fils,'Kilosort');
fils = fils(kp);


if ~isempty(fils) && length(fils)==1
    [a,b] = fileparts(fils{1});
    fet = readNPY(fils{1});
    
    clufil = [a filesep 'spike_clusters.npy'];
    
    tempfil = [a filesep 'spike_templates.npy'];
    
    copyfile(tempfil,clufil)
    
    clu = readNPY(tempfil);
    
    tsfil = [a filesep 'spike_times.npy'];
    ts = readNPY(tsfil);
    ts = ts/30000;
    uclu = unique(clu);
    maxClu = max(clu);
    for j = uclu'
        
        tst =  double(ts(clu==j));
        
        
        if mean(diff(tst)<reclusterISI)>recluserThres
            
            % recluster
            kp = any(squeeze(fet(clu==j,1,:))',2);
            
            fett = fet(clu==j,:,kp);
            fett = reshape(fett,size(fett,1),size(fett,2)*size(fett,3));
            
            %rescale for int64
            factor = 2^60;
            factor = int64(factor/max(abs(fett(:))));
            fet2 = int64(fett) * factor;
            
            
            
            fetname = fullfile(a, ['tmp.fet.1']);
            cluname = fullfile(a, ['tmp.clu.1']);
            
            
            
            %write temp feature (fet) file
            SaveFetIn(fetname,fet2);
            
            
            
            %recluster
            cmd = [program  ' ' fullfile(a,'tmp') ' 1'];
            cmd = [cmd ' -UseDistributional 0 -MaxClusters 20 -MaxPossibleClusters 20 -MinClusters 20 -UseFeatures ' sprintf('%d',ones(1,size(fet2,2))) ...
                ' -MaxIter 10000'];
       
            status = system(cmd);
            
            %relabel
            clut = load(cluname);
            clut = clut(2:end); % toss the first sample to match res/spk files
            clut = int32(maxClu)+int32(clut);
            
            clu(clu==j) = clut;
            writeNPY((clu), clufil)
            maxClu = max(clut);
        end
        
        
    end
    
    
% exclude noisy units
%[good,wf,shank] = sm_assign_noise(datfil);
%sm_MergeCluster(datfil,'good',good,'shank',shank,'wf',wf);
end




end


