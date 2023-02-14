function sm_reClusterKlustakwikCarc(dirN)

% find subdirectory

fils= getAllExtFiles(dirN,'npy',1);


kp = contains(fils,'pc_features');
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
    ts = double(ts)/30000;
    uclu = unique(clu);
    maxClu = max(clu);
    for j = uclu'
        
        tst =  double(ts(clu==j));
        
        
        if mean(diff(tst)<.002)>.01
            
            % recluster
            kp = any(squeeze(fet(clu==j,1,:))',2);
            
            fett = fet(clu==j,:,kp);
            fett = reshape(fett,size(fett,1),size(fett,2)*size(fett,3));
            
            factor = 2^60;
            factor = int64(factor/max(abs(fett(:))));
            fet2 = int64(fett) * factor;
            
            
            
            fetname = fullfile(a, ['tmp.fet.1']);
            cluname = fullfile(a, ['tmp.clu.1']);
            
            
            
            
            SaveFetIn(fetname,fet2);
            
            
            
            program = '/carc/scratch/projects/mckenzie2016183/code/software/KlustaKwik-1.5/KlustaKwik ';
            
            cmd = [program  fullfile(a,'tmp') ' 1'];
            cmd = [cmd ' -UseDistributional 0 -MaxPossibleClusters 20 -MinClusters 20'];
            
            status = system(cmd);
            
            clut = load(cluname);
            clut = clut(2:end); % toss the first sample to match res/spk files
            clut = int32(maxClu)+int32(clut);
            
            clu(clu==j) = clut;
            writeNPY((clu), clufil)
            maxClu = max(clut)
        end
        
        
    end
    
    
    
    assign_noise(clufil,ts)
  %  getAllZeroLag(clufil,ts)
end




end

function assign_noise(clufil,ts)
%%
[a,b] = fileparts(clufil);
fid = fopen([a filesep 'cluster_group.tsv'],'wt');
fprintf(fid, 'cluster_id	group\n');
clu = readNPY(clufil);


uclu = unique(clu);

maxT = double(max(ts));
ix = 1;
for i = uclu'
    tst = ts(clu==i);
    
    if length(tst)/maxT<.1 |  mean(diff(tst)<.001)>.01
        
        fprintf(fid, [num2str(i) '	noise\n']);
        
    end
    
    
end
fclose(fid);
%%
end



function getAllZeroLag(clufil,ts)

[a,b] = fileparts(clufil);

clu = readNPY(clufil);


uclu = unique(clu);

ix = 1;

CC = nan(length(uclu),length(uclu));
for i = 1:length(uclu)
    tsti = ts(clu==uclu(i));
    if ~(length(tsti)/maxT<.1 |  mean(diff(tsti)<.001)>.005)
        for j = i+1:length(uclu)
            tstj = ts(clu==uclu(j));
            if ~(length(tstj)/maxT<.1 |  mean(diff(tstj)<.001)>.005)
                tmp = CrossCorr(tsti,tstj,.001,3);
                CC(i,j) = tmp(2)/length(tsti);
            end
        end
    end
end
end

