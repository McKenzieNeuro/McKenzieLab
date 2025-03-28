
function [good,wf,shank] = sm_assign_noise(datfil,varargin)

%%







p = inputParser;


addParameter(p,'fs',30000,@isnumeric)
addParameter(p,'clufil',[],@isstr)
addParameter(p,'minRate',.1,@isnumeric)
addParameter(p,'noiseISI',.002,@isnumeric)
addParameter(p,'noiseThres',.005,@isnumeric)
addParameter(p,'wfCorrThres',.99,@isnumeric)
addParameter(p,'getWaveform',true,@islogical)


parse(p,varargin{:})


fs = p.Results.fs;
clufil = p.Results.clufil;
minRate = p.Results.minRate;
noiseISI = p.Results.noiseISI;
noiseThres = p.Results.noiseThres;
wfCorrThres = p.Results.wfCorrThres;
getWaveform = p.Results.getWaveform;
if isempty(clufil)
    
    [dirN] = fileparts(datfil);
    
    %find subdirectory
    fils= getAllExtFiles(dirN,'npy',1);
    
    kp = contains(fils,'spike_clusters') & contains(fils,'Kilosort_202');
    fils = fils(kp);
    
    
    if ~isempty(fils) && length(fils)==1
        clufil = fils{1};
    else
        error('must specify clufile')
    end
end

%%
[a,b] = fileparts(clufil);

if isempty(a)
    a = pwd;
end


% load data
tsfil = [a filesep 'spike_times.npy'];
ts = readNPY(tsfil);
ts = double(ts)/fs;



clu = readNPY(clufil);




if length(clu)~=length(ts)
    
    error('clu IDs do not match length of times')
end


%%
fid1 = fopen([a filesep 'cluster_deletion.tsv'],'wt');
fid = fopen([a filesep 'cluster_group.tsv'],'wt');
fprintf(fid, 'cluster_id	group\n');
fprintf(fid1, 'cluster_id	reason\n');
uclu = unique(clu);
nclu = length(uclu);
maxT = double(max(ts));

good = true(nclu,1);

ix = 1;
binSize = .00033;
conv_w = .5/binSize;  % 2ms window
alpha = .05;
% assign as noise due to ISI violations
for i = uclu'
    tst = ts(clu==i);
    cch = CrossCorr(tst,tst,binSize,101);
    
    
    pred = mean(cch(1:10));
    cch = cch(52:56);
    
    
    %if length(tst)/maxT<minRate ||  mean(diff(tst)<noiseISI)>noiseThres
    
    
    
    if any(cch>pred) || (length(tst)/maxT) < minRate
        fprintf(fid, [num2str(i) '	noise\n']);
        
        if any(cch>pred)
            
            fprintf(fid1, [num2str(i) '	noisy MUA\n']);
            
        else
            fprintf(fid1, [num2str(i) '	low Rate\n']);
        end
        good(ix) = false;
        isiVal(ix) =  mean(diff(tst)<noiseISI);
    end
    % end
    ix = ix+1;
end



if getWaveform
    [wf, shank,channelCorr] = sm_getWaveform(datfil,clu,ts,good);
    
    
    % assign as noise due to waveform abnormalities
    ix = 1;
    for i = uclu'
        
        if good(ix) && (channelCorr(ix) > wfCorrThres)
            fprintf(fid, [num2str(i) '	noise\n']);
            fprintf(fid1, [num2str(i) '	similar waveform\n']);
            
            good(ix) = false;
            i
        end
        ix = ix+1;
    end
    
    
end

fclose(fid);
fclose(fid1);
%%
end
