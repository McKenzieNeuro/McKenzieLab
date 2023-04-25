
function sm_MergeCluster(datfil,varargin)




p = inputParser;


addParameter(p,'fs',30000,@isnumeric)
addParameter(p,'clufil',[],@isstr)
addParameter(p,'minRate',.1,@isnumeric)
addParameter(p,'noiseISI',.002,@isnumeric)
addParameter(p,'noiseThres',.005,@isnumeric)
addParameter(p,'wfCorrThres',.5,@isnumeric)
addParameter(p,'wf',[],@iscell)
addParameter(p,'good',[],@islogical)
addParameter(p,'shank',[],@isnumeric)

parse(p,varargin{:})


fs = p.Results.fs;
clufil = p.Results.clufil;
minRate = p.Results.minRate;
noiseISI = p.Results.noiseISI;
noiseThres = p.Results.noiseThres;
wfCorrThres = p.Results.wfCorrThres;
wf = p.Results.wf;
good = p.Results.good;
shank = p.Results.shank;




if isempty(clufil)
    
    [dirN] = fileparts(datfil);
    
    %find subdirectory
    fils= getAllExtFiles(dirN,'npy',1);
    
    kp = contains(fils,'spike_clusters') & contains(fils,'Kilosort_2023');
    fils = fils(kp);
    
    
    if ~isempty(fils) && length(fils)==1
        clufil = fils{1};
    else
        error('must specify clufile')
    end
end


[phyDir,b] = fileparts(clufil);

if isempty(phyDir)
    phyDir = pwd;
end


% load data
tsfil = [phyDir filesep 'spike_times.npy'];
ts = readNPY(tsfil);
ts = double(ts)/fs;


copyfile(clufil,[clufil '.postKlust_premerge'])
clu = readNPY(clufil);




maxT = max(ts);
uclu = unique(clu);
nclu = length(uclu);
groupfil = fullfile(phyDir,'cluster_group.tsv');

if isempty(good)
    
    
    
    
    
    if exist(groupfil)
        cluster_group = tdfread(groupfil);
        badID = cluster_group.cluster_id(all(cluster_group.group =='noise',2));
         good = true(nclu,1);
         good(ismember(uclu,badID)) = false;
    else
        
        good = true(nclu,1);
    end
end

if isempty(wf)
    [wf, shank,channelCorr] = sm_getWaveform(datfil,clu,ts,good);
end




binSize = .0005; %.5ms
duration = .2; %200ms
conv_w = .002/binSize;  % 2ms window
alpha = 0.001;


%%


CC = nan(length(uclu),length(uclu));


for i = 1:length(uclu)
    
    %initialize mrg, indexed by uclu
    mrg{i} =[];
    if good(i)
        
        tsti = ts(clu==uclu(i));
        
        for j = i+1:length(uclu)
            if good(j) && shank(i)==shank(j)
                
                
                
                tstj = ts(clu==uclu(j));
                
                if length(tsti) > length(tstj)
                    ref = tsti;
                    target =tstj;
                else
                    ref = tstj;
                    target = tsti;
                    
                end
                
                
                cch = CrossCorr(ref,target,binSize,ceil(duration/binSize)+1);
                
                
                
                
                
                %delete other parts of the refractory period
                
                cch([200 202]) = [];
                
                
                [pvals,pred,qvals]=bz_cch_conv(cch,conv_w);
                
                hiBound=poissinv(1-alpha,pred);
                
                % collect some data
                
                confus(1) = cch(200)/hiBound(200);
                confus(2) =corr(wf{i}(:),wf{j}(:));
                
                
                
                %find spikes belonging to same unit
                if confus(1)>2 && confus(2)>.75
                    % check if merging would hurt the refractory period of
                    % high rate
                    
                    mrg{i} = [mrg{i} j];
                    
                    
                    
                end
                %                 if cch(200) > hiBound(200)
                %
                %                     subplot(1,3,1)
                %                     plot(cch)
                %                     subplot(1,3,2)
                %                     imagesc(wf{i}')
                %                     subplot(1,3,3)
                %                     imagesc(wf{j}')
                %                     waitforbuttonpress
                %                     close all
                %                 end
                %handle empty refractory period (ISI violation <
                
                % CC(i,j) = tmp(2)/length(tsti);
            end
        end
        
        
        
        
    end
end
%%

fid = fopen([phyDir filesep 'cluster_group.tsv'],'a+');
fid1 = fopen([phyDir filesep 'cluster_deletion.tsv'],'a+');
% now find all confusions and merge in order of rate

rate = [];
for j = 1:nclu
    
    rate(j) = sum(clu==uclu(j));
end

[~,b] = sort(rate,'descend');


mrg1 = mrg(b);
uclu1 = uclu(b);


clu1 = clu;
for i=1:nclu
    idx = b(i);
    if ~isempty(mrg1{i})
        
        prs = mrg1{i};
        otherprs = b(find(cellfun(@(a) any(ismember(a,prs)),mrg1)));
        otherprs = [otherprs cell2mat(mrg1(cellfun(@(a) any(ismember(a,prs)),mrg1)))];
        prs = unique([idx prs otherprs]);
        
        % check if all on sahnk shank
        
        sh = shank(prs);
        
        kp = sh == shank(idx);
        
        prs = prs(kp);
        
        % merge in order of highest rate
        [~,b1] = sort(rate(prs),'descend');
        
        prs = prs(b1);
        
        
        
        ref = uclu(prs(1));
        targets = prs(2:end);
        
        
        
        for j = 1:length(targets)
            tsti = ts(clu1==ref);
            tstj = ts(clu1==uclu(targets(j)));
            
            % merge spikes
            
         
     
            cch1 = CrossCorr(tsti,tsti,.0005,101);
            pred1 = mean(cch1(1:10));
            
           
            
            cch2 = CrossCorr(tstj,tstj,.0005,101);
            pred2 = mean(cch2(1:10));
            
            
            cch = CrossCorr(tsti,tstj,.0005,101);
            cch = cch(52:56);
            
            
            if ~(any(cch>pred1) | any(cch>pred2)) && corr(wf{prs(1)}(:),wf{targets(j)}(:))>.75
                
                % reassign good spikes that one have isi violations
                
                clu1(clu1==uclu(targets(j))) = ref;
                ixx = find(clu1==ref);
                %exclude ISI violations and set those spikes to noise
                
                ids = [ones(length(tsti),1);2*ones(length(tstj),1)];
                [mergets,bb]  = sort([tsti;tstj]);
                ids = ids(bb);
                
                
                %delete same spike 
                kp1 = [true;ids(1:end-1)~=ids(2:end)];
                kp = [true;diff([mergets])<.001] & kp1;
                
                
                clu1(ixx(kp)) = uclu(targets(j));
                fprintf(fid, [num2str(uclu(targets(j))) '	noise\n']);
                fprintf(fid1, [num2str(uclu(targets(j))) '	merged with ' num2str(ref) '\n']);
                %delete other merges
                mrg1{uclu1==uclu(targets(j))} = [];
                good(targets(j)) = false;
            else
                
                % if the lower rate unit has rate lower than thres2,
                % delete
                if rate(targets(j)) < .5
                    fprintf(fid, [num2str(uclu(targets(j))) '	noise\n']);
                    fprintf(fid1, [num2str(uclu(targets(j))) '	similar to ' num2str(ref) '\n']);
                    %delete other merges
                    mrg1{uclu1==uclu(targets(j))} = [];
                    good(targets(j)) = false;
                end
            end
        end
    end
end
%%
% 
% % merge all units with clean ISIs
% 
% 
% 
% 
% for i = 1:length(uclu)
%     mrg{i} =[];
%     if good(i)
%         
%         tsti = ts(clu1==uclu(i));
%         
%         for j = i+1:length(uclu)
%             
%          
%             if good(j) && shank(i)==shank(j)
%                 
%                 
%                 
%                 tstj = ts(clu1==uclu(j));
%                 
%                 if length(tsti) > length(tstj)
%                     ref = tsti;
%                     target =tstj;
%                 else
%                     ref = tstj;
%                     target = tsti;
%                     
%                 end
%                 
%                 
%                 cch1 = CrossCorr(tsti,tsti,.0005,101);
%                 pred1 = mean(cch1(1:10));
%                 
%             
%                 
%                 cch2 = CrossCorr(tstj,tstj,.0005,101);
%                 pred2 = mean(cch2(1:10));
%                 
%                 
%                 cch = CrossCorr(tsti,tstj,.0005,101);
%                 cch = cch(52:56);
%                 
%             
%                 
%                 
%                 
%                 
%                 
%                 if corr(wf{i}(:),wf{j}(:))>.75 &&  ~(any(cch>pred1) | any(cch>pred2))
%                     % check if merging would hurt the refractory period of
%                     % high rate
%                     
%                     mrg{i} = [mrg{i} j];
%                     
%                     
%                     
%                 end
%                 %                 if cch(200) > hiBound(200)
%                 %
%                 %                     subplot(1,3,1)
%                 %                     plot(cch)
%                 %                     subplot(1,3,2)
%                 %                     imagesc(wf{i}')
%                 %                     subplot(1,3,3)
%                 %                     imagesc(wf{j}')
%                 %                     waitforbuttonpress
%                 %                     close all
%                 %                 end
%                 %handle empty refractory period (ISI violation <
%                 
%                 % CC(i,j) = tmp(2)/length(tsti);
%             end
%         end
%         
%         
%         
%         
%     end
% end
% %%
% 
% 
% % now find all confusions and merge in order of rate
% 
% rate = [];
% for j = 1:nclu
%     
%     rate(j) = sum(clu==uclu(j));
% end
% 
% [~,b] = sort(rate,'descend');
% 
% 
% mrg1 = mrg(b);
% uclu1 = uclu(b);
% 
% 
% for i=1:nclu
%     idx = b(i);
%     if ~isempty(mrg1{i})
%         
%         prs = mrg1{i};
%         otherprs = b(find(cellfun(@(a) any(ismember(a,prs)),mrg1)));
%         otherprs = [otherprs cell2mat(mrg1(cellfun(@(a) any(ismember(a,prs)),mrg1)))];
%         prs = unique([idx prs otherprs]);
%         
%         % check if all on sahnk shank
%         
%         sh = shank(prs);
%         
%         kp = sh == shank(idx);
%         
%         prs = prs(kp);
%         
%         % merge in order of highest rate
%         [~,b1] = sort(rate(prs),'descend');
%         
%         prs = prs(b1);
%         
%         
%         
%         ref = uclu(prs(1));
%         targets = prs(2:end);
%         
%         
%         
%         for j = 1:length(targets)
%             tsti = ts(clu1==ref);
%             tstj = ts(clu1==uclu(targets(j)));
%             
%             % merge spikes
%             
%             [mergets,bb]  = sort([tsti;tstj]);
%             
%     
%             cch = CrossCorr(mergets,mergets,.0005,101);
%             pred = mean(cch(1:10));
%             
%             cch = cch(51:56);
%             
%             if ~(any(cch>pred))
%                 
%                 % reassign good spikes that one have isi violations
%                 
%                 clu1(clu1==uclu(targets(j))) = ref;
%                 fprintf(fid1, [num2str(uclu(targets(j))) '	merged with ' num2str(ref) '\n']);
%                 % fprintf(fid, [num2str(uclu(targets(j))) '	noise\n']);
%                 
%                 %delete other merges
%                 mrg1{uclu1==uclu(targets(j))} = [];
%             else
%                 
%                 % if the lower rate unit has rate lower than thres2,
%                 % delete
%                 if rate(targets(j)) < .5
%                     fprintf(fid, [num2str(uclu(targets(j))) '	noise\n']);
%                     fprintf(fid1, [num2str(uclu(targets(j))) '	similar to ' num2str(ref) '\n']);
%                     %delete other merges
%                     mrg1{uclu1==uclu(targets(j))} = [];
%                 end
%             end
%         end
%     end
% end

%%

writeNPY(clu1, clufil)
fclose(fid)
fclose(fid1)
end