fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);


%only keep files with NE2 and Novel Env
kp = cellfun(@any,regexp(fils,'Novel Env')) & (cellfun(@any,regexp(fils,'NE2')) | cellfun(@any,regexp(fils,'NE3'))) & contains(fils,'sessiondata.mat');


fils = fils(kp);

%%
clear datAll
ixx = 1;
for i = 1:length(fils)
    dirN = fileparts(fils{i});
    cd(dirN)
    load(fils{i})
    
    nSes = size(sessiondata.contextEntry,1);
    day1 = all(cell2mat(sessiondata.contextEntry(:,end))==0) & nSes==7;
    
    if day1
        data = sessiondata.neural.signal_DFoF;
        
        fs = sessiondata.neural.fs_neural;
        k = gaussian2Dfilter([ 1 fs*10 ],fs);
        data = nanconvn(data,k);
        
        siz = length(data);
        for j = 2:7
        evs = sessiondata.contextEntry{j,2};
         [ix,early,late,ts] = sm_getIndicesAroundEvent(evs,300,300,fs,siz);
         
         datAll{ixx,j-1} = data(ix);
        end
           ixx = ixx+1
    end
    
 
    
end
%%
close all
ok = cell2mat(datAll);
 ts = (1:1000:size(ok,2))/fs;
 figure
 plotMeanSEM(ts,ok(:,1:1000:end),'b')
 
 hold on
 plot([300:600:3600;300:600:3600],[-1.5 3],'k')
 plot([0 3600],[0 0 ],'k')
 set(gcf, 'Renderer', 'painters')