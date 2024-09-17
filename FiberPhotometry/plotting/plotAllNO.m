fils = getAllExtFiles('R:\DANEHippocampalResponse','mat',1);

kp = cellfun(@any,regexp(fils,'SOR')) & cellfun(@any,regexp(fils,'NE2'));

fils = fils(kp);




[dirs] = fileparts(fils);
% [dirs] = cellfun(@fileparts,fils,'UniformOutput',false);

dirs =  unique(dirs);
%%
subj =[];
for i = 1:length(dirs)
    [~,subjt] = fileparts(dirs{i});
    subj_split = split(subjt, "-");
    subjt = subj_split{1};
    subj = [subj;{subjt}];
end
usubjs = unique(subj);
%%
ix=  1;
clear   sesID subjID
obj_intr = [];obj_samp=[];
for i = 1:length(dirs)
    
    cd(dirs{i})
    if exist('sessiondata.mat') && exist('novelObject.mat')
        load('sessiondata.mat')
        %  v=load('novelObjectIntro.mat');
        %  vv= load('novelObject.mat');
        %  sessiondata.behavior.objectIntro=v.data;
        %  sessiondata.behavior.objectSample=vv.data;
        %  save('sessiondata.mat','sessiondata','-v7.3')
        
        fs = sessiondata.neural.fs_neural;
        
        k  = gaussian2Dfilter([fs*10 1],fs);
        signal_DFoF = nanconvn(sessiondata.neural.signal_DFoF, k');
        
        
        objectIntro = cell2mat(sessiondata.behavior.objectIntro(:,2));
        objectSample = cell2mat(sessiondata.behavior.objectSample(:,2));
        [ix_intro,early1,late1,ts_PETH1] = sm_getIndicesAroundEvent(objectIntro,150,150,fs,length(signal_DFoF));
        
        [ix_sample,early1,late1,ts_PETH2] = sm_getIndicesAroundEvent(objectSample,5,5,fs,length(signal_DFoF));
        
        obj_intr = [obj_intr;signal_DFoF(ix_intro)];
         obj_samp = [obj_samp;nanmean(signal_DFoF(ix_sample))];
         ts = (1:length(signal_DFoF))/fs;
        figure
         plot(ts,signal_DFoF)
         hold on
         plot([objectIntro objectIntro],[-2 3],'k')
         ylim([-2 3])
        ix = ix+1
        close all
    end
end

%%
close all
figure
plotMeanSEM(ts_PETH2(1:100:end),obj_samp(:,1:100:end),'k')
ylim([-.4 1.5])
figure
plotMeanSEM(ts_PETH1(1:100:end),obj_intr(:,1:100:end),'k')
ylim([-.4 1.5])