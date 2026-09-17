animals={'NE2h25','NE2h27','NE2h30'}
grand_PETH=[];
grand_fils=[];
for idx = 1:length(animals)
    
    animal=animals{idx}
    
    fils = getAllExtFiles(['R:\DANEHippocampalResponse\',animal],'Tbk',1);
    kp = contains(fils,'StimOn\');
    fils = fils(kp);
    
    all_PETH = [];
    all_stim_ts = [];

    for i = 1:length(fils)
        dirN = fileparts(fils{i});
       
        sesPath = [dirN,'\sessiondata.mat'];
        if ~exist(sesPath,"file")
            [signal_DFoF,ts_data2,fs,data] = sm_getSignal_DFoF(dirN,'baseline',[5 300]);
        else
            sessiondata = load(sesPath);
            signal_DFoF                = sessiondata.neural.signal_DFoF;
            ts_data                    = sessiondata.neural.ts_data;
            fs                         = sessiondata.neural.fs;
            data                       = sessiondata.neural.data;
        end
        
        novelT=data.epocs.Note.onset(strcmp(data.epocs.Note.notes,'Novel'));
        homeT=data.epocs.Note.onset(strcmp(data.epocs.Note.notes,'Home'));
        transfers=[novelT,homeT];

        k = gaussian2Dfilter([1 fs*10 ],fs);
        signal_DFoF_conv = nanconvn(signal_DFoF,k);
        
        stimT=data.epocs.Pe2_.onset;
        kp = abs(stimT - bestmatch(stimT,transfers))>100;
        stimT = stimT(kp);

        [ix,early,late,ts1] = sm_getIndicesAroundEvent(stimT,30,100,fs,length(signal_DFoF));
        all_stim_ts = [all_stim_ts;stimT];
        
        all_PETH = cat(1,all_PETH,signal_DFoF_conv(ix));
        
    end
    plotsem(ts1,all_PETH,'k',[],fs,'sem');
    title(sprintf([animal,' PETH OptoStim \n n=',num2str(size(all_PETH,1)),' stimulations \n n=',num2str(length(fils)),' sessions']));
    saveas(gcf,['R:\BHarvey\FPplots\',animal,'\',animal,'CA1Stim_PETH_Figure_nodrop.png']);
    
    grand_fils=cat(1,grand_fils,fils);
    grand_PETH=cat(1,grand_PETH,all_PETH);
end



plotsem(ts1,grand_PETH,'k',[],fs,'sem');
title(sprintf(['Average PETH OptoStim \n n=',num2str(size(grand_PETH,1)),' stimulations \n n=',num2str(length(grand_fils)),' sessions']));
saveas(gcf,['R:\BHarvey\FPplots\Average\Average_CA1Stim_PETH_Figure.png']);



%%



%%