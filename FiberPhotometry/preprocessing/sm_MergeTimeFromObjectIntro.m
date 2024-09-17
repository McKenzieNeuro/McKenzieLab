function sm_MergeTimeFromObjectIntro(dirName)


cd(dirName)

if exist([dirName filesep 'sessiondata.mat']) && exist([dirName filesep 'novelObjectIntro.mat'])
    
    load('sessiondata.mat','sessiondata')
    load('novelObjectIntro.mat','data')
    
       timeFrmCChange = nan(size(sessiondata.behavior.ts_video(:,1)));
    sessiondata.behavior.timeFrmEntry = nan(size(sessiondata.behavior.ts_video(:,1)));

    for i = 1:size(data,1)-1
        kp = sessiondata.behavior.ts_video>data{i,2} & sessiondata.behavior.ts_video<data{i+1,2};
        timeFrmCChange(kp) =   sessiondata.behavior.ts_video(kp) - data{i,2};
        
        
        
    end
    
    
    kp = sessiondata.behavior.ts_video>data{end,2};
    timeFrmCChange(kp) =   sessiondata.behavior.ts_video(kp) - data{end,2};
    sessiondata.behavior.timeFrmObj = timeFrmCChange;
    save('sessiondata.mat','sessiondata','-v7.3')
end

end




