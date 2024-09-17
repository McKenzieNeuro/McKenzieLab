function sm_MergeRewardTime(dirName)


cd(dirName)

if exist([dirName filesep 'sessiondata.mat'])
    
    load('sessiondata.mat','sessiondata')
    
    TDTdata = TDTbin2mat(dirName);
    
    if isfield(TDTdata.epocs,'U12_') &&isfield(TDTdata.epocs,'U11_')
        
        
        rew1 = TDTdata.epocs.U12_.onset;
        rew2 = TDTdata.epocs.U11_.onset;
        sessiondata.behavior.rewardTimeLeft = rew1;
        sessiondata.behavior.rewardTimeRight = rew2;
        
        
        save('sessiondata.mat','sessiondata','-v7.3')
    end
    
end




