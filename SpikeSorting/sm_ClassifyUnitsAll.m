fils = getAllExtFiles('R:\STDP\STDP4','npy',1);
fils = fils(contains(fils,'Kilosort_2023') & contains(fils,'pc_features'));

dirN = fileparts(fileparts(fils));
basename =cellfun(@bz_BasenameFromBasepath,dirN,'uni',0);
%%
for i = 1:length(dirN)
    
    fname = [dirN{i} filesep basename{i} '.dat'];
    
    if exist(fname)
        try
            % exclude noisy units
            [good,wf,shank] = sm_assign_noise(fname);
            sm_MergeCluster(datfil,'good',good,'shank',shank,'wf',wf);
        end
    end
end