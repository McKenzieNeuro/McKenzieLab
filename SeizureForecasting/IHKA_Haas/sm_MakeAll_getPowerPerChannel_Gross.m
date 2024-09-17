%find seizures with files that have been converted to edf

FilePath = 'E:\Dropbox\Data SamMcKenzie\Data_raw';

fils_h5 = getAllExtFiles(FilePath,'.h5',1);



% get unique session
fils_h5 = unique(cellfun(@(a,b) a(1:b(3)-1),fils_h5,regexp(fils_h5,'_'),'UniformOutput',false));

%%


masterDir = 'G:\data\IHKA_Haas';

for i = 1:length(fils_h5)
  
    [a,basename] = fileparts(fils_h5{i}); % define output directory
    dirOut = [masterDir filesep basename ];
    if ~exist(dirOut)
        
        for j = 1:3
            sm_getPowerPerChannel_Gross(fils_h5{i},j)
        end
    end
    i
end

