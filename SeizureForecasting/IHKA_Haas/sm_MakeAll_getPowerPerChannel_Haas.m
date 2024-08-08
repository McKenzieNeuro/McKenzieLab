%find seizures with files that have been converted to edf

FilePath = 'E:\Dropbox\Data SamMcKenzie\Data_raw';
FilePath = 'R:\IHKA_Haas\LFS';
fils_h5 = getAllExtFiles(FilePath,'.h5',1);



% get unique session
%fils_h5 = unique(cellfun(@(a,b) a(1:b(3)-1),fils_h5,regexp(fils_h5,'_'),'UniformOutput',false));

for i = 1:length(fils_h5)
    sl = regexp(fils_h5{i},'_');
fils_h5{i} = [fils_h5{i}(1:sl(4)) 'channel' fils_h5{i}(sl(5):end)]  ;  
end

%%


masterDir = 'G:\data\IHKA_Haas\LFS';

for i = 1:4%5:length(fils_h5)
  
    [a,basename] = fileparts(fils_h5{i}); % define output directory
    dirOut = [masterDir filesep basename ];
    
        
        for j = 1:2
            sm_getPowerPerChannel_Haas1(fils_h5{i},j)
        end
    
    i
end

