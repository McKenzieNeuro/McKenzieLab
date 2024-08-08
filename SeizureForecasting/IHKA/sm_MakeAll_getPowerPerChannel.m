%find seizures with files that have been converted to edf

FilePath = 'R:\IHKA_Scharfman\IHKA data';

fils_edf = getAllExtFiles(FilePath,'edf',1);
fils_txt = getAllExtFiles(FilePath,'txt',1);
[~,b_edf] = cellfun(@fileparts,fils_edf,'uni',0);
[~,b_txt] = cellfun(@fileparts,fils_txt,'uni',0);

goodFils = intersect(b_txt,b_edf);

seizure_fils = fils_txt(ismember(b_txt,goodFils));

edf_fils = fils_edf(ismember(b_edf,goodFils));



%%


masterDir = 'E:\data\IHKA';

for i = 74:length(fils_edf)
    
    % check if file exists
    fileOut = regexprep(fils_edf{i},' ','_');
    
    [a,basename] = fileparts(fileOut); % define output directory
    dirOut = [masterDir filesep basename ];
    if ~exist(dirOut)
        % loop over 4 channels
        for j = 1:4
            sm_getPowerPerChannel(fils_edf{i},j)
        end
    end
    i
end

