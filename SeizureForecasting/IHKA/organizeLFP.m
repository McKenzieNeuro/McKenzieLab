topDir = 'F:\data1\IHKA';

fils = getAllExtFiles(topDir,'dat',1);
dirN  = unique(fileparts(fils));

for i = 2:length(dirN)
    cd(dirN{i})
    fils = getAllExtFiles(dirN{i},'dat',0);
    if length(fils)==4
        [a,b] = fileparts(fils{1});
        
        outfil = ['R:\IHKA_Scharfman\IHKA data\mergedDat' filesep b(1:end-2) '_allCh.dat'];
        try
            sm_MergeDats(fils,outfil,1,41);
        catch
            disp(outfil)
        end
    end
end