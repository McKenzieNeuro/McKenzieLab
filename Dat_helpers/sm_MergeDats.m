function sm_MergeDats(fnameIn,fnameOut,channels,nbChan)

sizeInBytes = 2; % .
chunk = 1e5; % depends on the system

fInfo = dir(fnameIn{1});

nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(sizeInBytes*chunk*nbChan));
fidO = fopen(fnameOut,'w');

nFiles = length(fnameIn);
for i = 1:nFiles
    fidI(i) = fopen(fnameIn{i},'r');
end



for ii=1:nbChunks
    h=waitbar(ii/(nbChunks+1));
    
    
    dat = nan(nFiles,chunk);
    for jj = 1:nFiles
        
        
        tmp = fread(fidI(jj),nbChan*chunk,'int16');
        tmp = reshape(tmp,nbChan,[]);
        tmp = tmp(channels,:);
        dat(jj,:) = tmp;
        
    end
    
    
    fwrite(fidO,dat(:),'int16');
end

remainder = nBytes/(sizeInBytes) - chunk*nbChan;
if ~isempty(remainder)
    dat = nan(nFiles,floor(remainder/nbChan));
    for jj = 1:nFiles
       
        tmp = fread(fidI(jj),remainder,'int16');
        tmp = reshape(tmp,nbChan,[]);
        tmp = tmp(channels,:);
        dat(jj,:) = tmp;
        
    end
    
    
    fwrite(fidO,dat(:),'int16');
end
close(h);
for i = 1:nFiles
    fclose(fidI(i));
end
fclose(fidO);

end


