function sm_MergeDats(fnameIn,fnameOut)

sizeInBytes = 2; % .
chunk = 1e5; % depends on the system

fInfo = dir(fnameIn{1});

nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(sizeInBytes*chunk));
fidO = fopen(fnameOut,'w');

nFiles = length(fnameIn);
for i = 1:nFiles
    fidI(i) = fopen(fnameIn{i},'r');
end



for ii=1:nbChunks
    h=waitbar(ii/(nbChunks+1));
    
    
    dat = nan(nFiles,chunk);
    for jj = 1:nFiles
        dat(jj,:) = fread(fidI(jj),chunk,'int16');
    end
    
    
    fwrite(fidO,dat(:),'int16');
end

remainder = nBytes/(sizeInBytes) - nbChunks*chunk;
if ~isempty(remainder)
    dat = nan(nFiles,remainder);
    for jj = 1:nFiles
        dat(jj,:) = fread(fidI(jj),remainder,'int16');
    end
    
    
    fwrite(fidO,dat(:),'int16');
end
close(h);
for i = 1:nFiles
    fclose(fidI(i));
end
fclose(fidO);

end


