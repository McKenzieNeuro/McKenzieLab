function [ups,dwns]  = sm_getDigitalin(dirName,fname,fs,Nchan)


%[x,y] = Process_IntanDigitalChannels([dirName '/' basename '_digitalin.dat']);
fname = [dirName filesep fname];

if ~exist(fname)
    fname = [dirName '/digitalin.dat'];
end

fidI = fopen(fname,'r');


chunk = 1e5;
sizeInBytes = 2; % int16
nbChan = 1;
fInfo = dir(fname);
nBytes = fInfo.bytes;
nbChunks = floor(nBytes/(nbChan*sizeInBytes*chunk));

if isempty(Nchan)
Nchan = 16;
end
Nchan2 = Nchan+1;

pulses = cell(Nchan,1);
pulses2 = cell(Nchan,1);
for ii=1:nbChunks
    h=waitbar(ii/(nbChunks+1));
    x = fread(fidI,nbChan*chunk,'int16');
    
    
    digital_word2 = double(x);
    
    
    for k = 1:Nchan
        tester = (digital_word2 - 2^(Nchan-k))>=0;
        digital_word2 = digital_word2 - tester*2^(Nchan-k);
        
        
        pulses{Nchan2-k} = [pulses{Nchan2-k};find(diff(tester==1)== 1) + chunk*(ii-1)];
        pulses2{Nchan2-k} =[pulses2{Nchan2-k};find(diff(tester==1)== -1) + chunk*(ii-1)];
    end
 
end


remainder = nBytes/(sizeInBytes*nbChan) - nbChunks*chunk;
if ~isempty(remainder)
    x = fread(fidI,nbChan*remainder,'int16');
    digital_word2 = double(x);
    
  
    for k = 1:Nchan
        tester = (digital_word2 - 2^(Nchan-k))>=0;
        digital_word2 = digital_word2 - tester*2^(Nchan-k);
        
        
        pulses{Nchan2-k} = [pulses{Nchan2-k};find(diff(tester==1)== 1) + chunk*nbChunks];
        pulses2{Nchan2-k} =[ pulses2{Nchan2-k};find(diff(tester==1)== -1) + chunk*nbChunks];
    end
    
end



for ch = 1:Nchan
    ups{ch} = pulses{ch}/fs;
    dwns{ch}= pulses2{ch}/fs;
    
    if length(ups{ch}) ~= length(dwns{ch}) ||  ~all(dwns{ch}>ups{ch})
        
        warning('length mismatch, must manually fix')
    end
    
end


save([dirName filesep 'TTL_pulse.mat'],'ups','dwns')
disp('Saved TTL_pulse.mat')

end




