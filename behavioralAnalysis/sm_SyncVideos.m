function sm_SyncVideos(topDir)

%amplifier basename
keyword = 'amplifier_analogin_auxiliary_int16.dat';

%get amplifier data

fils = getAllExtFiles(topDir,'dat',1);
kp = contains(fils,keyword);
fils = fils(kp);



%exclude the merged in the top directory
subdirs = fileparts(fils);


kp = ~strcmp(topDir,subdirs);
fils = fils(kp);

subdirs = fileparts(fils);



%sort by time
tims = cellfun(@(a) str2num(cell2mat(regexp(a,'[0-9]','match'))),subdirs);
[~,b] = sort(tims);
fils = fils(b);


for i = 1:length(fils)
    
    xmlfil = strrep(fils{i},'dat','xml');
    
    if ~exist(xmlfil)
        xmlfil = [topDir filesep keyword(1:end-3) 'xml'];
    end
    xml = LoadXml(xmlfil);
    ok = dir(fils{i});
    siz = ok.bytes;
    dur(i) = siz/xml.SampleRate/xml.nChannels/2;
    
end

% now find videos

fils = getAllExtFiles(topDir,'mp4',1);
clear  outTable in
for i = 1:length(fils)
    disp(['Processing video ' num2str(i)])
    
    in{i} = sm_defineMaze(fils{i});
    
end


for i = 1:length(fils)
    disp(['Processing video ' num2str(i)])
    
    [whl,in,threshF,outTable{i}] = sm_ExtractLed2(fils{i},'in',in{i});
end

end