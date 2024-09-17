function [raw_data, time, FS, skip, fname, iOk] = read_smrx_data_KAsz(stim_metatable, row_to_extract, data_dir, channels)
%read in raw data once the proper files have been identified
%   Detailed explanation goes here
skip = 0;
raw_data = [];



%need to transform the trial number format in the excel file to the actual
%file name in the folder (i.e. add 'Trial '  or append '.smrx', etc.
if contains(stim_metatable.TrialNumber(row_to_extract),'KA92')
    fname = append('E:\HFS ANT92\',stim_metatable.TrialNumber(row_to_extract),'_ns3.smrx');
elseif max(strlength(stim_metatable.TrialNumber(row_to_extract)) == [5,6,9,13,16])  % || length(stim_metatable.TrialNumber(row_to_extract)) == 6 || length(stim_metatable.TrialNumber(row_to_extract)) == 9
    if contains(stim_metatable.TrialNumber(row_to_extract),'2020')
        fname = append(data_dir,'\',stim_metatable.TrialNumber(row_to_extract),'_ns3.smrx');
    else
        fname = append(data_dir,'\','Trial ',stim_metatable.TrialNumber(row_to_extract),'_ns3.smrx');
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%import the data

%open the file handle
fhand = CEDS64Open(convertStringsToChars(fname));      %should be 1 if loaded correctly, -1 otherwise
if fhand < 0
    error('Error loading file')
    skip = 1
end


%pull out the channel names by looping through all the channels
[iChans] = CEDS64MaxChan(fhand);
chanTitles = {};
for i = 1:iChans
    [~,chanTitles{i}] = CEDS64ChanTitle(fhand,i);
end



%find indices for 4 channels of interest by comparing channel titles in
%.smrx file with titles supplied by 'channels' variable
%data is stored by 'Cage' since we sometimes record from multiple animals
%at the same time; this information is located in the excel sheet
chInds = [];
for i = 1:length(channels)
    chName = ['Cage',num2str(stim_metatable.CageNumber(row_to_extract)),' ',channels{i}];
    chName2 = ['Cage',num2str(stim_metatable.CageNumber(row_to_extract)),'  ',channels{i}];
    chName3 = ['Cage ',num2str(stim_metatable.CageNumber(row_to_extract)),' ',channels{i}];
    if ~isempty(find(strcmp(chanTitles,chName)))
        chInds(i) = find(strcmp(chanTitles,chName),1,'first');
    elseif ~isempty(find(strcmp(chanTitles,chName2)))
        chInds(i) = find(strcmp(chanTitles,chName2),1,'first');
    elseif ~isempty(find(strcmp(chanTitles,chName3)))
        chInds(i) = find(strcmp(chanTitles,chName3),1,'first');
    end
end


%actually read in the raw data
for i = 1:length(chInds)
    maxTimeTicks = CEDS64ChanMaxTime(fhand,chInds(i))+1;
    [fRead, fVals, ~] = CEDS64ReadWaveF(fhand,chInds(i),1e9,0,maxTimeTicks);
    if fRead>0
        raw_data(i,:) = fVals;
    else
        error('Could not read data from channel')
    end
end



FS = 1/(CEDS64ChanDiv(fhand, 1)*CEDS64TimeBase(fhand));
time = [1:size(raw_data,2)]/FS;


[iOk] = CEDS64Close(fhand);      %should be 1 if loaded correctly, -1 otherwise

end

