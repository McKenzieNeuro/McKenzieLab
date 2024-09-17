function [rows_to_extract, stim_metatable] = read_metatable_KAsz(subject_IDs, sz_files)
%script to load in data from excel sheet
%   Detailed explanation goes here

%load meta data about stim trial organization
fname = "R:\IHKA_gross\KA Spike2 + Videos\Black Rock Recording Data Base Separate Updated_TEE_sm.xlsx";
stim_metatable =readtable(fname,'sheet',subject_IDs);
%convert TrialNumber to string, then use readcell and move strings over to
%first stim_metatable to complete the table
stim_metatable.TrialNumber = string(stim_metatable.TrialNumber);
%fname = "R:\IHKA_gross\KA Spike2 + Videos\Black Rock Recording Data Base Separate Updated_TEE.xlsx";

stim_temp = readcell(fname,'sheet',subject_IDs);
for i = 1:size(stim_metatable.TrialNumber,1)
    if ismissing(stim_metatable.TrialNumber(i)) && length(stim_temp{i+1,1}) > 4
        if strcmp(stim_temp{i+1,1}(1:4),'2020')
            stim_metatable.TrialNumber(i) = stim_temp{i+1,1};
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select rows to extract space

%only rows with the matching ID
name_inds = strcmp(stim_metatable.AnimalIdentification,subject_IDs);

%option to only load files which had a seizure
if sz_files
    sz_inds = ~cellfun(@isempty,cellfun(@str2num,stim_metatable.Seizures_sec_,'UniformOutput',false));
else
    sz_inds = ones(size(name_inds,1),1);
end

%remove rows which don't have converted files yet
no_file = (cellfun(@length,arrayfun(@num2str, stim_metatable.TrialNumber, 'UniformOutput', false)) > 3);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%determine which rows to load from above criteria 
options = [name_inds,sz_inds, no_file];
rows_to_extract = find(floor(mean(options,2)));

