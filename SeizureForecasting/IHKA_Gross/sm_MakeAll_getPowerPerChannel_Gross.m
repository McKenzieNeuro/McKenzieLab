%Read file for KA recordings

%This script has been adapted to read in all the data files for the animals
%in the recordings, using the excel datasheet.

%If you have questions, please contact me at teggers@emory.edu


clear all
close all
clc

%set random number generator to help isolate occasional breaks in the code
rng(123)


subject_IDs = {'KA5','KA6','KA8','KA12','KA13','KA14','KA11'};
%these are the subject IDs for all the animals on file
%subject_IDs = {'KA11'};  %for testing, only read in one animals data at a time


sz_files = true;                   %whether to load only the files with seizures (true) or all files (false)


%May need to change directory location based on
data_dir = 'R:\IHKA_gross\KA Spike2 + Videos';

cd('R:\IHKA_gross\KA Spike2 + Videos')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add paths to necessary folders

% Here's a url to the CED website where you can download the file (after
% giving them your name/email): https://ced.co.uk/upgrades/spike2matson
% I also included the files on the hard drive (CEDMATLAB folder), but you
% may have to download them separately to get them to work - I followed the
% instructions on the s64mat file and didn't have much issue. Let me know
% if you have trouble.

%load CED MATSON library for reading in .smrx files
cedpath = getenv('R:\Analysis\CEDMATLAB\CEDS64ML');        %this part may be where you need to edit;
cedpath = 'R:\Analysis\CEDMATLAB\CEDS64ML';
%you might be able to switch to the manual path location rather than use
%the 'getenv' function; see pg. 1-6 of the s64mat pdf for more info.
addpath(cedpath)
CEDS64LoadLib(cedpath)

%other libraries - again, may have to adjust when you start
%addpath(genpath('D:\file_for_Sam\'))


for k1 = 2:length(subject_IDs) %loops through all the animal IDs you provide at the top of the code
    clc
    fprintf('Starting %s... ',subject_IDs{k1});
    
    %load in excel file
    [rows_to_extract, stim_metatable] = read_metatable_KAsz(subject_IDs{k1}, sz_files);
    %this file should work without having to touch it - it decides which
    %rows from the excel file to load for processing. There are some
    %heuristic bits of code to make sure only the rows with proper data are
    %loaded, since there are some files missing here and there and some
    %files had unusual noise etc.
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop through each file and load the data (several channels from each
    %recording possible, i.e. L/R HPC, EEG, EMG)
    excel_length = [];
    file_length = [];
    
    %     I had this bit of code to remove some files from KA5 since I use the
    %     EMG data in my origianl sz detection, but this probably won't affect
    %     your approach so I commented it out; I also commented out anything
    %     related to KA92 - this animal received ANT stimulation, so the
    %     recordings are full of artifact
    
    %     %last 2 files of KA5 had unusually high noise in EMG, which throws off
    %     %seizure detection; ignore these trials for now
    %     if sz_files == false && strcmp(subject_IDs{k1}, 'KA5')
    %         rows_to_extract(end-1:end) = [];
    %
    %     elseif sz_files == false && strcmp(subject_IDs{k1}, 'KA92')
    %         %I only checked the files for the baseline/first stimulation, not
    %         %the second baseline/stim files
    %         rows_to_extract(26:end)=[];
    %     end
    
    %these are the times that have been manually marked for seizures; they
    %aren't all exactly the start time of a seizure, but are marked during
    %the seizure
    sz_manual{k1} = cellfun(@str2num,stim_metatable.Seizures_sec_(rows_to_extract),'UniformOutput',false);
    
    %now loop through each row that had been identified and load the data
    for k3 = 1:length(rows_to_extract)
        
        %determine which channels (EMG, L/R HPC, EEG) to load; you can set
        %this to only 'LHPC' if you want to only focus on the Left
        %Hippocampus
        channels = strsplit(stim_metatable.Channels{rows_to_extract(k3)},','); %{'EMG','LHPC','EEG','RHPC'};
        
        skip=0;
        %read in raw data
        [raw_data, time, FS, skip, fname, iOk] = read_smrx_data_KAsz(stim_metatable, rows_to_extract(k3), data_dir, channels);
        
        dat = stim_metatable(rows_to_extract(k3),:).TrialNumber;
        dirOut = ['G:\data\IHKA_gross' filesep subject_IDs{k1} filesep  subject_IDs{k1} '_' dat{1}];
        basename = [ subject_IDs{k1} '_' dat{1}];
        for ii = 1:size(raw_data,1)
            sm_getPowerPerChannel_Gross(raw_data(ii,:),num2str(ii),channels{ii},'Fs',FS,'dirOut',dirOut,'basename',basename)
        end
        
        %insert your seizure detection code here if you want to analyze one
        %file at a time
        % xxx
        
        
        
    end
    
end
