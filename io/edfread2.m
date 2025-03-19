function [hdr, dataMat] = edfread2(fname)
% Read European Data Format file into MATLAB
%
% [hdr, dataMat] = edfread2(fname)
%         Reads data from ALL RECORDS of file fname ('*.edf'). Header
%         information is returned in structure hdr, and the signals
%         (waveforms) are returned in dataMat, which is either a (samples x
%         channels) matrix, or a (1xchannels) cell array full of (samples x
%         1) arrays IF there are different numbers of samples in each
%         channel
%
% FORMAT SPEC: Source: http://www.edfplus.info/specs/edf.html SEE ALSO:
% http://www.dpmi.tu-graz.ac.at/~schloegl/matlab/eeg/edf_spec.htm
%
% The first 256 bytes of the header record specify the version number of
% this format, local patient and recording identification, time information
% about the recording, the number of data records and finally the number of
% signals (ns) in each data record. Then for each signal another 256 bytes
% follow in the header record, each specifying the type of signal (e.g.
% EEG, body temperature, etc.), amplitude calibration and the number of
% samples in each data record (from which the sampling frequency can be
% derived since the duration of a data record is also known). In this way,
% the format allows for different gains and sampling frequencies for each
% signal. The header record contains 256 + (ns * 256) bytes.
%
% Following the header record, each of the subsequent data records contains
% 'duration' seconds of 'ns' signals, with each signal being represented by
% the specified (in the header) number of samples. In order to reduce data
% size and adapt to commonly used software for acquisition, processing and
% graphical display of polygraphic signals, each sample value is
% represented as a 2-byte integer in 2's complement format. Figure 1 shows
% the detailed format of each data record.
%
% DATA SOURCE: Signals of various types (including the sample signal used
% below) are available from PHYSIONET: http://www.physionet.org/
%
%
% % EXAMPLE 1:
% % Read all waveforms/data associated with file 'ecgca998.edf':
% [header, recorddata] = edfRead('ecgca998.edf');
%
% Coded 8/27/09 by Brett Shoelson, PhD
% brett.shoelson@mathworks.com
% Copyright 2009 - 2012 MathWorks, Inc.
% 
% Modified 7/2/14 by Matt Fifer
% msfifer@gmail.com
% changed input options (removed some), changed output type to either
% samples x channels array or 1 x channels cell array of samples x 1 arrays
%
% Modifications:
% 7/2/14 Made major modifications to this function to make it a faster
% import, removed AssignToVariables and DesiredSignals options
%
% 5/6/13 Fixed a problem with a poorly subscripted variable. (Under certain
% conditions, data were being improperly written to the 'records' variable.
% Thanks to Hisham El Moaqet for reporting the problem and for sharing a
% file that helped me track it down.)
% 
% 5/22/13 Enabled import of a user-selected subset of signals. Thanks to
% Farid and Cindy for pointing out the deficiency. Also fixed the import of
% signals that had "bad" characters (spaces, etc) in their names.
% HEADER RECORD
% 8 ascii : version of this data format (0)
% 80 ascii : local patient identification
% 80 ascii : local recording identification
% 8 ascii : startdate of recording (dd.mm.yy)
% 8 ascii : starttime of recording (hh.mm.ss)
% 8 ascii : number of bytes in header record
% 44 ascii : reserved
% 8 ascii : number of data records (-1 if unknown)
% 8 ascii : duration of a data record, in seconds
% 4 ascii : number of signals (ns) in data record
% ns * 16 ascii : ns * label (e.g. EEG FpzCz or Body temp)
% ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode)
% ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC)
% ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
% ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
% ns * 8 ascii : ns * digital minimum (e.g. -2048)
% ns * 8 ascii : ns * digital maximum (e.g. 2047)
% ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
% ns * 8 ascii : ns * nr of samples in each data record
% ns * 32 ascii : ns * reserved
% DATA RECORD
% nr of samples[1] * integer : first signal in the data record
% nr of samples[2] * integer : second signal
% ..
% ..
% nr of samples[ns] * integer : last signal
if nargin > 5
    error('EDFREAD: Too many input arguments.');
end
if ~nargin
    error('EDFREAD: Requires at least one input argument (filename to read).');
end
[fid,msg] = fopen(fname,'r');
if fid == -1
    error(msg)
end
% HEADER
hdr.ver        = str2double(char(fread(fid,8)'));
hdr.patientID  = fread(fid,80,'*char')';
hdr.recordID   = fread(fid,80,'*char')';
hdr.startdate  = fread(fid,8,'*char')';% (dd.mm.yy)
% hdr.startdate  = datestr(datenum(fread(fid,8,'*char')','dd.mm.yy'), 29); %'yyyy-mm-dd' (ISO 8601)
hdr.starttime  = fread(fid,8,'*char')';% (hh.mm.ss)
% hdr.starttime  = datestr(datenum(fread(fid,8,'*char')','hh.mm.ss'), 13); %'HH:MM:SS' (ISO 8601)
hdr.bytes      = str2double(fread(fid,8,'*char')');
reserved       = fread(fid,44);
hdr.records    = str2double(fread(fid,8,'*char')');
hdr.duration   = str2double(fread(fid,8,'*char')');
% Number of signals
hdr.ns    = str2double(fread(fid,4,'*char')');
for ii = 1:hdr.ns
    hdr.label{ii} = regexprep(fread(fid,16,'*char')','\W','');
end
for ii = 1:hdr.ns
    hdr.transducer{ii} = fread(fid,80,'*char')';
end
% Physical dimension
for ii = 1:hdr.ns
    hdr.units{ii} = fread(fid,8,'*char')';
end
% Physical minimum
for ii = 1:hdr.ns
    hdr.physicalMin(ii) = str2double(fread(fid,8,'*char')');
end
% Physical maximum
for ii = 1:hdr.ns
    hdr.physicalMax(ii) = str2double(fread(fid,8,'*char')');
end
% Digital minimum
for ii = 1:hdr.ns
    hdr.digitalMin(ii) = str2double(fread(fid,8,'*char')');
end
% Digital maximum
for ii = 1:hdr.ns
    hdr.digitalMax(ii) = str2double(fread(fid,8,'*char')');
end
for ii = 1:hdr.ns
    hdr.prefilter{ii} = fread(fid,80,'*char')';
end
for ii = 1:hdr.ns
    hdr.samples(ii) = str2double(fread(fid,8,'*char')');
end
for ii = 1:hdr.ns
    reserved    = fread(fid,32,'*char')';
end
hdr.label = regexprep(hdr.label,'\W','');
hdr.units = regexprep(hdr.units,'\W','');
if nargout > 1
    % Scale data (linear scaling)
    scalefac = (hdr.physicalMax - hdr.physicalMin)./(hdr.digitalMax - hdr.digitalMin);
    dc = hdr.physicalMax - scalefac .* hdr.digitalMax;
    
    % RECORD DATA REQUESTED
    % multiple samples per record
    if ~all(hdr.samples == 1)
        
        % warn the user this could be a while
        warning('EDFREAD: More than one sample per record, might take a while');
        
        % initialize the output data cell and auxiliary var for writing to it
        dataMat = arrayfun(@(x) nan(x * hdr.records, 1), hdr.samples, 'UniformOutput', 0);
        writeIndices = arrayfun(@(x) (1:x)', hdr.samples, 'UniformOutput', 0);
        
        % create scaling and offsetting variables
        sampsPerRec = sum(hdr.samples);
        scalefac = cell2mat(arrayfun(@(x, y) ones(1, x) * y, hdr.samples, scalefac, 'UniformOutput', 0))';
        dc = cell2mat(arrayfun(@(x, y) ones(1, x) * y, hdr.samples, dc, 'UniformOutput', 0))';
        
        % initialize a waitbar for updating
        waitHandle = waitbar(0, 'Extracting data records...');
        
        % iterate through every record
        for recnum = 1:hdr.records
            
            % get all the data in this record
            tempData = mat2cell(fread(fid,sampsPerRec,'int16') .* scalefac + dc, hdr.samples, 1)';
            
            % write the data from each channel to a unique matrix cell
            for chans = 1:length(dataMat)
                dataMat{chans}(writeIndices{chans}) = tempData{chans};
                writeIndices{chans} = writeIndices{chans} + hdr.samples;
            end
            
            % update the waitbar
            waitbar(recnum/hdr.records, waitHandle);
        end
        
        % close the waitbar
        close(waitHandle);
        
        % if all the channels have the same number of samples, return a
        % matrix
        if length(dataMat{1}) == cellfun(@length, dataMat)
            dataMat = cell2mat(dataMat);
        end
        
    % one sample per record    
    else
        
        % create scaling and offsetting variables
        scalefac = repmat(scalefac, 1, hdr.records)';
        dc = repmat(dc, 1, hdr.records)';
        
        % read in all the data and reshape it into a matrix
        dataMat = fread(fid,hdr.records*hdr.ns,'int16') .* scalefac + dc;
        dataMat = reshape(dataMat, hdr.ns, hdr.records)'; % returns records (samples) x channels
    end
end
% close the file
fclose(fid);