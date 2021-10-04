%% Read TDT Fiber Photometry Data into UnfoldToolbox Format
%  This example reads fiber photometry data and converts it to 
% <https://www.unfoldtoolbox.org/index.html unfoldtoolbox.org> format.
%  Artifact removal, detrending, dFF calculation is also included.

%% Housekeeping
close all; clear all; clc;

BLOCK_PATH = 'C:\TDT\TDTExampleData\FiPho-180416';
GCAMP = 'x4654';
ISOS = 'x4054';
ART_TIME = 8; % time threshold below which we will discard, in seconds

%% Read TDT data into Matlab
data = TDTbin2mat(BLOCK_PATH);

%% Artifact removal
% There is often a large artifact on the onset of LEDs turning on
% Remove data below a set time ART_TIME

% Make a time array based on number of samples and sample freq of
% demodulated streams
time_arr = (1:length(data.streams.(GCAMP).data))/data.streams.(GCAMP).fs;

ind = find(time_arr > ART_TIME, 1); % find first index of when time crosses threshold
time_arr = time_arr(ind:end); % reformat vector to only include allowed time
data.streams.(GCAMP).data = data.streams.(GCAMP).data(ind:end);
data.streams.(ISOS).data = data.streams.(ISOS).data(ind:end);

%% Detrending and dFF
bls = polyfit(data.streams.(ISOS).data, data.streams.(GCAMP).data, 1);
Y_fit_all = bls(1) .* data.streams.(ISOS).data + bls(2);
Y_dF_all = data.streams.(GCAMP).data - Y_fit_all; %dF (units mV) is not dFF
dFF = 100*(Y_dF_all)./Y_fit_all;
std_dFF = std(double(dFF));

%% Put data and epoc events into UnfoldToolbox format
UF_DATA.filepath = BLOCK_PATH;
UF_DATA.data = dFF;
UF_DATA.pnts = size(dFF, 2);
UF_DATA.srate = data.streams.(GCAMP).fs;
UF_DATA.times = time_arr;
UF_DATA.event = [];
if isstruct(data.epocs)
    event_ind = 1;
    fff = fields(data.epocs);
    for n = 1:numel(fff)
        % adjust epocs based on artifact rejection time
        if ART_TIME > 0
            ind = data.epocs.(fff{n}).onset - ART_TIME > 0;
            data.epocs.(fff{n}).onset = data.epocs.(fff{n}).onset(ind) - ART_TIME;
            data.epocs.(fff{n}).offset = data.epocs.(fff{n}).offset(ind) - ART_TIME;
            data.epocs.(fff{n}).data = data.epocs.(fff{n}).data(ind);
            if isfield(data.epocs.(fff{n}), 'notes')
                data.epocs.(fff{n}).notes = data.epocs.(fff{n}).notes(ind);
            end
        end
        onsets = data.epocs.(fff{n}).onset;
        for j = 1:numel(onsets)
            UF_DATA.event(event_ind,1).type = fff{n};
            UF_DATA.event(event_ind,1).latency = onsets(j);
            UF_DATA.event(event_ind,1).value = data.epocs.(fff{n}).data(j);
            event_ind = event_ind + 1;
        end
    end
end
