function [h] = plot_sz_gross(featurefil_path,start_tim,varargin)


%INPUTS
% featurefil_path = path for each of the *.dat files that store the raw
% time series as well as the features for that channel (typically we use 41
% channels per file with channel 1 = time series, suffix = nCh_rec
% start_tim = time of event (s)
%
% OPTIONAL INPUTS
% nCh_dat = number of channels in *.dat file (typically 41 = 20 freq.*2
% nCh_rec = number of channels recorded
% (ph/amp) +1 (timeseries)
% fs = sampling rate of file (Hz)
% pre = duration before event to plot (s)
% duration = time to plot (s)

%% parse inputs
p = inputParser;


addParameter(p,'nCh_dat',41,@isnumeric)
addParameter(p,'nCh_rec',4,@isnumeric)
addParameter(p,'chan_2_load',1,@isnumeric)
addParameter(p,'fs',2000,@isnumeric)
addParameter(p,'pre',20,@isnumeric)
addParameter(p,'duration',100,@isnumeric)
addParameter(p,'basename',[],@isstr)
parse(p,varargin{:})

nCh_dat = p.Results.nCh_dat;
nCh_rec = p.Results.nCh_rec;
fs = p.Results.fs;
pre = p.Results.pre;
duration = p.Results.duration;
chan_2_load = p.Results.chan_2_load;
basename = p.Results.basename;

if isempty(basename)
    
    [~,basename] = fileparts(featurefil_path);
    
end
featurefil_path = [featurefil_path filesep basename];

%%
shift_factor =  2e4;
% loop over channel
h = figure;
for k = 1:nCh_rec
    
    fname = [featurefil_path '_' num2str(k) '.dat'];
    
    if exist(fname)
        dat = LoadBinary(fname,'frequency',fs,'start',start_tim-pre,...
            'duration',duration,'nchannels',nCh_dat,'channels',chan_2_load);
        ts = (1:length(dat))/fs - pre ;%+ start_tim ;
        plot(ts,double(dat) - shift_factor*(k-1),'k')
        hold on
        
    else
        error(['file: ' fname ' does not exist'])
    end
end

%calc spectrogram on channel 4 for each seizure with frequencies going from
%0.5Hz to the nyquist in 100 log spaced bins

%awt_freqlist


end