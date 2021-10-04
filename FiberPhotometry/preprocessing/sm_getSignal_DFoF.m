function [signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(dirName,varargin)
%% Loads corrected, normalized fiber photometry signal
% Neural data should be stored in TDT format and read with TDT's SDK TDTbin2mat
% neural data should be stored in dirName

% INPUTS
% dirName =directory with TDT data and evFile
% varargin =
%       photoBleachCorrection  = [linear, quadratic , exp1 (default = exp2)
%       skiptime = time into recording to skip (s)
%       endTime = time into recording to consider (s)
%       streams = name of TDT data stream
%       isosbestic = name of isosbestic data stream
%       returnedDataType = isosbestic correction (default = corrected)

% OUTPUTS
%  signal_DFoF = photobleach, isosbestic, zscore signal
%  ts_data = time stamps of data (s)
%  fs = sampling rate (Hz)
%
% USAGE
% sm_getSignal_DFoF(pwd,'photoBleachCorrection','quadtratic')
%
% DEPENDENCIES
% Clone Mckenzieneurolab github (https://github.com/McKenzieNeuro/McKenzieLab)
% TDT SDKs (https://www.tdt.com/support/matlab-sdk/)
%
% S. McKenzie, 08/23/2021
%%

warning off
p = inputParser;


addParameter(p,'photoBleachCorrection','exp2',@ischar) %method to correct for decaying signal
addParameter(p,'skiptime',1,@isnumeric) % how much time (s) to skip before calculating mean/std/decay
addParameter(p,'streams',{'x465A','x405A'},@iscell) %which LED streams are recorded?
addParameter(p,'isosbestic','x405A',@ischar) % which stream is the isosbestic?
addParameter(p,'endTime',[],@isnumeric) % when was the LED turned off (s)?
addParameter(p,'returnedDataType','corrected',@ischar) % do we want (corrected,signal, isosbestic,raw)

parse(p,varargin{:})


photoBleachCorrection = p.Results.photoBleachCorrection;
skiptime = p.Results.skiptime;
streams = p.Results.streams;
isosbestic = p.Results.isosbestic;
endTime = p.Results.endTime;
returnedDataType = p.Results.returnedDataType;


signal = setdiff(streams,isosbestic);
%%
%load data
data = TDTbin2mat(dirName);
% get sampling rate (assume it's fixed across data streams)
fs = data.streams.(streams{1}).fs;
ts_data = (1:length(data.streams.(signal{1}).data))/fs;

%%

if strmatch(returnedDataType,'raw_signal')
    signal_DFoF = data.streams.(signal{1}).data;
elseif strmatch(returnedDataType,'raw_iso')
    signal_DFoF = data.streams.(isosbestic).data;
    
else
    
    nSamples = length(data.streams.(signal{1}).data);
    if isempty(endTime)
        
        endTime = nSamples ;
        
    end
    
    
    %%
    %convert skip and end times into samples
    endTime = endTime*fs;
    
    
    if endTime>nSamples
        endTime = nSamples;
    end
    
    skiptime = skiptime*fs;
    
    if skiptime<=0
        skiptime = 1;
    end
    
    
    
    
    %%
    
    % subtract photobleach for each stream
    
    for i = 1:length(streams)
        
        if length(data.streams.(streams{i}).data)<endTime
            endTime1 = length(data.streams.(streams{i}).data);
        else
            endTime1 = endTime;
            
        end
        switch photoBleachCorrection
            
            case 'linear'
                a = polyfit(ts_data(skiptime:endTime1),data.streams.(streams{i}).data(skiptime:endTime1),1);
                signal_bleach_corr = polyval(a,ts_data);
            case 'quadratic'
                
                a = polyfit(ts_data(skiptime:endTime1),data.streams.(streams{i}).data(skiptime:endTime1),2);
                signal_bleach_corr = polyval(a,ts_data);
            case 'cubic'
                a = polyfit(ts_data(skiptime:endTime1),data.streams.(streams{i}).data(skiptime:endTime1),3);
                signal_bleach_corr = polyval(a,ts_data);
            case 'exp1'
                
                a = fit(ts_data(skiptime:endTime1)',data.streams.(streams{i}).data(skiptime:endTime1)','exp1');
                signal_bleach_corr = a(ts_data);
                signal_bleach_corr = signal_bleach_corr';
            case 'exp2'
                a = fit(ts_data(skiptime:endTime1)',data.streams.(streams{i}).data(skiptime:endTime1)','exp2');
                signal_bleach_corr = a(ts_data);
                signal_bleach_corr = signal_bleach_corr';
                
                
        end
        
        
        
        endTime1 = min(length( data.streams.(signal{1}).data(:)),length(data.streams.(isosbestic).data(:)));
        %detrend
        if strcmp(streams{i},signal)
            
            
            
            signal_photoBleachCorrected= data.streams.(signal{1}).data(1:endTime1) - signal_bleach_corr(1:endTime1);
        else
            
            isosbestic_photoBleachCorrected= data.streams.(isosbestic).data(1:endTime1) - signal_bleach_corr(1:endTime1);
        end
        
        
        
        
    end
    %%
    %rescale isosbestic
    
    %rescale
    a = polyfit(isosbestic_photoBleachCorrected(skiptime:endTime1),signal_photoBleachCorrected(skiptime:endTime1),1);
    
    isosbestic_rescaled = polyval(a,isosbestic_photoBleachCorrected);
    
    %subtract rescaled isobestic
    signal_isosbesticCorrect = signal_photoBleachCorrected - isosbestic_rescaled;
    
    %calc zscore
    
    
    signal_DFoF_u = mean(signal_isosbesticCorrect(skiptime:endTime1));
    signal_DFoF_std = std(signal_isosbesticCorrect(skiptime:endTime1));
    
    signal_DFoF = (signal_isosbesticCorrect(:)' - signal_DFoF_u) / signal_DFoF_std;
    
    
    switch returnedDataType
        
        
        case 'isosbestic'
            signal_DFoF_u = mean(isosbestic_rescaled(skiptime:endTime1));
            signal_DFoF_std = std(isosbestic_rescaled(skiptime:endTime1));
            
            signal_DFoF = (isosbestic_rescaled(:)' - signal_DFoF_u) / signal_DFoF_std;
            
            
        case 'signal'
            signal_DFoF_u = mean(signal_photoBleachCorrected(skiptime:endTime1));
            signal_DFoF_std = std(signal_photoBleachCorrected(skiptime:endTime1));
            
            signal_DFoF = (signal_photoBleachCorrected(:)' - signal_DFoF_u) / signal_DFoF_std;
            
        case 'corrected'
            signal_DFoF = signal_DFoF;
            
    end
end
end

