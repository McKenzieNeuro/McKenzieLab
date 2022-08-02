function [signal_DFoF,ts_data,ev_tims,ix,ts_PETH] = sm_PETH_DFoF(dirName,evFile,evType,varargin)
%% Plots the fiber photometry signal around events of interest
% Neural data should be stored in TDT format and read with TDT's SDK TDTbin2mat
% neural data should be stored in dirName
% Events are stored in a *.mat file under evFile.mat
% evFile.mat is a cell array with the first column = event name (str) and
% the second column = event time (seconds)

% INPUTS
% dirName =directory with TDT data and evFile
% evFile = *.mat file with events of interest (make using sm_labelTDTvideo)
% evType = string stored in evFile to index event times
% varargin =
%       photoBleachCorrection  = [linear, quadratic , etc] (default = quadratic)
%       savefig = [true,false] (default = false)
%       plotIntervals = [time before, time after] (s)
%       skiptime = time into recording to skip (s)
%       endTime = time into recording to consider (s)
%       streams = name of TDT data stream
%       isosbestic = name of isosbestic data stream
%       figureName = name of figure
%       returnedDataType = isosbestic correction (default = corrected)

% OUTPUTS
%  signal_DFoF = photobleach, isosbestic, zscore signal
%  ts_data = time stamps of data (s)
%  ev_tims = event times (s)
%
% USAGE
% sm_PETH_DFoF(pwd,'myEvents.mat','homeCage','plotIntervals',[25 25])
%
% DEPENDENCIES
% Clone Mckenzieneurolab github (https://github.com/McKenzieNeuro/McKenzieLab)
% TDT SDKs (https://www.tdt.com/support/matlab-sdk/)
%
% S. McKenzie, 08/23/2021


%%

% parse inputs


p = inputParser;


addParameter(p,'photoBleachCorrection','quadratic',@ischar)
addParameter(p,'savefig',false,@islogical)
addParameter(p,'skiptime',0,@isnumeric)
addParameter(p,'streams',{'x465A','x405A'},@iscell)
addParameter(p,'isosbestic','x405A',@ischar)
addParameter(p,'endTime',[],@isnumeric)
addParameter(p,'plotIntervals',[50 50],@isvector);
addParameter(p,'figureName','myPETH.fig',@ischar);
addParameter(p,'returnedDataType','corrected',@ischar);


parse(p,varargin{:})


photoBleachCorrection = p.Results.photoBleachCorrection;
savefig = p.Results.savefig;
skiptime = p.Results.skiptime;
streams = p.Results.streams;
isosbestic = p.Results.isosbestic;
endTime = p.Results.endTime;
plotIntervals = p.Results.plotIntervals;
figureName = p.Results.figureName;
returnedDataType = p.Results.returnedDataType;

%%
%load data
[signal_DFoF,ts_data,fs] = sm_getSignal_DFoF(dirName,'photoBleachCorrection',photoBleachCorrection,'skiptime',skiptime,'endTime',endTime, ...
    'isosbestic',isosbestic,'streams',streams,'returnedDataType',returnedDataType);


evs = load(evFile);


%%


%
% close all
% textLabel = {'Home','Home','Home','Home','Novel A','Novel B','Novel C'};
%
% signal_DFoF = nanconvn(signal_DFoF,k');
%
% signal_DFoF(1:5000) = nan;
% figure
%
% hold on
% plot([-100 max(ts_data)],[0 0],'k')
% plot([0 0; cell2mat(data(:,2)) cell2mat(data(:,2))]',[3*ones(size(data,1)+1,1) 4*ones(size(data,1)+1,1) ]','k')
% plot(ts_data,signal_DFoF,'r','linewidth',2)
%
% %h = text(cell2mat(data(:,2)),15*ones(size(data,1),1),data(:,1));
% %h = text([0 ;cell2mat(data(:,2))],4.5*ones(size(data,1)+1,1),textLabel);
% ylim([-5 8])
% xlim([-100 max(ts_data)])
%  set(h,'Rotation',45);
%
% xlabel ('Time (s)')
% ylabel ('Z-score')


%%
% get indices of events

for ii = 1:length(evType)
    
    evMatch = cellfun(@any,regexp(evs.data(:,1),evType{ii}));
    
    
    ev_tims{ii} = cell2mat(evs.data(evMatch,2));
    
    
    
    %data = TDTbin2mat(dirName);
    %ts_video = data.epocs.Cam1.onset;
    %ev_tims = ts_video(ev_tims);
    
    
    
    [ix{ii},early,late,ts_PETH] = sm_getIndicesAroundEvent(ev_tims{ii},plotIntervals(1),plotIntervals(2),fs,length(signal_DFoF));
    
    
    ix{ii} = ix{ii}(~(early |late),:); %exclude rows with time before/after recording
    
    
    
    
    if savefig
        figure
        
        
        d = signal_DFoF(ix{ii}); % sample with isosbestic substracted
        h = plot(ts_PETH,nanmean(d,1),'k');
        saveas(h,figureName);
        
    end
    
end


