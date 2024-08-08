function [spikes] = sm_LoadPhy(varargin)
% Load kilosort/phy clusters
%
% USAGE
%
%   [spikes] = bz_loadPhy(varargin)
%
% INPUTS
%
% basepath        -path to recording (where .dat are)
% kilosort_path   -path to kilosort* folder after phy curing (by default
%                   .\kilosort_\)
% getWaveforms    -logical (default=true) to load mean of raw waveform data
% saveMat         -logical (default=true) to save in buzcode format
% UID             -vector subset of UID's to load 
% fs              -scalar (default is taken from sessionInfo, otherwise
%                   30000). Sampling frequency. 
% nChannels       -scalar (default is taken from sessionInfo, otherwise
%                   32). Total number of channels.
% forceReload    -logical (default=false) to force loading from
%                   res/clu/spk files
% verbose        -logical (default=false)

%
% OUTPUTS
%
% spikes - cellinfo struct with the following fields
%   .sessionName    -name of recording file
%   .UID            -unique identifier for each neuron in a recording
%   .times          -cell array of timestamps (seconds) for each neuron
%   .spindices      -sorted vector of [spiketime UID], useful for 
%                           input to some functions and plotting rasters
%   .shankID        -shank ID that each neuron was recorded on
%   .amplitude:     amplitudes for cluster N
%   .maxWaveformCh  -channel # with largest amplitude spike for each neuron
%   .rawWaveform    -average waveform on maxWaveformCh (from raw .dat)
%   .filtWaveform   -average filtered waveform on maxWaveformCh
%   .region         -region ID for each neuron (especially important large scale, high density probes)
%   .numcells       -number of cells/UIDs
%
%  HISTORY:
%  9/2018  Manu Valero
%  10/2018 AntonioFR    
%  To do: Add a call to function for calculating cell features. 
%         Fix directory search for Linux

%% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'kilosort_path',[]); % probably this line only works in windows
%addParameter(p,'kilosort_path',pwd,@isstr);
addParameter(p,'getWaveforms',true,@islogical)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'UID',[],@isvector);
addParameter(p,'fs',30000,@isnumeric);
addParameter(p,'nChannels',32,@isnumeric);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'basename',[],@isstr);

parse(p,varargin{:});

basepath = p.Results.basepath;
kilosort_path = p.Results.kilosort_path;
getWave = p.Results.getWaveforms;
saveMat = p.Results.saveMat;
UID = p.Results.UID;
fs = p.Results.fs; % it will be overwritten if bz_getSessionInfo
nChannels = p.Results.nChannels; % it will be overwritten if bz_getSessionInfo
forceReload = p.Results.forceReload;
noPrompts = p.Results.noPrompts;
verbose = p.Results.verbose;

basename = p.Results.basename;


try [sessionInfo] = sm_getSessionInfo(basepath, 'noPrompts', noPrompts,'baseName',basename);
    fs = sessionInfo.rates.wideband;
    nChannels = sessionInfo.nChannels;
catch
    warning('SessionInfo file not found.');
end

if exist('sessionInfo') && exist([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'],'file') && forceReload == false
    disp('loading spikes from cellinfo file..')
    load([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'])
    
else
    
    if isempty(kilosort_path)
        fils = getAllExtFiles(basepath,'npy',1);
        kp = contains(fils,'spike_times');
        
        if sum(kp) ==1
        fils = fils{kp};
        kilosort_path = fileparts(fils);
        
        else
            error('kilosort_path must be specified')
        end
        
    end
    
    clufil = fullfile(kilosort_path, 'spike_clusters.npy');
    
    if ~exist(clufil)
       clufil = fullfile(kilosort_path, 'spike_templates.npy');
    end
    spike_cluster_index = readNPY(clufil);
    uclu = double(unique(spike_cluster_index));
    nclu = length(uclu);
    spike_times = readNPY(fullfile(kilosort_path, 'spike_times.npy'));
    % spike_amplitudes = readNPY(fullfile(kilosort_path, 'amplitudes.npy'));
    % spike_clusters = unique(spike_cluster_index);
    groupfil = fullfile(kilosort_path,'cluster_group.tsv');
    
    if exist(groupfil)
    cluster_group = tdfread(groupfil);
        cluster_group.group = [cluster_group.group;repmat('good ',sum(~ismember(uclu, cluster_group.cluster_id)),1)];
    cluster_group.cluster_id = [cluster_group.cluster_id ;uclu(~ismember(uclu, cluster_group.cluster_id))];

    
    else
        cluster_group.group = repmat('good',nclu,1);
        cluster_group.cluster_id = uclu;
    end
    try
        shanks = readNPY(fullfile(kilosort_path, 'channel_shanks.npy')); % done
    catch
        shanks = ones(nclu,1);
        warning('No shanks.npy file, assuming single shank!');
    end

    
    % if clusters not marked as good, default
    
    
    spikes = [];
    spikes.sessionName = sessionInfo.FileName;
    xml = LoadXml([basepath filesep sessionInfo.session.name '.xml']);
    jj = 1;
    for ii = 1:nclu
        if ~ (strcmpi(strtrim(cluster_group.group(ii,:)),'noise') | strcmpi(strtrim(cluster_group.group(ii,:)),'MUA'))
            ids = find(spike_cluster_index == cluster_group.cluster_id(ii)); % cluster id
            spikes.UID(jj) = cluster_group.cluster_id(ii);
            spikes.times{jj} = double(spike_times(ids))/fs; % cluster time
            spikes.ts{jj} = double(spike_times(ids)); % cluster time
            cluster_id = find(cluster_group.cluster_id == spikes.UID(jj));
          %  spikes.shankID(jj) = shanks(cluster_id);
            % spikes.amplitudes{jj} = double(spike_amplitudes(ids));
            jj = jj + 1;
        end
    end

    % get waveforms
    if getWave
    nPull = 1000; % number of spikes to pull out
    wfWin = 0.008; % Larger size of waveform windows for filterning
    filtFreq = 500;
    hpFilt = designfilt('highpassiir','FilterOrder',3, 'PassbandFrequency',filtFreq,'PassbandRipple',0.1, 'SampleRate',fs);
    % f = waitbar(0,'Getting waveforms...');
    wfWin = round((wfWin * fs)/2);
        for ii = 1 : size(spikes.times,2)
            spkTmp = spikes.ts{ii};
            if length(spkTmp) > nPull
                spkTmp = spkTmp(randperm(length(spkTmp)));
                spkTmp = spkTmp(1:nPull);
            end
            wf = [];
            for jj = 1 : length(spkTmp)
                if verbose
                    fprintf(' ** %3.i/%3.i for cluster %3.i/%3.i  \n',jj, length(spkTmp), ii, size(spikes.times,2));
                end
                 wf = cat(3,wf,bz_LoadBinary([basepath filesep sessionInfo.session.name '.dat'],'offset',spkTmp(jj) - (wfWin),...
                    'samples',(wfWin * 2)+1,'frequency',sessionInfo.rates.wideband,'nChannels',sessionInfo.nChannels));
            end
            
             tmp = mean(wf,3) - mean(wf(1,:,:),3);
            wf = mean(wf,3);
           
            [~,b] = min(min(tmp,[],1));
            shank = find(cellfun(@(a) any(a+1==b),{xml.AnatGrps.Channels}));
            spikes.shankID(ii) = shank;
            ch = xml.AnatGrps(shank).Channels;
            wf = wf(:,ch+1);
            %if isfield(sessionInfo,'badchannels')
            %    wf(:,ismember(sessionInfo.channels,sessionInfo.badchannels))=0;
            %end
            
            
            %asign ashank
            
             wfF = [];
            for jj = 1 : size(wf,2)          
                wfF(:,jj) = filtfilt(hpFilt,wf(:,jj) - mean(wf(:,jj)));
            end
            [~, maxCh] = max(abs(wfF(wfWin,:)));
            rawWaveform = detrend(wf(:,maxCh) - mean(wf(:,maxCh))); 
            filtWaveform = wfF(:,maxCh) - mean(wfF(:,maxCh));
            spikes.rawWaveform{ii} = rawWaveform(wfWin-(0.002*fs):wfWin+(0.002*fs)); % keep only +- 1ms of waveform
            spikes.filtWaveform{ii} = filtWaveform(wfWin-(0.002*fs):wfWin+(0.002*fs)); 
            spikes.maxWaveformCh(ii) = ch(maxCh);
            
            % figure;plot(spikes.filtWaveform{ii},'r');hold on;plot(spikes.rawWaveform{ii},'b');
            % waitbar(ii/size(spikes.times,2),f,'Pulling out waveforms...');
        end
        % close(f)
    end
end

% To match bz_GetSpikes
spikes.numcells = length(spikes.UID);
if ~isfield(spikes,'region') && isfield(spikes,'maxWaveformCh') && isfield(sessionInfo,'region')
    for cc = 1:spikes.numcells
        spikes.region{cc} = [sessionInfo.region{find(spikes.maxWaveformCh(cc)==sessionInfo.channels)} ''];
    end
end

% saveMat (only saving if no exclusions)
if saveMat
    save([basepath filesep sessionInfo.FileName '.spikes.cellinfo.mat'],'spikes');
end

% filter by UID input
if ~isempty(UID)
      [toRemove] = ~ismember(spikes.UID,UID);
    spikes.UID(toRemove) = [];
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.times{r} = [];
         % spikes.region{r} = [];
        end
    end
    spikes.times = removeEmptyCells(spikes.times);
    %spikes.region = removeEmptyCells(spikes.region);
    % spikes.cluID(toRemove) = [];
    spikes.shankID(toRemove) = [];
    
    if getWaveforms
    for r = 1:length(toRemove)
        if toRemove(r) == 1
         spikes.rawWaveform{r} = [];
        end
    end
    spikes.rawWaveform = removeEmptyCells(spikes.rawWaveform);
    spikes.maxWaveformCh(toRemove) = [];
    end
end

% Generate spindices matrics
for cc = 1:length(spikes.UID)
    groups{cc}=spikes.UID(cc).*ones(size(spikes.times{cc}));
end
if ~isempty(spikes.UID)
    alltimes = cat(1,spikes.times{:}); groups = cat(1,groups{:}); %from cell to array
    [alltimes,sortidx] = sort(alltimes); groups = groups(sortidx); %sort both
    spikes.spindices = [alltimes groups];
end

if ~isfield(spikes,'shankID')
    spikes.shankID = ones(spikes.numcells,1);


end
sm_WriteResClu(spikes,basepath)
% % Compute spike measures
% if ~isempty(spikes.UID) && getWave
%     for ii = 1:size(spikes.UID,2)
%         [~,tmp] = max(spikes.rawWaveform{ii}(1,41:end));
%         p2pWidth = tmp/fs; % peak (negative) to peak (second positive) duration 
%         [spikes.autocorr{ii}.xout,spikes.autocorr{ii}.r,spikes.autocorr{ii}.peakAutocorr] =...
%             autocorr_spikes(spikes.ts{ii},fs,26,1);
%     end 
% end