function [binnedPop,bin_times]=populationMatrix(spikes,secbefore,secafter,nbins,ts_onset,varargin)
%self is an object with the current epoch and all cells and ts
%secbefore the time before the event occurred
%sec after the time(s) after the event occurred
%binSize is size of each bin in sec
%row=cell %column=time % z=trial

normalize = false;
kp = true(length(spikes.times),1);
kernel =[];
for i = 1:2:length(varargin)
    
    switch varargin{i}
        case 'kernel'
            kernel = varargin{i+1};
        case 'kp'
            kp = varargin{i+1};
        case 'zscore'
            normalize = varargin{i+1};
            
    end
end

spikes.times = spikes.times(kp);
nCel = length(spikes.times);

spikes_str = [];

if isempty(ts_onset)
    binnedPop=nan(nCel,nbins,1);
    return
end

inbin=nan(length(ts_onset),nbins,nCel);
tot_duration = secbefore + secafter;
bin_duration = tot_duration/nbins;

bin_times=-secbefore:bin_duration:secafter;

ep=[ts_onset-secbefore ts_onset+secafter];

for u=1:nCel
    
    spks=spikes.times{u};
    
    
    adjustedSpikes=cellfun(@(a) spks(spks>a(1) & spks<a(2)),num2cell(ep,2) , 'uni',0);
    adjustedSpikes=cellfun(@(a,b) a-b, adjustedSpikes, num2cell(ts_onset) , 'uni',0);
    spkhist = cell2mat(cellfun(@(a) reshape(histoc(a,bin_times),1,[]),adjustedSpikes, 'UniformOutput', false));
    if ~isempty(kernel)
        spkhist=cell2mat(cellfun(@(a) imfilter(reshape(histoc(a,bin_times),1,[]),kernel,'same'),adjustedSpikes, 'UniformOutput', false));
    end
    % spkhist = spkhist(:,1:end-1)./bin_duration;%turn into rate
    t = spkhist(:,1:end-1)./bin_duration;
    
    
    
    %    inbin(:,:,u) =  inbin(:,:,u) + min(linearize( inbin(:,:,u)));
    
    
    if normalize
        inbin(:,:,u)= (t-nanmean(t(:)))/nanstd(t(:));%turn into z-score
    else
        
        inbin(:,:,u) = t;
    end
end


binnedPop=permute(inbin,[3 2 1]);
end
