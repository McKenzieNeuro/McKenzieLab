function conPETH = sm_make_PETH_per_context(spikes,pulseInfo,ev)
% spikes  = buzcode struct with spike info
%pulseInfo = struct with pulse info
%ev = struct with context transitions (  ev = LoadEvents('amplifier_analogin_auxiliary_int16.evt.ctx');

  [binnedPop,bin_times]=populationMatrix(spikes,.5,.5,100,pulseInfo.time(:,1),'zscore',false);

  maxT = max(cellfun(@max,spikes.times));

  
con_ep(:,1) = sort(ev.time(~contains(ev.description,'HC') & contains(ev.description,'in')));

con_ep(:,2) = sort(ev.time(~contains(ev.description,'HC') & contains(ev.description,'out')));
  HC_ep = excludeEpochs([0 maxT],con_ep);
  
  eps = [con_ep;HC_ep];
  
  eps = sortrows(eps);
  

  
  [~,b] = InIntervals(pulseInfo.time(:,1),eps);

  for i = 1:max(b)
      
      for j = 1:size(binnedPop,1)
          
          conPETH(j,:,i) = nanmean(binnedPop(j,:,b==i),3);
      end
  end
  
end