% save all feature files
topDir = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/bids';
fils = getAllExtFiles(topDir,'mat',1);

kp = contains(fils,'eeg_features.mat');

fils  = fils(kp);
master_feats = '/carc/scratch/projects/bshuttleworth2016391/mckenzie2016183/data/HarvardEEGDatabase/bids/allfeat.mat';
save(master_feats,'fils')
% 
% for i = 1:length(fils)
%     
%     v(i) = load(fils{i});
%     ok = whos('v');
%     if ok.bytes/1e9>4
%         save(master_feats,'v','-v7.3')
%         disp(['saved: ' master_feats])
%         break
%     end
%  
%     
% end



