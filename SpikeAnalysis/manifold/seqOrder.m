function [ord, rel, isMember] = seqOrder(P, kp, thr, relThresh)
% P         : nNeurons x nBins x nEvents (e.g., binnedPopRipple5z)
% kp        : events for this condition
% thr       : z threshold for "active" (replaces the 0.15 peak rule)
% relThresh : min split-half profile correlation to count as reliable
    ev   = find(kp);
    bins = 1:size(P,2);
    m    = nanmean(P(:,:,ev),3);                 % mean profile, nNeurons x nBins

    % split-half reliability of the temporal profile (avg a few splits to de-noise)
    rel = zeros(size(m,1),1); nRep = 10;
    for r = 1:nRep
        evp = ev(randperm(numel(ev))); h = floor(numel(evp)/2);
        m1 = nanmean(P(:,:,evp(1:h)),3); m2 = nanmean(P(:,:,evp(h+1:end)),3);
        for n = 1:size(m,1)
            rel(n) = rel(n) + corr(m1(n,:).', m2(n,:).', 'rows','pairwise')/nRep;
        end
    end

    % localized participation + center-of-mass order
    w        = max(m - thr, 0);                   % supra-threshold mass only
    hasField = sum(w,2) > 0;
    isMember = hasField & (rel >= relThresh);
    ord      = (w * bins.') ./ sum(w,2);          % weighted mean bin
%    ord(~isMember) = nan;                         % NaN -> pairwise does the intersection
end