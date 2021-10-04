function [x,y] = FixPosition(x,y,ts)
% root = root.FixPos;
%
% Takes all 0,0s and large jumps in position (greater than
% jitter_threshold=15 pixels/sample) that persist for less than
% max_allowed_flips (default = 5) and linearly interpolates the missing
% data. Smooths conservatively afterward, as well (convolution with a gaussian, standard
% deviation = 2 samples).
%
% andrew december 2009
% update andrew june 15 2011

import CMBHOME.Utils.*


    jitter_threshold = 10*10; % pixels in 10 cm; one-sample change in distance that qualifies as bad vector
    
 
    
    y(x==0) = NaN;
    x(x==0) = NaN;
    x(y==0) = NaN;
    y(y==0) = NaN;
      

        max_allowed_flips = 1; % samples
    
    
    flips = findOnsetsAndOffsets(isnan(x));
    
    flips(:,2) = flips(:,2)+1;
    
    flips = cat(1, 1, flips(:), find([0; sqrt(diff(x).^2 + diff(y).^2)]>jitter_threshold), length(x));
    
    flips = sort(unique(flips)); % indeces of NaNs or jumps
    
    flips = [flips(1:end-1), flips(2:end)];  % epochs formation
    
    flips(flips(:,2)-flips(:,1)>max_allowed_flips,:) = [];
    
    flips(:,2) = flips(:,2)-1; % adjust for diff shift
    
    flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps
   
    flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices
    
    x([flips{:}]) = []; % remove samples in ts and x
    y([flips{:}]) = []; % remove samples in ts and x
    ts1=ts;
    ts1([flips{:}]) = [];
    
    x = interp1(ts1, x, ts);
    y = interp1(ts1, y, ts);

    x = ndnanfilter(x, normpdf(-6:6, 0, 2)', [], 1, {}, {}, 1); % conv with gaussian and ignore NaNs
    y = ndnanfilter(y, normpdf(-6:6, 0, 2)', [], 1, {}, {}, 1);
    

end
