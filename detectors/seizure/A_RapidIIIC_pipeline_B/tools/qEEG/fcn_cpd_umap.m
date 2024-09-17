function [icp_, P_, flags1, flags2] = fcn_cpd_umap(sdata, nn) % DATA{iFile, 3} = {tmp.Sdata, tmp.stimes, tmp.sfreqs};
%         if size(sdata{1}, 2)~=nn
%             keyboard
%         end
        

        % total power %
        %P = mean(pow2db(cell2mat(sdata(:, 1))+eps), 1);
        
        thr_cp = .1;
        P = mean(pow2db(sdata+eps), 1);
        
        % 2. Smooth % 
        P_ = smooth(P, 10,'sgolay')';

        % 3 Clip at 1000dB %
        P_(P_>25) = 25;
        P_(P_<-10) = -10;

        [icp, ~] = findchangepts(P_, 'Statistic', 'mean','MinThreshold',thr_cp*var(P_));
        flags1 = zeros(nn, 1);
        flags1(icp) = 1;
        flags1(1) = 1; % both ends %
        flags1(end) = -1; % both ends %

        icp_ = unique([icp, 1, nn]);  % add 2 ends %
        icp_center = floor((icp_(1:end-1)+icp_(2:end))/2);

        flags2 = zeros(nn, 1);
        flags2(icp_center) = 1;

    end