function [isCPs, isCPcenters] = fcn_cpd(sdata, alpha_cpd)

    % 1. Mean power in dB
    P = mean(pow2db(cell2mat(sdata(:, 1))+eps), 1);

    % 2. smoothing 
    P = smooth(P, 10,'sgolay')';

    % 3. clipping at [-15 25dB]  
    P(P> 25) =  25;
    P(P<-15) = -15;

    % 4. CP detection with threshold on mean power 
    icp = findchangepts(P, 'Statistic', 'mean', 'MinThreshold', alpha_cpd*var(P));
    
    % 5. Outputs
    nn = size(sdata{1}, 2);
    
    isCPs = zeros(nn, 1);
    isCPs(icp) = 1;
    isCPs(1) = 1; isCPs(end) = -1;  % deal with both ends

    idx_cp = unique([icp, 1, nn]);      
    idx_cpcenter = floor((idx_cp(1:end-1)+idx_cp(2:end))/2);

    isCPcenters = zeros(nn, 1);
    isCPcenters(idx_cpcenter) = 1;
end
