function [feature, headers] = fcn_computeFeatures_powers_complex_v2(seg, ekg, Fs)
    
    % seg: bipolar %
    set0 = [min(seg, [], 2), max(seg, [], 2), mean(seg, 2), std(seg, [], 2),  median(seg, 2), ...
           (prctile(seg, 75, 2) - prctile(seg, 25, 2)), (prctile(seg, 97.5, 2) - prctile(seg, 2.5, 2)),...
            prctile(seg, 2.5, 2), prctile(seg, 25, 2), prctile(seg, 75, 2), prctile(seg, 97.5, 2)];
    

    % Features %
    P1 = bandpower(seg', Fs, [.5  4])';
    P2 = bandpower(seg', Fs, [4   7])';
    P3 = bandpower(seg', Fs, [8  15])';
    P4 = bandpower(seg', Fs, [16 31])';
    P5 = bandpower(seg', Fs, [32 64])';
    P  = (P1+P2+P3+P4+P5);

    % set1: powers %
    set1 = [pow2db(P1+eps), pow2db(P2+eps),  pow2db(P3+eps),  pow2db(P4+eps),  pow2db(P5+eps), pow2db(P+eps), ...
            P1./(P+eps), P2./(P+eps),  P3./(P+eps),  P4./(P+eps),  P5./(P+eps), ...
            P1./(P2+eps),  P1./(P3+eps), P1./(P4+eps), P1./(P5+eps)...
            P2./(P1+eps),  P2./(P3+eps), P2./(P4+eps), P2./(P5+eps), ...
            P3./(P1+eps),  P3./(P2+eps), P3./(P4+eps), P3./(P5+eps),...
            P4./(P1+eps),  P4./(P2+eps), P4./(P3+eps), P4./(P5+eps),...
            P5./(P1+eps),  P5./(P2+eps), P5./(P3+eps), P5./(P4+eps)];
        
    % Zero-crossings and line-length before z-score %
    entr =NaN(size(seg, 1), 1);
    ZCC = NaN(size(seg, 1), 1);
    LL = ZCC;
    for i = 1:size(seg, 1)
        x = seg(i, :);
        entr(i) = entropy(x);
        ZCC(i) = length(fcn_getZeroCrossings(x));
        LL(i)  = nanmean(abs(diff(x)));
    end
    set2 = [entr, ZCC, LL];
    
    seg_n = (seg - repmat(nanmean(seg, 2), 1, size(seg, 2)))./repmat(nanstd(seg, [], 2), 1, size(seg, 2));
    entr_n =NaN(size(seg, 1), 1);
    ZCC_n = NaN(size(seg, 1), 1);
    LL_n = NaN(size(seg, 1), 1);
    for i = 1:size(seg, 1)
        x = seg_n(i, :);
        entr_n(i) = entropy(x);
        ZCC_n(i) = length(fcn_getZeroCrossings(x));
        LL_n(i)  = nanmean(abs(diff(x)));
    end
    % featrue vec %
    set2 = [set2, entr_n, ZCC_n, LL_n];
    
    % xcorr - mean of top3 %
    Xcorr = NaN(size(seg, 1));
    for i = 1:(size(seg, 1)-1)
        x1 = seg(i, :);
        
        for j = (i+1):size(seg, 1)
           x2 = seg(j, :); 
           Xcorr(i, j) = max(abs(xcorr(x1, x2, 'coeff')));
           Xcorr(j, i) = Xcorr(i, j);
        end
    end
    
    set3 = NaN(size(seg, 1), 3);
    for i = 1:size(Xcorr, 1)
        x = Xcorr(i, :);
        
        [x1, ~] = sort(x(~isnan(x)), 'descend');
        [x2, ~] = sort(x(~isnan(x)), 'ascend');
        set3(i, 1) = mean(x1(1:min(length(x1),3)));
        set3(i, 2) = mean(x2(1:min(length(x2),3)));
        set3(i, 3) = max(abs(xcorr(seg(i, :), ekg, 'coeff')));
    end
    set3(isnan(set3)) = 0;

    feature =  [set0 set1 set2 set3];
    headers = {'V-min', 'V-max', 'V-mean', 'V-std', 'V-median', 'V-iqr', 'V-ci', 'V-p2.5', 'V-p25', 'V-p75', 'V-p97.5', ...
               'P-delta', 'P-theta', 'P-alpha', 'P-beta', 'P-gamma', 'P-total', ...
               'P-delta/total', 'P-theta/total', 'P-alpha/total', 'P-beta/total', 'P-gama/total',...
               'P-delta/theta', 'P-delta/alpha', 'P-delta/beta', 'P-delta/gamma', ...
               'P-theta/delta', 'P-theta/alpha', 'P-theta/beta', 'P-theta/gamma', ...
               'P-alpha/delta', 'P-alpha/theta', 'P-alpha/beta', 'P-alpha/gamma',...
               'P-beta/delta',  'P-beta/theta',  'P-beta/alpha', 'P-beta/gamma',...
               'P-gamma/delta', 'P-gamma/theta', 'P-gamma/alpha','P-gamma/beta',...
               'entropy', 'zero-crossing', 'line-length',...
               'entropy-norm', 'zero-crossing-norm', 'line-length-norm',...
               'xcorr-max', 'xcorr-min', 'xcorr-ekg'};
           
 
end