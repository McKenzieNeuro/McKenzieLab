function [hist_y, n_y,hist_ySEM]= avghist(x,y,edges,varargin)
%returns a histogram of the average values of y based upon binning the data
%in x
if isempty(x)
       hist_y = nan(length(edges),1)';
        hist_ySEM = nan(length(edges),1)';
        n_y = nan(length(edges),1)';
    return;
else
    if ~isempty(varargin)
        
        fnc = varargin{1};
    else
        fnc = @nanmean;
    end
    edges=edges(:);
%     if length(y)<length(x)
%         x=x(1:length(y));
%         
%     end
    y(isinf(y)) = nan;
    
    [n,bin]=histc(x,edges);
    
    if all(n==0)
       hist_y = nan(length(edges),1)';
        hist_ySEM = nan(length(edges),1)';
        n_y = nan(length(edges),1)';
        return;
    end
    bin=bin(:);
    y=y(:);
    kp = ~isnan(y) & ~isnan(bin) & bin~=0;
    if any(kp)
    hist_y=accumarray(bin(kp),y(kp),[length(edges) 1],fnc,nan)';
    hist_ySEM=accumarray(bin(kp),y(kp),[length(edges) 1],@SEM,nan)';
    n_y=accumarray(bin(kp),y(kp),[length(edges) 1],@(a) nnz(~isnan(a)),nan)';
    else
        hist_y = nan(length(edges),1)';
        hist_ySEM = nan(length(edges),1)';
        n_y = nan(length(edges),1)';
    end
        
end
end