function [x,b] = histoc(y,n)
    
    if ~isempty(y)
    [x,b] = histc(y,n);
    x = x(:);
    
    else
        x = zeros(length(n),1);
        b = zeros(length(n),1);
    end
end