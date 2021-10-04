
function Y = nanPad(X,len,varargin)
flip = false;
if ~isempty(varargin)
    flip = true;
end
Y=X;
if size(X,2)<=len
    Y=nan(size(X,1),len);
    
    if flip
        Y(:,end-size(X,2)+1:end)=X;
    else
        
        Y(:,1:size(X,2))=X;
    end
else
    %    Y=X(1:len);
    
end

if isempty(X)
    
    if flip
        Y = nan(len);
    else
        
        Y = nan(1,len);
    end
end
end