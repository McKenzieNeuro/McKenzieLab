function Y = nanmean1(X,varargin)
X(isinf(X)) = nan;
if size(X,1) == 1
   
  X=X(:);
end

if isempty(varargin)
    
    dim = 1;
else
    dim = varargin{1};
end

Y = nansum(X,dim)./nansum(~isnan(X),dim);

end