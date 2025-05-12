function y = nansum1(x,dim)
temp=x;
temp(isnan(x)) = 0;
if nargin == 1 % let sum figure out which dimension to work along
    y = nansum(temp);
    y(all(isnan(x),1))=nan;
else           % work along the explicitly given dimension
    y = nansum(x,dim);
    y(all(isnan(x),dim))=nan;
end
