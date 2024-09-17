function Y = nanprod(X,dim)
X1=X;
X1(isnan(X))=1;
Y=prod(X1,dim);
Y(all(isnan(X),dim))=nan;
end