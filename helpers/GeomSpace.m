function Y = GeomSpace(X,d,n)
%X = [min max]
% d = fraction to divide range(x)
% n = number of divisions

r = range(X);
s = r -(r./(d.^(1:n-2)));
Y = [X(1) s+X(1) X(2)];


end