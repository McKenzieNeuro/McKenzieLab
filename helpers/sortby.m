function [Y,b] =  sortby(X,Z)

[~,b] = sort(Z);
Y = X(b,:);

%Y = Y./repmat(max(Y,[],2),1,size(Y,2));
end

