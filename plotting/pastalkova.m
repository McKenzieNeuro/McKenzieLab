function [Y,b] = pastalkova(X)
[~,b]=max(X,[],2);
[~,b]=sort(b,'descend');
Y=X(b,:);

%Y = Y./repmat(max(Y,[],2),1,size(Y,2));

%Y = nanzscore(Y,[],2);
end