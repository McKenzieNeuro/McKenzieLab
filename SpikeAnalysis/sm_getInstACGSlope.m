function [inst_ACG,beta] = sm_getInstACGSlope(corrmat,nlag,nTr)


base = true(size(corrmat)) ;

inst_ACG = nan(size(corrmat,1),nlag);

for i = 1:size(corrmat,1)
    if i+nlag<=size(corrmat,2)
        tmp = corrmat(i,i+1:i+nlag);
    else
        tmp = corrmat(i,i+1:end);
        tmp = nanPad(tmp,nlag);
    end
    inst_ACG(i,:) = tmp;
    
end


% now fit decay
warning off
f = @(b,x) b(1).*exp(-b(2).*x) + b(3);
beta0 = [.2 .1 0];
beta = nan(size(inst_ACG,1),3);
for i = 1:size(inst_ACG,1)-nTr
    X = 1:nlag;
    Y = nanmean(inst_ACG(i:i+nTr-1,:));
    
    tbl = table(X(:),Y(:));
    try
        mdl = fitnlm(tbl,f,beta0);
        beta(i,:) = mdl.Coefficients.Estimate;
    end
end
end



