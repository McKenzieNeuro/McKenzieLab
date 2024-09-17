function D=fcn_distChiSq(X,Y)
    mm=size(X,1);nn=size(Y,1);
    mOnes=ones(1,mm);
    D=zeros(mm,nn);
    for j=1:nn
        yi=Y(j,:);yiRep=yi(mOnes,:);ss=yiRep+X;dd=yiRep-X;
        D(:,j)=sum(dd.^2./(ss+eps),2);
    end
    D=D/2;
end
