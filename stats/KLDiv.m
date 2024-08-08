function dist=KLDiv(P,Q)
%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1



if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end

% if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
%    error('the inputs contain non-finite values!') 
% end

% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./nansum(Q);
    P = P ./repmat((nansum(P,2)),[1 round(nansum(P,2))]);
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isinf(temp))=nan;% resolving the case when P(i)==0
    dist = nansum(temp,2);
    
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(nansum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(nansum(P,2),[1 size(P,2)]);
    temp =  P.*log(P./Q);
    temp(isinf(temp))=nan; % resolving the case when P(i)==0
    dist = nansum(temp,2);
end


