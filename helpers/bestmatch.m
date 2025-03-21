function [match, idx] = bestmatch(a,b,varargin)
%finds the match and index of the closest scalar to b in a

a1=a(:);
b1=b(:)';
b=b1(ones(length(a1),1),:);
a=a1(:,ones(length(b1),1));

if ~isempty(varargin)
[~, idx]=max((a-b),[],2);    
else
[~, idx]=min(abs(a-b),[],2);
end
match=b1(idx);

% idx(isnan(b1))=nan;
% match(isnan(b1))=nan;
% 
% idx(isinf(b1))=nan;
% match(isinf(b1))=nan;


match=match(:);

if isempty(match)
    match=nan;
    idx=nan;
end

end