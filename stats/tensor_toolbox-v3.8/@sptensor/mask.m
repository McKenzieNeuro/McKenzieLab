function vals = mask(X,W)
%MASK Extract values as specified by a mask tensor.
%
%   V = MASK(X,W) extracts the values in X that correspond to nonzero/known
%   values in the mask tensor W.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

% Error check
if any(size(W) > size(X))
    error('Mask cannot be bigger than the data tensor')
end

% Extract locations of nonzeros/values in W
wsubs = find(W);

% Find which values in the mask match values in X
[tf,loc] = ismember(wsubs,X.subs,'rows');

% Assemble the final array
nvals = size(wsubs,1);
if issparse(X)
    vals = zeros(nvals,1);
else
    vals = nan(nvals,1);
end
vals(tf) = X.vals(loc(tf));
