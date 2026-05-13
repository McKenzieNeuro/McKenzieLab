function B = full(A,fillval)
%FULL Convert a sptensor to a (dense) tensor.
%
%   B = FULL(A) converts a sptensor A to a (dense) tensor B. If A is
%   incomplete, then the missing values are NaN's.
%
%   B = FULL(A,0) fills the missing values with zero for an incomplete
%   tensor.
%
%   See also SPTENSOR, TENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



% Extract the order and size of A
siz = size(A);

% Handle the completely empty (no size) case
if isempty(siz)
    B = tensor;
    return;
end

% Create a dense zero or nan tensor B that is the same size as A
if issparse(A) || (nargin > 1 && fillval == 0)
    B = tensor(zeros([siz,1,1]),siz);
else
    B = tensor(nan([siz,1,1]),siz);
end

if isempty(A.subs)
    return;
end

% Extract the linear indices of entries in A
idx = tt_sub2ind(siz,A.subs);

% Copy the values of A into B using linear indices
B(idx) = A.vals;

