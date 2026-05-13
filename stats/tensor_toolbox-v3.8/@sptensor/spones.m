function X = spones(X)
%SPONES Replace sptensor elements with ones.
%
%   Y = SPONES(X) generates a tensor with the same structure as X,
%   but with ones for all values.
%
%   See also SPTENSOR, SPTENSOR/ONES.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


X.vals = ones(size(X.vals));
