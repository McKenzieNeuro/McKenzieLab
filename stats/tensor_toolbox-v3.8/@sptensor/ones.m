function t = ones(t)
%ONES Replace values of sptensor with ones.
%
%   S = ONES(T) generates a sptensor with the same sparsity
%   structure as T, but with ones in the nonzero/known value positions.
%
%   See also SPTENSOR, SPONES.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



t.vals = ones(size(t.vals));
