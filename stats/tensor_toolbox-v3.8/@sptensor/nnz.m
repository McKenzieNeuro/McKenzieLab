function a = nnz(t)
%NNZ Number of values in sptensor.
%
%   NNZ(T) is the number of nonzero elements in T if T is sparse;
%   otherwise it is the number of elements in T (if it is incomplete).
%
%   See also SPTENSOR, SPTENSOR/FIND.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if isempty(t.subs)
    a = 0;
else
    a = size(t.subs,1);
end
