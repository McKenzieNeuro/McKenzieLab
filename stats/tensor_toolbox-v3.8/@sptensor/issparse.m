function tf = issparse(X)
%ISSPARSE True if X is a sparse tensor (versus incomplete).
%
%   ISSPARSE(X) returns logical 1 (true) if X is a sparse tensor.
%
%   See also SPTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>
tf = strcmp(X.type,'sparse');
return