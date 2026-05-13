function tf = isincomplete(X)
%ISINCOMPLETE True if X is an incomplete tensor (versus sparse).
%
%   ISINCOMPLETE(X) returns logical 1 (true) if X is an incomplete tensor.
%
%   See also SPTENSOR
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>
tf = strcmp(X.type,'incomplete');
return