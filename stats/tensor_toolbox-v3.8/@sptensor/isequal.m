function z = isequal(x,y)
%ISEQUAL Compare sptensors for equality.
%
%   ISEQUAL(A,B) compares the sptensors A and B for equality.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Observations for sparse matrix case.
% The result of isequal(a,full(a)) is true!

%%
if ~isequal(x.size,y.size) 
    z = false;
elseif isa(x,'sptensor') && isa(y,'sptensor') 
    if ~isequal(x.type,y.type) 
        z = false;
    else
        diff = sptensor([x.subs; y.subs], [x.vals; -y.vals], x.size);
        z = isempty(diff.vals);
    end
elseif isa(y,'tensor')
    z = isequal(full(x),y);
else
    z = false;
end
