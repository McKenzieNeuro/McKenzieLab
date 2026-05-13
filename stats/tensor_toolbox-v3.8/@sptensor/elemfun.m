function a = elemfun(a,fun)
%ELEMFUN Manipulate the elements of a sptensor.
%
%   X = ELEMFUN(X,@FUN) modifies the elements of X according to the
%   function @FUN which should take and array and output an equally
%   sized array.
%
%   Examples
%   X = sptenrand([10,10,10],10);
%   X = elemfun(X,@sqrt) %<-- square root of every entry
%   X = elemfun(X, @(x) x+1) %<-- increase every entry by 1
%   X = elemfun(X, @(x) x ~= 0) %<-- change every nonzero to be 1
%
%   See also SPTENSOR, SPFUN.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>




if ~isa(a,'sptensor')
    error('First argument must be a sptensor.');
end

a.vals = fun(a.vals);

if isincomplete(a) % Nothing left to do for an incomplete tensor
    return
end

idx = find(a.vals);
if isempty(idx)
    a.vals = [];
    a.subs = [];
else
    a.vals = a.vals(idx);
    a.subs = a.subs(idx,:);
end
