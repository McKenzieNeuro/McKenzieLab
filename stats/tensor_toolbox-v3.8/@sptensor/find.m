function [subs,vals] = find(t)
%FIND Find subscripts of nonzero elements in a sptensor.
%
%   [SUBS,VALS] = FIND(T) returns the subscripts and corresponding
%   values of the nonzero (and known) elements of T.
%
%   Note that unlike the standard MATLAB find function for an array,
%   find does not return linear indices. Instead, it returns an M x N
%   array where M is the number of nonzero values and N = ndims(T).
%   Thus, I(k,:) specifies the subscript of value V(k).
%
%   If T is an incomplete spetensor, it returns the locations of all
%   nonzeros. If you want the locations of all known values, use T.subs 
%   and T.vals.
%
%   See also SPTENSOR, FIND.
% 
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


if isincomplete(t)
    idx = find(t.vals);
    subs = t.subs(idx,:);
    vals = t.vals(idx);
else
    subs = t.subs;
    vals = t.vals;
end