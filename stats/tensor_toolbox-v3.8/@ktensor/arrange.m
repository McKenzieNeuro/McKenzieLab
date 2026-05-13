function [X,p] = arrange(X,foo)
%ARRANGE Arranges the rank-1 components of a ktensor.
%
%   X = ARRANGE(X) normalizes the columns of the factor matrices and then
%   sorts the ktensor components by magnitude, greatest to least.
%
%   X = ARRANGE(X,N) absorbs the weights into the Nth factor matrix instead
%   of lambda. 
%
%   X = ARRANGE(X,P) rearranges the components of X according to the
%   permutation P. P should be a permutation of 1 to NCOMPONENTS(X).
%
%   [X,P] = ARRANGE(...) returns also the permutation of the components.
%
%   Examples
%   K = ktensor([3; 2], rand(4,2), rand(5,2), rand(3,2))
%   arrange(K) %<--Normalize and sort according to weight vector
%   arrange(K,[2, 1]) %<--Order components according to permutation
%
%   See also KTENSOR, NCOMPONENTS, NORMALIZE.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%% Just rearrange and return if second argument is a permutation
if exist('foo','var') && (length(foo) > 1)
    p = foo; % save the permutation
    X.lambda = X.lambda(foo);
    for i = 1 : ndims(X)
        X.u{i} = X.u{i}(:,foo);
    end   
    return;
end

%% Ensure that matrices are normalized
X = normalize(X);

%% Sort
[X.lambda, p] = sort(X.lambda, 1, 'descend');
for i = 1 : ndims(X)
    X.u{i} = X.u{i}(:,p);
end


%% Absorb the weight into one factor, if requested
if exist('foo','var')
    r = length(X.lambda);
    X.u{foo} = full(X.u{foo} * spdiags(X.lambda,0,r,r));
    X.lambda = ones(size(X.lambda));
end

