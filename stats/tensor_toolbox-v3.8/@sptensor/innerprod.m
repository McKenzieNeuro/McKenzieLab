function res = innerprod(X,Y)
%INNERPROD Efficient inner product with a sparse tensor.
%
%   R = INNERPROD(X,Y) efficiently computes the inner product between
%   two tensors X and Y.  If Y is a tensor or sptensor, the inner
%   product is computed directly and the computational complexity is
%   O(min(nnz(X),nnz(Y))). If Y is a ktensor or a ttensor, the
%   inner product method for that type of tensor is called.
%
%   See also SPTENSOR, KTENSOR/INNERPROD, TTENSOR/INNERPROD.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

if isa(X,'sptensor') && isincomplete(X)
    error('Cannot handle incomplete tensors');
end
if isa(Y,'sptensor') && isincomplete(Y)
    error('Cannot handle incomplete tensors');
end

% X is a sptensor
if nnz(X) == 0 %There are no nonzero terms in X.
    res = 0;
    return
end

if isa(Y,'sptensor')
    if ~isequal(size(X),size(Y))
        error('X and Y must be the same size.');
    end
    if nnz(Y) == 0 %There are no nonzero terms in Y.
        res = 0;
        return
    end
    if nnz(X) < nnz(Y);
        [SX,VX] = find(X);
        VY = extract(Y,SX);   %<-----VY = Y(SX);
    else
        [SY,VY] = find(Y);
        VX = extract(X,SY);   %<-----VX = X(SY);
    end
    res = VY'*VX;
    return;
    
elseif isa(Y,'tensor')
    if ~isequal(size(X),size(Y))
        error('X and Y must be the same size.');
    end
    [SX,VX] = find(X);
    VY = Y(SX,'extract');
    res = VY'*VX;
    return;
    
elseif isa(Y,'ktensor') || isa(Y,'ttensor')
    % Reverse arguments to call ktensor/ttensor implementation
    res = innerprod(Y,X);
    return;
    
else
    error(['Inner product not available for class ' class(Y)]);
        
end

