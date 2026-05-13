function z = eq(x,y)
%EQ Equal (==) for sptensors.
%
%   A == B compares the elements of A and B for equality. The arguments can
%   be a pair of sptensors, an sptensor and a tensor, or an sptensor and a
%   scalar.  Regardless, the result is always returned as a sptensor.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%% Observations for sparse matrix case.
% The result of a == 5 is sparse.
% The result of a == 0 is sparse.
% The result of a == full(a) is sparse.

%% Case 1: One argument is a scalar
if isscalar(y)
    if y == 0        
        z = ~x;
    else
        idx = (x.vals == y);
        z = sptensor(x.subs(idx,:),true,size(x));
    end   
    return;
end

% Call back with the arguments reversed.
if isscalar(x)
    z = eq(y,x);
    return;
end

%% Case 2: Both x and y are tensors of some sort
% Check that the sizes match
if ~isequal(x.size,y.size)
    error('Size mismatch');
end

% Case 2a: Two sptensors
if isa(x,'sptensor') && isa(y,'sptensor')

    % Find where their zeros/known values intersect
    if x.type ~= y.type
        zzerosubs = [];
    else
        xzerosubs = setdiff(allsubs(x),x.subs,'rows');
        yzerosubs = setdiff(allsubs(y),y.subs,'rows');
        zzerosubs = intersect(xzerosubs,yzerosubs,'rows');
    end

    % find where their nonzeros/known values intersect 
    [nzsubs,ix,iy] = intersect(x.subs,y.subs,'rows');
    znzsubs = nzsubs(x.vals(ix) == y.vals(iy),:);

    % Build z
    z = sptensor([zzerosubs;znzsubs],true,x.size);
    
    return;

end

% Case 2b: One dense tensor
if isa(y,'tensor')

    if issparse(x) == 1
        targetval = 0;
    else
        targetval = nan;
    end

    % Find where their zeros intersect
    yzerosubs = find(y == targetval);
    zzerosubs = yzerosubs(extract(x,yzerosubs) == 0,:);

    % Find where their nonzeros intersect 
    yvals = y(x.subs);
    znzsubs = x.subs(yvals == x.vals,:);
    
    % Build z
    z = sptensor([zzerosubs;znzsubs],true,x.size);
    
    return;
    
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');
