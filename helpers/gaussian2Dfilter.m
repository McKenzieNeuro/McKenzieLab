function h = gaussian2Dfilter(siz, sigma)
% GAUSSIAN2DFILTER Creates a 2-D Gaussian lowpass filter.
%   H = GAUSSIAN2DFILTER(HSIZE,SIGMA) returns a Gaussian lowpass filter of
%   size HSIZE with standard deviation SIGMA (positive). HSIZE can be a
%   vector specifying the number of rows and columns in H or a scalar, in
%   which case H is a square matrix. SIGMA can be a vector specifying the
%   standard deviation of the rows and columns separately, or a scalar, in
%   which case H will be rotationally symmetric. If SIGMA is a scalar, this
%   function behaves exactly the same as FSPECIAL called using the
%   'gaussian' input parameter.
%
%   The default HSIZE is [3 3], the default SIGMA is [0.5 0.5].
%
%   See also FSPECIAL, GAUSSIANNDFILTER

    if(exist('siz','var')~=1); siz = [3 3]; end
    if(exist('sigma','var')~=1); sigma = 0.5; end

    if(isscalar(siz)); siz = [siz siz]; end
    if(isscalar(sigma)); sigma = [sigma sigma]; end
    
    siz   = (siz-1)/2;

    [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
    arg   = -((x.*x)/(2*sigma(2)*sigma(2)) + (y.*y)/(2*sigma(1)*sigma(1)));

    h     = exp(arg);
    h(h<eps*max(h(:))) = 0;

    sumh = sum(h(:));
    if(sumh~=0); h  = h/sumh; end
end