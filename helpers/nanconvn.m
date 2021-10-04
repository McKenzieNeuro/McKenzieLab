function c = nanconvn(a, k, varargin)
% NANCONV Convolution in ND ignoring NaNs.
%   C = NANCONVN(A, K) convolves A and K, correcting for any NaN values
%   in the input vector A. The result is the same size as A (as though you
%   called 'conv', 'conv2', or 'convn' with the 'same' shape).
%
%   C = NANCONV(A, K, 'params', ...) specifies one or more of the following:
%     'edge'     - Apply edge correction to the output.
%     'noedge'   - Do not apply edge correction to the output (default).
%     'nanout'   - The result C should have NaNs in the same places as A.
%     'nonanout' - The result C should have ignored NaNs removed (default).
%                  Even with this option, C will have NaN values where the
%                  number of consecutive NaNs is too large to ignore.
%     'Nd'       - Treat the input vectors as ND matrices (default).
%     '1d'       - Treat the input vectors as 1D vectors.
%                  This option only matters if 'a' or 'k' is a row vector,
%                  and the other is a column vector. Otherwise, this
%                  option has no effect.
%
% See also conv, conv2, convn
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2012, Benjamin Kraus
% $Id$

% Process input arguments
for arg = 1:nargin-2
    switch lower(varargin{arg})
        case 'edge'; edge = true; % Apply edge correction
        case 'noedge'; edge = false; % Do not apply edge correction
        case {'same','full','valid'}; shape = varargin{arg}; % Specify shape
        case 'nanout'; nanout = true; % Include original NaNs in the output.
        case 'nonanout'; nanout = false; % Do not include NaNs in the output.
        case {'nd','isnd'}; is1D = false; % Treat the input as ND
        case {'1d','is1d'}; is1D = true; % Treat the input as 1D
    end
end

% Apply default options when necessary.
if(exist('edge','var')~=1); edge = true; end
if(exist('nanout','var')~=1); nanout = false; end
if(exist('is1D','var')~=1); is1D = false; end
if(exist('shape','var')~=1); shape = 'same';
elseif(~strcmp(shape,'same'))
    error([mfilename ':NotImplemented'],'Shape ''%s'' not implemented',shape);
end

% Get the size of 'a' for use later.
sza = size(a);

% If 1D, then convert them both to columns.
% This modification only matters if 'a' or 'k' is a row vector, and the
% other is a column vector. Otherwise, this argument has no effect.
if(is1D);
    if(~isvector(a) || ~isvector(k))
        error('MATLAB:conv:AorBNotVector','A and B must be vectors.');
    end
    a = a(:); k = k(:);
end

% Flat function for comparison.
o = ones(size(a));

% Flat function with NaNs for comparison.
on = ones(size(a));

% Find all the NaNs in the input.
n = isnan(a);

% Replace NaNs with zero, both in 'a' and 'on'.
a(n) = 0;
on(n) = 0;

% Check that the filter does not have NaNs.
if(any(isnan(k)));
    error([mfilename ':NaNinFilter'],'Filter (k) contains NaN values.');
end

% Calculate what a 'flat' function looks like after convolution.
if(any(n(:)) || edge)
    flat = convn(on,k,shape);
else flat = o;
end

% The line above will automatically include a correction for edge effects,
% so remove that correction if the user does not want it.
if(any(n(:)) && ~edge); flat = flat./convn(o,k,shape); end

% Do the actual convolution
c = convn(a,k,shape)./flat;

% If requested, replace output values with NaNs corresponding to input.
if(nanout); c(n) = NaN; end

% If 1D, convert back to the original shape.
if(is1D && sza(1) == 1); c = c.'; end

end