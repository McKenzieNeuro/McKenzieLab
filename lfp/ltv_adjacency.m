function A = ltv_adjacency(tsData)
% Ex.  A = ltv_adjacency(tsData)

%% NOTES
% - creates an adjacency matrix for the time series data: tsData
% - Uses matlab's least sqaures algorithm for solving a linear sysem ("\")
% - assumes x(0) was removed

% INPUTS:
% - tsData: time series (time x channels)

% OUTPUTS:
% - A: adjacency matrix for the time series

% REFERENCES: 
%   Li, Adam, Sara Inati, Kareem Zaghloul, and Sridevi Sarma. 2017. 
%       ìFragility in Epileptic Networks: The Epileptogenic Zone.î 
%       Proceedings of the American Control Conference: 2817ñ22.

% updated 03.29.21 - NF

%%

% --- SET PARAMS ---
T = size(tsData,1)-1; % number of time points -1
N = size(tsData,2); % number of channels

% --- BUILD H ---
H = kron(repmat(speye(N),T,1),ones(1,N)); % set up system of equations
r = repmat((1:T*N)',1,N); % rows
c = repmat(reshape(1:N^2,[],N)',(N*T)/N,1); % columns
H(sub2ind(size(H),r,c)) = repelem(tsData(1:T,:),N,1); % input data into H

% --- BUILD B ---
b = tsData(2:end,:)'; % use t(2) through end
b = b(:); % format
b = sparse(b);

% --- SOLVE FOR A ---
A = H\b; % vectorized adjacency matrix
A = reshape(A,N,N)'; % format
A = full(A);