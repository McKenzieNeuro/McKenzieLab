function [P,Uinit,output] = cp_nmu(X,R,varargin)
%CP_NMU Compute nonnegative CP with multiplicative updates.
%
%   P = CP_NMU(X,R) computes an estimate of the best rank-R CP
%   model of a tensor X with nonnegative constraints on the factors.
%   This version uses the Lee & Seung multiplicative updates from
%   their NMF algorithm.  The input X can be a tensor, sptensor,
%   ktensor, or ttensor. The result P is a ktensor. It is minimizing
%   the least squares objective function: norm(X - full(P))^2.
%
%   P = CP_NMU(X,R,'param',value,...) specify options:
%     'tol'       - Tolerance on difference in fit {1.0e-4}
%     'maxiters'  - Maximum number of iterations {500}
%     'dimorder'  - Order to loop through dimensions {1:ndims(A)}
%     'init'      - Initial guess [{'random'}|cell array|ktensor]
%     'printitn'  - Print fit every n iterations {1}
%     'trace'     - Return time and function value traces in output {false}
%
%   [P,U0] = CP_NMU(...) also returns the initial guess as a cell array.
%
%   [P,U0,out] = CP_NMU(...) also returns additional output that contains
%   the final fit and the number of iterations performed.
%
%   Examples:
%   X = sptenrand([5 4 3], 10);
%   P = cp_nmu(X,2);
%   P = cp_nmu(X,2,'dimorder',[3 2 1]);
%   U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
%   P = cp_nmu(X,2,'dimorder',[3 2 1],'init',U0);
%
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','cp_nmu_doc.html')))">Documentation page for CP-NMU</a>
%
%   See also KTENSOR, TENSOR, SPTENSOR, TTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>


%% Parse input parameters
N = ndims(X);
params = inputParser;
params.addParameter('tol',1e-4,@isscalar);
params.addParameter('maxiters',500,@(x) isscalar(x) & x > 0);
params.addParameter('dimorder',1:N,@(x) isequal(sort(x),1:N));
params.addParameter('init','random');
params.addParameter('printitn',1,@isscalar);
params.addParameter('trace',false,@islogical);
params.parse(varargin{:});

fitchangetol = params.Results.tol;
maxiters = params.Results.maxiters;
dimorder = params.Results.dimorder;
init = params.Results.init;
printitn = params.Results.printitn;
dotrace = params.Results.trace;
epsilon = 1e-12;  % Small number to protect against round-off error

%% Error checking on maxiters and dimorder
if maxiters < 0
    error('maxiters must be positive');
end
if ~isequal(1:N,sort(dimorder))
    error('dimorder must include all elements from 1 to ndims(X)');
end

%% Set up and error checking on initial guess for U.
if iscell(init) || isa(init,'ktensor')
    if isa(init,'ktensor')
        Uinit = init.U;
        init = 'ktensor';
    else
        Uinit = init;
        init = 'cell';
    end
    if numel(Uinit) ~= N
        error('init does not have %d cells',N);
    end
    for n = dimorder(1:end)
        if ~isequal(size(Uinit{n}),[size(X,n) R])
            error('init{%d} is the wrong size',n);
        end
    end
else
    if strcmp(init,'random')
        Uinit = cell(N,1);
        for n = dimorder(1:end)
            Uinit{n} = rand(size(X,n),R) + 0.1;
        end
    elseif strcmp(init,'nvecs') || strcmp(init,'eigs') 
        Uinit = cell(N,1);
        for n = dimorder(1:end)
            k = min(R,size(X,n)-2);
            fprintf('  Computing %d leading e-vectors for factor %d.\n',k,n);
            Uinit{n} = abs(nvecs(X,n,k));
            if (k < R)
              Uinit{n} = [Uinit{n} rand(size(X,n),R-k)]; 
            end
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;
normX = norm(X);
if dotrace
    fittrace = zeros(maxiters,1);
    timetrace = zeros(maxiters,1);
    itertime = tic;
end

if printitn>0
    fprintf('\nNonnegative CP:\n');
    fprintf('  tol = %g, ', fitchangetol);
    fprintf('  maxiters = %d\n', maxiters);
    if ~isequal(dimorder,1:N)
        fprintf('  dimorder = [%s]\n', num2str(dimorder));    
    end
    fprintf('  init = %s\n', init);
    if dotrace
        fprintf('  trace = true\n');
    end
end

%% Main Loop: Iterate until convergence
for iter = 1:maxiters
    fitold = fit;

    % Iterate over all N modes of the tensor
    for n = dimorder(1:end)

        % Compute the matrix of coefficients for linear system
        Y = ones(R,R);
        for i = [1:n-1,n+1:N]
            Y = Y .* (U{i}'*U{i});
        end
        Y = U{n} * Y;

        % Initialize matrix of unknowns
        Unew = U{n};

        % Calculate Unew = X_(n) * khatrirao(all U except n, 'r').
        tmp = mttkrp(X,U,n) + epsilon;

        % Update unknowns
        Unew = Unew .* tmp;
        Unew = Unew ./ (Y + epsilon);

        U{n} = Unew;
    end

    P = ktensor(U);
    normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
    fit = 1 - (normresidual / normX); %fraction explained by model
    fitchange = abs(fitold - fit);
    if dotrace
        fittrace(iter) = fit;
        timetrace(iter) = toc(itertime);
    end
    if mod(iter,printitn)==0
      fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
    end

    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end
end

%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);

if printitn>0
  normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
  fit = 1 - (normresidual / normX); %fraction explained by model
  fprintf(' Final fit = %e \n', fit);
end

output = struct;
output.params = params.Results;
output.iters = iter;
output.final_fit = fit;
if dotrace
    output.fit_trace = fittrace(1:iter);
    output.time_trace = timetrace(1:iter);
end

end
