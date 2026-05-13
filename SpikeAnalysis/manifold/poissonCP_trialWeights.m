function model = poissonCP_trialWeights(X, opts)
% POISSONCP_TRIALWEIGHTS
% Fit a Poisson CP decomposition to a 3D count tensor
%
% X    : N x T x M tensor of counts (neurons × time × trials)
% opts : struct with fields:
%   - K        : number of templates (default 4)
%   - maxIter  : maximum iterations (default 200)
%   - tol      : convergence tolerance (default 1e-4)
%   - lambdaC  : optional L1 penalty on trial weights (C) for sparsity (default 0)
%   - verbose  : true/false (default true)
%
% OUTPUT model struct with:
%   - A : N x K neuron/template factors
%   - B : T x K temporal factors
%   - C : M x K trial weights
%   - M : full CP tensor object
%   - info : convergence info from cp_apr
%   - opts : options used

%% ---------------- Defaults ----------------
if ~isfield(opts,'K'),       opts.K = 20;        end
if ~isfield(opts,'maxIter'), opts.maxIter = 200; end
if ~isfield(opts,'tol'),     opts.tol = 1e-4;   end
if ~isfield(opts,'lambda'), opts.lambda = zeros(size(X,1),1) ; end
if ~isfield(opts,'verbose'), opts.verbose = true; end

%% ---------------- Convert to tensor ----------------
if ~(isa(X,'tensor') | isa(X,'sptensor'))
    Xten = tensor(X);
else
    Xten = X;
end

%% ---------------- Fit Poisson CP decomposition ----------------
% Tensor Toolbox cp_apr uses alternating Poisson regression
% lambdaC is applied as L1 penalty on C (trial weights)
K = opts.K;

% options for cp_apr
cpOpts = struct;
cpOpts.maxiters = opts.maxIter;
%cpOpts.tol = opts.tol;
cpOpts.printitn = double(opts.verbose);
cpOpts.precompinds = true;
% Apply sparsity to trial mode if requested
%if opts.lambdaC > 0
%    cpOpts.lambda = {0, 0, opts.lambdaC}; % N,T,M modes
%end

cpOpts.lambda = opts.lambda;

% Run Poisson CP decomposition
[M, info] = sm_cp_apr(Xten, K, cpOpts);

%% ---------------- Extract factors ----------------
model.A = M.U{1};  % N x K neuron/template factors
model.B = M.U{2};  % T x K temporal factors
model.C = M.U{3};  % M x K trial weights
model.M = M;        % full tensor object
model.info = info;
model.opts = opts;

%% ---------------- Normalize for interpretability ----------------
% Normalize A and B, absorb scale into C
for k = 1:K
    s = norm(model.A(:,k),1) * norm(model.B(:,k),1);
    if s > 0
        model.A(:,k) = model.A(:,k) / norm(model.A(:,k),1);
        model.B(:,k) = model.B(:,k) / norm(model.B(:,k),1);
        model.C(:,k) = model.C(:,k) * s;
    end
end

if opts.verbose
    fprintf('Poisson CP decomposition finished. Templates: %d\n', K);
end
end
