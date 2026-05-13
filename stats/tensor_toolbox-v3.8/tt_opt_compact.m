function [xbest, fbest, info] = tt_opt_compact(xinit, fgh, varargin)
%TT_OPT_COMPACT Wrapper for the Compact Representation optimization method.
%
%   [X, F, INFO] = TT_OPT_COMPACT(X0, FGH, 'param', value, ...) is a wrapper
%   for the compact representation quasi-Newton methods from Brust [1]. 
%   The wrapper interfaces the Tensor Toolbox with compLS1, a
%   line-search implementation of compact quasi-Newton methods. Here X0 is in the
%   initial guess for the solution, FGH is a function handle to a function
%   that returns the function and gradient, and then there are optional
%   parameter-value pairs (see below).
%
%   Reference:
%   [1] Brust, J. J. (2024). Useful Compact Representations for Data-Fitting. 
%   arXiv preprint arXiv:2403.12206
%
%   For more optimzation algorithm choices and parameters, see
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','opt_options_doc.html')))">Tensor Toolbox Optimization Methods</a>, and
%   <a href="matlab:web(strcat('file://',fullfile(getfield(what('tensor_toolbox'),'path'),'doc','html','tt_opt_doc.html')))">Tensor Toolbox Optimization Methods for Developers</a>
%
%   See also TT_OPT_ADAM, TT_OPT_LBFGSB, TT_OPT_FMINUNC, TT_OPT_LBFGS.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

%%
setupTimer = tic;

%% Setup
if ~iscolumn(xinit)
    error('Initial guess must be a column vector');
end

n = length(xinit);

%% Algorithm Parameters
params = inputParser;
params.KeepUnmatched = true;
params.addParameter('printitn', 1); % DisplayIters
params.addParameter('m', 5); % M
params.addParameter('maxiters', 1000); % MaxIters
params.addParameter('subiters', 10) % MaxFuncEvals = maxiters * subiters
params.addParameter('ftol', 1e-10); % RelFuncTol
params.addParameter('gtol', 1e-5); % StopTol
params.addParameter('mdesc', 'CR solver (limited-memory) (https://github.com/johannesbrust/CR)');
params.addParameter('xdesc', []);
params.parse(varargin{:});

%% Options for printing
mdesc = params.Results.mdesc;
xdesc = params.Results.xdesc;

%% Setting optimization parameters
opts = params.Unmatched;
opts.l = params.Results.m;
opts.maxIt = params.Results.maxiters;
opts.maxItLS = 50*params.Results.subiters;
opts.tol = params.Results.gtol;

%% Other parameters (defaults used)
% opts.whichV = 's', 'g', 'ag', 'y'                         (limited-memory formulas)
% opts.whichL = 'armijoG', 'wolfeG', 'wolfeB', 'wolfeMT'    (search options)
% opts.nrmtol = 'p'                                         (norm(.,p) for norms)
% opts.c1     = 1e-4                                        (armijo condition)
% opts.c2     = 0.9                                         (wolfe condition)
% opts.stepMin = 1e-12                                      (min step)

printitn = params.Results.printitn;
if printitn == 0
    opts.print = 0;
else
    opts.printitn = printitn;
end

opts.store = true;

%% Welcome message
if printitn > 0
    fprintf('\n');
    fprintf('%s\n', mdesc);
    if isempty(xdesc)
        fprintf('Number of Variables: %d\n', n);
    else
        fprintf('%s\n', xdesc);
    end
    fprintf('Parameters: ');    
    fprintf('m=%d, ', opts.l);
    fprintf('gtol=%0.2g, ', opts.tol);
    fprintf('maxiters = %d ', opts.maxIt);
    fprintf('subiters = %d', opts.maxItLS);
    fprintf('\n');
    fprintf('\n');
    fprintf('Begin Main Loop\n');
end
setuptime = toc(setupTimer);

%% Check compact representation in path
if ~exist('compLS1.m','file')
    
    ttbpath = fileparts(which('tt_opt_compact.m'));
    crpath = fullfile(ttbpath,'libraries','compact','src');
    
    fprintf('*** Important Notice ***\n');
    fprintf('-> compLS1.m is not in your current path!\n')
    fprintf('-> Adding to your path: %s\n',crpath);
    fprintf('-> Save your path to avoid this notice in the future\n');
    fprintf('***\n');
    addpath(crpath,'-BEGIN');
    fprintf('***\n');
end

%% Run optimization

[xbest,fbest,outs] = compLS1(xinit, fgh, opts);


%% Save stuff
info.params = params.Results;
info.f_trace = outs.fs;
info.gnorm_trace = outs.gs;
info.time_trace = [];
info.setup_time = setuptime;
info.opt_time = outs.time;
info.f_final = fbest;
info.iters = outs.it;
info.subiters = outs.nf;

% parse exit condition
if outs.ex == 1
    info.exit_condition = 'Convergence: Norm tolerance';
else
    info.exit_condition = 'Non-convergence: Maximum iterations';
end

%% Goodbye message
if printitn > 0
    fprintf('End Main Loop\n');
    fprintf('\n');
    fprintf('Final f: %10.4e\n', info.f_final);
    fprintf('Setup time: %.2g seconds\n', info.setup_time);
    fprintf('Optimization time: %.2g seconds\n', info.opt_time);
    fprintf('Iterations: %d\n', info.iters);
    fprintf('Total iterations: %d\n', info.subiters);
    fprintf('Exit condition: %s\n', info.exit_condition);
end


info.out = outs;

