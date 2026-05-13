%% Nonnegative CP Decomposition using Multiplicative Updates
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_nmu_doc.html">CP-NMU</a>
% </p>
% </html>
%
% The function |cp_nmu| computes a nonnegative CP factorization using 
% multiplicative updates. It does this using Lee & Seung multiplicative 
% updates for nonnegative matrix/tensor factorization with a least squares
% objective function. The input X can be a tensor,
% sptensor, ktensor, or ttensor. The output CP model is a |ktensor|.
%
% CP-NMU is a simple method for nonnegative CP decomposition and
% is not guaranteed to converge to a stationary 
% point. Use <cp_opt_doc.html CP-OPT> 
% with nonnegativity constraints for a more robust method.
%
% Additionally, we have only implemented here the method for 
% the least squares objective function.
% Use <cp_apr_doc.html CP-APR> for KL Divergence.
%
% The CP-NMU method (for matrices) is described in:
%
% * D. D. Lee and H. S. Seung.
%   Algorithms for Non-negative Matrix Factorization.
%   Advances in Neural Information Processing Systems, 13, 2001.
%   <https://proceedings.neurips.cc/paper_files/paper/2000/hash/f9d1152547c0bde01830b7e8bd60024c-Abstract.html>
%

%% Set up a sample problem
% We create a nonnegative tensor with a known solution and add a small amount of noise.

rng('default') % For reproducibility
sz = [30 25 20];
R = 3;
A = cell(3,1);
for n = 1:3
    A{n} = rand(sz(n), R);
end
lambda = rand(R,1);
S = ktensor(lambda, A);
X = full(S) + 0.01 * tensor(rand(sz));
M_true = S;


%% Run CP-NMU with default options
% The default is random initialization and 500 iterations.
%
% *NOTE* The random initialization uses random values in [0.1,1.1] to avoid
% zeros, which cause problems with the multiplicative updates.

rng('default') % For reproducibility
M = cp_nmu(X, R);
score_val = score(M, M_true);
fprintf('Score (run with default options): %.4f\n', score_val);


%% Run with custom initialization
% You can provide a cell array as the initial guess.

rng('default') % For reproducibility
Uinit = cell(3,1);
for n = 1:3
    Uinit{n} = rand(sz(n), R);
end
M1 = cp_nmu(X, R, 'init', Uinit, 'maxiters', 100);
score_val = score(M1, M_true);
fprintf('Score (run with custom init): %.4f\n', score_val);

%% Run with ktensor initialization
% Alternatively, you can pass in a ktensor as the initial guess.
% Here, I'm using the same factors as the custom initialization 
% and the result will be the same.

kt_init = ktensor(Uinit);
M2 = cp_nmu(X, R, 'init', kt_init, 'maxiters', 100);
score_val = score(M2, M_true);
fprintf('Score (run with ktensor init): %.4f\n', score_val);

fprintf('Solutions equal? %d\n', isequal(M1, M2));

%% Trace option
% You can enable tracing to get fit and time per iteration.

rng('default') % For reproducibility
[M, ~, out] = cp_nmu(X, R, 'trace', true, 'maxiters', 50);
score_val = score(M, M_true);
fprintf('Score (run with default options): %.4f\n', score_val);
figure(1); clf;
plot(out.time_trace, out.fit_trace, '-o');
xlabel('Time (seconds)'); ylabel('Fit');
title('CP-NMU Fit Trace');
grid on;


%% Varying options
% You can change the order of updates or the tolerance or the maximum number of iterations.

rng('default') % For reproducibility
[M, U0, out] = cp_nmu(X, R, 'dimorder', [3 2 1], 'tol', 1e-5, 'maxiters', 100);
score_val = score(M, M_true);
fprintf('dimorder/tol score: %.4f\n', score_val);

