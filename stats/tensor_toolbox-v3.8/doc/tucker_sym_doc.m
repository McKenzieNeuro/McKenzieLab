%% Symmetric Tucker Decomposition
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a>
% &#62;&#62; <a href="tucker.html">Tucker Decompositions</a>
% &#62;&#62; <a href="tucker_sym_doc.html">Symmetric Tucker</a>
% </p>
% </html>
%
% The function |tucker_sym| computes the best symmetric rank-(R,R,...,R)
% approximation of a symmetric tensor S, according to the specified rank R.
% The result returned in T is a ttensor with all factors equal. S must be a
% tensor (not a symtensor) and must be exactly symmetric.
%
% The method is based on the algorithm described in:
% * Phillip A. Regalia, Monotonically Convergent Algorithms for Symmetric Tensor Approximation,
%   Linear Algebra and its Applications 438(2):875-890, 2013, http://dx.doi.org/10.1016/j.laa.2011.10.033


%% Create a symmetric tensor of size 4x4x4
rng(0); %<-- Set seed for reproducibility
S = tensor(@randn,[4,4,4]);
S = symmetrize(S);

%% Compute a symmetric Tucker decomposition with rank 2 with default parameters
% Except for the rank, all parameters are optional. We specify here that
% we want to print the fit every 10 iterations. The usual default is every iteration.
% By default, the initial guess is random. Since this is a nonconvex problem,
% we will run the algorithm multiple times to find the best solution.
best_err = inf;
best_T = [];
errs = zeros(3,1); % Store errors for the first 3 examples
for i = 1:10
	rng(i); % Set different seed each time
	T_tmp = tucker_sym(S,2,'printitn',10);
	err = norm(S - full(T_tmp)) / norm(S);
	if err < best_err
		best_err = err;
		best_T = T_tmp;
	end
end
errs(1) = best_err;
fprintf('\n\n*** Best relative error over 10 runs: %g ***\n', best_err);
disp(best_T)

%% We do the same test again, but provide user-specified initial guesses
best_err = inf;
best_T = [];
for i = 1:10
	rng(i); % Set different seed each time
	X0 = randn(4,2); % Initial guess for factor matrix
	T_tmp = tucker_sym(S,2,'init',X0,'printitn',10);
	err = norm(S - full(T_tmp)) / norm(S);
	if err < best_err
		best_err = err;
		best_T = T_tmp;
	end
end
errs(2) = best_err;
fprintf('\n\n*** Best relative error over 10 runs (user init): %g ***\n', best_err);
disp(best_T)

%% We do the same test again, but use the n-vecs initialization method
% This initialization is more expensive but generally works very well.
% This initialization method is deterministic, so we don't need to set 
% the random seed.  
T = tucker_sym(S,2,'init','nvecs','printitn',10);
errs(3) = norm(S - full(T)) / norm(S);

%% Comparison of errors with different initialization methods
% The first two methods are repeated runs with random initializations,
% while the third method is a single run with the n-vecs initialization.
% Generally, the n-vecs initialization is better than a random initialization,
% but it is not guaranteed to find the best solution. Repeated runs with
% random initializations may find a better solution than a single run with n-vecs.
fprintf('\nSummary of best relative errors for each initialization method:\n');
fprintf('  Default random init (best of 10): %g\n', errs(1));
fprintf('  User random init (best of 10):    %g\n', errs(2));
fprintf('  N-vecs init (single run):         %g\n', errs(3));

%% Example: Symmetric tensor with known rank-2 decomposition
% If we know the rank, then |tucker_sym| should compute an exact decomposition
% with an error of zero (or very close to zero due to numerical errors).

% Create a symmetric factor matrix
A = rand(4,2);
% Construct a symmetric core tensor (2x2x2)
G = tensor(rand(2,2,2));
% Build the symmetric Tucker tensor
S_known = ttensor(G, {A, A, A});
% Symmetrize to ensure exact symmetry
S_known = symmetrize(full(S_known));
% Compute the symmetric Tucker decomposition
T_est = tucker_sym(S_known, 2);
disp('Relative error between original and estimated tensor:');
disp(norm(S_known - full(T_est)) / norm(S_known));