%% Alternating Poisson Regression for fitting CP to sparse count data
%
% <html>
% <p class="navigate">
% &#62;&#62; <a href="index.html">Tensor Toolbox</a> 
% &#62;&#62; <a href="cp.html">CP Decompositions</a> 
% &#62;&#62; <a href="cp_apr_doc.html">CP-APR</a>
% </p>
% </html>
%
% References: 
% 
% * E. C. Chi, T. G. Kolda, On Tensors, Sparsity, and Nonnegative Factorizations,
%   SIAM J. Matrix Analysis and Applications, 33:1272-1299, 2012, 
%   <https://doi.org/10.1137/110859063>
% * S. Hansen, T. Plantenga and T. G. Kolda, Newton-Based Optimization
%   for Kullback-Leibler Nonnegative Tensor Factorizations, 
%   Optimization Methods and Software, 30(5):955-979, 2015, 
%   <http://dx.doi.org/10.1080/10556788.2015.1009977>
%
%% Set up a sample problem
% We follow the general procedure for creating a problem outlined by Chi and Kolda (2012).
% This creates a sparse count tensor with a known solution. The solution is
% a CP decomposition with a few large entries in each column of the factor
% matrices. The solution is normalized and sorted by component size in 
% descending order.

rng('default') %<- Setting random seed for reproducibility of this script

% Pick the size and rank
sz = [100 80 60];
R = 5;

% Generate factor matrices with a few large entries in each column; this
% will be the basis of our soln.
A = cell(3,1);
for n = 1:length(sz)
    A{n} = rand(sz(n), R);
    for r = 1:R
        p = randperm(sz(n));
        nbig = round( (1/R)*sz(n) );
        A{n}(p(1:nbig),r) = 100 * A{n}(p(1:nbig),r);
    end
end
lambda = rand(R,1);
S = ktensor(lambda, A);
S = normalize(S,'sort',1);

% Create sparse test problem based on provided solution. 
nz = prod(sz) * .05;
info = create_problem('Soln', S, 'Sparse_Generation', nz);

% Extract data and solution
X = info.Data;
M_true = info.Soln;

%% Call CP-APR
% Alternating Poisson Regression (APR) is a method for fitting a CP
% decomposition to sparse count data. It is a nonnegative method that
% minimizes the Kullback-Leibler divergence between the data and the
% model. The method is implemented in the |cp_apr| function, which is
% similar to the |cp_als| function, but uses a different objective
% function and optimization method. The |cp_apr| function is designed to
% handle sparse count data and is particularly useful for fitting
% nonnegative CP decompositions.
%
% The |cp_apr| function is a wrapper that calls one of three specific
% algorithms, selected by the 'alg' parameter:
%
% * |'pqnr'| (Default): Row subproblems are solved by Projected Quasi-Newton
%   with L-BFGS. This method generally offers a good balance of speed and
%   robustness and is suitable for a wide range of problems. It approximates
%   the Hessian using gradient information from previous iterations.
%   It is based on the work by Hansen, Plantenga, and Kolda (2015).
% * |'pdnr'|: Row subproblems are solved by Projected Damped Newton. This
%   method uses the exact Hessian for the row subproblems, which can lead
%   to higher accuracy per iteration but may be more computationally
%   intensive, especially for large R, as it involves forming and solving
%   an R x R system at each inner iteration for each row.
%   It is based on the work by Hansen, Plantenga, and Kolda (2015).
% * |'mu'|: Multiplicative Update. This is a simpler algorithm, often with
%   cheaper iterations. 
%   It can be slower to converge to high accuracy compared to Newton-based
%   methods but can be effective for very large, sparse problems or for
%   obtaining an initial guess quickly. It is based on the work by
%   Chi & Kolda (2012).
%
% The following example uses the default 'pqnr' algorithm.
% We run the method 5 times and keep the best solution (lowest objective
% value).

% Compute a solution using the default 'pqnr' algorithm
fprintf('--- Running CP-APR with PQNR (default) ---\n');
rng('default')
best_obj = inf;
for trial = 1:5
    [M_tmp,~,output_tmp] = cp_apr(X, R, 'printitn', 2);
    if output_tmp.obj < best_obj
        best_obj = output_tmp.obj;
        M_pqnr = M_tmp;
        output_pqnr = output_tmp;
    end
end

% Score the solution (a score of 1 is perfect)
factor_match_score_pqnr = score(M_pqnr, M_true, 'greedy', true);

%% Example using the 'pdnr' algorithm
% Here, we explicitly select the 'pdnr' algorithm. We also reduce the
% maximum number of iterations for this demonstration.
% We run the method 5 times and keep the best solution (lowest objective
% value).

fprintf('--- Running CP-APR with PDNR ---\n');
rng('default')
best_obj_pdnr = inf;
for trial = 1:5
    [M_tmp,~,output_tmp] = cp_apr(X, R, 'alg', 'pdnr', 'printitn', 2);
    if output_tmp.obj < best_obj_pdnr
        best_obj_pdnr = output_tmp.obj;
        M_pdnr = M_tmp;
        output_pdnr = output_tmp;
    end
end

% Score the solution
factor_match_score_pdnr = score(M_pdnr, M_true, 'greedy', true);

%% Example using the 'mu' algorithm
% This example demonstrates the 'mu' algorithm. We can also set parameters
% specific to 'mu', like 'kappa'. Because this is faster per iteration by
% requires more iterations, we only print the results every 10 iterations.
% We run the method 5 times and keep the best solution (lowest objective
% value).

fprintf('--- Running CP-APR with MU ---\n');
rng('default')

best_obj_mu = inf;
for trial = 1:5
    [M_tmp,~,output_tmp] = cp_apr(X, R, 'alg', 'mu', 'printitn', 10, 'maxiters', 200, 'kappa', 50);
    if output_tmp.obj < best_obj_mu
        best_obj_mu = output_tmp.obj;
        M_mu = M_tmp;
        output_mu = output_tmp;
    end
end
% Score the solution
factor_match_score_mu = score(M_mu, M_true, 'greedy', true);

%% Comparing results
% We can see that all methods can find a reasonable solution, though
% convergence speed and final accuracy might differ. The 'pqnr' and 'pdnr'
% methods are generally more sophisticated and may converge to a better
% solution or faster in terms of outer iterations, while 'mu' iterations
% are typically cheaper.
%
% For this particular problem and random initialization:
fprintf('Factor Match Score (PQNR): %.4f\n', factor_match_score_pqnr);
fprintf('Factor Match Score (PDNR): %.4f\n', factor_match_score_pdnr);
fprintf('Factor Match Score (MU):   %.4f\n', factor_match_score_mu);

%% Visualize the results
% The function values are only computed when printing the output to the
% screen. The legend indicates the method and, in parentheses, how
% frequently the objective value was recorded. 
figure(1); clf;
tf = ~isnan(output_pqnr.fnVals);
plot(output_pqnr.times(tf), 1-output_pqnr.fnVals(tf),'-+','DisplayName','PQNR (2)');

hold on;
tf = ~isnan(output_pdnr.fnVals);
plot(output_pdnr.times(tf), 1-output_pdnr.fnVals(tf),'-.o', 'DisplayName','PDNR (2)');

tf = ~isnan(output_mu.fnVals);
plot(output_mu.times(tf), 1-output_mu.fnVals(tf),'--*', 'DisplayName','MU (10)');
hold off;

title('CP-APR');
xlabel('Time (seconds)');
ylabel('Function Value');
legend('Location','SouthEast')