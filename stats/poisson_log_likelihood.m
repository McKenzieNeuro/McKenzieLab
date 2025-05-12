function logL = poisson_log_likelihood(y, lambda,delta_t)
% Ensure inputs are column vectors
y = y(:);
lambda = lambda(:);

% Avoid log(0) issues by adding a small epsilon
epsilon = 1e-10;
lambda_per_bin = lambda * delta_t;       % Convert rate to expected count
lambda_safe = lambda_per_bin + epsilon;

% Compute log-likelihood
logL = sum(y .* log(lambda_safe) - lambda_per_bin - gammaln(y + 1));
end
