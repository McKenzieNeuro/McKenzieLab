function p = permtest(x, y, nperm)
%PERMTEST  Unpaired sign-flip / label-permutation test on the difference of means.
%
%   p = permtest(x, y, nperm)
%
%   Two-sided p for mean(y) - mean(x), by randomly reassigning the pooled
%   values to two groups of the original sizes, nperm times. NaNs are
%   dropped first. Default nperm = 1e4 if omitted.
%
%   This is what was called 'local_permtest' inside runSyncAnalysis.m --
%   local functions at the bottom of a script are scoped to that script's
%   run and are not callable afterward from the command line. Pulled out
%   here so it's usable standalone.

if nargin < 3, nperm = 1e4; end
x = x(:); x = x(~isnan(x));
y = y(:); y = y(~isnan(y));

obs = mean(y) - mean(x);
pool = [x; y]; nx = numel(x);
null = zeros(nperm,1);
for k = 1:nperm
    s = pool(randperm(numel(pool)));
    null(k) = mean(s(nx+1:end)) - mean(s(1:nx));
end
p = (1 + sum(abs(null) >= abs(obs))) / (1 + nperm);
end
