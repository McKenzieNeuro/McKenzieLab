THIS SOFTWARE AND ANY ACCOMPANYING DOCUMENTATION IS RELEASED "AS IS."  THE U.S. GOVERNMENT MAKES NO WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, CONCERNING THIS SOFTWARE AND ANY ACCOMPANYING DOCUMENTATION, INCLUDING, WITHOUT LIMITATION, ANY WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT WILL THE U.S. GOVERNMENT BE LIABLE FOR ANY DAMAGES, INCLUDING ANY LOST PROFITS, LOST SAVINGS OR OTHER INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE, OR INABILITY TO USE, THIS SOFTWARE OR ANY ACCOMPANYING DOCUMENTATION, EVEN IF INFORMED IN ADVANCE OF THE POSSIBILITY OF SUCH DAMAGES.

==================================================

Title

Code for Accurate Confidence Intervals for Proportion, Rate, and Difference Estimations

==================================================

Change Log:
	2003-02-18 Original
	2003-03-06 Minor corrections to prop_ci.m and README.txt 
	2003-12-23 Minor changes to handle NaN inputs more gracefully
	2004-02-03 Minor correct to ensuring integer inputs in prop_diff_ci.m
	2004-02-10 Added stubs for vestigial functions left over from development.  This prevents compiler warnings.

==================================================

Use

There are three functions each for rate and proportion estimation problems.  These six functions use a few dozen other functions that are included as separate files, but the documentation assumes that only these six functions are used directly by a user.

The six functions are:
	prop_ci.m		rate_ci.m
	prop_diff.m		rate_diff.m
	prop_diff_ci.m		rate_diff_ci.m

MatLab function calls:
	function ci = prop_ci(x,n,alpha,[method, [verbose]]) 
	function ci = rate_ci(x,A,alpha,[method,[verbose]]) 
	function y = prop_diff(x1,n1,x2,n2,delta)
	function y = rate_diff(x1,A1,x2,A2,delta)
	function ci = prop_diff_ci(x1,n1,x2,n2,alpha,[method, [verbose]]) 
	function ci = rate_diff_ci(x1,A1,x2,A2,alpha,[method,[verbose]]) 

Inputs
	*** All inputs are scalars ***
	x - 	number of positive outcomes (numerator of the proportion estimate), 0 <= x <= n 
	n - 	number of test cases (denominator of the proportion estimate), 1 <= n <= 10^5 
	A - 	Area of test (denominator of the rate estimate), 10^-5 < A < 10^5 
	alpha -	1 - confidence level (e.g., for 0.95 confidence, alpha = 0.05), 10^-4 <= alpha <= 1-10^-4 
	method-	optional, default is 2 for prop_ci and rate_ci, default is 3 for prop_diff_ci and rate_diff_ci 
		Method 1 provides one-sided CIs 
		Method 2 provides two-sided minimum length CIs (not available for diff_ci's)
		Method 3 provides two-sided balanced length CIs 
		Method 4 provides two-sided balanced tail CIs 
		Method 5 is the so called "exact" Clopper-Pearson approach.  (not available for rate_diff_ci's)
		Method 6 is the Normal approximation.  (not available for rate_diff_ci's)
		Note: Methods 1 through 4 provided here are accurate (alpha is within 0.00005 of requested value) for most meaningful values of n (0, 10^5], x [0, n], and alpha (1e-4,1).  These are the recommended methods. 
	verbose - optional, default is 0 (not verbose) 
	delta - the difference of interest (proportion or rate)
	
Outputs
	ci - 1 x 3 vector where 
		c(1) - p_hat = x/n, r_hat = x/A, delta_p_hat = x1/n1 - x2/n2, or delta_r_hat = n1/A1 - n2/A2 
		c(2) - lower bound on confidence interval (ignore for one-sided methods when interested in upper bound) 
		c(3) - upper bound on confidence interval (ignore for one-sided methods when interested in lower bound) 
		units of rate CI are the same as input A (e.g., square kilometers)
	y - 	Pr (p1 - p2 >= delta), where p1_hat = x1/n1 and p2_hat = x2/n2 or
		Pr ( r1 - r2 >= delta), where r1_hat = x1/A1 and r2_hat = x2/A2
	Additional outputs for verbose = 1, see references 

==================================================

Installation

Current version available at MatLab Central File Exchange (http://www.mathworks.com/matlabcentral/fileexchange) in the Statistics Category.  Users may want to check for updates.

Requirements: MatLab Release 13 with the Statistics Toolbox.  Older releases may work, especially if you replace all quadl calls with quad8 calls, but only limited testing has been performed with releases 11 and 12.  Note also that all testing was performed on a Windows PC (Windows 98 and 2000).

Installation Procedure: unzip, copy source code to a local directory and include that directory in your MatLab path.

Test Cases:
There are six test functions that may help verify proper installation
	test_prop_ci.m		test_rate_ci.m
	test_prop_diff_ci.m	test_rate_diff_ci.m
	test_prop_diff.m	test_rate_diff.m
All of the test programs run without arguments and produce a .csv data file (in the current working directory).  The test programs should run without error (various warnings are normally produced) and produce results not significantly different from the "reference" data file provided with the code. Note that the tests can take several hours to run, with the diff_ci tests taking most of the time.

==================================================

Theory Summary

We determine the posterior distribution of a parameter (Pd or FAR) conditioned on the measured values (no. of detections, no. of test targets, no. of false alarms, area searched).  The posterior distribution is computed assuming a binomial distribution for Pd (this is a rock solid assumption) and a Poisson distribution for FAR (this is a pretty good assumption).  "Diffuse" priors are assumed, i.e., we assume that we don't have any idea about the metric prior to getting measurements.  This is a good assumption.  Finally, we use Bayes theorem to get the posterior.  We then search for upper and lower bounds around the estimated values (Pd_hat = detects/tested, FAR_hat = false alarms/area) such that the integral of the posterior from the lower bound to the upper bound is equal to 1-alpha.  I believe you can refer to these a "Bayesian confidence intervals" and you can say they have the following meaning:  The probability that the observed values were produced by actual parameters that fall outside the CI is equal to alpha.  Under most circumstances, Bayesian CIs are the same as classical CI's, particularly for diffuse priors.

==================================================

Contact

Tim Ross, Air Force Research Lab, t.ross@ieee.org

==================================================

References

	Ross, T. D., "Confidence Intervals for ATR Performance Metrics", Proceedings of SPIE, Algor. for SAR Imagery VIII, 2001, Paper 4382-41
	Ross, T. D., "Accurate Confidence Intervals for Binomial Proportions and Poisson Rate Estimation", Computers in Biology and Medicine, Vol.33, Issue 6, pp.509-531, 2003.

