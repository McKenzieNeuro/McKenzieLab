% THIS SOFTWARE AND ANY ACCOMPANYING DOCUMENTATION IS RELEASED "AS IS."  THE U.S. GOVERNMENT MAKES NO WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, CONCERNING THIS SOFTWARE AND ANY ACCOMPANYING DOCUMENTATION, INCLUDING, WITHOUT LIMITATION, ANY WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  IN NO EVENT WILL THE U.S. GOVERNMENT BE LIABLE FOR ANY DAMAGES, INCLUDING ANY LOST PROFITS, LOST SAVINGS OR OTHER INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE, OR INABILITY TO USE, THIS SOFTWARE OR ANY ACCOMPANYING DOCUMENTATION, EVEN IF INFORMED IN ADVANCE OF THE POSSIBILITY OF SUCH DAMAGES.
%
% file: get_r_ci_cs1_1s.m
% computes the COMPASE standard confidence intervals for a rate estimate
%  - assumes x is the number of occurrances of an event in area A
%  - finds confidence interval around the estimate for rate R
%  - with confidence level 1 - alpha
%  - method used is based on a normal approximation or numerical integration
%  - finds both the lower and upper one-sided intervals
%
% This version is always accurate, but is slower than CS2
% for small x and small alpha.
%
% 010129 tdr created from get_r_ci_cs1.m for one-sided CIs
% 010412 tdr turned off the na thing.

function ci = get_r_ci_cs1_1s(x,A,alpha)

if nargin < 3, 
    error('Requires three input arguments (x,A,alpha)');
end

%na_ok = r_ci_na_ok(x,alpha);
na_ok = 0;

if na_ok,
   ci = get_r_ci_na_1s(x,A,alpha);
else
   ci = get_r_ci_ibp_1s(x,A,alpha);
end;
