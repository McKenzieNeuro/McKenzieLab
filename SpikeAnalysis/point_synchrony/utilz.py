import numpy as np
import logging



NEG_INF = -float('inf')
EPS = np.finfo(float).eps / 2

def _next_power_of_two(n):
    return 2 ** int(np.ceil(np.log2(n)))


def pairwise_convolution_lengths(a, b):
    true = a + b - 1
    return true, _next_power_of_two(true)


def log_min_pos(log_pmf):
    return log_pmf[log_pmf > NEG_INF].min()

def log_sum(log_u):
    """Compute `log(sum(exp(log_u)))`"""
    if len(log_u) == 0:
        return NEG_INF
    maxi = np.argmax(log_u)
    max = log_u[maxi]
    if max == NEG_INF:
        return max
    else:
        exp = log_u - max
        np.exp(exp, out = exp)
        return np.log1p(np.sum(exp[:maxi]) + np.sum(exp[maxi + 1:])) + max

      
def shift(log_pmf, theta):
    shifted = log_pmf + theta * np.arange(len(log_pmf))
    log_mgf = log_sum(shifted)
    shifted -= log_mgf
    return shifted, log_mgf


def unshift(convolved, theta, *mgfs):
    c = convolved - theta * np.arange(len(convolved))
    for (mgf, multiplicity) in mgfs:
        c += multiplicity * mgf
    return c
    



def log_dynamic_range_shifted(log_pmf, theta):
    # this is equivalent to log_dynamic_range(shift(log_pmf,
    # theta)[0]), but is more efficient, as it avoids an unnecessary
    # log_sum computation.
    shifted = log_pmf + np.arange(float(len(log_pmf))) * theta
    lo = log_min_pos(shifted)
    hi = log_sum(shifted * 2.0) / 2
    return hi - lo