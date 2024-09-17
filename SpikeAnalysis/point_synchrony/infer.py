import numpy as np
from scipy.stats import binom, norm
from scipy.fftpack import fft, ifft
from scipy.optimize import brentq # Find the root of the function mu
import logging
import utilz
NEG_INF = -float('inf')
EPS = np.finfo(float).eps / 2

def synapse(Rspk,Tspk,Delta,delta,tau,confidence_level,M,return_exact,return_normal):
    results = {}
    params = preprocess_spike_trains(Rspk,Tspk,Delta,delta,tau)
    params['alpha_ovr_2'] = (1-confidence_level)/2
    params['M_div'] = M

    ptem = np.append(params['q_syn'],params['q_asyn'])
    mutem = np.sum(ptem)
    sigmatem = np.sqrt(np.sum(ptem*(1-ptem)))
    ncdf = norm(loc=mutem, scale=sigmatem);
    pval = ncdf.cdf(params['cardS'])
    logging.warning(pval)
    
    if pval > (1-params['alpha_ovr_2']) or pval < params['alpha_ovr_2']:
        # Get chernoff bounds
        if pval > (1-params['alpha_ovr_2']):
            R_0=params['cardS']
            L_0=0
            is_pos_h = True
            is_left_tail = True
            params['chernoff_high_pos'] = binary_search_pivot_for_bound_chernoff(params,R_0,L_0,is_pos_h,is_left_tail)

            R_0=params['cardS']
            L_0=0
            is_pos_h = True
            is_left_tail = False
            params['chernoff_low_pos'] = binary_search_pivot_for_bound_chernoff(params,R_0,L_0,is_pos_h,is_left_tail)

        elif pval < params['alpha_ovr_2']:
            R_0=0
            L_0=-params['cardE']
            is_pos_h=False
            is_left_tail= True
            params['chernoff_low_neg'] = binary_search_pivot_for_bound_chernoff(params,R_0,L_0,is_pos_h,is_left_tail)

            R_0=0
            L_0=-params['cardE']
            is_pos_h=False
            is_left_tail= False
            params['chernoff_high_neg'] = binary_search_pivot_for_bound_chernoff(params,R_0,L_0,is_pos_h,is_left_tail)


        # Normal approximation CI's
        if return_normal:
            if pval > (1-params['alpha_ovr_2']):
                params['is_pos_h'] = True
                params['is_left_tail'] = True
                params['ci_normal_high_pos'] = binary_search_pivot_for_bound_normal(params)


                params['is_pos_h'] = True
                params['is_left_tail'] = False
                params['ci_normal_low_pos'] = binary_search_pivot_for_bound_normal(params)

            elif pval < params['alpha_ovr_2']:
                params['is_pos_h'] = False
                params['is_left_tail'] = True
                params['ci_normal_low_neg'] = binary_search_pivot_for_bound_normal(params)


                params['is_pos_h'] = False
                params['is_left_tail'] = False
                params['ci_normal_high_neg'] = binary_search_pivot_for_bound_normal(params)   

        # Exact CI's
        if return_exact:
            if pval > (1-params['alpha_ovr_2']):
                params['is_pos_h'] = True
                params['is_left_tail'] = True
                params['ci_exact_high_pos'] = binary_search_pivot_for_bound_exact_modified(params) # ,xx1,yy1


                params['is_pos_h'] = True
                params['is_left_tail'] = False
                params['ci_exact_low_pos'] = binary_search_pivot_for_bound_exact_modified(params) #,xx0,yy0

            elif pval < params['alpha_ovr_2']:
                params['is_pos_h'] = False
                params['is_left_tail'] = True
                params['ci_exact_low_neg'] = binary_search_pivot_for_bound_exact_modified(params)


                params['is_pos_h'] = False
                params['is_left_tail'] = False
                params['ci_exact_high_neg'] = binary_search_pivot_for_bound_exact_modified(params)
                 
        results['Delta'] = params['Delta']
        results['delta'] = params['delta']
        results['tau'] = params['tau']
        results['pnt_estimate'] = params['pnt_estimate']
        results['alpha_ovr_2'] = params['alpha_ovr_2']
        results['M_div'] = params['M_div']
        
        if pval > (1-params['alpha_ovr_2']):
            results['chernoff_high_pos'] = params['chernoff_high_pos']
            results['chernoff_low_pos'] = params['chernoff_low_pos']
            results['classification'] = 'e'

        else:
            results['chernoff_high_neg'] = params['chernoff_high_neg']
            results['chernoff_low_neg'] = params['chernoff_low_neg']
            results['classification'] = 'i'

           
            
        if return_normal:
            if pval > (1-params['alpha_ovr_2']):
                results['ci_normal_high_pos'] = params['ci_normal_high_pos']
                results['ci_normal_low_pos'] = params['ci_normal_low_pos']
                results['classification'] = 'e'
            else:
                results['ci_normal_high_neg'] = params['ci_normal_high_neg']
                results['ci_normal_low_neg'] = params['ci_normal_low_neg']
                results['classification'] = 'i'

        if return_exact:
            if pval > (1-params['alpha_ovr_2']):
                results['ci_exact_high_pos'] = params['ci_exact_high_pos']
                results['ci_exact_low_pos'] = params['ci_exact_low_pos']
                results['classification'] = 'e'
            else:
                results['ci_exact_high_neg'] = params['ci_exact_high_neg']
                results['ci_exact_low_neg'] = params['ci_exact_low_neg']
                results['classification'] = 'i'

    else:
        results['Result:'] = 'judged unconnected with pvalue = ' + str(pval)
        results['pnt_estimate'] = params['pnt_estimate']
        results['classification'] = 'u'
    return results
    
    
def preprocess_spike_trains(R,T,Delta,delta,tau): # input in units of samples
    params = {} # Define params that will aggregate and flow through functions
    params['Delta'] = Delta
    params['delta'] = delta
    params['tau'] = tau
    
    T = T-tau
    #syn_rng = np.arange(-delta/2,delta/2+1,1) # ZERO CENTERED SINCE I'M SHIFTING THE TARGET TRAIN
    syn_rng = np.arange(0,delta+1,1).astype('int') # add as a parameter in the generation stage
    Rrepeat = np.repeat(R, len(syn_rng))
    syn_regions = Rrepeat + np.tile(syn_rng, len(R))
    syn_region_label = np.repeat(np.arange(0,len(R),1), len(syn_rng))
    SR = np.unique(syn_regions)

    synT = np.intersect1d(T,SR) # Synchronous target spikes
    nosynT = np.setdiff1d(T,SR) # non-synchronous target spikes
    fulll = (np.max(np.append(T,SR))// Delta).astype(int)
    fup = (np.max(np.append(T,SR))).astype(int)
    lab_syn_region_no_target = np.setdiff1d(syn_region_label,syn_region_label[np.isin(syn_regions,T)]) # These syn region labels have no target spike at all
    non_synch_ref = np.unique(Rrepeat[np.isin(syn_region_label,lab_syn_region_no_target)]) # \vect G

    bin_edges = np.arange(0,np.max(fup)+Delta,Delta)
    valid_prob = np.histogram(SR, bins=bin_edges)[0]/Delta
    NkS = np.histogram(synT, bins=bin_edges)[0]
    NkT = np.histogram(T, bins=bin_edges)[0]
    
    params['pnt_estimate'] = np.nansum((NkS-valid_prob*NkT)/(1-valid_prob))
    
    if np.sum(valid_prob==1) > 0:
        logging.warning('Delta and delta are too close relative the presynaptic firing rate in these data, leading to a violation of Timescale separation. The point estimate ignored these cases to be well-defined, and the confidence intervals are still defined, but you might consider a different choice of Delta and delta to be further apart.')
    
    SRbi = (SR[SR>0] // Delta).astype(int)
    max_value = (np.max(SRbi) + 2*Delta).astype('int')
    valid_prob = np.bincount(SRbi, minlength=fulll+1)/Delta
    synTbin = (synT // Delta).astype(int)
    nosynTbin = (nosynT // Delta).astype(int)
    ghostbin = (non_synch_ref // Delta).astype(int)

    
    
    # possible to get errors below!: need to change valid_prob minlength argument of bincount
    q_syn = valid_prob[synTbin]
    q_syn = q_syn[q_syn>0]
    q_syn = np.sort(q_syn)
    q_asyn = valid_prob[nosynTbin]
    q_asyn = q_asyn[q_asyn>0] # target spikes alone in the whole window have zero contribution to variability
    q_ghost = valid_prob[ghostbin]
    q_ghost = np.sort(q_ghost)

    cardS = len(synT)
    cardT = cardS + len(nosynT)
    cardE = len(q_ghost)
    
    params['cardS'] = len(synT)
    params['cardT'] = cardS + len(nosynT)
    params['cardE'] = len(q_ghost)
    
    params['q_syn'] = q_syn
    params['q_asyn'] = q_asyn
    params['q_ghost'] = q_ghost
    
    return params


def chernoff_bound_left_tail(tt,p_vect):
    mm = np.sum(p_vect)
    val = np.exp(-1*mm + tt + tt*np.log(mm/tt))
    if mm < tt:
        val = 1 - val
    return val

def chernoff_bound_right_tail(tt,p_vect):
    mm = np.sum(p_vect)
    val = np.exp(-1*mm + tt + tt*np.log(mm/tt))
    if mm > tt:
        val = 1 - val
    return val


def prepare_convolution_power(ptem,M_div):
    # compute for this tail, the unique bernoulli's
    up,ucnts = np.unique(ptem,return_counts=True)
    ix = np.argsort(ucnts) # this sort is probably not necessary... it was just mentally helpful?
    ucnts = ucnts[ix]
    up = up[ix]

    # choose a M_div (subdivide the unique's into a blob that can easily obtained via M_div fold convolution)
    #M_div = 2**5 # Number of blobs to initiate FFT with (the computation time is a parabola and 2**5 seems optimal in large example)
    size_blobs_before_fft = np.floor(np.sum(ucnts)/M_div) # approximate
    proportion_counts = ucnts/np.sum(ucnts) # the approximate proportion of each unique binom in the data
    freq_unique_bern_per_blob = np.floor(size_blobs_before_fft*proportion_counts) # apply M_div fold convolution to this blob
    temp = (freq_unique_bern_per_blob*M_div)
    left_over_bernoullis_per_unique = ucnts - temp # append these with direct convolution post hoc (and append them with those obtained via more hypotheses


    # Construct the fundemental blob unit of FFT you can compute as an "iid slice" of the total
    tempn = freq_unique_bern_per_blob[freq_unique_bern_per_blob!=0]
    tempp = up[freq_unique_bern_per_blob!=0]
    temp_blob = binom(tempn[0],tempp[0]).pmf(np.arange(0,tempn[0]+1,1))
    for k in range(1,len(tempn)):
        temp_append =  binom(tempn[k],tempp[k]).pmf(np.arange(0,tempn[k]+1,1))
        temp_blob = np.convolve(temp_blob,temp_append)
    
    temp_blob_log = np.zeros_like(temp_blob)
    temp_blob_log[temp_blob==0] = NEG_INF
    temp_blob_log[temp_blob!=0] = np.log(temp_blob[temp_blob!=0])
    # prepare the leftovers
    tempnn = left_over_bernoullis_per_unique
    temppp = up
    #add_bern0 = []
    #for k in range(len(left_over_bernoullis_per_unique)):
    #    add_bern0.append(temppp[k]*np.ones(int(tempnn[k])))
    #add_bern0 = np.concatenate(add_bern0)
    
    temp_bern = binom(tempnn[0],temppp[0]).pmf(np.arange(0,tempnn[0]+1,1))
    for k in range(1,len(tempnn)):
        temp_append =  binom(tempnn[k],temppp[k]).pmf(np.arange(0,tempnn[k]+1,1))
        temp_bern = np.convolve(temp_bern,temp_append)
    

    return temp_blob_log, temp_bern




def binary_search_pivot_for_bound_chernoff(params,R_0,L_0,is_pos_h,is_left_tail):
    if is_pos_h:
        p =  params['q_syn']
        pp = params['q_asyn']
    else:
        p =  params['q_ghost'] 
        pp = np.append(params['q_syn'],params['q_asyn'] )  
    L = L_0
    R = R_0
    kk = 0
    while L < R:
        h = (L + R) // 2
        if is_pos_h:
            ncnv = params['cardS'] - h # positive theta
        else: 
            ncnv = np.abs(h) # negative theta
        if is_left_tail: # left tail
            ptem = np.append(p[:ncnv],pp)
            pval = chernoff_bound_left_tail(params['cardS']-h,ptem)  
            if pval > params['alpha_ovr_2']:
                L = h + 1
            else:
                R = h   
        else: # right tail
            ptem = np.append(np.flipud(p)[:ncnv],pp)  #flipud is not necessary but easier for me mentally, consider changing later to -ncnv:
            pval = chernoff_bound_right_tail(params['cardS']-h,ptem)  
            if pval > params['alpha_ovr_2']:
                R = h 
            else:
                L = h + 1        
        kk += 1
    logging.warn(str(kk))
    if (is_left_tail == True):
        ci = L - 1 
        
    if  (is_left_tail == False):
        ci = L
    return ci



def binary_search_pivot_for_bound_normal(params):
    if params['is_pos_h']:
        ncnv_hi = params['cardS'] - params['chernoff_low_pos']
        ncnv_low = params['cardS'] - params['chernoff_high_pos']
        p =  params['q_syn']
        pp = params['q_asyn']
        L = 0
        R = params['chernoff_high_pos']-params['chernoff_low_pos']
    else:
        ncnv_hi = np.abs(params['chernoff_high_neg'])
        ncnv_low = np.abs(params['chernoff_low_neg'])
        p =  params['q_ghost']
        pp = np.append(params['q_asyn'],params['q_syn'])
        L = -1*(ncnv_hi-ncnv_low)
        R = 0

    if params['is_left_tail']:
        p_2 = p
    else:
        p_2 = np.flipud(p)
    # Pre-prepare for exact convolution
    params['pivot'] = p_2[ncnv_low:ncnv_hi] # Anything above ncnv_hi are gauranteed to be injected cause chernoff says we need at least that much: so we don't need to convolve them
    params['to_blob'] = np.hstack([p_2[:ncnv_low],pp]) # on the other hand we know for sure everything below ncnv_low will be convolved, so use a fast method to exploit that.
    #params['blob_log'], params['add_bern'] = prepare_convolution_power(params['to_blob'],params['M_div'])
    while L < R:
        h = (L + R) // 2
        if params['is_pos_h']:
            ncnv = len(params['pivot']) - h # positive theta
        else: 
            ncnv = np.abs(h) # negative theta

        params['temp_pivot'] = params['pivot'][:ncnv]
        if params['is_pos_h']:
            thta_temp = (params['chernoff_low_pos'] + h)
            params['cgfO'] = params['cardS'] - thta_temp
        else:
            thta_temp = params['chernoff_low_neg'] + h 
            params['cgfO'] = params['cardS'] - thta_temp    
        #pval = keich_pval_all(params)
        ptem = np.append(params['to_blob'],params['temp_pivot'])
        if params['is_left_tail']: # left tail
            pval = normal_approximation(ptem,params['cgfO'],is_cdf=True)
            if pval > params['alpha_ovr_2']:
                L = h + 1
            else:
                R = h   
        else: # right tail
            pval = normal_approximation(ptem,params['cgfO'],is_cdf=False)
            if pval > params['alpha_ovr_2']:
                R = h 
            else:
                L = h + 1  
 
        
    if (params['is_left_tail'] == True):
        ci = L - 1 

    if  (params['is_left_tail'] == False):
        ci = L
    
    if params['is_pos_h']:
        ci = params['chernoff_low_pos'] + h
    else:
        ci = params['chernoff_low_neg'] + h  
    
    return ci


def normal_approximation(ptem,obs,is_cdf):
    mutem = np.sum(ptem)
    sigmatem = np.sqrt(np.sum(ptem*(1-ptem)))
    normal = norm(loc=mutem, scale=sigmatem);
    if is_cdf:
        pval = normal.cdf(obs)
    else:
        pval = 1-normal.cdf(obs)
    return pval

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



def binary_search_pivot_for_bound_exact_modified(params):
    if params['is_pos_h']:
        ncnv_hi = params['cardS'] - params['chernoff_low_pos']
        ncnv_low = params['cardS'] - params['chernoff_high_pos']
        p =  params['q_syn']
        pp = params['q_asyn']
        L = 0
        R = params['chernoff_high_pos']-params['chernoff_low_pos']
    else:
        ncnv_hi = np.abs(params['chernoff_high_neg'])
        ncnv_low = np.abs(params['chernoff_low_neg'])
        p =  params['q_ghost']
        pp = np.append(params['q_asyn'],params['q_syn'])
        L = -1*(ncnv_hi-ncnv_low)
        R = 0

    if params['is_left_tail']:
        p_2 = p
    else:
        p_2 = np.flipud(p)
    # Pre-prepare for exact convolution
    params['pivot'] = p_2[ncnv_low:ncnv_hi] # Anything above ncnv_hi are gauranteed to be injected cause chernoff says we need at least that much: so we don't need to convolve them
    params['to_blob'] = np.hstack([p_2[:ncnv_low],pp]) # on the other hand we know for sure everything below ncnv_low will be convolved, so use a fast method to exploit that.
    params['blob_log'], params['add_bern'] = prepare_convolution_power(params['to_blob'],params['M_div'])

    while L < R:
        h = (L + R) // 2
        if params['is_pos_h']:
            ncnv = len(params['pivot']) - h # positive theta
        else: 
            ncnv = np.abs(h) # negative theta

        params['temp_pivot'] = params['pivot'][:ncnv]
        if params['is_pos_h']:
            params['cgfO'] = params['cardS'] - (params['chernoff_low_pos'] + h)
        else:
            params['cgfO'] = params['cardS'] - params['chernoff_low_neg'] - h    
        
        pval = keich_pval_modified(params)
        if params['is_left_tail']: # left tail
            if pval > params['alpha_ovr_2']:
                L = h + 1
            else:
                R = h   
        else: # right tail
            if pval > params['alpha_ovr_2']:
                R = h 
            else:
                L = h + 1  
    
    if (params['is_left_tail'] == True):
        ci = L - 1 

    if  (params['is_left_tail'] == False):
        ci = L
    
    if params['is_pos_h']:
        ci = params['chernoff_low_pos'] + h
    else:
        ci = params['chernoff_low_neg'] + h  
    
    return ci





def keich_pval_modified(params):
    if params['is_left_tail']: # Compute the random variables -X
        cgfP = 1-np.hstack([params['to_blob'],params['temp_pivot']])
        log_blob = np.flipud(params['blob_log'])
        add_bern = 1-params['temp_pivot']
        pre_comp = np.flipud(params['add_bern'])
        Lfold = params['M_div']
        cgfO = (1+len(cgfP)) - params['cgfO']
    else:
        cgfP = np.hstack([params['to_blob'],params['temp_pivot']])
        log_blob = params['blob_log']
        add_bern = params['temp_pivot']
        pre_comp = params['add_bern']
        Lfold = params['M_div']
        cgfO = params['cgfO']
   
    for k in range(len(add_bern)):
        pre_comp = np.convolve(pre_comp,np.array([1-add_bern[k],add_bern[k]]))
    all_else_log = np.log(pre_comp)
    
    #tempp,tempn = np.unique(add_bern,return_counts=True)
    #for k in range(len(tempn)):
    #    temp_append =  binom(tempn[k],tempp[k]).pmf(np.arange(0,tempn[k]+1,1))
    #    temp_blob = np.convolve(temp_blob,temp_append)
    
    
    def Lfoldshift(s): # cumulant generating function set equal to observed value
        return np.sum((np.exp(s)*cgfP)/((1-cgfP)+np.exp(s)*cgfP)) - cgfO

    shat = brentq(Lfoldshift,-50,50,rtol=10**-15,xtol=10**-200) # find root of Lfoldshift
    
    log_blob_shifted, log_blob_cgf = utilz.shift(log_blob, shat) # shift the log of the big blob distribution
    blob_shifted = np.exp(log_blob_shifted)
    blob_shifted = blob_shifted/np.sum(blob_shifted)

    
    all_else_log_shifted, all_else_cgf = utilz.shift(all_else_log, shat) # shift the log of the big blob distribution
    all_else_shifted = np.exp(all_else_log_shifted)
    all_else_shifted = all_else_shifted/np.sum(all_else_shifted)
    
    # Raw FFT with no shift
    nspt = (Lfold)*(len(blob_shifted)-1) + 1 + (len(all_else_shifted) - 1)# the size of all pmfs you've convolving without zero and then adding one zero
    K = np.ceil(np.log2(nspt))
    dimension_of_DFT_operator = int(2**K)
    F1 = np.zeros(dimension_of_DFT_operator)
    F1[:len(blob_shifted)] = blob_shifted
    
    F2 = np.zeros(dimension_of_DFT_operator)
    F2[:len(all_else_shifted)] = all_else_shifted
    
    Fpmf = np.fft.fft(F2)*np.fft.fft(F1)**(Lfold)
    convolved = np.abs(np.fft.ifft(Fpmf))
    resultlog = np.log(convolved) - shat*np.arange(len(convolved)) + log_blob_cgf*(Lfold) + all_else_cgf
    
    resultlog = resultlog[:nspt]
    resultlog[resultlog>1] = 0
    fft_result = np.exp(resultlog)

    tail_prob =  np.sum(fft_result[cgfO:])
    return tail_prob # the precise p-value



