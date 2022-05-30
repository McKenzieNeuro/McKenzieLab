import numpy as np
from scipy.signal import coherence
from itertools import combinations as combin


# the feature funcitons
# this is a utility script that contains all the functions the compute features
# each function in this file takes a 2d numpy np.ndarray as input
# and returns a float as output: the float, or a 1darray of features
# In the case of returning a 1d array, those features are computed for each 
# channel individually, so the 1darray has length n_channels
# In shape is (len_window,n_channels)


"""
Each feature is a function that takes an ordered list of N_RAW_CHAN 
windows: the index of the window is the index of the raw electrode 
channel. The feature function returns a dictionary with key-value pairs
    key   : str   = name of sub-feature
    value : float = value of the feature evaluated
The keys of this returned dictionary are to be used a column names for
the data-frame in which we save our features.

Example 1. Mean Power raw. 
Takes four, 5-second windows of data
For each window (all windows are same point in time but correspond to 
different electrodes recordings)
    Compute the mean power of only the raw channel (i.e. not wavelets),
    Save this do dictionary under key=f"mean_power_ch{n_raw_chan}"

Example 2. Correlation feature.


"""

"""
Convention:
    ch_raw is the raw channel index
    chan_freq is the binary (wavelet freq) file channel index
"""

def _select_freq_chans(windows,chan_freq_idxs) -> (np.ndarray,np.ndarray):
    """Helper for all features. 

    Selects frequency (wavelet binary) channels of interest and formats 
    chan_freq_idxs as a 1d numpy ndarray.
    """
    # If None, select all frequency channels
    windows = np.asarray(windows)
    if chan_freq_idxs is None: 
        chan_freq_idxs = np.arange(windows[0].shape[1])
    elif type(chan_freq_idxs) == type(1): # int
        chan_freq_idxs = np.asarray([chan_freq_idxs])
        windows = windows[:,:,chan_freq_idxs] 
    else:
        chan_freq_idxs = np.asarray(chan_freq_idxs)
        windows = windows[:,:,chan_freq_idxs] 
    return windows,chan_freq_idxs

def _feats_as_dict_from_2d_array(
        feats_all : np.ndarray,
        chan_freq_idxs : np.ndarray
        ):
    """Helper for single-channel features.

    Turns a 2d feats array with dims = (n_raw_chans , n_freq_chans)
    into our standard dictionary format that can easily be added to a
    pandas dataframe.
    """
    feats_dict = {}
    for ch_raw,feats_ch_raw in enumerate(feats_all):
        for chb,feat in zip(chan_freq_idxs,feats_ch_raw):
            feats_dict[f"mean_power_chraw_{ch_raw}_chfreq_{chb}"] = feat
    return feats_dict


### Features on individual channels

def mean_power(
        windows : list or np.ndarray,
        chan_freq_idxs : np.ndarray or int = None
        ) -> dict:
    """Get mean power across each channel.

    Parameters
    ----------
    windows : list
        A list of windows. Dims = (n_raw , n_samples , n_freq_chan)

    chan_freq_idxs : ndarray or int = None
        The 1d list of binary file indices we are interested in, e.g.
        if we're interested in all wavelet channels (aka binary file 
        channels) pass None. If we're interested in the first two we 
        can pass [0,1] or np.array([0,1]). If we're only interested in
        the raw data we can simply pass an int, 0 in this case. 

    Returns
    -------
    dict
        A dictionary where the key,value pairs are `name of feature`
        and `computed feature float`. All keys are strings, all values
        are floats.
    """
    # Select binary (aka wavelet/freq) file channels
    windows,chan_freq_idxs = _select_freq_chans(windows,chan_freq_idxs)
    # Compute the features along samples axis
    mp_all = np.mean(windows,axis=1)
    # Save and return, formatted as dict
    mp_feats = _feats_as_dict_from_2d_array(mp_all,chan_freq_idxs)
    return mp_feats


def var(
        windows         : list or np.ndarray,
        chan_freq_idxs  : np.ndarray or int = None
        ) -> dict:
    """Get the variance of the signal"""
    windows,chan_freq_idxs = _select_windows(windows,chan_freq_idxs)
    var_all = np.var(windows,axis=1)
    var_feats = _feats_as_dict_from_2d_array(mp_all,chan_freq_idxs)
    return var_feats

    
### Features comparing two channels

def coherence(
        windows         : list or np.ndarray,
        chan_freq_idxs  : np.ndarray or int = None,
        fs             : float or int # = 2000
        ) -> dict:
    """Cohere each pair of signals"""
    windows,chan_freq_idxs = _select_windows(windows,chan_freq_idxs)
    raw_chans = list(range(len(windows)))    
    coherence_dict = {}
    for (chx,chy),(x,y) in zip(combin(raw_chans),combin(windows)):
        # There parameters defined here are similar to the default 
        # MatLab options
        sample_freqs,coh_xy = coherence(x,y,fs=fs,window='hann',
                nperseg256,noverlap=None,nfft=None,
                detrend='constant',axis=0)
        # NB: FFT and therefore scipy.signal.coherence samples 
        # frequencies linearly, not logarithmically

        # TODO: continute here
        
        # TODO: verify that scipy's coherence outputs the dims you're
        # expecting. 






















