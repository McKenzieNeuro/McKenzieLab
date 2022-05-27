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


This approach is the best I could come up with. 

"""



# example, though I think useless because zscore normalizes this stuff
def rms(window:np.ndarray) -> np.ndarray:
    """Returns Root Mean Squared of each channel."""
    return np.sqrt(np.power(a,2).mean(axis=0))





