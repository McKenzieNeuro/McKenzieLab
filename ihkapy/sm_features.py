# the feature funcitons
# this is a utility script that contains all the functions the compute features
# each function in this file takes a 2d numpy np.ndarray as input
# and returns a float as output: the float, or a 1darray of features
# In the case of returning a 1d array, those features are computed for each 
# channel individually, so the 1darray has length n_channels
# In shape is (len_window,n_channels)

# example, though I think useless because zscore normalizes this stuff
def rms(window:np.ndarray) -> np.ndarray:
    """Returns Root Mean Squared of each channel."""
    return np.sqrt(np.power(a,2).mean(axis=0))





