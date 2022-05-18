

"""
MatLab code is structures as follows:
    sm_PredictIHKA_getAllFeatures.m 

"""

# pure function
def _get_x_pct_time_of_interval(
        start_time : float      # in seconds
        end_time : float        # in seconds
        window_length : float   # in seconds
        pct : float             # proportion of times to sample
        ) -> list:
    """Returns an array of start times"""
    assert end_time > start_time
    assert pct >= 0.0 and pct <= 1.0
    # The number of windows that fit in the interval
    n_windows = int((end_time - start_time) // window_length) 
    # List of all start times of max number of non-overlapping windows
    # that fit in the interval.
    window_start_times = np.linspace(start_time,end_time-window_length,n_windows)
    if pct == 1.0: return all_window_start_times
    # Choose and sort a random sample according to pct
    n_select = int(np.ceil(n_windows * pct)) # num win to select
    window_start_times = np.choice(window_start_times,n_select,replace=False)
    window_start_times.sort()
    return window_start_times


# pure function
def _get_random_sample_window_start_times(
        seizure_start_times : list,
        seizure_end_times : list,
        session_length : float,
        window_length : float
        ) -> list:
    """Generates random window samples from seizure metadata.

    Parameters
    ----------
    seizure_start_times : list
        List of floats. The start and end time of each seizure, in seconds.  
    seizure_end_times : list
        List of floats. The end time of each seizure, in seconds.
    session_length : float
        Length of the whole session, in seconds. (order of 1 day ~500_000 seconds)
    window_length : float
        Length of each window, in seconds. (order of 5s)

    Returns
    -------
    list
        List of start times for windows to be sampled.
    """
    assert len(seizure_start_times) == len(seizure_end_times)
    for i,j in zip(seizure_start_times,seizure_end_times): assert i<j

    # TODO: Q: Do we need to preserve metadata about features we compute? Labels?
    # I think so... It's hard to know because all that is done with indices in MLab
    

    return # dummy


def calc_features_seizure(
        fname,
        start_times,
        ) -> np.ndarray:
    """Compute all features in all channels for one seizure."""
    # Get list of features that we want to use from a features file
    # each feature should be a function with standard in/out types
    # a segment (window) is the input, the output is a float
    windows = # TODO: 3d np array shape = (n_chan,n_windows,window_length)
    feature_functions = # TODO: this is a list of functions
    features = # TODO: 2d np array shape = (n_windows,n_feats_per_window)



    # Read data params
    # Do some checks on the data and metadata params
    # call _calc_features
    return # the output of _calc_features
    


# Why seperate these? 
# It's to make a distinction between pure and non-pure functions
# As much as possible we want our functions to be pure
# When our functions are not pure, because, for instance, they need
# to use user defined hyper-parameters, it is desirable to create a
# pure function that does all the heavy lifting, and then create a non-pure
# wrapper for it. Pure functions are easier to understand for people
# read the code, and it's good for the core functionality to be written
# in pure functions.
def _calc_features(
        fname_channels,
        window_start_times,
        fs : float = 2000,
        # bunch of options
        ) -> np.ndarray:
    """
    fname_chan_bin_paths : list
        Ordered list of paths to the filenames of the binary .dat files
        for each channel. NB: each channel has it's own binary.
    window_start_times : list
        List of start times of all windows we will sample.
        Same as the 'tim' cell array in MatLab
    fs : float
        Sample frequency in Hz. Defaults to 2000Hz
    
    

    Returns
    -------
    np.ndarray
        2d array shape (n_windows,n_features)
    """
    # Check that all the files specified exist
    for path in fname_chan_bin_paths: assert os.path.exists(path), "File Not Found"
    # TODO: For now we implement as it is implement as is in MatLab
    #       but it may be an idea to change the channel input from a list
    #       to a dictionary, so that we can encode more information about
    #       the channel if need be. May be useful for training. 

    return 



"""
PSEUDO-CODE

Based on pre-ictal and post-ictal hyper-parameters, seizure start and end times,
and the boundaries of the recordings, randomly generate list (or dict for labels?)
of window-start times to sample and compute features for.

"""






