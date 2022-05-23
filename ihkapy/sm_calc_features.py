import numpy as np

"""
MatLab code is structures as follows:
    sm_PredictIHKA_getAllFeatures.m 

"""

# pure function
def _get_x_pct_time_of_interval(
        start_time : float,     # in seconds
        end_time : float,       # in seconds
        window_length : float,  # in seconds
        pct : float             # proportion of times to sample
        ) -> np.ndarray:
    """Get a randomly sampled 1d array of start times

    Parameters
    ----------
    start_time : float
        The beginning timestamp, in seconds, of the interval we sample from.
    end_time : float
        The end timestamp, in seconds, of the interval we sample from.
    window_length : float
        The time, in seconds of our sample windows.
    pct : float
        The proportion of windows to select, must be between 0 and 1. 

    Returns
    -------
    np.ndarray
        A 1d numpy array of start times for the windows. 
        Together with the window length, these fully define the windows
        of interest to us that we would like to sample from. 
    """
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
    window_start_times = np.random.choice(window_start_times,n_select,replace=False)
    window_start_times.sort()
    return window_start_times


# # pure function
# def _get_random_sample_window_start_times(
#         seizure_start_times : list,
#         seizure_end_times : list,
#         session_length : float,
#         window_length : float
#         ) -> list:
#     """Generates random window samples from seizure metadata for single seizure.
# 
#     Parameters
#     ----------
#     seizure_start_times : list
#         List of floats. The start and end time of each seizure, in seconds.  
#     seizure_end_times : list
#         List of floats. The end time of each seizure, in seconds.
#     session_length : float
#         Length of the whole session, in seconds. (order of 1 day ~500_000 seconds)
#     window_length : float
#         Length of each window, in seconds. (order of 5s)
# 
#     Returns
#     -------
#     list
#         List of start times for windows to be sampled.
#     """
#     assert len(seizure_start_times) == len(seizure_end_times)
#     for i,j in zip(seizure_start_times,seizure_end_times): assert i<j
# 
#     # TODO: Q: Do we need to preserve metadata about features we compute? Labels?
#     # I think so... It's hard to know because all that is done with indices in MLab
# 
#     return # dummy


def get_feats_window(
        window : np.ndarray,
        featurelist : list) -> np.ndarray:
    """Computes all features specified in featurelist."""
    # featurelist is a list of functions. 
    # Open features.py, (it's a script that contains all our feature functions)

    # If featurelist == "all", select all the features defined in features.py

    # Check that all our features in featurelist correspond to names of 
    # feature functions in features.py

    # 
    return 

# Helper for calc_features()
def _get_bin_absolute_intervals(
        seizure_start_time      : float,
        seizure_end_time        : float,
        preictal_bins           : list,
        postictal_delay         : int,
        bin_names               : list,
        all_start_times         : list,
        all_end_times           : list,
        total_recoreding_time   : float
        ):
    """Returns dictionary of all valid bin intervals for single seizure.

    Parameters
    ----------
    start_seiz : float or int
        The start time of the seizure in question, in seconds. (Relative
        to the start of the recording.)

    end_seiz : float or int
        The end time of the seizure in questionm, in seconds. (Relative
        to the start of the recording.)

    preictal_bins : list
        The list found in Options.toml config file. It's a list of floats
        that specify, in seconds, the pre-ictal bin timestamps.

    postictal_delay : int or float
        Specifies post-ictal period after the end of a seizure in seconds
        (e.g. 600)

    bin_names : list
        A list of bin names (strings). Must have length of preictal_bins + 2. 
        This is because we have one single intra-ictal and one single post-
        ictal bin. 

    all_start_times : list
        A list of the start times of every seizure (in seconds from start 
        of recording). Needed to check validity of intervals.

    all_end_times : list 
        A list of all end times of every seizure (in seconds from start of
        recording). Needed to check validity of intervals. 

    total_recording_time : float,
        Total number of seconds from start to end of the sessions.
    
    Returns
    -------
    dict
        The keys are names of time bins (e.g. pre_ictal_bin1, post_ictal, etc.)
        If the corresponding bin for the seizure is valid, the value will be a 
        tuple (bin_start_time,bin_end_time). If the bin is not valid, the value
        will be None. A siezure's time bin is considered valid if all of these
        conditions are satisfied:
            1. The bin does not over-lap with any other seizure.
            2. The bin starts after the beginning of the recording. 
            3. The bin ends before the end of the recording. 
    """
    # Some tests and checks
    assert len(preictal_bins) + 2 == len(bin_names), "Args have incompatible length"
    assert len(start_times) == len(end_times), "Logic err, smthng is terribly wrong"

    bin_name_intraictal = bin_names[-2]     # string
    bin_name_postictal = bin_names[-1]      # string
    bin_names_preictal = bin_names[:-2]     # list of strings
    
    # Init bin_intervals (the dict we will return)
    # NB: The seizure interval is always valid
    bin_intervals = {bin_name_intraictal : (start_seiz,end_seiz)} 

    # The absolute time of pre-ictal intervals, we iterate through them
    pre_bins_abs = [start_time - i for i in preictal_bins] 
    # Important, don't forget this final implied pre-ictal slice 
    # (the start of seizure)
    pre_bins_abs.append(start_time) 
    for bin_name,start_bin,end_bin in zip(bin_names_preictal,pre_bins_abs[:-1],pre_bins_abs[1:]):
        # TODO: compare with MatLab
        #       Notably verify if there are missing checks

        # Check whether another other seizure overlaps with our bin
        bin_intervals[bin_name_intraictal] = (start_bin,end_bin)
        for start_seiz_2,end_seiz_2 in zip(all_start_times,all_end_times):
            # TODO: check MatLab confirm that postictal delay term necessary
            if start_bin < end_seiz_2 + postictal_delay 
                and start_seiz > start_seiz_2 
                and start_bin > 0:
                # The pre-ictal bin interval is not valid
                # Re-define it as None and break out of for-loop
                bin_intervals[bin_name_intraictal] = None
                break 

    # If the post-ictal interval is valid, define it.
    # Post-ictal is valid iff that there is no second seizure right 
    # after the one in question. 
    end_postictal = end_seiz + postictal_delay
    bin_intervals[bin_name_postictal] = (end_seiz , end_postictal)
    # Check if invalid: redefine to None
    if end_postictal > total_recording_time: 
        bin_intervals[bin_name_postictal] = None
    for start_seiz_2 in all_start_times:
        if end_postictal < start_seiz_2 and end_seiz > start_seiz_2:
            bin_intervals[bin_name_postictal] = None
    return bin_intervals


def _get_total_recording_time(filepath,fs,num_bin_chan=41,precision="int16"):
    """Determines the total recording time from binary files."""
    fsize_bytes = os.path.getsize(filepath)
    bytes_per_sample = np.dtype(precision).itemsize
    n_samples_per_chan = fsize_bytes / bytes_per_sample / num_bin_chan
    assert n_samples_per_chan == int(n_samples_per_chan) , "Logic error, possibly num_bin_chan is incorrect: this the number of channels saved in the binary file."
    total_duration_seconds = n_samples_per_chan / fs # total duration in seconds
    return total_duration_seconds


def calc_features(
        binary_basename_list : list = [],
        fileio : dict,
        data_ops : dict,
        feature_ops : dict,
        ) -> np.ndarray:
    """Compute all features in all channels for multiple sessions.

    Several noteworthy assumptions are made: It is assumed that the
    binary data file and the metadata text file share the same basename.
    (Asside: metadata txt file is parsed by parse_metadata in metadata_io.py)
    
    Parameters
    ----------

    Returns
    -------
    dict
        The keys are the bin names to which the features belong, and the
        values are pandas dataframes with columns { sessions_basename , 
        start_time , feat_01_name , feat_02_name , ... }
        
    dict
        Options that where used to generate the data. This could be a 
        serialized copy of the Options.toml file. 
    """
    # feature_functions = # TODO: this is a list of functions
    # features = # TODO: 2d np array shape = (n_windows,n_feats_per_window)

    # Unpack relevant parameters
    BIN_NAMES       = feature_ops["BIN_NAMES"]
    PREICTAL_BINS   = feature_ops["PREICTAL"]["BINS"] # in (negative) seconds
    POSTICTAL_DELAY = feature_ops["POSTICTAL"]["DELAY"] # in seconds
    PREICTAL_PCT    = feature_ops["PREICTAL"]["PCT"]
    # unpack other stuff
    bins = {} # TODO: define this, keys=bin_names, values=relative_intervals
              # relative interval is a tuple (rel start, rel end)
              # should be negative for pre-ictal

    data = {b:[] for b in BIN_NAMES} # TODO: turn list into pandas dataframe

    for session_basename in binary_basename_list:
        # TODO: implement get_seizure_start_end_times()
        start_times,end_times = get_seizure_start_end_times(session_basename,fileio)
        total_recording_time = _get_total_recording_time(basename,precision) # in seconds
        for start_time,end_time in zip(start_times,end_times):
            # TODO: define absolute start and end intervals for each bin
            # TODO: test _get_bin_absolute_intervals()
            #       returns a dic with keys=BIN_NAMES,values=(strtbin,endbin) (in s)
            bins_absolute_intervals = _get_bin_absolute_intervals(
                    seizure_start_time  = start_time,
                    seizure_end_time    = end_time,
                    preictal_bins       = PREICTAL_BINS, 
                    posictal_delay      = POSTICTAL_DELAY,
                    all_start_times     = start_times,
                    all_end_times       = end_times
                    total_recording_time = # TODO get this one
                    )
        
            for bin_name,interval in bins_absolute_intervals:
                # TODO: get the windows xPctTims thingy
                # TODO: compute all features for this interval
                # TODO: add a row to the data[bin_name] pandas data-frame

    # Read data params
    # Do some checks on the data and metadata params
    # call _calc_features
    return # the output of _calc_features



# # in pure functions.
# def _calc_features(
#         fname_channels,
#         window_start_times,
#         fs : float = 2000,
#         # bunch of options
#         ) -> np.ndarray:
#     """
#     fname_chan_bin_paths : list
#         Ordered list of paths to the filenames of the binary .dat files
#         for each channel. NB: each channel has it's own binary.
#     window_start_times : list
#         List of start times of all windows we will sample.
#         Same as the 'tim' cell array in MatLab
#     fs : float
#         Sample frequency in Hz. Defaults to 2000Hz
#     
#     
# 
#     Returns
#     -------
#     np.ndarray
#         2d array shape (n_windows,n_features)
#     """
#     # Check that all the files specified exist
#     for path in fname_chan_bin_paths: assert os.path.exists(path), "File Not Found"
#     # TODO: For now we implement as it is implement as is in MatLab
#     #       but it may be an idea to change the channel input from a list
#     #       to a dictionary, so that we can encode more information about
#     #       the channel if need be. May be useful for training. 
# 
#     return 


"""
PSEUDO-CODE

Based on pre-ictal and post-ictal hyper-parameters, seizure start and end times,
and the boundaries of the recordings, randomly generate list (or dict for labels?)
of window-start times to sample and compute features for.
"""


if __name__=="__main__":
    ### TEST _get_x_pct_time_of_interval()
    arr = _get_x_pct_time_of_interval(
            start_time    = 5.0,
            end_time      = 152.6,
            window_length = 1.0,
            pct           = 0.05)
    print(arr)




