import numpy as np
from ihkapy.fileio.utils import get_all_valid_session_basenames,check_session_basenames_are_valid
from ihkapy.fileio.metadata_io import get_seizure_start_end_times
from scipy.signal import coherence


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
        start_seiz              : float,
        end_seiz                : float,
        preictal_bins           : list,
        postictal_delay         : int,
        bin_names               : list,
        all_start_times         : list,
        all_end_times           : list,
        total_session_time      : float
        ):
    """Returns dictionary of all valid bin intervals for single seizure.

    The function is named 'absolute' because it measures the absolute 
    number of seconds from the begining of the session (aka recording/file)


    Parameters
    ----------
    start_seiz : float or int
        The start time of the seizure in question, in seconds. (Relative
        to the start of the session.)

    end_seiz : float or int
        The end time of the seizure in questionm, in seconds. (Relative
        to the start of the session.)

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
        of session). Needed to check validity of intervals.

    all_end_times : list 
        A list of all end times of every seizure (in seconds from start of
        session). Needed to check validity of intervals. 

    total_session_time : float,
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
            2. The bin starts after the beginning of the session. 
            3. The bin ends before the end of the session. 
    """
    # Some tests and checks
    assert len(preictal_bins) + 2 == len(bin_names), "Args have incompatible length."
    assert len(all_start_times) == len(all_end_times), "Logic err, start and end times lists must match in length."

    bin_name_intraictal = bin_names[-2]     # string
    bin_name_postictal = bin_names[-1]      # string
    bin_names_preictal = bin_names[:-2]     # list of strings
    
    # Init bin_intervals (the dict we will return)
    # NB: The seizure interval is always valid
    bin_intervals = {bin_name_intraictal : (start_seiz,end_seiz)} 

    ### PRE-ICTAL
    # The absolute time of pre-ictal intervals, (we will iterate through them)
    pre_bins_abs = [start_seiz - i for i in preictal_bins] 
    # Important, not to forget this final implied pre-ictal slice 
    # (the start of seizure)
    pre_bins_abs.append(start_seiz) 
    for bin_name_preictal,start_bin,end_bin in zip(bin_names_preictal,pre_bins_abs[:-1],pre_bins_abs[1:]):
        # TODO: compare with MatLab
        #       Notably verify if there are missing checks

        bin_intervals[bin_name_preictal] = (start_bin,end_bin)
        # Check whether another (different) seizure overlaps with our bin
        # if so, define it as none
        for start_seiz_2,end_seiz_2 in zip(all_start_times,all_end_times):
            # TODO: check MatLab confirm that postictal delay term necessary
            condition0 = start_bin < end_seiz_2 + postictal_delay 
            condition1 = start_seiz > start_seiz_2
            condition2 = start_bin < 0
            if (condition0 and condition1) or condition2:
                # The pre-ictal bin interval is not valid
                # Re-define it as None and break out of for-loop
                bin_intervals[bin_name_preictal] = None
                break 

    ### POST-ICTAL
    # If the post-ictal interval is valid, define it.
    # Post-ictal is valid iff that there is no second seizure right 
    # after the one in question. 
    end_postictal = end_seiz + postictal_delay
    bin_intervals[bin_name_postictal] = (end_seiz , end_postictal)
    # Check if invalid: redefine to None
    if end_postictal > total_session_time: 
        bin_intervals[bin_name_postictal] = None
    for start_seiz_2 in all_start_times:
        if end_postictal > start_seiz_2 and end_seiz < start_seiz_2:
            bin_intervals[bin_name_postictal] = None
    return bin_intervals


def _get_total_session_time(filepath,fs=2000,num_bin_chan=41,precision="int16"):
    """Determines the total session time from binary files."""
    fsize_bytes = os.path.getsize(filepath)
    bytes_per_sample = np.dtype(precision).itemsize
    n_samples_per_chan = fsize_bytes / bytes_per_sample / num_bin_chan
    assert n_samples_per_chan == int(n_samples_per_chan) , "Logic error, possibly num_bin_chan is incorrect: this the number of channels saved in the binary file."
    total_duration_in_seconds = n_samples_per_chan / fs # total duration in seconds
    return total_duration_in_seconds


# TODO: once completed, refactor this method
def calc_features(
        fileio : dict,
        data_ops : dict,
        feature_ops : dict,
        session_basenames_list : list or str = "all",
        ) -> np.ndarray:
    """Compute all features in all channels for multiple sessions.

    Several noteworthy assumptions are made: It is assumed that the
    binary data file and the metadata text file share the same basename.
    (Asside: metadata txt file is parsed by parse_metadata in metadata_io.py)
    
    Parameters
    ----------
    session_basenames_list : list or str
        If the list is empty, it will default to all valid edf files in the 
        directory. 

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

    # Unpack fileio params
    RAW_DATA_PATH = fileio["RAW_DATA_PATH"]

    # Unpack data params
    SCALE_PHASE     = feature_ops["SCALE_PHASE"]
    SCALE_POWER     = feature_ops["SCALE_POWER"]

    # Unpack relevant parameters, params.feature in Options.toml
    BIN_NAMES       = feature_ops["BIN_NAMES"]
    PREICTAL_BINS   = feature_ops["PREICTAL"]["BINS"] # in (negative) seconds
    POSTICTAL_DELAY = feature_ops["POSTICTAL"]["DELAY"] # in seconds
    PREICTAL_PCT    = feature_ops["PREICTAL"]["PCT"]
    INTRAICTAL_PCT  = feature_ops["INTRAICTAL"]["PCT"]
    POSTICTAL_PCT   = feature_ops["POSTICTAL"]["PCT"]
    PCT = PREICTAL_PCT + [INTRAICTAL_PCT] + [POSTICTAL_PCT] # Concatenation
    PCT_DIC = {b:pct for b,pct in zip(BIN_NAMES,PCT)} # a handy pct dictionary 
    DUR_FEAT        = feature_ops["DUR_FEAT"]
    amp_idx         = feature_ops["AMP_IDX"]
    amp_idx = np.array(amp_idx) * 2 + 1 
    ph_idx          = feature_ops["PH_IDX"]
    ph_idx  = np.array(ph_idx) *  2 + 2 


    ### BEGIN, TEMPORARY HACK
    # TODO: implement this properly so that it uses Options.toml instead
    # of hard coding so that user doesn't have to edit code. 
    # The ideal situation is to have a codebase robust enough that, most
    # of the time, the researcher doesn't have to worry about the code,
    # and when they do, they only have to write a new feature function in
    # the features module, and then add it's name to a list in the Options.toml
    # config file, and thats all—the code should work. 

    # Define all the feature functions
    # Features on each channel
    def mean_power(window) -> np.ndarray:
        """Returns mean across all Power channels. 
        Warning! Power channel indices hard-coded."""
        return np.mean(window[:,np.arange(1,40,2)],axis=0)/SCALE_POWER
    # Features on each pair of channels n*(n-1)/2
    def cohere_two_windows(window1,window2) -> np.ndarray:
        """Returns coherence across select power channels"""
        w1,w2 = window1[amp_idx],window2[amp_idx]
        coherence(w1, w2, fs=FS, window='hann', nperseg=256, noverlap=None, nfft=None, detrend='constant', axis=0)

        return # Dummy TODO: complete this method and continue


    
    # This is a dictionary with key=feature_name, value=feature_function
    features_one_channel = {
            "mean_power":mean_power
            }
    features_two_channel = {
            "coherence":cohere_two_windows
            }
    

    ### END  , TEMPORARY HACK


    if session_basenames_list == "all": 
        session_basenames_list = get_all_valid_session_basenames(RAW_DATA_PATH)
    elif type(session_basenames_list)==type([]):
        assert check_session_basenames_are_valid(session_basenames_list,RAW_DATA_PATH)
    else: raise Exception()
        

    # Init Pandas DataFrame
    df_cols = {ft:[] for ft in FEATS}
    df_cols["session_basename"] = []
    df_cols["bin_name"] = [] # the name of the bin, one of BIN_NAMES
    for ft in feature_list: df_cols[f"{ft}"] = [] # TODO, this line
    df = pd.DataFrame(dataframe_init)

    for session_basename in session_basenames_list:
        # Retrieve seizure times from metadata
        session_metadata_path = os.path.join(RAW_DATA_PATH, session_basename+".txt")
        start_times,end_times = get_seizure_start_end_times(session_metadata_path)
        total_session_time = _get_total_session_time(basename,precision) # in secs
        for start_time,end_time in zip(start_times,end_times):
            # _get_bin_absolute_intervals() returns a dic with 
            # keys=BIN_NAMES,values=(strtbin,endbin) (in s). 
            # If the intervals of a bin are not valid, either because 
            # they start before or end after a file or because they 
            # overlap with another seizure: the value at bin_name 
            # is None.
            bins_absolute_intervals = _get_bin_absolute_intervals(
                    seizure_start_time   = start_time,
                    seizure_end_time     = end_time,
                    preictal_bins        = PREICTAL_BINS, 
                    posictal_delay       = POSTICTAL_DELAY,
                    all_start_times      = start_times,
                    all_end_times        = end_times,
                    total_session_time   = total_session_time
                    )
        
            for bin_name in BIN_NAMES:
                interval = bins_absolute_intervals[bin_name]
                pct      = PCT_DIC[bin_name]
                if interval is not None:
                    start_time,end_time = interval
                    print("dummy print statement")
                    window_starts = _get_x_pct_time_of_interval(
                            start_time      = start_time,
                            end_time        = end_time,
                            window_length   = DUR_FEAT,
                            pct             = pct
                            )
                    ### Pseudocode
                    # For each start time
                    # Init a list of windows
                    # Load each channel into a window, add each window to the list
                    # Compute all features
                    # Add a row to our features data-frame

                    # TODO: compute all features for each pct window in interval
                    # TODO: add a row to the data[bin_name] pandas data-frame
                    #   row includes session_basename, the start times of the windows
                    #   see the jupyter notebook for inspiration on how to add row

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
and the boundaries of the sessions, randomly generate list (or dict for labels?)
of window-start times to sample and compute features for.
"""


if __name__=="__main__":
    ### TEST _get_x_pct_time_of_interval()
    arr = _get_x_pct_time_of_interval(
            start_time    = 5.0,
            end_time      = 152.6,
            window_length = 1.0,
            pct           = 0.05)
    print("Test _get_x_pct_time_of_interval(), array returned:")
    print(arr)
    print()


    ### TEST _get_bin_absolute_intervals()  
    # Test 1
    bin_abs_intervals = _get_bin_absolute_intervals(
            start_seiz           = 100,
            end_seiz             = 130,
            preictal_bins        = [50 , 25 , 10],
            postictal_delay      = 60,
            bin_names            = ["pre1","pre2","pre3","intra","post"],
            all_start_times      = [100,500,1000],
            all_end_times        = [130,550,1100],
            total_session_time = 2000
            )
    assert bin_abs_intervals["pre1"] == (50,75)
    assert bin_abs_intervals["pre2"] == (75,90)
    assert bin_abs_intervals["pre3"] == (90,100)
    assert bin_abs_intervals["intra"] == (100,130)
    assert bin_abs_intervals["post"] == (130,190)

    # Test 2
    bin_abs_intervals = _get_bin_absolute_intervals(
            start_seiz           = 100,
            end_seiz             = 130,
            preictal_bins        = [50 , 25 , 10],
            postictal_delay      = 10,
            bin_names            = ["pre1","pre2","pre3","intra","post"],
            all_start_times      = [50,100,135],
            all_end_times        = [60,130,170],
            total_session_time = 2000 
            )
    assert bin_abs_intervals["pre1"] == None    # Overlaps with previous post-ictal
    assert bin_abs_intervals["pre2"] == (75,90)
    assert bin_abs_intervals["pre3"] == (90,100)
    assert bin_abs_intervals["intra"] == (100,130)
    assert bin_abs_intervals["post"] == None        # Overlaps with next seizure

    # Test 3
    bin_abs_intervals = _get_bin_absolute_intervals(
            start_seiz           = 15,
            end_seiz             = 100,
            preictal_bins        = [50 , 25 , 10],
            postictal_delay      = 60,
            bin_names            = ["pre1","pre2","pre3","intra","post"],
            all_start_times      = [15],
            all_end_times        = [100],
            total_session_time = 150 
            )
    assert bin_abs_intervals["pre1"] == None        # Before file start
    assert bin_abs_intervals["pre2"] == None        # Before file start
    assert bin_abs_intervals["pre3"] == (5,15)      # Valid
    assert bin_abs_intervals["intra"] == (15,100)   # Valid
    assert bin_abs_intervals["post"] == None        # Ends after end of file

    # Not every single edge-case is tested... (low priority TODO)
    print("Tests All Passed: _get_bin_absolute_intervals()")








