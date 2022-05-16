

def calc_features(
        fname,
        start_times,
        
        ) -> np.ndarray:
    # Get list of features that we want to use from a features file
    # each feature should be a function with standard in/out types
    # a segment (window) is the input, the output is a float

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
        fname,
        start_times,
        # bunch of options
        ) -> np.ndarray:
    """

    Returns
    -------
    np.ndarray
        2d array shape (n_windows,n_features)
    """
    return xyz

