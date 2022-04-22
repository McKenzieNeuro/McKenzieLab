import os
import numpy as np


def load_binary(
        file_name: str,
        sample_rate: int,
        start: int,
        duration: int,):
    """Load data from a multiplexed binary file.

    Reading a subset of data can be done in two different manners: 
    either by specifying start time and duration (more intuitive), or
    by indicating the position and size of the subset in terms of 
    number of samples per channel (more accurate)

    Parameters
    ----------
    file_name : str
        Path to a .dat binary file
    sample_rate : int or float
        Sample rate in Hz, (aka fs, frequency, sr is the MNE convention) 
    start : int or float
        Position to start reading in seconds
    duration : int or float
        Duration to read in seconds, (can be Inf)
    offset : int 
        Position to start reading (in samples per channel)
    samples : int
        Number of samples (per channel) to read, (can be Inf)
    n_channels : int
        Number of data channels in the file
    channels_to_read : str or list
        Indices of channels to read from (default = 'all').
    precision : str, optional
        Sample precision (default = 'int16').
    skip : int
        Number of bytes to skip after each value is read (default = 0).
        

    """
    # Tests
    assert filename[-4:] == ".dat" , "load_binary can only be called on binary files with the .dat extension."
    return # dummy


def _load_binary() -> np.ndarray: 
    """Helper for load_binary; this is the method that does all the work.

    Parameters
    ----------
    file_name : str
        ...
    n_samples : int
        Number of samples per channel
    n_chan : int
        Number of channels. 
    """
    # TODO: email John D. Long jlong29@gmail.com or Michaël Zugaro about this
    # they are the authors of the matlab script upon which this script is based
    # I don't understand memory allocation stuff well enough to understand
    # why this max_samples_per_chunk monkey business is required
    MAX_SAMPLES_PER_CHUNK = 10000 
    n_samples = n_samples_per_chan * n_chan
    with open(file_name , "rb") as file:
        # Rem.  data_offset: uint = 
        #           start_time * sample_rate * n_chan * bytes_per_sample
        # Rem.  bytes_per_sample = np.dtype(precision).itemsize
        file.seek(data_offset)
        if n_samples <= MAX_SAMPLES_PER_CHUNK:
            data = _load_chunk(file,n_chan,n_samples,precision)
        else:
            def ciel(x): return int(x + 1) 
            # Preallocate memory
            data = np.zeros((n_chan , n_samples))


    return arr



# %LoadBinary - Load data from a multiplexed binary file.
# %
# %  Reading a subset of the data can be done in two different manners: either
# %  by specifying start time and duration (more intuitive), or by indicating
# %  the position and size of the subset in terms of number of samples per
# %  channel (more accurate).
# %
# %  USAGE
# %
# %    data = LoadBinary(filename,<options>)
# %
# %    filename       file to read
# %    <options>      optional list of property-value pairs (see table below)
# %
# %    =========================================================================
# %     Properties    Values
# %    -------------------------------------------------------------------------
# %     'frequency'   sampling rate (in Hz, default = 20kHz)
# %     'start'       position to start reading (in s, default = 0)
# %     'duration'    duration to read (in s, default = Inf)
# %     'offset'      position to start reading (in samples per channel,
# %                   default = 0)
# %     'samples'     number of samples (per channel) to read (default = Inf)
# %     'nChannels'   number of data channels in the file (default = 1)
# %     'channels'    channels to read (default = all) (base 1)
# %     'precision'   sample precision (default = 'int16')
# %     'skip'        number of bytes to skip after each value is read
# %                   (default = 0)
# %    =========================================================================
# 
# % Copyright (C) 2004-2011 by Michaël Zugaro
# %
# % This program is free software; you can redistribute it and/or modify
# % it under the terms of the GNU General Public License as published by
# % the Free Software Foundation; either version 3 of the License, or
# % (at your option) any later version.
# 
# % 03/20/2014 Modified by John D. Long to use only built-in Matlab 8.1
# % functions. Contact: jlong29@gmail.com



def _load_chunk(
        file,
        n_chan : int,
        n_samples : int,
        precision : type) -> np.ndarray:
    """Helper. Loads a chunk of size (n_chan,n_samples) from buffered reader f.
    
    This is really just a MatLab-esque wrapper for numpy.fromfile().

    Parameters
    ----------
    file : an io file buffered reader object
        The binary file that you are reading. The python built-in type
        that you get from open(file_name , "rb")
    n_chan : int
        The number of channels. 
    n_samples : int
        The number of units (measurements) in the sample 
    precision : type (or a str representation of a valid type)
        The precision of the binary data, 
        e.g. numpy.int16 or "int16" are both valid
    
    Returns
    -------
    numpy.ndarray
        2D array of shape (n_chan , n_samples)
        If the binary file contains an ascending sequence [0,1,2,3,...]
        then calling _load_chunk with n_chan = 2 and n_samples = 3 will
        result in the following array [[0, 1], [2, 3], [4, 5]]
    """
    arr = np.fromfile(
        file,
        dtype = precision,
        count = n_chan * n_samples).reshape((n_samples,n_chan))
    return arr



if __name__=="__main__":
    # TESTS
    import array
    
    # Test _load_chunk
    print("\nTesting _load_chunk()...",end="\t")
    # Write a binary file
    arrin = array.array("h" , np.arange(50))
    with open(".testdat.dat","wb") as file:
        arrin.tofile(file)

    # Read the binary file we just created with _load_chunk
    with open(".testdat.dat","rb") as file:
        arrout = _load_chunk(file, n_chan=2, n_samples=5, precision="int16")

    # Assert that result must be equal to matlab version
    assert (arrout == np.array([[0,1],[2,3],[4,5],[6,7],[8,9]])).all()
    print("Passed")

















