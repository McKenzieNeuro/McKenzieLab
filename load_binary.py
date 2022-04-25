import os
import numpy as np


# TODO
# Write Tests
# Move tests into seperate file
def load_binary(
        file_path : str,
        n_chan : int = 1,
        sample_rate : int = None,
        offset_time : float = None,
        duration_time : float = None,
        offset_size : int = None,
        duration_size : int = None,
        channels : list = [],
        precision : type = "int16"):
    """Load data from a multiplexed binary file.

    Reading a subset of data can be done in two different manners: 
    either by specifying start time ("offset_time") and duration ("duration_time") 
    (more intuitive), or by indicating the position ("offset_size") and size of 
    the subset in terms of number of samples per channel ("duration_size") 
    (more accurate). The script will complain if both 'time' and 'size'
    arguments are provided. 

    Parameters
    ----------
    file_path : str
        Path to a .dat binary file
    n_chan : int
        Number of data channels in the file (defaults to 1)
    sample_rate : int or float
        Sample rate in Hz, (aka fs, frequency, sr is the MNE convention) 
        Defaults to None, if none, must specify offset_size and duration_size
    offset_time : int or float or None
        Position to start reading in seconds, (aka start_time) (defaults to None)
    duration_time : int or float or None
        Duration to read in seconds, (defaults to Inf)
    offset_size : int or None
        Position to start reading in samples (per channel) (defaults to None)
    duration_size : int or None
        Duration to read in number of samples (per channel) (defaults to None)
    channels : str or list
        Indices of channels to read from (default = 'all').
    precision : str, optional
        Sample precision (default = 'int16').

    Returns
    -------

    """
    # Make sure the intput is correct
    assert n_chan == int(nchan)
    assert n_chan >= 1
    print(f"{n_chan} channel(s) in this binary file")
    assert os.path.exists(file_path) , f"{file_path} appears not to exist."
    assert sample_rate > 0 , f"Sample rate must be positive {sample_rate}"
    if channels: 
        assert len(channels) <= n_chan , "Too many channels passed"
        assert len(set(channels)) == len(channels) , "Repeating channels"
        for chan in channels: 
            assert chan < n_chan and chan >= 0 , "Channel out of range"
            assert int(chan) == chan , "Wrong type, must be int"
    else: channels = [i for i in range(n_chan)]

    # Either all four args are none -> read whole file xor:
    #     offset_time,duration_time xor offset_size,duration_size
    #     are both None (not just Falsy!)
    if sample_rate == None: assert (offset_time,duration_time)==(None,)*2
    if (offset_time,duration_time,offset_size,duration_size)==(None,)*4:
        offset_size = 0
        duration_size = np.inf
    elif (offset_time,duration_time) == (None,)*2:
        if offset_size == None: offset_size = 0
        if duration_size == None: duration_size = np.inf
    elif (offset_size,duration_size) == (None,)*2:
        assert sample_rate
        offset_size = 0
        duration_size = np.inf
        if offset_time: 
            offset_size = int(offset_time * sample_rate + 0.5)
        if duration_time: 
            duration_size = int(offset_time * sample_rate + 0.5)
    else:
        raise Exception("Invalid Argument Combination!\\
                \nYou cannot specify both size-like and a time-like arguments \\
                for the duration and offset.")
    assert offset_size >= 0 and int(offset) == offset , f"Bad offset {offset_size}"
    assert duration_size > 0 , f"Non-positive duration size {duration_size}"

        

    # Figure out what the data offset is in bytes
    bytes_per_sample = np.dtype(precision).itemsize
    fsize_bytes = os.path.getsize(file_path)        # file size in num of bytes
    fsize_samples = fsize_bytes // bytes_per_sample # file size in num of samples
    assert fsize_bytes / bytes_per_sample == fsize_samples
    fsize_samples_tail = fsize_samples - offset_size

    # Make sure duration_size is compatible with file size and offset
    if duration_size == np.inf:
        duration_size = fsize_samples_tail // n_chan
        assert fsize_samples_tail / n_chan == duration_size , f"Incompatability of parameters with shape of file. Either n_chan={nchan} is incorrect or your file {file_path} is corrupted."
    else: 
        assert duration_size * n_chan <= fsize_samples_tail , f"Duration size ={duration_size} and offset={offset_size} exceed the end of the file {file_name}"


    data_offset = offset_size * n_chan * bytes_per_sample
    n_samples = duration_size * bytes_per_sample
    
    return _load_binary(file_path,n_chan,n_samples,precision,data_offset)


def _load_binary(
        file_path : str,
        n_chan : int,
        n_samples : int,
        precision : type,
        data_offset : int = 0) -> np.ndarray: 
    """Helper for load_binary; this is the method that contains the logic.

    Parameters
    ----------
    file_path : str
        Path to binary file with multiplexed data.
    n_chan : int
        The number of channels. 
    n_samples : int
        The number of units (measurements) in the sample 
    precision : type (or a str representation of a valid type)
        The precision of the binary data, 
        e.g. numpy.int16 or "int16" are both valid
    data_offset : int
        Exact index of starting time.
    """
    # TODO: email John D. Long jlong29@gmail.com or Michaël Zugaro 
    # about this they are the authors of the matlab script upon which 
    # this script is based I don't understand memory allocation stuff 
    # well enough to understand why this max_samples_per_chunk monkey 
    # business is required
    MAX_SAMPLES_PER_CHUNK = 10000 
    n_samples = n_samples_per_chan * n_chan
    with open(file_path , "rb") as file:
        # Rem.  data_offset: uint = 
        #           start_time * sample_rate * n_chan * bytes_per_sample
        # Rem.  bytes_per_sample = np.dtype(precision).itemsize
        file.seek(data_offset)
        if n_samples <= MAX_SAMPLES_PER_CHUNK:
            data = _load_chunk(file,n_chan,n_samples,precision)
        else:
            # Preallocate memory
            data = np.zeros((n_samples , n_chan) , dtype=precision)

            # Read all chunks
            n_samples_per_chunk = int(MAX_SAMPLES_PER_CHUNK / n_chan) * n_chan
            n_chunks = n_samples // n_samples_per_chunk 
            if not n_chunks: m=0 # extreme rare case, required for assertion
            for j in range(n_chunks):
                d =  _load_chunk(file,n_chan,n_samples,precision)
                m,_ = d.shape
                data[j*m:(j+1)*m , :] = d
            # If data size not multiple of chunk size, read remainder
            remainder = n_samples - n_chunks * n_samples_per_chunk
            if remainder:
                d = _load_chunk(file,n_chan,remainder//n_chan,precision)
                m_rem,_ = d.shape[0]
                assert m_remi # sanity check: logically m_rem cannot be zero
                assert n_chunks*m == data.shape[0] - m_rem # sanity check
                data[-m_rem: , :] = d
    return data



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
        that you get from open(file_path , "rb")
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
    d = np.fromfile(
        file,
        dtype = precision,
        count = n_chan * n_samples).reshape((n_samples,n_chan))
    assert (n_samples,n_chan) == d.size , f"Incompatible size ({n_chan},{n_samples} == {d.shape})"
    return d



if __name__=="__main__":
    # TESTS
    import array
    
    ### Test _load_chunk
    print("\nTesting _load_chunk()...",end="\t")
    # Write a binary file
    arrin = array.array("h" , np.arange(50))
    with open("temp_test.dat","wb") as file:
        arrin.tofile(file)

    # Read the binary file we just created with _load_chunk
    with open("temp_test.dat","rb") as file:
        arrout = _load_chunk(file, n_chan=2, n_samples=5, precision="int16")

    # Assert that result must be equal to matlab version
    assert (arrout == np.array([[0,1],[2,3],[4,5],[6,7],[8,9]])).all()

    # Remove temp binary file
    os.remove("temp_test.dat")
    print("Passed")

    ### Test load_binary
    print("\nTesting load_binary()...",end="\t")

    # Write a binary file
    arrin = array.array("h" , np.arange(50))
    with open("temp_test.dat","wb") as file:
        arrin.tofile(file)

    def load(fname="temp_test.dat" , **kwargs): # helper
        with open(fname,"rb") as file:
            arrout = load_binary(fname,**kwargs)
        return arrout
    arrout = load(n_chan=2)
    assert (arrout == np.arange(50).reshape(25,2)).all()
    arrout = load(n_chan=2,offset_size=2,duration_size=2)
    assert (arrout == np.arange(4,8)).reshape(2,2)).all()
    













