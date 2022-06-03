# Copyright (C) 2004-2011 by Michaël Zugaro
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# 22/04/2022 Modified and translated to Python (3.9.5) 
#            by Stephen Fay. Contact: dcxstephen@gmail.com 
# 03/20/2014 Modified by John D. Long to use only built-in Matlab 8.1
#            functions. Contact: jlong29@gmail.com

import os                        # I/O
import numpy as np               # Scientific computing
import logging                   # Debug
from tqdm import tqdm            # Progress Bar
from contextlib import ExitStack # Context manager for opening many files at once

# # Init logger and set the logging level
# logging.basicConfig(level=logging.DEBUG)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG) # DEBUG < INFO < WARNING < ERROR < CRITICAL

# Constant, used in _load_binary (=> it's parent load_binary too) and merge_dats
MAX_SAMPLES_PER_CHUNK = 10000 

# Helper, used in other modules too (don't repeat yourself principle) 
def get_n_samples_from_dur_fs(dur,fs):
    return int(dur * fs + 0.5) 

def load_binary_multiple_segments(
        file_path : str,
        n_chan          : int = 1,
        sample_rate     : int = None,
        offset_times    : list = [], 
        duration_time   : float or None = None,
        offset_sizes    : list = [],
        duration_size   : int or None = None,
        channels        : list = [],
        precision       : type = "int16"
        ) -> np.ndarray:
    """Load many segments of data from multiplexed binary file.

    Either provide a list of offset times and a duration time in seconds
    XOR provide a list of offset sizes and a duration size for the window
    in number of samples. 

    Parameters
    ----------
    file_path : str
        Path to a .dat binary file
    n_chan : int
        Number of data channels in the file (defaults to 1)
    sample_rate : int or float
        Sample rate in Hz, (aka fs, frequency, sr is the MNE convention)
        Defaults to None, if none, must specify offset_size and duration_size
    offset_times : list or np.ndarray
        Positions to start reading in seconds, (aka start_time), (defaults to empty)
    duration_time : float or None = None 
        Duration to read in seconds (per channel) (defaults to None)
    offset_sizes : list or np.ndarray
        Positions to start reading in num of samples, defaults to empty.
    duration_size : int or None
        Duration to read in number of samples (per channel) (defaults to None)
    channels : list 
        Indices of channels to read from, defaults to empty and uses all chs.
    precision : str
        Sample precision, defaults to 'int16'.

    Returns
    -------
    numpy.ndarray
        A 3d array containg the segments' data.
        Shape = (n_segments , n_samples , n_binary_channels)
    """
    # If required, convert time to n_samples (aka sizes) 
    if list(offset_times): # falsy
        assert duration_time is not None, "Duration time must be specified"
        assert duration_time > 0 , "Duration time must be specified"
        assert not offset_sizes, "Cannot specify both times and sizes" 
        assert not duration_size, "Cannot specify both times and sizes"
        offset_sizes = [get_n_samples_from_dur_fs(dt,sample_rate) for dt in offset_times]
        duration_size = get_n_samples_from_dur_fs(duration_time,sample_rate)
    assert list(offset_sizes)
    assert duration_size > 0
    if not channels: channels = [i for i in range(n_chan)]
    # TODO: check whether they are integer values? Prob not necessary
    # TODO: check channels are valid ints in valid range

    n_segments = len(offset_sizes) # the number of segments aka windows
    # Allocate space in memory
    segments_data = np.zeros((n_segments, duration_size, len(channels)),dtype=precision) 
    for idx,offset_size in enumerate(offset_sizes):
        segments_data[idx,:,:] = load_binary(
                file_path,
                n_chan,
                sample_rate,
                offset_size=offset_size,
                duration_size=duration_size,
                channels=channels,
                precision=precision)
    return segments_data

def load_binary(
        file_path : str,
        n_chan : int = 1,
        sample_rate : int = None,
        offset_time : float = None,
        duration_time : float = None,
        offset_size : int = None,
        duration_size : int = None,
        channels : list = [],
        precision : type = "int16") -> np.ndarray:
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
    channels : list or None
        Indices of channels to read from, defaults to None, if None uses all chs. 
    precision : str
        Sample precision, defaults to 'int16'.

    Returns
    -------
    numpy.ndarray
        A 2d array containg the specified segment's data. (1d if only one chan)
    """
    # Checks to make sure the intput is correct
    assert n_chan == int(n_chan)
    assert n_chan >= 1
    logger.debug(f"{n_chan} channel(s) in this binary file") 
    assert os.path.exists(file_path) , f"{file_path} appears not to exist."
    if sample_rate is not None: assert sample_rate > 0 , f"Sample rate must be positive {sample_rate}"
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
            offset_size = get_n_samples_from_dur_fs(offset_time,sample_rate)
        if duration_time: 
            duration_size = get_n_samples_from_dur_fs(duration_time,sample_rate)
    else:
        raise Exception("Invalid Argument Combination!\nYou cannot specify both size-like and a time-like arguments for the duration and offset.")
    assert offset_size >= 0 and int(offset_size) == offset_size , f"Bad offset {offset_size}"
    assert duration_size > 0 , f"Non-positive duration size {duration_size}"

        
    # Figure out what the data offset is in bytes
    bytes_per_sample = np.dtype(precision).itemsize
    fsize_bytes = os.path.getsize(file_path)        # file size in num of bytes
    fsize_samples = fsize_bytes // bytes_per_sample # file size in num of samples
    assert fsize_bytes / bytes_per_sample == fsize_samples
    fsize_samples_tail = fsize_samples - offset_size

    # Make sure duration_size is compatible with file size and offset
    if duration_size == np.inf:
        logger.info("duration_size is np.inf")
        duration_size = fsize_samples_tail // n_chan
        assert fsize_samples_tail / n_chan == duration_size , f"Incompatability of parameters with shape of file. Either n_chan={nchan} is incorrect or your file {file_path} is corrupted."
    else: 
        assert duration_size * n_chan <= fsize_samples_tail , f"Duration size ={duration_size} and offset={offset_size} exceed the end of the file {file_name}"


    data_offset = offset_size * n_chan * bytes_per_sample
    n_samples = duration_size # number of samples per channel

    return _load_binary(file_path,n_chan,n_samples,precision,data_offset)[:,channels]


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
        The number of units (samples/measurements) per channel 
    precision : type (or a str representation of a valid type)
        The precision of the binary data, 
        e.g. numpy.int16 or "int16" are both valid
    data_offset : int
        Exact index of starting time.

    Returns
    -------
    np.ndarray
        the loaded segment of size (n_samples , n_chan)
    """
    # TODO: email John D. Long jlong29@gmail.com or Michaël Zugaro 
    # about this they are the authors of the matlab script upon which 
    # this script is based I don't understand memory allocation stuff 
    # well enough to understand why this max_samples_per_chunk monkey 
    # business is required
    total_n_samples = n_samples * n_chan 
    with open(file_path , "rb") as file:
        # Rem.  data_offset: uint = 
        #           start_time * sample_rate * n_chan * bytes_per_sample
        # Rem.  bytes_per_sample = np.dtype(precision).itemsize
        file.seek(data_offset)
        if total_n_samples <= MAX_SAMPLES_PER_CHUNK:
            data = _load_chunk(file,n_chan,n_samples,precision)
        else:
            # Preallocate memory
            data = np.zeros((n_samples , n_chan) , dtype=precision)

            # Read all chunks
            n_samples_per_chunk = MAX_SAMPLES_PER_CHUNK // n_chan * n_chan
            n_chunks = n_samples // n_samples_per_chunk 
            if not n_chunks: m=0 # extreme rare case, required define m for assertion
            for j in range(n_chunks):
                d =  _load_chunk(file,n_chan,n_samples_per_chunk,precision)
                m,_ = d.shape
                data[j*m:(j+1)*m , :] = d
            # If data size not multiple of chunk size, read remainder
            remainder = n_samples - n_chunks * n_samples_per_chunk
            if remainder:
                d = _load_chunk(file,n_chan,remainder,precision)
                m_rem,_ = d.shape
                assert m_rem # sanity check: logically m_rem cannot be zero
                assert n_chunks*m == data.shape[0] - m_rem # sanity check
                data[-m_rem: , :] = d
    return data


def merge_dats(
        fpaths_in: list,
        dir_out: str,
        fname_out: str,
        precision: str = "int16"
        ):
    """Merges all binary files fnames from the directory dir_in. 

    Returns nothing (void). 

    Parameters
    ----------
    fpaths_in : list
        The ordered list of binary file paths (names) we are merging. 
    dir_out : str
        The directory we want to save the output to. 
    fname_out : str
        The name of the output file we are saving in dir_out
        (including the extension, e.g. '.bin' or '.dat')
    precision : str (optional, defaults to "int16")
        The precision of the data stored in our binary files e.g. "int16"
    """

    assert os.path.exists(dir_out)
    # Assert that all the binary files exist and have equal num of bytes
    # Also, get the size of all the files, in bytes
    size_in_bytes = _assert_all_files_same_size(fpaths_in)
    fpath_out = os.path.join(dir_out,fname_out)
    
    # Define loading parameters
    n_files = len(fpaths_in) # Equal to number of channels in the output file
    n_samples_per_chunk = MAX_SAMPLES_PER_CHUNK // n_files * n_files
    bytes_per_sample = np.dtype(precision).itemsize
    assert size_in_bytes % bytes_per_sample == 0 # Sanity check
    n_samples = size_in_bytes // bytes_per_sample # Number of samples in each file
    # n_chunks = num of full chunks we need to load (there will be a remainder)
    chunk_size = MAX_SAMPLES_PER_CHUNK
    n_chunks = n_samples // chunk_size
    remainder_chunksize = n_samples % chunk_size # In n of samples
    
    logger.info("Started merging files...")
    with ExitStack() as stack, open(fpath_out,"wb") as f_out:
        files = [stack.enter_context(open(fpath,"rb")) for fpath in fpaths_in]

        d_buffer = np.zeros([chunk_size,n_files],dtype=precision) # data buffer, load into f_out
        for _ in tqdm(range(n_chunks)): # tqdm is a progress bar
            # Load a chunk from each of the files we are merging into memory
            for idx,f in enumerate(files):
                d_buffer[:,idx] = np.squeeze(_load_chunk(f,1,chunk_size,precision))
            # Combine the chunks and write them to file
            f_out.write(bytearray(d_buffer.flatten().tobytes())) # TODO: make sure this is saving things at same precision

        # Add the left over chunk
        if remainder_chunksize:
            d_buffer = np.zeros([remainder_chunksize,n_files],dtype=precision)
            for idx,f in enumerate(files):
                d_buffer[:,idx] = np.squeeze(_load_chunk(f,1,remainder_chunksize,precision))
                # Verify that we truely have reached the end of the file
                assert not f.read(1), "Logic Error! Wrongly calculated file size."
            # Combine the chunks and write them to file
            f_out.write(bytearray(d_buffer.flatten().tobytes()))
    logger.info("...Done merging files.")
    return 


# Helper, make sure all files contain the same number of bytes
def _assert_all_files_same_size(filepaths:list):
    if len(filepaths) == 0:
        logger.debug("Zero files provided")
        return None
    os.path.exists(filepaths[0])
    size = os.path.getsize(filepaths[0]) # they should all be same size
    for fpath in filepaths[1:]:
        assert os.path.exists(fpath) # Make sure fname exists
        logger.debug(f"size: {size}")
        logger.debug(f"getsize fpath: {os.path.getsize(fpath)}")
        logger.debug(f"fpath = '{fpath}'")
        assert size == os.path.getsize(fpath)
    logger.debug(f"All {len(filepaths)} files passed are of size={size} bytes")
    return size


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
    assert (n_samples,n_chan) == d.shape , f"Incompatible size (({n_samples},{n_chan}) == {d.shape})"
    return d



if __name__=="__main__":
    # TESTS
    import array
    
    ### Test _load_chunk
    logger.debug("Testing _load_chunk()...")
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
    logger.debug("_load_chunk() passed all tests.")


    ### Test load_binary
    logger.debug("Testing load_binary()...")
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
    assert (arrout == np.arange(4,8).reshape(2,2)).all()
    arrout = load(n_chan=5,sample_rate=1,offset_time=5,duration_time=3)
    assert (arrout == np.arange(25,40).reshape(3,5)).all()
    # Remove temp binary file
    os.remove("temp_test.dat")
    logger.debug("load_binary() passed all tests.")


    ### Test load_binary_multiple_segments
    logger.debug("Testing load_binary_multiple_segments()...")
    # Write a binary file
    arrin = array.array("h", np.arange(500))
    with open("temp_test.dat","wb") as file:
        arrin.tofile(file)
    def load_segs(fname="temp_test.dat" , **kwargs): # helper
        with open(fname,"rb") as file:
            arrout = load_binary_multiple_segments(fname, **kwargs)
        return arrout
    arrout = load_segs(n_chan=2,sample_rate=2,offset_times=[10,20,30],duration_time=5,channels=[0,1],precision="int16")
    assert (arrout.shape == np.array([3,10,2])).all()
    assert (arrout[0,:,:] == np.arange(40,60).reshape(10,2)).all()
    assert (arrout[2,:,:] == np.arange(120,140).reshape(10,2)).all()
    # Same array load, therefore same tests, but with different args
    arrout = load_segs(n_chan=2,sample_rate=2,offset_sizes=[20,40,60],duration_size=10,precision="int16") 
    assert (arrout.shape == np.array([3,10,2])).all()
    assert (arrout[0,:,:] == np.arange(40,60).reshape(10,2)).all()
    assert (arrout[2,:,:] == np.arange(120,140).reshape(10,2)).all()
    logger.debug("load_binary_multiple_segments() all tests Passed")
    
    
    ### Test merge_dats() 
    # by nature, merge_dats() contains some very dense code 
    logger.debug("Testing merge_dats()...")
    # TEST 1
    # Create some binary files
    arr1 = np.asarray([0,2,4,6,8,10,12,14])
    arr2 = np.asarray([1,3,5,7,9,11,13,15])
    arr1.astype('int16').tofile("./arr1.dat")
    arr2.astype('int16').tofile("./arr2.dat")
    # Merge them
    merge_dats(
            fpaths_in=["./arr1.dat","./arr2.dat"],
            dir_out="./",
            fname_out="arrs1and2.dat",
            precision="int16")
    # Examine outfile
    merged = np.fromfile("./arrs1and2.dat",dtype="int16")
    expected_merged = np.arange(16,dtype="int16") 
    print(f"Merged dats: {merged}")
    print(f"Expected: {expected_merged}")
    for i,j in zip(merged,expected_merged): assert i==j
    # TEST 2
    # Create some binary files
    arr1 = np.arange(0,2*MAX_SAMPLES_PER_CHUNK,2) % 32768 # max val int16
    arr2 = np.arange(1,2*MAX_SAMPLES_PER_CHUNK,2) % 32768 # max val int16
    arr1.astype('int16').tofile("./arr1.dat")
    arr2.astype('int16').tofile("./arr2.dat")
    # Merge them
    merge_dats(
            fpaths_in=["./arr1.dat","./arr2.dat"],
            dir_out="./",
            fname_out="arrs1and2.dat",
            precision="int16")
    merged = np.fromfile("./arrs1and2.dat",dtype="int16")
    expected_merged = np.arange(0,2*MAX_SAMPLES_PER_CHUNK) % 32768
    # Examine outfile
    for i,j in zip(merged,expected_merged): assert i==j
    # Clean up
    os.remove("./arr1.dat")
    os.remove("./arr2.dat")
    os.remove("./arrs1and2.dat")


    # TODO: Write test for _load_binary()

    




