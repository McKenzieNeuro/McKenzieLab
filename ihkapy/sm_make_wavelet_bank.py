"""
Author statement
----------------
This file is based off of sm_getPowerPerChannel.m written by Samuel McKenzie, 
and awt_freqlist.m and Maureen Clerc, Christian Benar, october 2007 
Translated and adapted to python in May 2022 by Stephen Fay dcxstephen@gmail.com
"""


from ihkapy.fileio.binary_io import merge_dats # local dependency
import toml                     # Parameters/config file Options.toml
import os                       # I/O
import shutil                   # I/O
from tqdm import tqdm           # Progressbar
import logging                  # For debugging and following code
import warnings                 # Bulletproof code
import re                       # Regexp library, to bulletproof code
import pyedflib                 # Read from edf files | this is "terrible in MatLab", look into it, apparently it loads everything into RAM
import numpy as np              # Array manipulation, Scientific computing
from numpy.fft import fft, ifft # Signal processing
from scipy.stats import zscore  # Signal processing


# Init logger and set the logging level
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# For variables containing strings of with absolute path, we explicitly 
# include the word "path" in the variable name. For those with relative 
# or leaf paths, we do not put "path" in the name. 

def _load_fileio_and_data_ops(options_path="./Options.toml"):
    """Load the 'fileio' and 'data' dictionaries from options config file

    Helper for make_wavelet_bank() 
    Parses the .toml config file at `options_path` into a dictionary. 
    Returns the sub-dictionaries at the files "fileio" and "params.data"

    Parameters
    ----------
    options_path : str
        Path to the config file

    Returns
    -------
    dict
        A dictionary containing fileio params, i.e. paths to data files.
    dict
        A dictionary containing data parameters required for data manip.
    """

    warnings.warn("Change this relative path once package is configured properly. \nWe need a more reliable way of accessing the options.toml config file")
    with open(options_path,"r") as f:
        ops_string = f.read()
    ops = toml.loads(ops_string)
    fileio = ops["fileio"]
    data_ops = ops["params"]["data"]
    return fileio,data_ops

# TODO: implement for Lusin and Sombrero wavelets too
# Note, Lusin wasn't implemented in Matlab, and pipeline only used Gabor
# 
# compute_wavelet_gabor corresponds to awt_freqlist in the MatLab code
# 
# %  Maureen Clerc, Christian Benar, october 2007
# %  modified from awt from wavelab toolbox
# 
# % History of changes
# % 1/11/2007: (cgb) psi_array: output in complex 
# % 3/06/2008: (cgb) init of psi_array to size of wt
#   3/05/2022: SF translated awt_freqlist to 

def compute_wavelet_gabor(
        signal: np.ndarray,
        fs: int or float,
        freqs: list or float,
        xi: int = 5 # only needed for Gabor
        ) -> np.ndarray: 
    """Computes one or multiple wavelet transforms of the input signal.

    Follows awt_freqlist.m from the buzzcode repository.

    Parameters
    ----------
    signal : np.ndarray
        The input signal. Only accepts 1D signals. 
    fs : int or float
        The sampling frequency. 
    freqs : list or float
        The frequency or list of frequencies to compute. 
    xi : int
        The number of oscillations parameter, only needed for Gabor wavelet.

    Returns
    -------
    np.ndarray
        A numpy array of dim (len(freqs),len(signal))
    """
    # Make sure all types are correct
    if isinstance(freqs, float) or isinstance(freqs, int): freqs = [freqs]
    freqs = np.asarray(freqs)
    signal = np.asarray(signal)
    assert fs > 0 and (isinstance(fs, float) or isinstance(fs, int))
    # assert wavelet_type.lower() in ("gabor","lusin","sombrero")
    assert signal.ndim == 1, "Must be single dim signal" 
    # TODO: implement multi-dim and remove above assertion
    # Note, not crucial because we don't (yet) use that in pipeline

    (len_sig,) = signal.shape
    sigma2 = 1
    omega = np.concatenate((np.arange(0,len_sig//2+1) , np.arange(-((len_sig+1)//2)+1,0))) * fs / len_sig
    # omega *= fs / len_sig

    # TODO: If you feel like diving in and making this more rigorous, 
    # we know that it can't be more than half the nyquist and 
    # on the lower bound the window must be at least double the 
    # lowest frequency (sometimes advise 5x the lowest)
    tolerance = 0.5
    mincenterfreq = 2*tolerance*np.sqrt(sigma2)*fs*xi / len_sig
    maxcenterfreq = fs*xi/(xi+tolerance/np.sqrt(sigma2)) # Shouldn't this be divided by two because of aliasing? 
    logger.debug(f"fs = {fs}")
    logger.debug(f"freqs = {freqs}")
    logger.debug(f"\n\tLowest freq = {min(freqs)}\n\tHighest freq = {max(freqs)}")
    logger.debug(f"\n\tmincenterfreq = {mincenterfreq}\n\tmaxcenterfreq = {maxcenterfreq}")

    s_arr = xi / freqs
    minscale = xi / maxcenterfreq
    maxscale = xi / mincenterfreq
    # reject frequencies that are outside the given scale
    assert ((s_arr >= minscale) & (s_arr <= maxscale)).all() , "Invalid frequencies"
    # TODO: Turn the above assert into a warning

    n_freqs = len(freqs)
    # np.complex64 is numpy's coarsest complex numpy type
    wt = np.zeros((len_sig,n_freqs),dtype=np.complex64) 
    
    for idx,s in enumerate(s_arr):
        freq = (s * omega - xi)
        psi = np.power(4*np.pi*sigma2,0.25) * np.sqrt(s) * np.exp(-sigma2/2 * freq*freq)
        wt[:,idx] = ifft(fft(signal) * psi)

    return np.squeeze(wt) # turns 2d into 1d IFF single freq 

# Helper
def _contains_all(directory:str,*args) -> bool:
    """Determines whether the directory contains all of the file strings provided"""
    all_files = os.listdir(directory)
    for filestr in args:
        if filestr not in all_files: return False
    return True


# Helper, test to make sure our cache folder is not corrupted
def _assert_all_ext_type_match_regexp( 
        directory: str,
        extension: str,
        regexp_base: str):
    for fname in os.listdir(directory):
        base,ext = os.path.splitext()
        if ext == extension: assert bool(re.search(regexp_base,base))
    logger.debug(f"Test passed: all '{ext}' files in {directory} match the regexp:\n{regexp_base}")
    return 


def make_wavelet_bank(
        edf_fname:str,
        fileio:dict,
        data_ops:dict): 
    """Computes and saves a wavelet decomposition of each channel. 

    Uses user defined options from Options.toml (options_filepath) 
    file to compute the Gabor wavelet decomposition of the raw signals 
    in the provided edf file (edf_fname). 

    - Reads edf raw signal specified by edf_fname (and fileio params)
    - Iterates through each channel, computing wavelet convolutions
        for frequencies in a range specified by data_ops
    - Saves output binaries, one binary file for each hardware channel
        but all the frequencies are saved according to the below order

    Save order convention array flattening convention: 
    - Read 'sn' as 'sample number n'
    - A is for Amplitude (=Power), and P is for Phase
    - K is the index of the last frequency (= num of freqs - 1)
    [raw_s0,freq00_A_s0,freq00_P_s0,freq01_A_s0,freq01_P_s0,...,freqk_A_s0,
    freqK_P_s0,raw_s1,freq00_A_s1,freq00_P_s1,...,freqK_A_s1,freqK_P_s1,...
    ...
    raw_sn,freq00_A_sn,freq00_P_sn,freq01_A_s0,...,freqK_P_sn]

    Note: it is important the above convention is respected because this is
    how the binary_io tools read the files. It's the same convention as the 
    MatLab suit. 
 
    Parameters
    ----------

    edf_fname
        The name of the '.edf' raw data file. We look for all edf files 
        in fileio.RAW_DATA_PATH from Options.toml

    fileio : dict
        The fileio parameters defined in the Options.toml config file.

    data_ops : dict
        Data parameters from the Options.toml config file. 
        The Options.toml config file contains user defined parameters, 
        including fileio params, data and data-processing params, ML
        model hyper-params, and training params. 
    """
    assert len(edf_fname.split("."))==2, f"There can be no dots in the base file name, {edf_fname}"
    basename,ext = edf_fname.split(".") 
    assert ext == "edf", f"Incorrect file format, expected .edf, got {edf_fname}"

    # Unpack File IO constants
    RAW_DATA_PATH         = fileio["RAW_DATA_PATH"]
    WAVELET_BINARIES_PATH = fileio["WAVELET_BINARIES_PATH"]
    # Unpack Data config constants
    FS           = data_ops["FS"]
    NUM_FREQ     = data_ops["NUM_FREQ"]
    LOW_FREQ     = data_ops["LOW_FREQ"]
    HIGH_FREQ    = data_ops["HIGH_FREQ"]
    SPACING      = data_ops["SPACING"]
    ZSCORE_POWER = data_ops["ZSCORE_POWER"] # bool
    SCALE_PHASE  = data_ops["SCALE_PHASE"]
    SCALE_POWER  = data_ops["SCALE_POWER"]
    N_CHAN_RAW   = data_ops["N_CHAN_RAW"]
    TS_FEATURES  = data_ops["TS_FEATURES"]
    # Define features space (based on num, lo, high-freq)
    FREQS = np.power(10, np.linspace(np.log10(LOW_FREQ), np.log10(HIGH_FREQ), NUM_FREQ))

    ### IO Checks
    # Check edf file and wavelet binaries directory exist
    edf_path = os.path.join(RAW_DATA_PATH , edf_fname)
    assert os.path.exists(edf_path), f"Invalid edf file path. Make sure edf file exists and is inside of the {RAW_DATA_PATH} directory."
    assert os.path.exists(WAVELET_BINARIES_PATH), f"Invalid path {WAVELET_BINARIES_PATH}\nCheck your configuration file Options.toml"
    # Check if the cache folder for binaries exists, delete it
    cache_dir_path = os.path.join(WAVELET_BINARIES_PATH, "cache")
    if os.path.exists(cache_dir_path):
        logger.info("Deleting cache")
        shutil.rmtree(cache_dir_path)
    else: logger.debug("No cache, proceeding")

    # Read edf file and loop through each channel one at a time
    for channel in range(N_CHAN_RAW):
        # TODO: simplify and sacrifice some of the checks in this untidy loop
        #       in the name of tidyness and readability and the good lord
        os.mkdir(cache_dir_path) # Create the cache directory
        sig = None
        with pyedflib.EdfReader(edf_path) as f:
            assert N_CHAN_RAW == f.signals_in_file, f"N_CHAN_RAW={N_CHAN_RAW} incompatible with detected number of channels in file ={f.signals_in_file}"
            sig = f.readSignal(channel).astype("int16")
            _rms = np.sqrt(np.power(sig,2).mean())
            logger.debug(f"Sample of quantized sig {sig[:10]}")
            logger.debug(f"Rood Mean Square raw quantized signal = {_rms:.3f}")
        assert sig.shape==(f.getNSamples()[0],) # make sure exists and right shape

        # Save raw channel data as .dat binary
        cached_bin_fname_raw = f"{basename}_ch_{str(channel).zfill(3)}_raw.dat"
        logger.debug(f"cache_dir_path = '{cache_dir_path}'")
        sig.tofile(os.path.join(cache_dir_path, cached_bin_fname_raw)) # TODO: This has yet to be tested

        print(f"Computing {NUM_FREQ} Gabor wavelet convolutions for channel {channel}.")
        # Loop through each frequency, convolve with Gabor wavelet
        for f_idx,freq in tqdm(enumerate(FREQS)):
            # Define cache filepath for this channel & freq
            cached_bin_fname_power = f"{basename}_ch_{str(channel).zfill(3)}_freqidx_{str(f_idx).zfill(2)}_A.dat" # A for amplitude
            cached_bin_fname_phase = f"{basename}_ch_{str(channel).zfill(3)}_freqidx_{str(f_idx).zfill(2)}_P.dat" # P for phase
            # Check if this frequency has already been computed & cached
            if _contains_all(cache_dir_path,cached_bin_fname_power,cached_bin_fname_phase):
                continue 
            
            # Convolve signal with the the wavelet (see awt_freqlist)
            wt = compute_wavelet_gabor(signal=sig,fs=FS,freqs=[freq])

            # Conditionally Zscore, re-scale, and save the power
            if "wavelet_power" in TS_FEATURES:
                wt_power = np.abs(wt) # deep copy
                if ZSCORE_POWER==True:
                    wt_power = zscore(np.abs(wt))
                # add comment intoml, SCALE_POWER shld be smaller if no zscore
                logger.debug(f"wt_power.dtype={wt_power.dtype}")
                wt_power = (wt_power * SCALE_POWER).astype("int16")
                wt_power.tofile(os.path.join(cache_dir_path, cached_bin_fname_power), format="int16")

            # Conditionally re-scale and save the phase
            if "wavelet_phase" in TS_FEATURES:
                wt_phase = (np.arctan(np.real(wt) / np.imag(wt)) * SCALE_PHASE).astype("int16")
                
                wt_phase.tofile(os.path.join(cache_dir_path, cached_bin_fname_phase), format="in16")
        logger.info("Finished computing wavelet transforms for channel={channel}.")

        # Check all the dat files in the cache match with the regex—–
        # this is a test to make sure our cached folder is not corrupted
        # Filenames must not contain any dots!
        # Regexp must match all the Amplicude (=Power) and Phase transforms
        # as well as the raw .dat binary. 
        regex = f"^{basename}_ch_{str(channel).zfill(3)}_(freqidx_\d\d|raw)(_A|_P|)$"
        _assert_all_ext_type_match_regexp( 
                directory = cache_dir_path,
                extension = "dat",
                regexp_base = regex)

        ### Merge all the cached frequency .dat files into a single one.
        # Sort the cached binaries, this will put lower f_idx first and 
        # alternate A before P, exact same order they are created
        # All of the files in the cache have the same channel number (checked above)
        sorted_cache_binaries = [i for i in os.listdir(cache_dir_path) if i[-4:]==".dat"]
        sorted_cache_binaries.sort() # The are already be sorted, this line
                                     # is mainly a crutch for the reader
        scb_paths = [os.path.join(cache_dir_path,i) for i in sorted_cache_binaries]
        merge_dats(
                fpaths_in = scb_paths,
                dir_out = WAVELET_BINARIES_PATH,
                fname_out = f"{basename}_ch_{str(channel).zfill(3)}.dat")

        # Delete the cached single-channel binary files
        shutil.rmtree(cache_dir_path) 
    return


def make_wavelet_bank_all(options_path="Options.toml"):
    # Unpack parameters and user-defined constants
    fileio,data_ops = _load_fileio_and_data_ops(options_path)
    # Bind the consts we use to local vars for readability
    RAW_DATA_PATH = fileio["RAW_DATA_PATH"]

    edf_files = [i for i in os.listdir(RAW_DATA_PATH) if os.path.splitext(i)[1]==".edf"] 
    print(f"Make wavelet bank all, convolving {len(edf_files)} edf files.")
    for edf_fname in tqdm(edf_files):
        logger.info(f"Running make_wavelet_bank on {edf_fname}")
        # Convolve with wavelets and write binary files
        make_wavelet_bank(edf_fname, fileio, data_ops)

    return

if __name__ == "__main__":
    # Only when you run this file directly
    make_wavelet_bank_all(options_path="Options.toml")
    




