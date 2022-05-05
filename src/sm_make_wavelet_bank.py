"""
Author statement
----------------
This file is based off of sm_getPowerPerChannel.m written by Samuel McKenzie, 
and awt_freqlist.m and Maureen Clerc, Christian Benar, october 2007 
Translated and adapted to python in May 2022 by Stephen Fay dcxstephen@gmail.com
"""


from fileio.binary_io import merge_dats # local dependency
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
        A dictionary containing paths to data files.
    dict
        A dictionary containing data parameters needed in data
        manipulation 
    """

    warnings.warning("Change this relative path once package is configured properly. \nWe need a more reliable way of accessing the options.toml config file")
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
        try:
            base,ext = fname.split(".")
        except ValueError:
            # enforce that there can be no dots in files with 
            # the 'extension' extension
            assert fname.split(".")[-1] != extension, "ValueError, all {extension} files must have no dots in the basename!" 
        else:
            if ext==extension:
                assert bool(re.search(regexp,base))
    logger.debug(f"Test passed: all '{ext}' files in {directory} match the regexp:\n{regexp_base}")
    return 



def make_wavelet_bank(edf_fname,options_filepath):
    """Computes and saves a wavelet decomposition of each channel. 

    Uses user defined options from Options.toml (options_filepath) 
    file to compute the Gabor wavelet decomposition of the raw signals 
    in the provided edf file (edf_fname). 
 
    Parameters
    ----------

    edf_fname
        The name of the '.edf' raw data file. We look for this in 
        fileio.RAW_DATA_PATH from Options.toml

    options_filepath
        The path to our Options.toml config file. This file contains
        user defined parameters, including fileio params, data and 
        data-processing params, ML model hyper-params, and training 
        params.
    """
    assert len(edf_fname.split("."))==2, f"There can be no dots in the base file name, {edf_fname}"
    basename,ext = edf_fname.split(".") 
    assert ext == "edf", f"Incorrect file format, expected .edf, got {edf_fname}"
    assert options_path[-5:] == ".toml", f"Incorect file format, expected .toml extension, got {options_path}"

    ### UNPACK parameters and user-defined constants
    fileio,data_ops = _load_fileio_and_data_ops(options_path)
    # Unpack File IO constants
    RAW_DATA_PATH = fileio["RAW_DATA_PATH"]
    WAVELET_BINARIES_PATH = fileio["WAVELET_BINARIES_PATH"]
    # Unpack Data config constants
    FS = data_ops["FS"]
    NUM_FREQ = data_ops["NUM_FREQ"]
    LOW_FREQ = data_ops["LOW_FREQ"]
    HIGH_FREQ = data_ops["HIGH_FREQ"]
    SPACING = data_ops["SPACING"]
    ZSCORE_POWER = data_ops["ZSCORE_POWER"] # bool
    SCALE_PHASE = data_ops["SCALE_PHASE"]
    SCALE_POWER = data_ops["SCALE_POWER"]
    TS_FEATURES = data_ops["TS_FEATURES"]

    # Check edf file exists
    edf_path = os.path.join(RAW_DATA_PATH , edf_fname)
    assert os.path.exists(edf_path), f"Invalid edf file path. Make sure edf file exists and is inside of the {RAW_DATA_PATH} directory."

    # Check wavelet binaries path directory exists
    assert os.path.exists(WAVELET_BINARIES_PATH), f"Invalid path {WAVELET_BINARIES_PATH}\nCheck your configuration file Options.toml"

    # Check if the cache folder for binaries exists
    cache_dir_path = os.path.join(WAVELET_BINARIES_PATH, "cache")
    if not os.path.exists(cache_dir_path): 
        logger.info(f"No cache directory, creating one at\n{cache_dir_path}")
        os.mkdir(cache_dir_path)
    else:
        # Check if cached binaries are right, if not delete all and start over
        cached_files = [i for i in os.listdir(cache_dir_path) if i[-4:]==".dat"]
        for i in cached_files:
            # if any of the cached dat files don't match the base name, remove them
            if i[:len(basename)] != basename:
                shutil.rmtree(cache_dir_path)
                os.mkdir(cache_dir_path)
                break

    # Define features space (based on num, lo, high-freq)
    FREQS = np.power(10, np.linspace(np.log10(LOW_FREQ), np.log10(HIGH_FREQ), NUM_FREQ))

    # Read edf file and loop through each channel one at a time
    for channel in range(N_CHAN_RAW):
        sig = None
        with pyedflib.EdfReader(edf_path) as f:
            assert N_CHAN_RAW == f.signals_in_file, f"N_CHAN_RAW={N_CHAN_RAW} incompatible with detected number of channels in file ={f.signals_in_file}"
            sig = f.readSignal(i)
        assert sig.shape==(f.getNSamples()[0],) # make sure exists and right shape

        # Save raw channel data as .dat binary
        sig.tofile(os.path.join(cache_dir_path, cached_bin_fname_power), format="int16") # TODO: This has yet to be tested

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
                # add comment intoml, SCALE_POWER shld be different if no zscore
                wt_power *= SCALE_POWER 
                wt_power.tofile(os.path.join(cache_dir_path, cached_bin_fname_power), format="int16")

            # Conditionally re-scale and save the phase
            if "wavelet_phase" in TS_FEATURES:
                wt_phase = np.arctan(np.real(wt) / np.imag(wt)) * SCALE_PHASE
                wt_phase.tofile(os.path.join(cache_dir_path, cached_bin_fname_phase), format="in16")
        logger.info("Finished computing wavelet transforms for channel={channel}.")

        # Check all the dat files in the cache match with the regex—–
        # this is a test to make sure our cached folder is not corrupted
        # Filenames must not contain any dots!
        regex = f"^{basename}_ch_{str(channel).zfill(3)}_freqidx_\d\d_(A|P)$"
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
        shutil.rmtree(cache_dir_path) # TODO: This has yet to be tested
    return



if __name__ == "__main__":
    # Init logger and set the logging level
    logger = logging.getLogger(__name__)
    console_logger = logging.StreamHandler()
    console_logger.setLevel(logging.DEBUG)
    logger.addHandler(console_logger) # DEBUG < INFO < WARNING < ERROR < CRITICAL

    #[test]
    def test_compute_wavelet_gabor(plot=False):
        signal = np.random.normal(0,1000,np.power(2,16))
        fs = 16.0 # Hz
        freqs = np.power(10,np.linspace(np.log10(0.5),np.log10(8),5))
        # logger.debug(f"Test compute_wavelet_gabor\nfreqs = {freqs}")
        xi = 5
        wt = compute_wavelet_gabor(signal,fs,freqs,xi)
        logger.debug(f"\n\twt.shape={wt.shape}\n\tInput signal shape={signal.shape}\n\tLen freqs={len(freqs)}")
        assert len(signal),len(freqs) == wt.shape
        if plot==True:
            logging.getLogger('matplotlib').setLevel(logging.WARNING)
            import matplotlib.pyplot as plt
            # plt.figure(figsize=(12,8))
            plt.subplots(2,1,figsize=(12,7))
            plt.suptitle(f"Wavelet Transforms of Gaussian Random Noise (sigma=1000)\nSample Frequency = {fs}")

            plt.subplot(2,1,1)
            plt.title("REAL")
            labels = [f"freq = {i:.2f}" for i in freqs]
            plt.plot(np.real(wt[100:200,:]),".",alpha=0.5,label=labels)
            plt.plot(np.real(wt[100:200,:]),"--",color="k",linewidth=0.5,alpha=0.1)
            plt.legend()

            plt.subplot(2,1,2)
            plt.title("IMAGINARY")
            labels = [f"freq = {i:.2f}" for i in freqs]
            plt.plot(np.imag(wt[100:200,:]),".",alpha=0.5,label=labels)
            plt.plot(np.imag(wt[100:200,:]),"--",color="k",linewidth=0.5,alpha=0.1)
            plt.legend()

            plt.show(block=True)
    test_compute_wavelet_gabor(plot=True) # set plot=True to display plots in test 
    # logger.info("TEST PASSED: compute_wavelet_gabor()") # strange logging bug, see https://stackoverflow.com/questions/72127312/python-logger-setlevel-bug 
    print("TEST PASSED: compute_wavelet_gabor()")



    def test_make_wavelet_bank():
        # TODO: implement this test
        return
    test_make_wavelet_bank()
    # uncomment below once implemented
    # logger.info("TEST PASSED: make_wavelet_bank()") 





