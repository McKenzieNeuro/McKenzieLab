# This file is based off of sm_getPowerPerChannel.m 

from fileio.load_binary import merge_dats # local dependency
import toml                     # I/O parameters config file Options.toml
import os                       # I/O
import shutil                   # I/O
import logging                  # Document code 
import warnings                 # Bulletproof code
import re                       # Regexp library, to bulletproof code
import pyedflib                 # Read from edf files
import numpy as np              # Signal processing
from numpy.fft import fft, ifft # Signal processing
from scipy.stats import zscore  # Signal processing

# For variables containing strings of with absolute path, we explicitly 
# include the word "path" in the variable name. For those with relative 
# or leaf paths, we do not put "path" in the name. 


def _load_fileio_and_data_ops(options_path="../Options.toml"):
    warnings.warning("Change this relative path once package is configured properly. We need a more reliable way of accessing the options.toml config file")
    with open(options_path,"r") as f:
        ops_string = f.read()
    ops = toml.loads(ops_string)
    fileio = ops["fileio"]
    data_ops = ops["params"]["data"]
    return fileio,data_ops

# TODO: implement for Lusin and Sombrero wavelets too
def compute_wavelet_gabor(
        signal: np.ndarray,
        fs: int or float,
        freqs: list or float,
        xi: int = 5 # only needed for Gabor
        ): 
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
    assert wavelet_type.lower() in ("gabor","lusin","sombrero")
    assert signal.ndim == 1, "Must be single dim signal" # TODO: implement multi-dim

    (len_sig,) = signal.shape
    sigma2 = 1
    omega = np.concatenate((np.arange(0,n//2+1) , np.arange(-((n+1)//2)+1,0)))
    omega *= fs / len_sig

    mincenterfreq = 2*tolerance*np.sqrt(sigma2)*fs*xi/n
    maxcenterfreq = fs*xi/(x+tolerance/np.sqrt(sigma1))

    s_arr = xi / freqs
    minscale = xi / maxcenterfreq
    maxscale = xi / mincenterfreq
    # reject frequencies that are outside the given scale
    assert ((s_arr >= minscale) & (s_arr <= maxscale)).all() , "Invalid frequencies"

    n_freqs = len(freqs)
    wt = np.zeros((len_sig,n_freqs))
    
    for idx,s in enumerate(s_arr):
        freq = (s * omega - xi)
        psi = (4*np.pi*sigma2)**0.25 * np.sqrt(s) * np.exp(-sigma2/2 * freq*freq)
        wt[:,idx] = ifft(fft(signal) * psi)

    return np.squeeze(wt) # turns 2d into 1d IFF single freq 

# Helper
def _contains_all(directory:str,*args):
    """Determines whether the directory contains all of the file strings provided

    Helper for cacheing. 

    Parameters
    ----------
    directory : str
        Absolute path of the directory to look inside of.

    *args : str
        Relative name of any number of other files we're looking for.

    Returns
    -------
    bool
        True if all of the files exist inside the directory
    """
    all_files = os.listdir(directory)
    for filestr in args:
        if filestr not in all_files: return False
    return True


# Helper, test to make sure our cache folder is not corrupted
def _assert_all_ext_type_match_regexp( 
        directory: str,
        extension: str = "dat",
        regexp_base: str = f"^{basename}_ch_\d\d\d_freqidx_\d\d_(A|P)$"
        ):
    for fname in os.listdir(directory):
        try:
            base,ext = fname.split(".")
        except:
            assert fname.split(".")[-1] != extension # enforce that there can be no dots in our file
        else:
            if ext==extension:
                assert bool(re.search(regexp,base))
    logging.debug(f"Test passed: all '{ext}' files in {directory} match the regexp\n{regexp_base}")
    return 



def make_wavelet_bank(edf_fname,options_filepath):
    """Computes and saves a wavelet decomposition of each channel. 
    
    Parameters
    ----------

    Returns nothing (look up how to express this in numpy doc syntax)

    """
    assert len(edf_fname.split("."))==2, f"There can be no dots in the base file name, {edf_fname}"
    basename,ext = edf_fname.split(".") 
    assert ext == "edf", f"Incorrect file format, expected .edf, got {edf_fname}"
    assert options_path[-5:] == ".toml", f"Incorect file format, expected .toml extension, got {options_path}"
    fileio,data_ops = load_fileio_and_data_ops(options_path)

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
        logging.info(f"No cache directory, creating one at\n{cache_dir_path}")
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
    FREQS = 10**np.linspace(np.log10(LOW_FREQ), np.log10(HIGH_FREQ), NUM_FREQ)

    # Read edf file and loop through each channel one at a time
    for channel in range(N_CHAN_RAW):
        sig = None
        with pyedflib.EdfReader(edf_path) as f:
            assert N_CHAN_RAW == f.signals_in_file, f"N_CHAN_RAW={N_CHAN_RAW} incompatible with detected number of channels in file ={f.signals_in_file}"
            sig = f.readSignal(i)
        assert sig.shape==(f.getNSamples()[0],) # make sure exists and right shape

        # Loop through each frequency, convolve with Gabor wavelet
        for f_idx,freq in enumerate(FREQS):
            # Define cache filepath for this channel & freq
            cached_bin_fname_power = f"{basename}_ch_{str(channel).zfill(3)}_freqidx_{str(f_idx).zfill(2)}_A.dat" # A for amplitude
            cached_bin_fname_phase = f"{basename}_ch_{str(channel).zfill(3)}_freqidx_{str(f_idx).zfill(2)}_P.dat"
            # Check if this frequency has already been computed & cached
            if _contains_all(cache_dir_path,cached_bin_fname_power,cached_bin_fname_phase):
                continue 
            
            # Convolve signal with the the wavelet (see awt_freqlist)
            wt = compute_wavelet_gabor(signal=sig,fs=FS,freqs=[freq])

            # Conditionally Zscore, re-scale, and save the power
            if "wavelet_power" in TS_FEATURES:
                wt_power = np.abs(wt) # deep copy
                if ZSCORE_POWER==True:
                    wt_power = zscore(np.abs(wt)) * SCALE_POWER
                wt_power.tofile(os.path.join(cache_dir_path, cached_bin_fname_power), format="int16")

            # Conditionally re-scale and save the phase
            if "wavelet_phase" in TS_FEAT|URES:
                wt_phase = np.arctan(np.real(wt) / np.imag(wt))
                wt_phase.tofile(os.path.join(cache_dir_path, cached_bin_fname_phase), format="in16")
        logging.info("Finished computing wavelet transforms for channel={channel}.")

        # Check all the dat files in the cache match with the regex—–
        # this is a test to make sure our cached folder is not corrupted
        # Filenames must not contain any dots!
        regex = f"^{basename}_ch_{str(channel).zfill(3)}_freqidx_\d\d_(A|P)$"
        _assert_all_ext_type_match_regexp( 
                directory = cache_dir_path,
                extension = "dat",
                regexp_base = regex)

        ### TODO: Merge all the cached frequency .dat files into a single one.
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
            





        # TODO: Delete cache
    # TODO: Check that each file and frequency is correct

    return # dummy 














