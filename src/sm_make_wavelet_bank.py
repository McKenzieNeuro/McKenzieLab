# This file is based off of sm_getPowerPerChannel.m 

import toml
import warnings
import os
import logging


def _load_fileio_and_data_ops(options_path="../Options.toml"):
    warnings.warning("Change this relative path once package is configured properly. We need a more reliable way of accessing the options.toml config file")
    with open(options_path,"r") as f:
        ops_string = f.read()
    ops = toml.loads(ops_string)
    fileio = ops["fileio"]
    data_ops = ops["params"]["data"]
    return fileio,data_ops

# def get_wavelet_bank(
#         path_to_edf_file,
#         n_chan,
#         now,
#         list_,
#         all,
#         the,
#         options,
#         its,
#         unfortunately,
#         the,
#         least,
#         of_evils): # get power per channel
#     return

def wavelet_convolve(signal, fs, freq,...):
    return # dummy

def make_wavelet_bank(edf_fname,options_filepath):
    """Computes and saves a wavelet decomposition of each channel. 
    
    Parameters
    ----------

    Returns nothing (look up how to express this in numpy doc syntax)

    """
    assert len(edf_fname.split("."))==2, f"There can be no dots in the base file name, {edf_fname}"
    basename,ext = edf_fname.split(".") 
    assert edf_fname[-4:] == ".edf", f"Incorrect file format, expected .edf, got {edf_fname}"
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
    ZSCORE_POWER = data_ops["ZSCORE_POWER"]
    SCALE_PHASE = data_ops["SCALE_PHASE"]
    SCALE_POWER = data_ops["SCALE_POWER"]

    # Check edf file exists
    edf_path = os.path.join(RAW_DATA_PATH , edf_fname)
    assert os.path.exists(edf_path), f"Invalid edf file path. Make sure edf file exists and is inside of the {RAW_DATA_PATH} directory."
    # Check wavelet binaries path directory exists
    assert os.path.exists(WAVELET_BINARIES_PATH), f"Invalid path {WAVELET_BINARIES_PATH}\nCheck your configuration file Options.toml"
    cache_dir = os.path.join(WAVELET_BINARIES_PATH, "cache")
    # Check if the cache folder for binaries exists
    if not os.path.exists(cache_dir): 
        logging.info(f"No cache directory, creating one at\n{cache_dir}")
        os.mkdir(cache_dir)
    else:
        cached_files = 

    # TODO: Check if the binaries are already there
    # TODO: Check, and conditionally clear cache dire (for temp storage of convs)
    # TODO: Define features space (based on num, lo, high-freq)
    # TODO: Read edf file and loop through each channel one at a time
        # TODO: Loop through each frequency
            # TODO: Define cache filepath for this channel & freq
            # TODO: Check if this frequency has already been computed & cached
            # TODO: Convolve signal with the the wavelet (see awt_freqlist)
            # TODO: Conditionally Zscore, re-scale, and save the power
            # TODO: Conditionally re-scale and save the phase
        # TODO: Merge all the cached frequency .dat files into a single one.
        # TODO: Delete cache
    # TODO: Check that each file and frequency is correct

    return # dummy 

