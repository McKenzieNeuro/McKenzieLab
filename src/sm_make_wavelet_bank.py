# This file is based off of sm_getPowerPerChannel.m 

import toml
import warnings


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

    return # dummy 

