import os

# fileio utility methods

def get_all_valid_session_basenames(dir_path):
    """Returns list of all base-filenames with edf and txt pair."""
    valid_basenames = []
    for fname in os.listdir(dir_path):
        basename,ext = os.path.splitext(fname)
        if ext == ".edf" and basename+".txt" in os.listdir(dir_path):
            valid_basenames.append(basename)
    return valid_basenames

# # Dead code, or potentially useful?
# def check_session_basenames_are_valid(basenames,dir_path):
#     """Raises error if edf and txt files not found in dir_path"""
#     if not basenames: print("Warning: no basenames supplied, True by default")
#     for basename in basenames:
#         if basename+".txt" not in os.listdir(dir_path) \
#                 or basename+".edf" not in os.listdir(dir_path):
#             return False
#     return True

### FORMAT PATHS
# All the properly formatted naming conventions
def fmt_binary_chan_raw_path(
        wavelet_binaries_path:str,
        session_basename:str,
        raw_ch_idx:int):
    """Format multiplexed wavelet bank binary file."""
    return os.path.join(wavelet_binaries_path,f"{session_basename}_ch_{str(raw_ch_idx).zfill(3)}.dat")

def fmt_binary_cache_wav_path(
        cache_path:str, # path to directory of cache
        session_basename:str,
        raw_ch_idx:int,
        amp_phase:str,
        freq_ch_idx:int=None):
    """
    Format path to single channel wavelet binary temporarily cached
    while computing transform bank.
    """
    assert amp_phase in ("RAW","AMP-PHASE"), "The Amplitude / Phase must be denoted AMP or PHASE"
    if amp_phase == "RAW":
        fname_raw = f"{session_basename}_ch_{str(raw_ch_idx).zfill(3)}_0RAW.dat"
        return os.path.join(cache_path,fname_raw)
        # Above, the zero infront of RAW is intentional, it's for sorting
        # "0..." < "freqidx..."
    else:
        fname_amp = f"{session_basename}_ch_{str(raw_ch_idx).zfill(3)}_freqidx_{str(freq_ch_idx).zfill(2)}_AMP.dat"
        fname_phase = f"{session_basename}_ch_{str(raw_ch_idx).zfill(3)}_freqidx_{str(freq_ch_idx).zfill(2)}_PHASE.dat"
        return os.path.join(cache_path,fname_amp),os.path.join(cache_path,fname_phase) 

def fmt_features_df_csv_path(
        features_path:str,
        basename:str):
    """Filepath for csv of features DataFrame of a session."""
    return os.path.jion(features_path,f"{basename}.csv")



# TESTS
if __name__ == "__main__":
    # Test fmt_binary_chan_raw_path
    wav_path = "/path/to/wavelet/"
    sb = "sb"
    rch_idx = 5
    assert fmt_binary_chan_raw_path(wav_path,sb,rch_idx) == "/path/to/wavelet/sb_ch_005.dat"
    print("Test passed: fmt_binary_chan_raw_path()")

    # Test fmt_binary_cache_wav_path()
    cache_path = "path/to/cache"
    assert fmt_binary_cache_wav_path(cache_path,sb,4,6,"P") == f"path/to/cache/sb_ch_004_freqidx_06_P.dat"
    print("Test passed: fmt_binary_cache_wav_path()")

    # Test 




