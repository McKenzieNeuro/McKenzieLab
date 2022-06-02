import os
import warnings

# fileio utility methods

def get_all_valid_session_basenames(dir_path):
    """Returns list of all base-filenames with edf and txt pair."""
    valid_basenames = []
    for fname in os.listdir(dir_path):
        basename,ext = os.path.splitext(fname)
        if ext == ".edf" and basename+".txt" in os.listdir(dir_path):
            valid_basenames.append(basename)
    return valid_basenames

def check_session_basenames_are_valid(basenames,dir_path):
    """Raises error if edf and txt files not found in dir_path"""
    if not basenames: print("Warning: no basenames supplied, True by default")
    for basename in basenames:
        if basename+".txt" not in dir_path or basename+".edf" not in dir_path:
            return False
    return True
