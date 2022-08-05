# Script to convert .edf files to the CSV format

"""
Run this file with one command line argument: the path to a edf file, it
will parse the data into an mne.raw object, then save that to a csv
file in the same directory.
"""


import numpy as np
import mne
import sys 

edf_path = sys.argv[1] # csv path passed as argument to the script
head,edf_fname = os.path.split(edf_path)
basename,_ = os.path.splitext(edf_fname)
csv_path = os.path.join(head,f"{basename}.csv")


edf = mne.io.read_raw_edf(edf_path)
header = ','.join(edf.ch_names)
np.savetxt(csv_path, edf.get_data().T, delimiter=',', header=header)

