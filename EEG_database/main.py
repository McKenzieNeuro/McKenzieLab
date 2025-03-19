import pyedflib
import argparse
from datetime import datetime, timedelta
from pyedflib import highlevel
import numpy as np

def remove_sensitive_data_and_adjust_time(edf_file_path, output_file_path):
    # Open the original EDF file

    signals, signal_headers, header = pyedflib.highlevel.read_edf(edf_file_path, digital=True)
    startdate = header['startdate']
    startdate = startdate - timedelta(days=7)
    header['annotations'] = []
    header['patientname'] = ''
    header['birthdate'] = ''
    header['admin_code'] = ''
    header['patientcode'] = ''
    header['technician'] = ''
    header['recording_additional'] = ''
    header['startdate'] = startdate
    #for i, signal in enumerate(signals):

    #    signal_min = signal_headers[i]['physical_min']
    #    signal_max = signal_headers[i]['physical_max']

        # Clip the signal based on the individual physical_min and physical_max
    #    signals[i] = np.clip(signal, signal_min, signal_max)

    pyedflib.highlevel.write_edf(output_file_path, signals, signal_headers, header, digital=True)

    # to_remove = ['patientname', 'birthdate', 'admin_code','startdate','patientcode','technician','recording_additional']
    # new_values = ["xxx", '', 'test', startdate, '', '', '']
    #highlevel.anonymize_edf(edf_file_path, to_remove=to_remove, new_values=new_values, verify=False, verbose=False)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(
        description="Modify an EDF file by removing sensitive data and adjusting start time.")
    parser.add_argument('input_file', type=str, help="Path to the input EDF file.")
    parser.add_argument('output_file', type=str, help="Path to the output EDF file.")

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with the provided input and output file paths
    remove_sensitive_data_and_adjust_time(args.input_file, args.output_file)


if __name__ == '__main__':
    main()