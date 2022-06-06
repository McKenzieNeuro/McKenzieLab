import numpy as np
import logging
import pandas as pd

def parse_metadata(
        session_metadata_path
        ) -> (dict,dict):
    """Opens and parses the '.edf' seizure recording metadata '.txt' file.

    Parameters
    ----------
    session_metadata_path
        The path of the file where the metadata is stored.
        e.g. "/Users/steve/Documents/code/unm/data/AC75a-5 DOB 072519_TS_2020-03-23_17_30_04.txt"

    Returns
    -------
    dict 
        Frontmatter information. The key-value pairs contain information 
        such as the experimentor, the subject's (or animal's) id. 

    dict 
        Seizure intervals. The keys are column headers, the values are
        numpy arrays. Column headers are typically 
        ['Number','Start Time','End Time','Time From Start','Channel','Annotation']
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO) # DEBUG < INFO < WARNING < ERROR < CRITICAL

    frontmatter = {}
    seizure_intervals = {}

    with open(session_metadata_path,"r") as f:
        def fmt(line): return "".join(line.split("\n")).split("\t")

        # read the frontmatter into a dictionary
        line = fmt(f.readline())
        while line != [""]:
            assert len(line) == 2 , "Frontmatter format error, must contain key-value pairs seperated by tabs."
            key,value = line
            frontmatter[key] = value
            line = fmt(f.readline())

        # skip lines whitespace lines
        while line == [""]:
            line = fmt(f.readline())

        # Read seizure interval data into list
        columns = line # Number, Start Time, End Time, Time From Start, Channel, Annotation
        len_col = len(columns)
        line = fmt(f.readline())
        table = []
        while line !=[""]:
            assert len(line) == len_col , "Length of row data incompatible with number of columns."
            table.append(line)
            line = fmt(f.readline())
        print("Finished reading annotation file.")

        # transpose and put into columns dictionary
        table = np.asarray(table).T
        seizure_metadata = {colname:data for colname,data in zip(columns,table)}
    return frontmatter,seizure_metadata


def get_seizure_start_end_times(session_metadata_path : str) -> (list, list):
    """Returns the start and end times (in s) specified by the file metadata"""
    _,seizure_metadata = parse_metadata(session_metadata_path)
    # Put the metadata into a pandas dataframe, and select the columns
    df = pd.DataFrame(seizure_metadata)
    start_times = df[df['Annotation'] == "Seizure starts"]["Time From Start"]
    end_times = df[df['Annotation'] == "Seizure ends"]["Time From Start"]
    start_times = [float(i) for i in start_times]
    end_times = [float(i) for i in end_times]

    # TODO: depending on exactly how the annotation (aka metadata)  
    # files are specified and what formats are used, we may want to use
    # a more elaborate regex for selecting start and end of seizures

    # Verify there is no formatting error
    assert len(start_times) == len(end_times), "Incompatible seizure start/end times"
    for i,j in zip(start_times,end_times):
        assert i < j, "Seizure cannot end before it starts!"

    return start_times,end_times








