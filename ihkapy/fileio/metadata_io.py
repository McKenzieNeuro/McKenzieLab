import numpy as np
import logging

def parse_metadata(
        metadata_filepath
        ):
    """Opens and parses the '.edf' seizure recording metadata '.txt' file.

    Parameters
    ----------
    metadata_filepath
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

    with open("metadata.txt","r") as f:
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

        # read seizure interval data into list
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
        table_transpose = np.asarray(table).T
        seizure_intervals = {colname:data for colname,data in zip(columns,table_transpose)}
    return frontmatter,seizure_intervals
