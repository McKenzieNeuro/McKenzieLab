# utility methods for loading options from Options.toml
import os
import toml                     # Parameters/config file Options.toml
import warnings


def load_fio_ops_and_data_ops(options_path="./Options.toml"):
    """Load the 'fio_ops' and 'data_ops' dictionaries from options config file

    Parses the .toml config file at `options_path` into a dictionary. 
    Returns the sub-dictionaries at the files "fio" and "params.data"

    Parameters
    ----------
    options_path : str
        Path to the config file

    Returns
    -------
    dict
        A dictionary containing fio params, i.e. paths to data files.
    dict
        A dictionary containing data parameters required for data manip.
    """

    warnings.warn("Change this relative path once package is configured properly. \nWe need a more reliable way of accessing the options.toml config file")
    ops = load_ops_as_dict(options_path)
    fio_ops = ops["fio"]
    data_ops = ops["params"]["data"]
    return fio_ops,data_ops

def load_ops_as_dict(options_path="./Options.toml"):
    """Returns the Options.toml file as a dictionary

    Parses the .toml config file at `options_path` into a dictionary
    and returns it.

    Parameters
    ----------
    options_path : str
        Path to the config file
    """
    with open(options_path,"r") as f:
        ops_string = f.read()
    ops = toml.loads(ops_string)
    return ops

