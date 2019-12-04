"""
Utils to be shared between scripts
"""

import os
import sys
import pandas as pd
import pyreadr



def read_rds(filename, location="data"):
    """Read an RDS-format gene or transcript expression matrix into Pandas.
    Location can be data, scratch, or results.
    Returns a pandas dataframe"""
    filepath = get_filepath(filename, location)
    raw_df = pyreadr.read_r(filepath)[None]
    if raw_df.isnull().values.any():
        raise ValueError("NaN's were found in the data matrix.")
    return raw_df.set_index(raw_df.columns[0], drop=True)

def write_rds(df, filename, location="scratch"):
    """Write a pandas dataframe to RDS file in the dir specified.
    Valid locations are scratch or results."""
    filepath = get_filepath(filename, location)
    pyreadr.write_rds(filepath, df)
    # Ensure that file can be read, because pyreadr.write_rds fails silently when volume is mounted read-only
    with open(filepath, "rb"):
        pass

def get_filepath(filename, location):
    """Given location of data, scratch, or results,
    construct path to file in that location"""
    valid_locations={"data":"../../data",
        "scratch":"../../scratch",
        "results": "results"
    }
    if location not in valid_locations.keys():
        raise ValueError("location must be one of {}".format(" ".join(valid_locations.keys())))
    return os.path.join(valid_locations[location], filename)
