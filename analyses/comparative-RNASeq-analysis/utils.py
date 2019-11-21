"""
Utils to be shared between scripts
"""

import sys
import pandas as pd
import pyreadr

def read_rds(filepath):
    """Read an RDS file into Pandas.
    Takes path to file; returns pandas dataframe"""
    result = pyreadr.read_r(filepath)
    return result[None]
