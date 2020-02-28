import pandas as pd 
import numpy as np


sample_file = pd.read_csv("caseid_and_160samplesnames", sep="\t", header=None, names=["UUID", "sample"])
BED_filenames = pd.read_csv("UUID_andBED.txt", sep="\t", header=None, names=["UUID", "BED"])

final = sample_file.merge(BED_filenames, how='outer', on='UUID')

final.to_csv("UUID_sample_BED.txt", header=None, index=None)

