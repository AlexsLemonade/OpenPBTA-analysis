#!/usr/bin/env python
# Author Teja Koganti (D3B)

# This script takes in the final EPN table with all the columns (in results folder  results/EPN_all.tsv), 
#	sets some rules and thresholds so samples can be  subgrouped
#	and adds the sub group name in the last "subgroup" column

## python 03-subgrouping_samples.py 
#	--final_table results/EPN_all_data.tsv 
#	--subgroup_table results/EPN_all_data_withsubgroup.tsv 

## Importing modules to be used
import argparse
import os
import sys
import numpy as np 
import pandas as pd

from matplotlib import pyplot as plt
parser = argparse.ArgumentParser()
parser.add_argument('--final_table', required = True,
                    help = 'path to the notebook with all EPN columns')
parser.add_argument('--subgroup_table', required = True,
		    help = 'output table  with subgroup')
args = parser.parse_args()

# REading in input table  
EPN_final = pd.read_csv(args.final_table, sep="\t")
# Adding a column with empty strings so subgroup name can be added
EPN_final["subgroup"] = ""

# File to write all the data to
outfile = open(args.subgroup_table, "w")

# This function takes in every row from EPN_final, 
#	subgroup name to be returned(if values are higher than threshold),
#       column values from EPN_final table  that need to be assessed for sub-grouping
#If values higher than threshold and last column empty, this will return the subgroupname
#	If there are other subgroup names already, it  appends the current subgroup name and returns
def subgroup_func(row, subgroupname, column_values):
    current_subgroup = row["subgroup"]
    for columnname, value in column_values:
        if row[columnname] > value:
            if current_subgroup == '':
                current_subgroup = subgroupname
            elif subgroupname in current_subgroup.split(","):
                pass
            else:
                current_subgroup = current_subgroup + "," + subgroupname
    return(current_subgroup)

# This function will look at consensus_focal_CN_CDKN2 and return
# ST-EPN-RELA subgroup name, if column value is "loss"
def subgroup_CDKN2Aloss(row):
     current_subgroup = row["subgroup"]
     if row["consensus_focal_CN_CDKN2"] == "loss":
             if current_subgroup == '': 
                 current_subgroup = "ST_EPN_RELA"
             elif "ST_EPN_RELA" in current_subgroup.split(","):
                 pass
             else:
                 current_subgroup = current_subgroup + "," + "ST_EPN_RELA"
     return(current_subgroup)

# This functions will look at GISTIC_focal_CN_CDKN2A and return
# ST-EPN-RELA subgroup name, if column value is less than 0
def subgroup_GISTIC_CDKN2A(row):
     current_subgroup = row["subgroup"]
     if row["GISTIC_focal_CN_CDKN2A"] < 0.0:
             if current_subgroup == '':
                 current_subgroup = "ST_EPN_RELA"
             elif "ST_EPN_RELA" in current_subgroup.split(","):
                 pass
             else:
                 current_subgroup = current_subgroup + "," + "ST_EPN_RELA"
     return(current_subgroup)


#### Looking for  ST_EPN_RELA sub-group samples
# PTEN--TAS2R1 fusion column is all zeros
# Creating  a list with tuples 
# Each tuple containin g column header  and threshold value
st_epn_rela_tests = [("C11orf95--RELA", 0), 
                     ("LTBP3--RELA", 0),
		     ("PTEN--TAS2R1",  0),
                     ("9p_loss", 0),
                     ("9q_loss", 0),
		     ("L1CAM_expr_zscore",3)]
# Calling function subgroup_func to  set  the values for last column "subgroup"		     	
EPN_final["subgroup"] = EPN_final.apply(subgroup_func,
                                        axis=1,
                                        subgroupname="ST_EPN_RELA",
                                        column_values=st_epn_rela_tests)
# Checking for CDKN2A loss 
EPN_final["subgroup"] = EPN_final.apply(subgroup_CDKN2Aloss, axis=1)
EPN_final["subgroup"] = EPN_final.apply(subgroup_GISTIC_CDKN2A, axis=1)

#### Looking for ST_EPN_YAP1 sub-group samples
st_epn_yap1_tests = [("C11orf95--YAP1", 0),
		     ("YAP1--MAMLD1", 0),
		     ("C11orf95--MAML2", 0),
		     ("YAP1--FAM118B", 0),
		     ("11q_loss",  0),
		     ("11q_gain", 0),
		     ("ARL4D_expr_zscore", 3), 
		     ("CLDN1_expr_zscore", 3)]  

EPN_final["subgroup"] = EPN_final.apply(subgroup_func,
                                        axis=1,
                                        subgroupname="ST_EPN_YAP1",
                                        column_values=st_epn_yap1_tests)

####  Looking for PT_EPN_A sub-group samples
pt_epn_a_tests = [("1q_loss", 0),
		  ("TKTL1_expr_zscore", 3),   
		  ("CXorf67_expr_zscore", 3)] 

EPN_final["subgroup"] = EPN_final.apply(subgroup_func,
                                         axis=1,
                                         subgroupname="PT_EPN_A",
                                         column_values=pt_epn_a_tests)


#### Looking for PT_EPN_B  sub-group samples
pt_epn_b_tests = [("6q_loss", 0),
		  ("6p_loss",  0),
                  ("IFT46_expr_zscore", 3),
                  ("GPBP1_expr_zscore", 3)] 

EPN_final["subgroup"] = EPN_final.apply(subgroup_func,
                                          axis=1,
                                          subgroupname="PT_EPN_B",
                                          column_values=pt_epn_b_tests)

EPN_final.to_csv(outfile, sep="\t", header=True, index=False)
outfile.close()

