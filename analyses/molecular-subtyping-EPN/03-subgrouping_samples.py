#!/usr/bin/env python
# Author Teja Koganti (D3B)

# This script takes in the final EPN table with all the columns (in results folder  results/EPN_all.tsv), 
#	sets some rules and thresholds so samples can be  subgrouped
#	and adds the sub group name in the last "subgroup" column

## python 03-subgrouping_samples.py 
#	--final_table results/EPN_all_data.tsv 
#	--subgroup_table results/EPN_all_data_withsubgroup.tsv 
#	--temp_folder_name subgrouping_samples   ## This folder will keep all the histograms to view  the outliers but can  be deleted  after review


## Importing modules to be used
import pandas as pd
import numpy as np 
import sys
import argparse
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--final_table', required = True,
                    help = 'path to the notebook with all EPN columns')
parser.add_argument('--subgroup_table', required = True,
		    help = 'output table  with subgroup')
parser.add_argument('--temp_folder_name', required = True,
                    help = 'outout table  with subgroup')
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

## There is no set thresholds for gene expression outliers, so creating histograms for each of them.  
## If gene expression Z-score is very high(outliers), that gene might be overexpressed in that sample and
## helps the sample to be subgrouped accordingly??


#### Looking for  ST_EPN_RELA sub-group samples
# PTEN--TAS2R1 fusion column is all zeros
# Creating  a list with tuples 
# Each tuple containin g column header  and threshold value
st_epn_rela_tests = [("C11orf95--RELA", 0), 
                     ("LTBP3--RELA", 0),
                     ("9p_loss", 0),
                     ("9q_loss", 0),
		     ("L1CAM_expr_zscore",3)]
# Calling function subgroup_func to  set  the values for last column "subgroup"		     	
EPN_final["subgroup"] = EPN_final.apply(subgroup_func,
                                        axis=1,
                                        subgroupname="ST_EPN_RELA",
                                        column_values=st_epn_rela_tests)

# Cannot find clear outliers from the histogram plots for RELA Z-scores 
# One Z-score above 3 
L1CAM_fig = plt.figure(figsize=(20,10))
#L1CAM_fig = plt.figure()
plt.hist(list(EPN_final["L1CAM_expr_zscore"])) 
plt.suptitle('L1CAM expression')
L1CAM_fig.savefig(args.temp_folder_name+"/L1CAM_hist.png", dpi=L1CAM_fig.dpi)




#### Looking for ST_EPN_YAP1 sub-group samples
st_epn_yap1_tests = [("C11orf95--YAP1", 0),
		     ("YAP1--MAMLD1", 0),
		     ("C11orf95--MAML2", 0),
		     ("YAP1--FAM118B", 0),
		     ("11q_loss",  0),
		     ("11q_gain", 0),
		     ("ARL4D_expr_zscore", 1.7), ## hist plot included to show threshold
		     ("CLDN1_expr_zscore", 1)]   ## hist plot included to show threshold

EPN_final["subgroup"] = EPN_final.apply(subgroup_func,
                                        axis=1,
                                        subgroupname="ST_EPN_YAP1",
                                        column_values=st_epn_yap1_tests)

# Two sample higher than 1.7, look like outliers based on histogram
ARL4D_fig = plt.figure()
plt.hist(list(EPN_final["ARL4D_expr_zscore"])) 
plt.suptitle('ARL4D expression')
ARL4D_fig.savefig(args.temp_folder_name+"/ARL4D_hist.png", dpi=ARL4D_fig.dpi)
# One sample higher than 1.0, another sample at 0.28 also look like an outlier but a cutoff value of 0.25 maybe too small?? 
CLDN1_fig = plt.figure()
plt.hist(list(EPN_final["CLDN1_expr_zscore"])) 
plt.suptitle('CLDN1 expression')
CLDN1_fig.savefig(args.temp_folder_name+"/CLDN1_hist.png", dpi=CLDN1_fig.dpi)




####  Looking for PT_EPN_A sub-group samples
pt_epn_a_tests = [("1q_loss", 0),
		  ("TKTL1_expr_zscore", 2),   #hist plot for threshold included
		  ("CXorf67_expr_zscore", 2)] #hist plot fot threshold included

EPN_final["subgroup"] = EPN_final.apply(subgroup_func,
                                         axis=1,
                                         subgroupname="PT_EPN_A",
                                         column_values=pt_epn_a_tests)

# Anything above 2 Z-score looks like an outlier
TKTL1_fig = plt.figure()
plt.hist(list(EPN_final["TKTL1_expr_zscore"])) 
plt.suptitle('TKTL1 expression')
TKTL1_fig.savefig(args.temp_folder_name+"/TKTL1_hist.png", dpi=TKTL1_fig.dpi)
# Any score above 2.0? Two samples at 2.03 and 2.07
CXorf67_fig = plt.figure()
plt.hist(list(EPN_final["CXorf67_expr_zscore"])) 
plt.suptitle('CXorf67 expression')
CXorf67_fig.savefig(args.temp_folder_name+"/CXorf67_hist.png", dpi=CXorf67_fig.dpi)



#### Looking for PT_EPN_B  sub-group samples
pt_epn_b_tests = [("6q_loss", 0),
		  ("6p_loss",  0),
                  ("IFT46_expr_zscore", 3),   #hist plot for threshold included
                  ("GPBP1_expr_zscore", 2)] #hist plot fot threshold included


# Two sample with Z-scores above 3
IFT46_fig = plt.figure()
plt.hist(list(EPN_final["IFT46_expr_zscore"])) 
plt.suptitle('IFT46 expression')
IFT46_fig.savefig(args.temp_folder_name+"/IFT46_hist.png", dpi=IFT46_fig.dpi)
# Look at GPB1 expression Z-scores  
GPBP1_fig = plt.figure()
plt.hist(list(EPN_final["GPBP1_expr_zscore"])) 
plt.suptitle('GPBP1 expression')
IFT46_fig.savefig(args.temp_folder_name+"/GPBP1_hist.png", dpi=GPBP1_fig.dpi)

EPN_final.to_csv(outfile, sep="\t", header=True, index=False)
outfile.close()


